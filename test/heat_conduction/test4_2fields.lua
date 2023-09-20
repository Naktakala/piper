nodes = {}
N = 100
xmin = 0.0
L = 1.0
nodes = {}

dx = L / N
for i = 0, N do
  nodes[i + 1] = xmin + dx * i
end
meshgen1 = chi_mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
chi_mesh.MeshGenerator.Execute(meshgen1)


system1 = chi_math.FEMKernelSystem.Create
({
  fields = {
    chi_physics.FieldFunctionGridBased.Create({
      name = "Te",
      sdm_type = "PWLC"
    }),
    chi_physics.FieldFunctionGridBased.Create({
      name = "T2",
      sdm_type = "PWLC"
    })
  },
  kernels = {
    { type = hcm.ThermalConductionKernel.type, var="Te", k = 16.0 },
    { type = chi_math.SinkSourceFEMKernel.type, var="Te", value = 100.0e2 },
    { type = hcm.ThermalConductionKernel.type, var="T2", k = 1.0 },
    { type = chi_math.SinkSourceFEMKernel.type, var="T2", value = 100.0e2 },
  },
  bcs = {
    {
      type = chi_math.FEMDirichletBC.type,
      boundaries = { "XMIN", "YMIN", "YMAX" },
      var="Te"
    },
    {
      type = hcm.ConvectiveHeatFluxBC.type,
      boundaries = { "XMAX" },
      T_bulk = 100.0,
      convection_coefficient = 10000.0,
      var="Te"
    },
    {
      type = chi_math.FEMDirichletBC.type,
      boundaries = { "XMIN", "XMAX", "YMIN", "YMAX" },
      var="T2"
    },
  },
  --verbosity = 2
})

phys1 = chi_math.SteadyNonLinearExecutioner.Create
({
  name = "phys1",
  system = system1,
  solver_params =
  {
    nl_method = "PJFNK",
    l_rel_tol = 1.0e-5,
  }
})

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

--############################################### PostProcessors
chi.CellVolumeIntegralPostProcessor.Create
({
  name = "avgval",
  field_function = "Te",
  compute_volume_average = true
})
chi.ExecutePostProcessors({"avgval"})

if (master_export == nil) then
  chiExportMultiFieldFunctionToVTK({ chiGetFieldFunctionHandleByName("Te"),
                                     chiGetFieldFunctionHandleByName("T2") }, "test4_2fields")
end