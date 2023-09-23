znodes = {}
N = 100
zmin = 0.0
L = 1.0
znodes = {}

dx = L / N
for i = 0, N do
  znodes[i + 1] = zmin + dx * i
end
meshgen1 = chi_mesh.OrthogonalMeshGenerator.Create({ node_sets = { znodes, znodes } })
chi_mesh.MeshGenerator.Execute(meshgen1)


system1 = chi_math.KernelSystem.Create
({
  fields = {
    chi_physics.FieldFunctionGridBased.Create({
      name = "T",
      sdm_type = "PWLC"
    })
  },
  kernels = {
    { type = hcm.ThermalConductionKernel.type, var="T", k = 16.0 },
    { type = chi_math.SinkSourceFEMKernel.type, var="T", value = 100.0e2 }
  },
  bcs = {
    {
      type = chi_math.FEMDirichletBC.type,
      boundaries = { "XMIN", "YMIN", "YMAX" },
      var="T"
    },
    {
      type = hcm.ConvectiveHeatFluxBC.type,
      boundaries = { "XMAX" },
      T_bulk = 100.0,
      convection_coefficient = 10000.0,
      var="T"
    }
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
  field_function = "T",
  compute_volume_average = true
})
chi.ExecutePostProcessors({"avgval"})

if (master_export == nil) then
  chiExportMultiFieldFunctionToVTK({ "T" }, "test3")
end