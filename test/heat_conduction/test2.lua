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

hcsystem = hcm.HeatConductionSystem.Create({
  kernels = {
    { type = hcm.ThermalConductionKernel.type, k = 16.0 },
    --chi_math.SinkSourceFEMKernel.Create({ value = 100.0e2 })
  },
  bcs = {
    { type = chi_math.FEMDirichletBC.type, boundaries = { "XMIN" }, bc_value = -1.0 },
    { type = chi_math.FEMDirichletBC.type, boundaries = { "XMAX" }, bc_value = 10.0 },
    { type = chi_math.FEMDirichletBC.type, boundaries = { "YMIN" }, bc_value = 20.0 },
    { type = chi_math.FEMDirichletBC.type, boundaries = { "YMAX" }, bc_value = 30.0 }
  }
})

phys1 = hcm.HCSteadyExecutor.Create({
  conduction_system = hcsystem,
  solver_params =
  {
    nl_method = "PJFNK",
    l_rel_tol = 1.0e-5
    --nl_method = "NEWTON"
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
  chiExportMultiFieldFunctionToVTK({ "T" }, "test2")
end