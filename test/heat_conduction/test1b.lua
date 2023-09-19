nodes = {}
N = 10
xmin = 0.0
L = 1.0
nodes = {}

dx = L / N
for i = 0, N do
  nodes[i + 1] = xmin + dx * i
end
meshgen1 = chi_mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes } })
chi_mesh.MeshGenerator.Execute(meshgen1)

hcsystem = hcm.HeatConductionSystem.Create({
  kernels = {
    { type = hcm.ThermalConductionKernel.type, k = 16.0 },
    { type = chi_math.SinkSourceFEMKernel.type, value = 100.0e2 }
  },
  bcs = {
    { type = chi_math.FEMDirichletBC.type,  boundaries = { "ZMAX" }, bc_value = 1.0 },
    { type = chi_math.FEMDirichletBC.type,  boundaries = { "ZMIN" }, bc_value = 1.0 }
  }
})

phys1 = hcm.HCSteadyExecutor.Create({
  conduction_system = hcsystem,
  solver_params =
  {
    --l_max_its = 10
    l_rel_tol = 1.0e-5,
    nl_method = "JFNK"
  }
})

chi.AggregateNodalValuePostProcessor.Create
({
  name = "maxval",
  field_function = "T",
  operation = "max",
  print_on = {"ProgramExecuted"}
})

chi.PostProcessorPrinterSetOptions
(
  {
    print_scalar_time_history = false
  }
)

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

if (master_export == nil) then
  chiExportMultiFieldFunctionToVTK({ "T" }, "test1")
end