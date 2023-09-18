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
    { type = hcm.ThermalConductionKernel.type, k = 16.0 }
  },
  bcs = {
    chi_math.FEMDirichletBC.Create({ boundaries = { "ZMAX" }, bc_value = -1.0 }),
    chi_math.FEMDirichletBC.Create({ boundaries = { "ZMIN" }, bc_value = 10.0 })
  }
})

phys1 = hcm.HCSteadyExecutor.Create({
  conduction_system = hcsystem,
  solver_params =
  {
    --l_max_its = 10
    nl_method = "PJFNK"
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

chiExportMultiFieldFunctionToVTK({ "T" }, "test1")