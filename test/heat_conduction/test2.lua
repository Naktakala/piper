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
    hcm.ThermalConductionKernel.Create({ k = 16.0 }),
    --chi_math.SinkSourceFEMKernel.Create({ value = 100.0e2 })
  },
  bcs = {
    chi_math.FEMDirichletBC.Create({ boundaries = { "XMIN" }, bc_value = -1.0 }),
    chi_math.FEMDirichletBC.Create({ boundaries = { "XMAX" }, bc_value = 10.0 }),
    chi_math.FEMDirichletBC.Create({ boundaries = { "YMIN" }, bc_value = 20.0 }),
    chi_math.FEMDirichletBC.Create({ boundaries = { "YMAX" }, bc_value = 30.0 })
  }
})

phys1 = hcm.HeatConductionSteadyStateExecutor.Create({
  conduction_system = hcsystem,
  solver_params =
  {
    nl_method = "PJFNK"
    --nl_method = "NEWTON"
  }
})

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

chiExportMultiFieldFunctionToVTK({ "T" }, "test2")
