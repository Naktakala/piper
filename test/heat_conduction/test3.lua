nodes = {}
N = 1000
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
    { type = chi_math.SinkSourceFEMKernel.type, value = 100.0e2 }
  },
  bcs = {
    chi_math.FEMDirichletBC.Create
    ({
      boundaries = { "XMIN", "YMIN", "YMAX" }
    }),
    hcm.ConvectiveHeatFluxBC.Create
    ({
      boundaries = { "XMAX" },
      T_bulk = 100.0,
      convection_coefficient = 10000.0
    })
  }
})

phys1 = hcm.HCSteadyExecutor.Create({
  conduction_system = hcsystem,
  solver_params =
  {
    nl_method = "PJFNK",
    --nl_method = "NEWTON",
    l_rel_tol = 1.0e-5,
    --pc_options =
    --{
    --  pc_type = "hypre",
    --  pc_hypre_type = "boomeramg",
    --  pc_hypre_boomeramg_coarsen_type = "HMIS"
    --}
  }
})

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

chiExportMultiFieldFunctionToVTK({ "T" }, "test3")
