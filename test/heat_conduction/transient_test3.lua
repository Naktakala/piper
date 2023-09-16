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
    hcm.ThermalConductionKernel.Create({ k = 16.0 }),
    chi_math.SinkSourceFEMKernel.Create({ value = 100.0e2 }),
    chi_math.FEMTimeDerivativeKernel.Create({})
  },
  bcs = {
    chi_math.FEMDirichletBC.Create
    ({
      boundaries = { "XMIN", "YMIN", "YMAX", "XMAX" }
    }),
    --hcm.ConvectiveHeatFluxBC.Create
    --({
    --  boundaries = { "XMAX" },
    --  T_bulk = 100.0,
    --  convection_coefficient = 10000.0
    --})
  }
})

phys1 = hcm.HCTransientExecutor.Create({
  conduction_system = hcsystem,
  solver_params =
  {
    nl_method = "PJFNK",
    l_rel_tol = 1.0e-5,
    nl_rel_tol = 1.0e-8,
    nl_abs_tol = 1.0e-8
  },
  time_controls =
  {
    dt = 0.001,
    end_time = 1.0
  },
  time_integrator = chi_math.ImplicitEulerTimeIntegrator.Create({})
})

chiSolverInitialize(phys1)
chiSolverStep(phys1)
chiSolverAdvance(phys1)
chiSolverStep(phys1)
chiSolverAdvance(phys1)

--phys1 = hcm.HCSteadyExecutor.Create({
--  conduction_system = hcsystem,
--  solver_params =
--  {
--    nl_method = "PJFNK",
--    --nl_method = "NEWTON",
--    l_rel_tol = 1.0e-5
--  }
--})
--
--chiSolverInitialize(phys1)
--chiSolverExecute(phys1)
--
chiExportMultiFieldFunctionToVTK({ "T" }, "transient_test3")
