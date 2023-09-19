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
    {type = hcm.ThermalConductionKernel.type, k = 16.0 },
    {type = chi_math.SinkSourceFEMKernel.type, value = 100.0e2 },
    {type = chi_math.FEMTimeDerivativeKernel.type,}
  },
  bcs = {
    {
      type = chi_math.FEMDirichletBC.type,
      boundaries = { "XMIN", "YMIN", "YMAX", "XMAX" }
    },
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

chi.AggregateNodalValuePostProcessor.Create
({
  name = "maxval",
  field_function = "T",
  operation = "max",
  print_on = {"ProgramExecuted"}
})

chi.PostProcessorPrinterSetOptions({ print_scalar_time_history = false })

chiSolverInitialize(phys1)
for k=1,2 do
  chiSolverStep(phys1)
  chiSolverAdvance(phys1)
end


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
if (master_export == nil) then
  chiExportMultiFieldFunctionToVTK({ "T" }, "transient_test3")
end