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
      sdm_type = "LagrangeC",
    })
  },
  kernels = {
    {type = hcm.ThermalConductionKernel.type, var="T", k = 16.0 },
    {type = chi_math.SinkSourceFEMKernel.type, var="T", value = 100.0e2 },
    {type = chi_math.FEMTimeDerivativeKernel.type, var="T",}
  },
  bcs = {
    --{
    --  type = chi_math.FEMDirichletBC.type,
    --  boundaries = { "XMIN", "YMIN", "YMAX", "XMAX" }
    --},
    {
      type = hcm.ConvectiveHeatFluxBC.type,
      boundaries = { "XMIN", "YMIN", "YMAX", "XMAX" },
      T_bulk = 100.0,
      convection_coefficient = 10000.0,
      var="T"
    }
  },
  --verbosity = 2
})

phys1 = chi_math.TransientNonLinearExecutioner.Create
({
  name = "phys1",
  system = system1,
  solver_params =
  {
    nl_method = "PJFNK",
    --nl_method = "NEWTON",
    l_rel_tol = 1.0e-5,
  },
  dt = 0.001,
  end_time = 1.0,
  max_time_steps = 100,

  print_footer = false,
  print_nl_residual = false,
  print_l_residual = false
})

chi.AggregateNodalValuePostProcessor.Create
({
  name = "maxval",
  field_function = "T",
  operation = "max",
  print_on = {"ProgramExecuted"}
})

chi.PostProcessorPrinterSetOptions
({
  print_scalar_time_history = false
})

chiSolverInitialize(phys1)
chiSolverExecute(phys1)
--for k=1,100 do
--  chiSolverStep(phys1)
--  chiSolverAdvance(phys1)
--end

if (master_export == nil) then
  chiExportMultiFieldFunctionToVTK({ "T" }, "transient_test3b")
end

chiLogPrintTimingGraph()