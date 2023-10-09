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

k_function = chi_math.functions.PiecewiseLinear1D.Create
({
  x_values = {0.0,200.0,250.0,300.0,350.0,400.0,450.0,500.0,550.0,600.0,650.0,700.0},
  y_values = {0.01,4.63,10.10,12.68,13.96,14.60,14.89,14.98,14.95,14.84,14.69,14.52}
  --x_values = {0.0, 1000.0},
  --y_values = {12.0,12.0}
})

cp_function = chi_math.functions.PiecewiseLinear1D.Create
({
  x_values = {0.0,  100.0,  300.0, 800.0},
  y_values = {1.0, 0.18e3, 0.32e3, 0.65e3}
  --x_values = {0.0, 1000.0},
  --y_values = {500.0, 500.0}
})

mat_props = chi.MaterialPropertiesData.Create
({
  properties =
  {
    hcm.ThermalConductivity.Create({name="k", value_function = k_function}),
    chi.ConstantMaterialProperty.Create({name="density", scalar_value = 5586.0}),
    hcm.HeatCapacity.Create({name="heat_capacity", value_function = cp_function})
  }
})

system1 = chi_math.KernelSystem.Create
({
  material_properties = mat_props,
  fields = {
    chi_physics.FieldFunctionGridBased.Create({
      name = "T",
      sdm_type = "LagrangeC",
      initial_value = 300.0
    })
  },
  kernels = {
    {type = hcm.ThermalConductionKernel2.type, var="T"},
    {type = chi_math.SinkSourceFEMKernel.type, var="T", value = 100.0e4 },
    {type = hcm.ThermalConductionTimeDerivative.type, var="T",}
  },
  bcs = {
    {
      type = chi_math.FEMDirichletBC.type,
      boundaries = { "XMIN", "YMIN", "YMAX", "XMAX" },
      var="T",
      bc_value = 293.15
    },
  },
  --time_integrator = chi_math.CrankNicolsonTimeIntegrator.Create({})
  --verbosity = 2
})

phys1 = chi_math.TransientNonLinearExecutioner.Create
--phys1 = chi_math.SteadyNonLinearExecutioner.Create
({
  name = "phys1",
  system = system1,
  solver_params =
  {
    --nl_method = "PJFNK",
    nl_method = "NEWTON",
    l_rel_tol = 1.0e-5,
    --pc_options =
    --{
    --  pc_type = "hypre",
    --  pc_hypre_type = "boomeramg",
    --  pc_hypre_boomeramg_coarsen_type = "HMIS"
    --}
  },
  dt = 20.0,
  end_time = 20000.0,
  --max_time_steps = 20
})

chi.AggregateNodalValuePostProcessor.Create
({
  name = "maxval",
  field_function = "T",
  operation = "max",
  --print_on = {"ProgramExecuted"}
})

--chi.PostProcessorPrinterSetOptions({ print_scalar_time_history = false })

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

if (master_export == nil) then
  chiExportMultiFieldFunctionToVTK({ "T" }, "transient_test4_hc")
end

chiLogPrintTimingGraph()