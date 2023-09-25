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

func1 = chi_math.functions.PiecewiseLinear1D.Create
({
  x_values = {0.0, 50.0, 100.0},
  y_values = {0.01, 20.0, 30.0}
})

mat_props = chi.MaterialPropertiesData.Create
({
  properties =
  {
    --chi.ConstantMaterialProperty.Create({name = "k", scalar_value = 10.0}),
    --hcm.ThermalConductivity.Create({name="k", constant_value = 8.0})
    hcm.ThermalConductivity.Create({name="k", value_function = func1})
  }
})

tbulk = chi_physics.FieldFunctionGridBased.Create({
  name = "T_bulk",
  sdm_type = "FiniteVolume",
  pwl_allow_lagrange = true,
  initial_value = 50.0
})

system1 = chi_math.KernelSystem.Create
({
  material_properties = mat_props,
  fields = {
    chi_physics.FieldFunctionGridBased.Create({
      name = "Te",
      sdm_type = "PWLC",
      pwl_allow_lagrange = true
    }),
    chi_physics.FieldFunctionGridBased.Create({
      name = "Tb",
      sdm_type = "PWLC",
      pwl_allow_lagrange = true
    })
  },
  kernels = {
    --{ type = hcm.ThermalConductionKernel.type, var="Te", k=10.0 },
    { type = hcm.ThermalConductionKernel2.type, var="Te" },
    { type = chi_math.SinkSourceFEMKernel.type, var="Te", value = 100.0e2 },
    { type = hcm.ThermalConductionKernel.type, var="Tb", k = 1.0 },
    --{ type = hcm.ThermalConductionKernel2.type, var="Tb" },
    { type = chi_math.SinkSourceFEMKernel.type, var="Tb", value = 10.0e2 },
  },
  bcs = {
    {
      type = chi_math.FEMDirichletBC.type,
      boundaries = { "XMIN", "YMIN", "YMAX", "XMAX" },
      var="Te"
    },
    --{
    --  type = chi_math.FEMDirichletBC.type,
    --  boundaries = { "XMIN", "YMIN", "YMAX", "XMAX" },
    --  var="Tb"
    --},
    {
      type = hcm.ConvectiveHeatFluxBC.type,
      boundaries = { "XMIN", "YMIN", "YMAX" },
      T_bulk = 80.0,
      convection_coefficient = 10000.0,
      var="Tb"
    },
    --{
    --  type = hcm.ConvectiveHeatFluxBC.type,
    --  boundaries = { "XMAX" },
    --  T_bulk = 100.0,
    --  convection_coefficient = 10000.0,
    --  var="Tb"
    --},
    {
      type = hcm.CoupledConvectiveHeatFluxBC.type,
      boundaries = { "XMAX" },
      T_bulk = "T_bulk",
      convection_coefficient = 10000.0,
      var="Tb"
    },
  },
  --verbosity = 2,
  output_filename_base = "test4_2fields"
})

phys1 = chi_math.SteadyNonLinearExecutioner.Create
({
  name = "phys1",
  system = system1,
  solver_params =
  {
    nl_method = "PJFNK",
    l_rel_tol = 1.0e-5,
    --nl_max_its = 1
  },
  print_timing_info = true
})

--############################################### PostProcessors
chi.AggregateNodalValuePostProcessor.Create
({
  name = "maxval",
  field_function = "Te",
  --compute_volume_average = true
  operation = "max"
})

chi.PostProcessorPrinterSetOptions({ print_scalar_time_history = false })

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

