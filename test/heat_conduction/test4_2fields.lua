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

chiVolumeMesherSetMatIDToAll(0)
lvSW = chi_mesh.RPPLogicalVolume.Create({ xmin=0.0, xmax=L/2.0, ymin=0.0, ymax=L/2.0, infz=true})
chiVolumeMesherSetProperty(MATID_FROMLOGICAL, lvSW, 1)

kSW_function = chi_math.functions.PiecewiseLinear1D.Create
({
  x_values = {0.0, 50.0, 100.0},
  y_values = {0.01, 20.0, 30.0}
})

kOther_function = chi_math.functions.PiecewiseLinear1D.Create
({
  x_values = {0.0, 100.0},
  y_values = {10.0, 10.0}
})

cp_function = chi_math.functions.PiecewiseLinear1D.Create
({
  x_values = {0.0, 50.0, 100.0},
  y_values = {5000.0, 4750.0, 4500.0}
})

mat_props = chi.MaterialPropertiesData.Create
({
  properties =
  {
    hcm.ThermalConductivity.Create({name="kSW", value_function = kSW_function}),
    hcm.ThermalConductivity.Create({name="kOther", value_function = kOther_function}),
    chi.ConstantMaterialProperty.Create({name="density", scalar_value = 1000.0}),
    hcm.HeatCapacity.Create({name="Cp", value_function = cp_function})
  }
})

tbulk = chi_physics.FieldFunctionGridBased.Create({
  name = "T_bulk",
  sdm_type = "FV",
  initial_value = 50.0
})

system1 = chi_math.KernelSystem.Create
({
  material_properties = mat_props,
  fields =
  {
    chi_physics.FieldFunctionGridBased.Create({ name = "Te", sdm_type = "LagrangeC" }),
    chi_physics.FieldFunctionGridBased.Create({ name = "Tb", sdm_type = "LagrangeC" })
  },
  kernels =
  {
    { type = hcm.ThermalConductionKernel2.type, var="Te", k_property_name = "kSW", mat_ids = {1} },
    { type = hcm.ThermalConductionKernel2.type, var="Te", k_property_name = "kOther", mat_ids = {0} },
    { type = chi_math.SinkSourceFEMKernel.type, var="Te", value = 100.0e2, mat_ids={1} },
    { type = hcm.ThermalConductionKernel2.type, var="Tb", k_property_name = "kOther" },
    { type = chi_math.SinkSourceFEMKernel.type, var="Tb", value = 10.0e2 },
  },
  bcs =
  {
    { type = chi_math.FEMDirichletBC.type, boundaries = { "XMIN", "YMIN", "YMAX", "XMAX" }, var="Te" },
    { type = chi_math.FEMDirichletBC.type, boundaries = { "XMIN", "YMIN", "YMAX" }, var="Tb" },
    { type = hcm.CoupledConvectiveHeatFluxBC.type, boundaries = { "XMAX" },
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
    --nl_method = "PJFNK",
    nl_method = "NEWTON",
    l_rel_tol = 1.0e-5,
    --nl_max_its = 1
  }
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

chiLogPrintTimingGraph()