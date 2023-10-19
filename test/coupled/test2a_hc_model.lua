--==================================================== TRIGA Fuel element
--
clad_thickness = 0.0508e-2
clad_RO = 0.5 * 3.5839e-2
zrh_RO = clad_RO - clad_thickness
zrh_RI = 0.5 * 0.635e-2
zrh_A = math.pi * (zrh_RO^2 - zrh_RI^2)
zrh_L = 3 * 12.7e-2

zrh_V = zrh_A * zrh_L
one_MW_peak_element_power = 11.0e3
one_MW_power_density = one_MW_peak_element_power/zrh_V
if (reactor_power == nil) then reactor_power = 1.0e6 end

--==================================================== Define fuel mesh
hc_mesh_handler = chiMeshHandlerCreate()
rnodes = {}
znodes = {}

Nz = 5
zmin = 0.0

dz = zrh_L / Nz
for i = 0, Nz do
  znodes[i + 1] = zmin + dz * i
end

-- 0.635cm Zirconium rod
-- 0.0508 cm SS cladding
-- 3.5839 cm Outer diameter
rnodes = { 0.0 }
rc = 1
function AddRLayer(Ro, Nlayers)
  Ri = rnodes[rc]
  dr = (Ro - Ri) / Nlayers
  r = Ri
  for _ = 1, Nlayers do
    rnodes[rc + 1] = r + dr
    rc = rc + 1
    r = r + dr
  end
end

AddRLayer(zrh_RI, 2)
AddRLayer(zrh_RO, 3)
AddRLayer(zrh_RO + clad_thickness, 1)

meshgen1 = chi_mesh.OrthogonalMeshGenerator.Create({ node_sets = { rnodes, znodes } })
chi_mesh.MeshGenerator.Execute(meshgen1)

function SetMatIDs(r, _, _, _)
  if (r > zrh_RO) then
    return 2
  elseif (r > zrh_RI) then
    return 1
  else
    return 0
  end
end

chiVolumeMesherSetProperty(MATID_FROM_LUA_FUNCTION, "SetMatIDs")

--==================================================== Define interface fields
--===================================== Heat generation
field_heat_gen = chi_physics.FieldFunctionGridBased.Create({
  name = "p_density",
  sdm_type = "FV",
  initial_value = 50.0,
  coordinate_system = "rz"
})

function HeatField(vals)
  r = vals[1]
  if (r > zrh_RI and r < zrh_RO) then
    return { (reactor_power / 1.0e6) * one_MW_power_density }
  else
    return { 0.0 };
  end
end
func_heatfield = chi_math.functions.LuaDimAToDimB.Create
({
  input_dimension = 4,
  output_dimension = 1,
  lua_function_name = "HeatField"
})
field_op_heat = chi_physics.field_operations.MultiFieldOperation.Create
({
  result_field_handle = field_heat_gen,
  function_handle = func_heatfield
})
chiFieldOperationExecute(field_op_heat)

--===================================== Fluid coupling
field_Tbulk = chi_physics.FieldFunctionGridBased.Create({
  name = "Tbulk",
  sdm_type = "FV",
  initial_value = 300.0,
  coordinate_system = "rz"
})

field_hcoeff = chi_physics.FieldFunctionGridBased.Create({
  name = "HTC",
  sdm_type = "FV",
  initial_value = 10000.0,
  coordinate_system = "rz"
})

--==================================================== Define heat conduction
ZrH_k_function = chi_math.functions.PiecewiseLinear1D.Create({
  x_values = { 0.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0 },
  y_values = { 0.01, 4.63, 10.10, 12.68, 13.96, 14.60, 14.89, 14.98, 14.95, 14.84, 14.69, 14.52 }
})

ZrH_cp_function = chi_math.functions.PiecewiseLinear1D.Create({
  x_values = { 0.0, 100.0, 300.0, 800.0 },
  y_values = { 1.0, 0.18e3, 0.32e3, 0.65e3 }
})

SS304_k_function = chi_math.functions.PiecewiseLinear1D.Create({
  x_values = {200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1600.0},
  y_values = {13.5,   17.0,  19.8,  22.5,   25.2,   28.0,   32.5}
})

SS304_cp_function = chi_math.functions.PiecewiseLinear1D.Create({
  x_values = {200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1600.0},
  y_values = {400.0, 520.0, 560.0, 580.0,  610.0,  640.0,  690.0}
})


mat_props = chi.MaterialPropertiesData.Create({
  properties = {
    hcm.ThermalConductivity.Create({ name = "k", value_function = ZrH_k_function, mat_ids = {0,1} }),
    chi.ConstantMaterialProperty.Create({ name = "density", scalar_value = 5586.0, mat_ids = {0,1} }),
    hcm.HeatCapacity.Create({ name = "heat_capacity", value_function = ZrH_cp_function, mat_ids = {0,1} }),

    hcm.ThermalConductivity.Create({ name = "k", value_function = SS304_k_function, mat_ids = {2} }),
    chi.ConstantMaterialProperty.Create({ name = "density", scalar_value = 8018.0, mat_ids = {2} }),
    hcm.HeatCapacity.Create({ name = "heat_capacity", value_function = SS304_cp_function, mat_ids = {2} })
  }
})

hc_system = chi_math.KernelSystem.Create({
  material_properties = mat_props,
  fields = {
    chi_physics.FieldFunctionGridBased.Create({
      name = "THC",
      sdm_type = "LagrangeC",
      initial_value = 300.0,
      coordinate_system = "rz"
    })
  },
  kernels = {
    { type = hcm.CoupledHeatGeneration.type, var = "THC", e_gen = "p_density", mat_ids={1} },

    { type = hcm.ThermalConductionKernel2.type, var = "THC", mat_ids={0,1}},
    { type = hcm.ThermalConductionTimeDerivative.type, var = "THC", mat_ids={0,1} },
    --
    { type = hcm.ThermalConductionKernel2.type, var = "THC", mat_ids={2}},
    { type = hcm.ThermalConductionTimeDerivative.type, var = "THC", mat_ids={2} }
  },
  bcs = {
    {
      type = hcm.CoupledConvectiveHeatFluxBC.type,
      boundaries = { "XMAX" },
      T_bulk = "Tbulk",
      convection_coefficient = "HTC",
      var = "THC"
    }
  },
  time_integrator = chi_math.CrankNicolsonTimeIntegrator.Create({}),
  --verbosity = 2
  --output_filename_base = "transient_cylTest1"
})

do_Tsurf = chi.derived_object.LayeredBoundaryAverage.Create
({
  name = "do_Tsurf",
  nodes = znodes,
  field_function = "THC",
  direction = {0.0,0.0,1.0},
  boundaries = {"XMAX"},
  parent_forward_direction = {0.0,0.0,1.0},
  parent_up_direction = {0.0,-1.0,0.0},
})

hc_model = chi_math.TransientNonLinearExecutioner.Create
--hc_model = chi_math.SteadyNonLinearExecutioner.Create
({
  name = "hc_model",
  system = hc_system,
  solver_params = {
    nl_method = "NEWTON",
    l_rel_tol = 1.0e-5,
  },

  dt = 0.001,
  end_time = 1.0,

  print_header = false,
  print_footer = false,
  print_nl_residual = false,
  print_l_residual = false,
})

chi.AggregateNodalValuePostProcessor.Create({
  name = "maxval",
  field_function = "THC",
  operation = "max",
  print_on = { "ProgramExecuted" }
})
chi.CellVolumeIntegralPostProcessor.Create({
  name = "element_power",
  field_function = "p_density",
  print_on = { "ProgramExecuted" }
})
chi.CellVolumeIntegralPostProcessor.Create({
  name = "element_avg_temperature",
  field_function = "THC",
  compute_volume_average = true,
  solvername_filter = "hc_model",
  execute_on = { "SolverExecuted" },
  print_on = { "ProgramExecuted" }
})
chi.AggregateNodalValuePostProcessor.Create({
  name = "element_max_temperature",
  field_function = "THC",
  operation = "max",
  solvername_filter = "hc_model",
  execute_on = { "SolverExecuted" },
  print_on = { "ProgramExecuted" }
})

chiSolverInitialize(hc_model)
if (main == nil) then
  chiSolverExecute(hc_model)
else
  hc_ss_model = chi_math.SteadyNonLinearExecutioner.Create
  ({
    name = "hc_model",
    system = hc_system,
    solver_params = {
      nl_method = "NEWTON",
      l_rel_tol = 1.0e-5,
    },
  })
  chiSolverInitialize(hc_ss_model)
  chiSolverExecute(hc_ss_model)
end

--chiLogPrintTimingGraph()

--if (master_export == nil) then
--  chiExportMultiFieldFunctionToVTK({ "T", "p_density"}, "ZTCoupling_test1a_HC")
--end