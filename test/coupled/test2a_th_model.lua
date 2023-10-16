hc_mesh_handler = chiMeshHandlerCreate()

-- fuel_pitch x = 1.66 in = 4.2164e-2 m
-- fuel_pitch y = 1.51 in = 3.8354e-2 m
-- fuel OD = 3.5839e-2 m
-- Area = 6.0837E-04 m2
-- Dh = 2.1613E-02 m
-- wetted_perimeter = 1.1259E-01 m

-- ============================================== Components
c = 1
components = {}
components[c] = piper.BoundaryComponent.Create({ name = "bcompA" })
c = c + 1

nodes = {}
Np = 5
L = 15 * 2.54e-2
--L = 3.0
nodes[1] = 0.0
for k = 1, Np do
  components[c] = piper.SingleVolume.Create({
    name = "Pipe" .. tostring(k),
    Dh = 2.1613E-02,
    A = 6.0837E-04,
    length = L/Np,
    orientation = { polar = 0.0 }
  })
  nodes[k+1] = k*L/Np
  c = c + 1
end

components[c] = piper.BoundaryComponent.Create({ name = "bcompB" })
c = c + 1
-- ============================================== Junctions

junctions = {}
for j = 0, Np do
  if (j == 0) then
    name = "jA"
    from = { "bcompA", "to_or_from" }
    to = { "Pipe1", "inlet" }
  elseif (j == Np) then
    name = "jB"
    from = { "Pipe" .. tostring(Np), "outlet" }
    to = { "bcompB", "to_or_from" }
  else
    name = "j" .. tostring(j)
    from = { "Pipe" .. tostring(j), "outlet" }
    to = { "Pipe" .. tostring(j + 1), "inlet" }
  end
  junctions[j] = piper.SingleJunction.Create({
    name = name,
    from = from,
    to = to
  })
end

-- ============================================== Simulation setup
pipe_system = piper.Piper.Create({
  name = "ShortPiper",
  components = components,
  connections = junctions,
  --print_nodalization = true,
  datum = { 0.0, 0.0, 0.0 }
})

field_Tsurf = chi_physics.FieldFunctionGridBased.Create({
  name = "T_surface",
  sdm_type = "FV",
  initial_value = 273.15 + 40.0
})

heat_generation_models =
{
  {
    type = piper.ConstantHeatGeneration.type,
    value = 7.5069E+07
  }
}

heat_flux_models =
{
  {
    type = piper.ConvectionHeatFlux.type,
    --surface_temperature_value = 350.0,
    surface_temperature_field = "T_surface",
    wetted_perimeter = 1.1259E-01
  }
}

th_model = piper.LiquidPhysics.Create({
  name = "LiquidPhysics",
  pipe_system = pipe_system,
  fluid_name = "Water",
  initializer = { type = "StaticGravity", state = { T = 300.0, p = 100.0e3 } },
  dt = 0.1,
  end_time = 1.0,
  max_time_steps = -1,
  model_parameters = {
    --{comp_name = "j1", form_loss_forward = 5.0},

    { comp_name = "Pipe1",  heat_flux_models = heat_flux_models },
    { comp_name = "Pipe2",  heat_flux_models = heat_flux_models },
    { comp_name = "Pipe3",  heat_flux_models = heat_flux_models },
    { comp_name = "Pipe4",  heat_flux_models = heat_flux_models },
    { comp_name = "Pipe5",  heat_flux_models = heat_flux_models },
    { comp_name = "Pipe6",  heat_flux_models = heat_flux_models },
    { comp_name = "Pipe7",  heat_flux_models = heat_flux_models },
    { comp_name = "Pipe8",  heat_flux_models = heat_flux_models },
    { comp_name = "Pipe9",  heat_flux_models = heat_flux_models },
    { comp_name = "Pipe10", heat_flux_models = heat_flux_models },

    --{ comp_name = "Pipe1",  heat_generation_models = heat_generation_models },
    --{ comp_name = "Pipe2",  heat_generation_models = heat_generation_models },
    --{ comp_name = "Pipe3",  heat_generation_models = heat_generation_models },
    --{ comp_name = "Pipe4",  heat_generation_models = heat_generation_models },
    --{ comp_name = "Pipe5",  heat_generation_models = heat_generation_models },
    --{ comp_name = "Pipe6",  heat_generation_models = heat_generation_models },
    --{ comp_name = "Pipe7",  heat_generation_models = heat_generation_models },
    --{ comp_name = "Pipe8",  heat_generation_models = heat_generation_models },
    --{ comp_name = "Pipe9",  heat_generation_models = heat_generation_models },
    --{ comp_name = "Pipe10", heat_generation_models = heat_generation_models },
  },
  print_header = false,
})

--========================================== Derived objects
do_tbulk = chi.derived_object.LayeredVolumeAverage.Create
({
  name = "do_tbulk",
  nodes = nodes,
  field_function = "T",
  direction = {0.0, 1.0, 0.0},
  parent_forward_direction = {0.0,0.0,-1.0},
  parent_up_direction = {0.0,1.0,0.0},
})
do_hcoeff = chi.derived_object.LayeredVolumeAverage.Create
({
  name = "do_hcoeff",
  nodes = nodes,
  field_function = "hcoeff",
  direction = {0.0, 1.0, 0.0},
  parent_forward_direction = {0.0,0.0,-1.0},
  parent_up_direction = {0.0,1.0,0.0},
})

--========================================== Post-Processors
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "p_Pipe1",
  solver = th_model,
  print_on = { "ProgramExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe1", var_name = "p" },
})
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "p_Pipe2",
  solver = th_model,
  print_on = { "ProgramExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe2", var_name = "p" },
})
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "p_Pipe3",
  solver = th_model,
  print_on = { "ProgramExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe3", var_name = "p" },
})
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "u_Pipe1",
  solver = th_model,
  print_on = { "ProgramExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe1", var_name = "u" },
})
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "u_Pipe5",
  solver = th_model,
  print_on = { "ProgramExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe5", var_name = "u" },
})
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "u_Pipe5",
  solver = th_model,
  print_on = { "ProgramExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe5", var_name = "u" },
})
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "hcoeff_Pipe5",
  solver = th_model,
  print_on = { "ProgramExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe5", var_name = "hcoeff" },
})
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "Re_Pipe5",
  solver = th_model,
  print_on = { "ProgramExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe5", var_name = "Re" },
})
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "T_Pipe5",
  solver = th_model,
  print_on = { "ProgramExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe5", var_name = "T" },
})

chiLog(LOG_0, "Before init")
chiSolverInitialize(th_model)
--chiSolverExecute(th_model)

chiSolverSetProperties(th_model, {end_time=2.0})

--chiSolverExecute(th_model)
if (main == nil) then
  chiSolverExecute(th_model)
end

--chiLogPrintTimingGraph()



