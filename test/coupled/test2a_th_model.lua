hc_mesh_handler = chiMeshHandlerCreate()

-- fuel_pitch x = 1.66 in = 4.2164e-2 m
-- fuel_pitch y = 1.51 in = 3.8354e-2 m
-- fuel OD = 3.5839e-2 m
-- Area = 6.0837E-04 m2
-- Dh = 6.7900E-02 m 
-- wetted_perimeter = 3.5839E-02 m

-- ============================================== Components
c = 1
components = {}
components[c] = piper.BoundaryComponent.Create({ name = "bcompA" })
c = c + 1

Np = 10
L = 15 * 2.54e-2
--L = 3.0
for k = 1, Np do
  components[c] = piper.SingleVolume.Create({
    name = "Pipe" .. tostring(k),
    Dh = 6.7900E-02,
    A = 6.0837E-04,
    length = L/Np,
    orientation = { polar = 0.0 }
  })
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
  datum = { 0.1, 0.1, 0.1 }
})

th_model = piper.LiquidPhysics.Create({
  name = "LiquidPhysics",
  pipe_system = pipe_system,
  fluid_name = "Water",
  initializer = { type = "StaticGravity", state = { T = 273.15 + 40.0, p = 100.0e3 } },
  dt = 0.1,
  end_time = 1.0,
  max_time_steps = -1,
  model_parameters = {
    --{comp_name = "j1", form_loss_forward = 5.0},
    { comp_name = "Pipe1", volumetric_heat_generation = 0.0 },
    { comp_name = "Pipe2", volumetric_heat_generation = 0.0 },
    { comp_name = "Pipe3", volumetric_heat_generation = 0.0 },
    { comp_name = "Pipe4", volumetric_heat_generation = 0.0 },
    { comp_name = "Pipe5", volumetric_heat_generation = 0.0 },
    { comp_name = "Pipe6", volumetric_heat_generation = 0.0 },
    { comp_name = "Pipe7", volumetric_heat_generation = 0.0 },
    { comp_name = "Pipe8", volumetric_heat_generation = 0.0 },
    { comp_name = "Pipe9", volumetric_heat_generation = 0.0 },
    { comp_name = "Pipe10", volumetric_heat_generation = 0.0 },
  },
  print_header = false,
})

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
  name = "u_Pipe10",
  solver = th_model,
  print_on = { "ProgramExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe10", var_name = "u" },
})
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "hcoeff_Pipe10",
  solver = th_model,
  print_on = { "ProgramExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe10", var_name = "hcoeff" },
})
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "Re_Pipe10",
  solver = th_model,
  print_on = { "ProgramExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe10", var_name = "Re" },
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



