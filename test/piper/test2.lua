chiLog(LOG_0, "Test: Single-phase physics")

-- ============================================== Components
c = 1
components = {}
components[c] = piper.BoundaryComponent.Create({ name = "bcompA" })
c = c + 1

Np = 10
--L = 15 * 2.54e-2
L = 3.0
for k = 1, Np do
  components[c] = piper.SingleVolume.Create({
    name = "Pipe" .. tostring(k),
    Dh = 0.1,
    A = 0.25 * math.pi * 0.1 ^ 2,
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

phys1 = piper.LiquidPhysics.Create({
  name = "LiquidPhysics",
  pipe_system = pipe_system,
  fluid_name = "Water",
  initializer = { type = "StaticGravity", state = { T = 273.15 + 40.0, p = 100.0e3 } },
  dt = 0.1,
  end_time = 1.0,
  max_time_steps = -1,
  model_parameters = {
    --{comp_name = "j1", form_loss_forward = 5.0},
    { comp_name = "Pipe1", volumetric_heat_generation = 1.0e6 },
    { comp_name = "Pipe2", volumetric_heat_generation = 1.0e6 },
    { comp_name = "Pipe3", volumetric_heat_generation = 1.0e6 },
    { comp_name = "Pipe4", volumetric_heat_generation = 1.0e6 },
    { comp_name = "Pipe5", volumetric_heat_generation = 1.0e6 },
    { comp_name = "Pipe6", volumetric_heat_generation = 1.0e6 },
    { comp_name = "Pipe7", volumetric_heat_generation = 1.0e6 },
    { comp_name = "Pipe8", volumetric_heat_generation = 1.0e6 },
    { comp_name = "Pipe9", volumetric_heat_generation = 1.0e6 },
    { comp_name = "Pipe10", volumetric_heat_generation = 1.0e6 },
  }
})

pp1 = chi.SolverInfoPostProcessor.Create({
  name = "p_Pipe1", solver = phys1, print_on = { "SolverExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe1", var_name = "p" },
})
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "p_Pipe2", solver = phys1, print_on = { "SolverExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe2", var_name = "p" },
})
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "p_Pipe3", solver = phys1, print_on = { "SolverExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe3", var_name = "p" },
})
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "u_Pipe1", solver = phys1, print_on = { "SolverExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe1", var_name = "u" },
})
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "u_Pipe5", solver = phys1, print_on = { "SolverExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe5", var_name = "u" },
})
pp1 = chi.SolverInfoPostProcessor.Create({
  name = "u_Pipe10", solver = phys1, print_on = { "SolverExecuted" };
  info = { name = "ComponentVariable", component_name = "Pipe10", var_name = "u" },
})

chiLog(LOG_0, "Before init")
chiSolverInitialize(phys1)
chiSolverExecute(phys1)

chiSolverSetProperties(phys1, {end_time=2.0})

chiSolverExecute(phys1)

chiLogPrintTimingGraph()



