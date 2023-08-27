chiLog(LOG_0, "Test: Single-phase physics")

-- ============================================== Components
bcompA = piper.BoundaryComponent.Create({name="bcompA"})

pipe1 =
piper.SingleVolume.Create
({
  name="Pipe1",
  Dh = 0.1,
  A = 0.25 * math.pi * 0.1^2,
  length = 1.0,
  orientation={polar=0.0}
})

pipe2 =
piper.SingleVolume.Create({
  name="Pipe2",
  Dh = 0.1,
  A = 0.25 * math.pi * 0.1^2,
  length = 1.0,
  orientation={polar=0.0}
})

pipe3 =
piper.SingleVolume.Create
({
  name="Pipe3",
  Dh = 0.1,
  A = 0.25 * math.pi * 0.1^2,
  length = 1.0,
  orientation={polar=0.0}
})

bcompB = piper.BoundaryComponent.Create({name="bcompB"})

-- ============================================== Junctions

j0 = piper.SingleJunction.Create({ name="jA", from={"bcompA","to_or_from"}, to={"Pipe1","inlet"} })
j1 = piper.SingleJunction.Create({ name="j1", from={"Pipe1","outlet"}, to={"Pipe2","inlet"} })
j2 = piper.SingleJunction.Create({ name="j2", from={"Pipe2","outlet"}, to={"Pipe3", "inlet"} })
j3 = piper.SingleJunction.Create({ name="jB", from={"Pipe3","outlet"}, to={"bcompB","to_or_from"} })

-- ============================================== Simulation setup
phys1 = piper.Piper.Create
({
  name="ShortPiper",
  components = {bcompA, pipe1, pipe2, pipe3, bcompB},
  connections = {j0, j1, j2, j3},
  --print_nodalization = true,
  datum = {0.1, 0.1, 0.1},
  physics = piper.LiquidPhysics.Create
  ({
    fluid_name = "Water",
    initializer = { type = "StaticGravity", state={T=273.15+40.0, p=100.0e3} },
    dt = 0.1,
    end_time = 100.0,
    max_time_steps = -1,
    model_parameters =
    {
      {comp_name = "j1", form_loss_forward = 5.0},
      {comp_name = "Pipe1", volumetric_heat_generation = 1.0e6},
      {comp_name = "Pipe2", volumetric_heat_generation = 1.0e6},
      {comp_name = "Pipe3", volumetric_heat_generation = 1.0e6}
    }
  })
})
chiLog(LOG_0, "Before init")
chiSolverInitialize(phys1)
chiSolverExecute(phys1)

--for k=1,2000 do
--  print("*** "..tostring(k).." t="..tostring(k*0.1))
--  chiSolverStep(phys1)
--end

