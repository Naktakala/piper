rnodes = {}
znodes = {}

Nz = 100
zmin = 0.0
L = 100

rnodes = {0.0, 0.25, 0.5, 0.75, 1.0, 1.27, 1.35}

dz = L / Nz
for i = 0, Nz do
  znodes[i + 1] = zmin + dz * i
end
meshgen1 = chi_mesh.OrthogonalMeshGenerator.Create({ node_sets = { rnodes, znodes } })
chi_mesh.MeshGenerator.Execute(meshgen1)

system1 = chi_math.KernelSystem.Create
({
  fields = {
    chi_physics.FieldFunctionGridBased.Create({
      name = "T",
      sdm_type = "PWLC",
      initial_value = 300.0
    })
  },
  kernels = {
    {type = hcm.ThermalConductionKernel.type, var="T", k = 16.0/1e2 },
    {type = chi_math.SinkSourceFEMKernel.type, var="T", value = 100.0e2/1e2 },
    {type = chi_math.FEMTimeDerivativeKernel.type, var="T",}
  },
  bcs = {
    {
      type = hcm.ConvectiveHeatFluxBC.type,
      boundaries = { "XMAX" },
      T_bulk = 300.0,
      convection_coefficient = 10000.0/1e2,
      var="T"
    }
  },
  --time_integrator = chi_math.CrankNicolsonTimeIntegrator.Create({})
  --verbosity = 2
  output_filename_base = "transient_cylTest1"
})

phys1 = chi_math.TransientNonLinearExecutioner.Create
--phys1 = chi_math.SteadyNonLinearExecutioner.Create
({
  name = "phys1",
  system = system1,
  solver_params =
  {
    nl_method = "PJFNK",
    l_rel_tol = 1.0e-5,
  },
  --time_controls =
  --{
  --  dt = 0.1,
  --  end_time = 1.0
  --},
  timestep_controller = chi_physics.AdaptiveTimeStepController.Create
  ({
    dt = 0.5,
    end_time = 60.0,
    --iteration_window_size = 100
    dt_max = 1.0,
    max_time_steps = 10000
  })
})

chi.AggregateNodalValuePostProcessor.Create
({
  name = "maxval",
  field_function = "T",
  operation = "max",
  --print_on = {"ProgramExecuted"}
})

chi.PostProcessorPrinterSetOptions
({
  print_scalar_time_history = false,
  csv_filename = "transient_cylTest1.csv"
})

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

