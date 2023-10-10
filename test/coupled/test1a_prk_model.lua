prk_system = prk2.PRK2System.Create
({
  initial_source = 0.0,
  time_integrator = chi_math.CrankNicolsonTimeIntegrator.Create({})
})

prk_model = chi_math.TransientNonLinearExecutioner.Create
({
  name = "prk_model",
  system = prk_system,
  solver_params =
  {
    nl_method = "LINEAR",
    --nl_method = "PJFNK",
    l_method = "preonly",
    l_rel_tol = 1.0e-5,
    pc_options =
    {
      pc_type = "lu"
    }
  },
  dt = 0.0001,
  end_time = 1.0,
  max_time_steps = -1,

  print_header = false,
  print_footer = false,
  print_nl_residual = false,
  print_l_residual = false,
  --print_timing_info = true
})

chi.SolverInfoPostProcessor.Create
({
  name = "period(s)",
  solver = prk_model,
  info = {name = "period"},
  print_on = {"ProgramExecuted" }
})
chi.SolverInfoPostProcessor.Create
({
  name = "population",
  solver = prk_model,
  info = {name = "population"},
  print_on = {"ProgramExecuted" }
})

chiSolverInitialize(prk_model)
if (main == nil) then
  chiSolverExecute(prk_model)
end