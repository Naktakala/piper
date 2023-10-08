system1 = prk2.PRK2System.Create
({
  initial_source = 0.0,
  time_integrator = chi_math.CrankNicolsonTimeIntegrator.Create({})
})

phys1 = chi_math.TransientNonLinearExecutioner.Create
({
  name = "prk_system",
  system = system1,
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
  dt = 0.01,
  end_time = 1.0,
  max_time_steps = -1,

  print_header = false,
  print_footer = false,
  print_nl_residual = false,
  print_l_residual = false,
  --print_timing_info = true
})

step_latched = false
function StepChange(event_name, params)
  if (params.time >= 0.03 and not step_latched) then
    chiPRK2SetProperty(system1, {rho=0.21})
    step_latched = true
  end
end

chi.LuaEventHook.Create
({
  name = "rho_change",
  lua_function_name = "StepChange",
  execute_on = { "SolverAdvanced" }
})

chi.SolverInfoPostProcessor.Create
({
  name = "period(s)",
  solver = phys1,
  info = {name = "period"},
  print_on = {"ProgramExecuted" }
})
chi.SolverInfoPostProcessor.Create
({
  name = "population",
  solver = phys1,
  info = {name = "population"},
  print_on = {"ProgramExecuted" }
})

chiSolverInitialize(phys1)
chiSolverExecute(phys1)