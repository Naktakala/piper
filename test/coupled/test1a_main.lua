main = true
reactor_power = 1.0e4
dofile("test1a_prk_model.lua")
dofile("test1a_hc_model.lua")

chi.ExecutePostProcessors({ "element_avg_temperature" })
init_temp = chi.PostProcessorGetValue("element_avg_temperature")

population_power_factor = reactor_power
rho_feedback = 0.0
rho_added = 0.0

step_latched = false
function StepChange(event_name, params)
  if (params.solver_name ~= "prk_model") then return end
  if (params.time >= 0.2 and not step_latched) then
    rho_added = 3.25
    step_latched = true
  end
end

chi.LuaEventHook.Create
({
  name = "rho_change",
  lua_function_name = "StepChange",
  execute_on = { "SolverAdvanced" }
})

dt = 0.005
for k=1,600 do
  rho = rho_feedback + rho_added
  chiPRK2SetProperty(prk_system, {rho=rho})
  chiSolverSetProperties(prk_model, {end_time=k*dt})
  chiSolverExecute(prk_model)

  reactor_power = chi.PostProcessorGetValue("population") * population_power_factor

  chiFieldOperationExecute(field_op_heat)

  chiSolverSetProperties(hc_model, {end_time=k*dt})
  chiSolverExecute(hc_model)

  avg_temp = chi.PostProcessorGetValue("element_avg_temperature")

  rho_feedback = (-1.0/300.0) *(avg_temp - init_temp)

  chi.PrintPostProcessors({ "population", "element_avg_temperature", "period(s)"})
  chiLog(LOG_0, "reactor_power = "..string.format("%.3g", reactor_power/1e6)..
  string.format(" time=%.2f", k*dt))

end

chi.PostProcessorPrinterSetOptions({ csv_filename = "test1a.csv" })