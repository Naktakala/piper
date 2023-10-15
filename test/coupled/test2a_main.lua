main = true
dofile("test2a_prk_model.lua")
dofile("test2a_hc_model.lua")
dofile("test2a_th_model.lua")

server = chi.NetworkServer.Create
({
  port_number = 49468,
  --verbose = true
})

chi.ExecutePostProcessors({ "element_avg_temperature" })
init_temp = chi.PostProcessorGetValue("element_avg_temperature")

population_power_factor = reactor_power
rho_feedback = 0.0
rho_added = 0.0
last_rho_added = rho_added

time = chiProgramTime()
display_old_time = time
display_intvl = 2.0
sixty_hertz_old_time = time
dt_sixty_hertz = 1.0/60.0

while (true) do
  time = chiProgramTime()
  if (time >= (sixty_hertz_old_time + dt_sixty_hertz)) then

    delta_t_required = time - sixty_hertz_old_time
    num_steps_required = math.floor(delta_t_required/0.001)
    local delta_t = delta_t_required/num_steps_required

    for t=1,num_steps_required do
      next_time = sixty_hertz_old_time + t*delta_t

      rho = rho_feedback + rho_added
      chiPRK2SetProperty(prk_system, {rho=rho})
      chiSolverSetProperties(prk_model, {end_time=next_time})
      chiSolverExecute(prk_model)

      reactor_power = chi.PostProcessorGetValue("population") * population_power_factor

      chiFieldOperationExecute(field_op_heat)

      chiSolverSetProperties(hc_model, {end_time=next_time})
      chiSolverExecute(hc_model)

      avg_temp = chi.PostProcessorGetValue("element_avg_temperature")

      rho_feedback = (-1.0/300.0) *(avg_temp - init_temp)
    end
    next_time = time + dt_sixty_hertz

    chiSolverSetProperties(th_model, {end_time=next_time,
                                      dt=dt_sixty_hertz})
    chiSolverExecute(th_model)

    sixty_hertz_old_time = time
  end

  if (time >= (display_old_time + display_intvl)) then
    display_old_time = time

    chi.PrintPostProcessors({
      "power",
      "element_avg_temperature",
      "period(s)",
      "p_Pipe1",
      "u_Pipe10",
      "hcoeff_Pipe10",
      "Re_Pipe10"
    })
    chiLog(LOG_0, "reactor_power = "..string.format("%.3g", reactor_power/1e6)..
      string.format(" time=%.2f program_time=%.2f", time, chiProgramTime()))
  end

  chiNetworkServerSynchronize(server)
end

chi.PostProcessorPrinterSetOptions({ csv_filename = "test1a.csv" })