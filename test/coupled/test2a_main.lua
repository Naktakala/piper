main = true
dofile("test2a_prk_model.lua")
dofile("test2a_hc_model.lua")
dofile("test2a_th_model.lua")

transfer_Tbulk = chi_physics.field_operations.DerivedObjectAssign.Create({
  derived_object = do_tbulk,
  field_function = "Tbulk"
})
transfer_hcoeff = chi_physics.field_operations.DerivedObjectAssign.Create({
  derived_object = do_hcoeff,
  field_function = "HTC"
})
transfer_Tsurf = chi_physics.field_operations.DerivedObjectAssign.Create({
  derived_object = do_Tsurf,
  field_function = "T_surface"
})

chiFieldOperationExecute(transfer_Tbulk)
chiFieldOperationExecute(transfer_hcoeff)
chiFieldOperationExecute(transfer_Tsurf)

ffs_timestamp1 = { "THC", "Tbulk", "HTC" }
chiExportMultiFieldFunctionToVTK(ffs_timestamp1, "ZTimeStamp1")
ffs_timestamp2 = { "T", "hcoeff", "T_surface" }
chiExportMultiFieldFunctionToVTK(ffs_timestamp2, "ZTimeStamp2")

server = chi.NetworkServer.Create({
  port_number = 49468,
  --verbose = true
})
--os.exit()

chi.ExecutePostProcessors({ "element_avg_temperature" })
init_temp = chi.PostProcessorGetValue("element_avg_temperature")

population_power_factor = reactor_power
rho_feedback = 0.0
rho_added = 0.0
last_rho_added = rho_added

real_time = chiProgramTime()
display_old_time = real_time
display_intvl = 2.0
sixty_hertz_old_time = real_time
thirty_hertz_old_time = real_time
dt_sixty_hertz = 1.0 / 60.0
dt_thirty_hertz = 1.0 / 60.0
sim_time60Hz = real_time
sim_time30Hz = real_time

alive = true
while (alive) do
  real_time = chiProgramTime()

  if (real_time >= (sixty_hertz_old_time + dt_sixty_hertz)) then
    sim_time60Hz = sim_time60Hz + dt_sixty_hertz

    num_steps_required = math.floor(dt_sixty_hertz / 0.001)
    local delta_t = dt_sixty_hertz / num_steps_required

    --print(num_steps_required)

    for t = 1, num_steps_required do
      next_time = sixty_hertz_old_time + t * delta_t

      rho = rho_feedback + rho_added
      chiPRK2SetProperty(prk_system, { rho = rho })
      chiSolverSetProperties(prk_model, { end_time = next_time })
      chiSolverExecute(prk_model)

      reactor_power = chi.PostProcessorGetValue("population") * population_power_factor

      chiFieldOperationExecute(field_op_heat)

      chiSolverSetProperties(hc_model, { end_time = next_time })
      chiSolverExecute(hc_model)

      avg_temp = chi.PostProcessorGetValue("element_avg_temperature")

      rho_feedback = (-1.0 / 300.0) * (avg_temp - init_temp)
    end

    sixty_hertz_old_time = real_time
  end

  if (real_time >= (thirty_hertz_old_time + dt_thirty_hertz)) then
    sim_time30Hz = sim_time30Hz + dt_thirty_hertz

    chiFieldOperationExecute(transfer_Tsurf)

    chiSolverSetProperties(th_model, { end_time = sim_time30Hz,
                                       dt = dt_thirty_hertz })
    chiSolverExecute(th_model)

    chiFieldOperationExecute(transfer_Tbulk)
    chiFieldOperationExecute(transfer_hcoeff)

    thirty_hertz_old_time = real_time
  end

  if (real_time >= (display_old_time + display_intvl)) then
    display_old_time = real_time

    chi.PrintPostProcessors({
      "power",
      "element_power",
      "element_avg_temperature",
      "element_max_temperature",
      "period(s)",
      "p_Pipe1",
      "u_Pipe5",
      "hcoeff_Pipe5",
      "Re_Pipe5",
      "T_Pipe5"
    })
    chiLog(LOG_0, "reactor_power = " .. string.format("%.3g", reactor_power / 1e6) ..
      string.format(" real_time=%.2f program_time=%.2f", real_time, chiProgramTime()) ..
    string.format(" sim_time30Hz=%.2f sim_time60Hz=%.2f", sim_time30Hz, sim_time60Hz))
  end

  chiNetworkServerSynchronize(server)
end

chiNetworkServerShutdown(server)
chi.PostProcessorPrinterSetOptions({ csv_filename = "test2a.csv" })

