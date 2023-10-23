main = true
dofile("test2b_prk_model.lua")
dofile("test2b_hc_model.lua")
dofile("test2b_th_model.lua")

population_power_factor = reactor_power
rho_feedback = 0.0
rho_added = 0.0

real_time = chiProgramTime()
display_old_time = real_time
display_intvl = 2.0
sixty_hertz_old_time = real_time
thirty_hertz_old_time = real_time
dt_sixty_hertz = 1.0 / 60.0
dt_thirty_hertz = 1.0 / 60.0
sim_time60Hz = real_time
sim_time30Hz = real_time

function ComputeRho() return rho_feedback + rho_added end
rho = ComputeRho()

chi.PostProcessorPrinterSetOptions({
  print_scalar_time_history = false,
  csv_filename = "test2a.csv"
})

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

ffs_timestamp1 = { "THC", "Tbulk", "HTC", "p_density" }
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

alive = true -- Variable normally set to false by HTTP client
while (alive) do
  real_time = chiProgramTime()

  if (real_time >= (sixty_hertz_old_time + dt_sixty_hertz)) then
    sim_time60Hz = sim_time60Hz + dt_sixty_hertz
    sixty_hertz_old_time = real_time

    chiTimingSectionBegin("PRKHC")

    rho = ComputeRho()
    local prkhc_dt = 0.005
    if (rho > 0.95) then prkhc_dt = 0.0005 end
    num_steps_required = math.floor(dt_sixty_hertz / prkhc_dt)
    local delta_t = dt_sixty_hertz / num_steps_required

    -- Substepping the PRK+HC
    for t = 1, num_steps_required do
      next_time = sixty_hertz_old_time + t * delta_t

      rho = ComputeRho()
      chiPRK2SetProperty(prk_system, { rho = rho })
      chiSolverSetProperties(prk_model, { end_time = next_time, dt = prkhc_dt })
      chiSolverExecute(prk_model)

      reactor_power = chi.PostProcessorGetValue("population") * population_power_factor

      chiFieldOperationExecute(field_op_heat)

      chiSolverSetProperties(hc_model, { end_time = next_time, dt = prkhc_dt })
      chiSolverExecute(hc_model)

      avg_temp = chi.PostProcessorGetValue("element_avg_temperature")

      rho_feedback = (-2.0 / 300.0) * (avg_temp - init_temp)
    end

    chiTimingSectionEnd("PRKHC")
  end

  if (real_time >= (thirty_hertz_old_time + dt_thirty_hertz)) then
    sim_time30Hz = sim_time30Hz + dt_thirty_hertz
    thirty_hertz_old_time = real_time

    chiTimingSectionBegin("TH")
    chiFieldOperationExecute(transfer_Tsurf)

    chiSolverSetProperties(th_model, { end_time = sim_time30Hz,
                                       dt = dt_thirty_hertz })
    chiSolverExecute(th_model)

    chiFieldOperationExecute(transfer_Tbulk)
    chiFieldOperationExecute(transfer_hcoeff)
    chiTimingSectionEnd("TH")
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

    t1_count, t1_total, t1_avg = chiTimingSectionValues("PRKHC")
    chiLog(LOG_0, "t_PRKHC: " .. tostring(t1_count) .. " " .. string.format("T=%10.3g ms Avg=%.3g ms", t1_total, t1_avg))
    chiTimingSectionReset("PRKHC")

    t2_count, t2_total, t2_avg = chiTimingSectionValues("TH")
    chiLog(LOG_0, "t_TH   : " .. tostring(t2_count) .. " " .. string.format("T=%10.3g ms Avg=%.3g ms", t2_total, t2_avg))
    chiTimingSectionReset("TH")
  end

  chiNetworkServerSynchronize(server)
end

chiNetworkServerShutdown(server)

