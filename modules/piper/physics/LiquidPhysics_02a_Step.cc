#include "LiquidPhysics.h"

#include "piper/Piper.h"
#include "piper/utils/CoolPropInterface.h"
#include "piper/components/HardwareComponent.h"

#include "physics/TimeSteppers/TimeStepper.h"

#include "math/chi_math.h"

#include "chi_log.h"
#include "utils/chi_timer.h"

#include <functional>

namespace piper
{

void LiquidPhysics::Step()
{
  const double start_time = Chi::program_timer.GetTime();
  ChiLogicalErrorIf(Chi::mpi.process_count != 1, "ONLY SERIAL!!!");

  const size_t tag_jnc_time = Chi::log.GetExistingRepeatingEventTag("JNC_TIME");
  const size_t tag_vol_time = Chi::log.GetExistingRepeatingEventTag("VOL_TIME");
  const size_t tag_slv_time = Chi::log.GetExistingRepeatingEventTag("SLV_TIME");

  const std::vector<size_t> tags = {
    Chi::log.GetExistingRepeatingEventTag("TIMING0"),
    Chi::log.GetExistingRepeatingEventTag("TIMING1"),
    Chi::log.GetExistingRepeatingEventTag("TIMING2"),
    Chi::log.GetExistingRepeatingEventTag("TIMING3"),
    Chi::log.GetExistingRepeatingEventTag("TIMING4")};

  const double projected_end_time =
    timestepper_->Time() + timestepper_->TimeStepSize();
  //if (projected_end_time > timestepper_->EndTime()) dt_ = (projected_end_time - Time());

  const auto& pipe_system = *pipe_system_ptr_;
  const auto& vol_comp_ids = pipe_system.VolumeComponentIDs();
  const auto& jnc_comp_ids = pipe_system.JunctionComponentIDs();

  //============================================= Assemble momentum eqs
  Chi::log.LogEvent(tag_jnc_time, chi::ChiLog::EventType::EVENT_BEGIN);
  for (const size_t jnc_comp_id : jnc_comp_ids)
  {
    auto& jnc_model = GetComponentLiquidModel(jnc_comp_id);
    jnc_model.AssembleEquations();
  } // for junction component
  Chi::log.LogEvent(tag_jnc_time, chi::ChiLog::EventType::EVENT_END);

  //============================================= Assemble conservation of
  //                                              energy, mass, eos
  Chi::log.LogEvent(tag_vol_time, chi::ChiLog::EventType::EVENT_BEGIN);
  for (const auto& vol_comp_id : vol_comp_ids)
  {
    auto& vol_model = GetComponentLiquidModel(vol_comp_id);
    vol_model.AssembleEquations();
  } // for volume components
  Chi::log.LogEvent(tag_vol_time, chi::ChiLog::EventType::EVENT_END);

  //============================================= Init pressure system
  Chi::log.LogEvent(tag_slv_time, chi::ChiLog::EventType::EVENT_BEGIN);
  const size_t num_volumetric_components = vol_comp_ids.size();

  MatDbl A(num_volumetric_components, VecDbl(num_volumetric_components, 0.0));
  VecDbl p(num_volumetric_components, 0.0);
  VecDbl b(num_volumetric_components, 0.0);

  std::map<size_t, size_t> vol_comp_id_2_row_i_map;
  {
    size_t i = 0;
    for (const auto& vol_comp_id : vol_comp_ids)
      vol_comp_id_2_row_i_map[vol_comp_id] = i++;
  }

  //============================================= Assemble pressure system
  for (const auto& vol_comp_id : vol_comp_ids)
  {
    auto& vol_model = GetComponentLiquidModel(vol_comp_id);

    const size_t ir = vol_comp_id_2_row_i_map.at(vol_comp_id);

    const auto& eq_COE3 = vol_model.GetEquationCoefficients(2);

    b[ir] = eq_COE3.rhs_;

    for (auto [j, jnc_j_model, connection_j] : vol_model.Connections())
    {
      auto& jnc_j_liq_model = GetComponentLiquidModel(jnc_j_model.GetID());
      const auto& eq_cons_of_mom_j = jnc_j_liq_model.GetEquationCoefficients(0);

      const double a_j_COE3 = eq_COE3.coeff_sets_[0][j];

      b[ir] -= a_j_COE3 * eq_cons_of_mom_j.rhs_;

      const size_t JI = eq_cons_of_mom_j.coeff_sets_[0].size();
      for (size_t i = 0; i < JI; ++i)
      {
        const size_t comp_id = eq_cons_of_mom_j.id_maps_[0][i];
        const auto it =
          std::find(vol_comp_ids.begin(), vol_comp_ids.end(), comp_id);
        if (it != vol_comp_ids.end())
        {
          const size_t adj_comp_id = *it;
          const size_t ir_adj = vol_comp_id_2_row_i_map[adj_comp_id];

          A[ir][ir_adj] += a_j_COE3 * eq_cons_of_mom_j.coeff_sets_[0][i];
        }
        else
        {
          const auto& bndry_comp_model = component_models_[comp_id];
          const double bndry_p = bndry_comp_model->VarOld("p");
          b[ir] -= a_j_COE3 * eq_cons_of_mom_j.coeff_sets_[0][i] * bndry_p;
        }
      } // for i
    }   // for j
  }     // for volume components

  p = b;
  chi_math::GaussElimination(A, p, static_cast<int>(num_volumetric_components));
  Chi::log.LogEvent(tag_slv_time, chi::ChiLog::EventType::EVENT_END);
  Chi::log.LogEvent(tags[0], chi::ChiLog::EventType::EVENT_BEGIN);

  Chi::log.LogEvent(tags[0], chi::ChiLog::EventType::EVENT_END);

  Chi::log.LogEvent(tags[1], chi::ChiLog::EventType::EVENT_BEGIN);
  // Update pressure
  for (const auto& [vol_comp_id, ir] : vol_comp_id_2_row_i_map)
  {
    auto& vol_comp_model = *component_models_.at(vol_comp_id);
    vol_comp_model.VarNew("p") = p[vol_comp_id_2_row_i_map.at(vol_comp_id)];
  }
  Chi::log.LogEvent(tags[1], chi::ChiLog::EventType::EVENT_END);
  Chi::log.LogEvent(tags[2], chi::ChiLog::EventType::EVENT_BEGIN);
  // Update junction velocity
  for (const size_t junc_comp_id : jnc_comp_ids)
  {
    auto& jnc_model = GetComponentLiquidModel(junc_comp_id);
    const auto& jnc_connections = jnc_model.ConnectionPoints();
    const size_t I = jnc_connections.size();

    const auto& eq_COM = jnc_model.GetEquationCoefficients(0);

    double u_j_tp1 = eq_COM.rhs_;
    for (size_t i = 0; i < I; ++i)
    {
      const size_t comp_id = eq_COM.id_maps_[0][i];
      const double coeff = eq_COM.coeff_sets_[0][i];

      const auto& vol_comp_model = component_models_.at(comp_id);
      u_j_tp1 += coeff * vol_comp_model->VarNew("p");
    }

    jnc_model.VarNew("u") = u_j_tp1;
  }
  Chi::log.LogEvent(tags[2], chi::ChiLog::EventType::EVENT_END);

  Chi::log.LogEvent(tags[3], chi::ChiLog::EventType::EVENT_BEGIN);

  // Update density and internal energy + other state variables
  for (const auto& [vol_comp_id, ir] : vol_comp_id_2_row_i_map)
  {
    auto& vol_model = GetComponentLiquidModel(vol_comp_id);
    const size_t J = vol_model.Connections().size();

    const double Ai = vol_model.Area();

    const auto& eq_COMass = vol_model.GetEquationCoefficients(0);

    double rho_i_tp1 = eq_COMass.rhs_;
    double avg_flowrate = 0.0;
    for (auto [j, jnc_j_model, connection_j] : vol_model.Connections())
    {
      const double A_j = jnc_j_model.Area();

      const bool junction_j_outgoing = jnc_j_model.IsOutgoingRelToConPoint(
        connection_j.connected_comp_connection_point_id_);

      const double u_j_old = jnc_j_model.VarNew("u");

      rho_i_tp1 += eq_COMass.coeff_sets_[0][j] * u_j_old;

      const double ujj = junction_j_outgoing ? u_j_old : -u_j_old;
      const double Axj = j == 0 ? -A_j : A_j;

      avg_flowrate += Axj * ujj;
    }
    vol_model.VarNew("rho") = rho_i_tp1;
    avg_flowrate /= double(J);

    vol_model.VarNew("u") = avg_flowrate / Ai;

    const auto& eq_EOS = vol_model.GetEquationCoefficients(1);
    double rho_e_i_tp1 = eq_EOS.coeff_sets_[0][0] * rho_i_tp1 + eq_EOS.rhs_;

    const double e_i_tp1 = rho_e_i_tp1 / rho_i_tp1;
    vol_model.VarNew("e") = e_i_tp1;

    const std::vector<StateVal> state_specs = {{"p", vol_model.VarNew("p")},
                                               {"e", e_i_tp1}};

    Chi::log.LogEvent(tags[4], chi::ChiLog::EventType::EVENT_BEGIN);
    const auto state =
      EvaluateState({"rho", "e", "T", "p", "k", "mu"}, state_specs);
    Chi::log.LogEvent(tags[4], chi::ChiLog::EventType::EVENT_END);

    const std::vector<std::string> var_names = {
      "rho", "e", "T", "p", "k", "mu"};

    for (const auto& var_name : var_names)
      vol_model.VarNew(var_name) = state.at(var_name);

    // Chi::log.Log() << "T=" << vol_model.VarNew("T")
    //                << " rho=" << vol_model.VarNew("rho")
    //                << " p=" << vol_model.VarNew("p")
    //                << " e=" << vol_model.VarNew("e")
    //                << " mu=" << vol_model.VarNew("mu")
    //                << " u=" << vol_model.VarNew("u");

    // Compute Reynold's number
    const double rho = rho_i_tp1;
    const double u = vol_model.VarNew("u");
    const double Dh = vol_model.HydraulicDiameter();
    const double mu = vol_model.VarNew("mu");

    const double Re = std::fabs(rho * u * Dh / mu);

    vol_model.VarNew("Re") = Re;
  }
  Chi::log.LogEvent(tags[3], chi::ChiLog::EventType::EVENT_END);

  step_time_ = Chi::program_timer.GetTime() - start_time;
}

} // namespace piper