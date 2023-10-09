#include "LiquidPhysics.h"

#include "piper/Piper.h"
#include "piper/components/HardwareComponent.h"

#include "physics/TimeSteppers/TimeStepper.h"

#include "math/chi_math.h"

#include "mesh/Cell/cell.h"

#include "chi_log.h"
#include "utils/chi_timer.h"

#include <functional>

#define scint64_t static_cast<int64_t>

namespace piper
{

void LiquidPhysics::Step()
{
  const double start_time = Chi::program_timer.GetTime();

  Chi::log.Log() << "Solver \"" + TextName() + "\" " +
                      timestepper_->StringTimeInfo(false);

  auto& t_solve = Chi::log.CreateOrGetTimingBlock("LiquidPhysics::Step");
  auto& t_junctions = Chi::log.CreateOrGetTimingBlock("Junction assembly time",
                                                      "LiquidPhysics::Step");
  auto& t_volumes = Chi::log.CreateOrGetTimingBlock("Volume assembly time",
                                                    "LiquidPhysics::Step");
  auto& t_psystem = Chi::log.CreateOrGetTimingBlock("Pressure system solve",
                                                    "LiquidPhysics::Step");
  t_solve.TimeSectionBegin();

  const auto& pipe_system = *pipe_system_ptr_;
  const auto& vol_comp_ids = pipe_system.VolumeComponentIDs();
  const auto& jnc_comp_ids = pipe_system.JunctionComponentIDs();

  //============================================= Assemble momentum eqs
  t_junctions.TimeSectionBegin();
  for (const size_t jnc_comp_id : jnc_comp_ids)
  {
    auto& jnc_model = GetComponentLiquidModel(jnc_comp_id);
    jnc_model.AssembleEquations();
  } // for junction component
  t_junctions.TimeSectionEnd();

  //============================================= Assemble conservation of
  //                                              energy, mass, eos
  t_volumes.TimeSectionBegin();
  for (const auto& vol_comp_id : vol_comp_ids)
  {
    auto& vol_model = GetComponentLiquidModel(vol_comp_id);
    const auto& cell = *vol_model.GetCellPtr();
    if (cell.partition_id_ != Chi::mpi.location_id) continue;
    vol_model.AssembleEquations();
  } // for volume components
  t_volumes.TimeSectionEnd();

  //============================================= Init pressure system
  t_psystem.TimeSectionBegin();
  const size_t num_volumetric_components = vol_comp_ids.size();

  MatZeroEntries(A_);
  VecSet(b_, 0.0);
  for (const auto& vol_comp_id : vol_comp_ids)
  {
    auto& vol_model = GetComponentLiquidModel(vol_comp_id);
    const auto& cell = *vol_model.GetCellPtr();
    if (cell.partition_id_ != Chi::mpi.location_id) continue;

    const int64_t ir = scint64_t(cell.global_id_);

    const auto& eq_COE3 = vol_model.GetEquationCoefficients(2);

    VecSetValue(b_, ir, eq_COE3.rhs_, ADD_VALUES);

    for (auto [j, jnc_j_model, connection_j] : vol_model.Connections())
    {
      auto& jnc_j_liq_model = GetComponentLiquidModel(jnc_j_model.GetID());
      const auto& eq_cons_of_mom_j = jnc_j_liq_model.GetEquationCoefficients(0);

      const double a_j_COE3 = eq_COE3.coeff_sets_[0][j];

      VecSetValue(b_, ir, -a_j_COE3 * eq_cons_of_mom_j.rhs_, ADD_VALUES);

      const size_t JI = eq_cons_of_mom_j.coeff_sets_[0].size();
      for (size_t i = 0; i < JI; ++i)
      {
        const size_t comp_id = eq_cons_of_mom_j.id_maps_[0][i];
        const auto it =
          std::find(vol_comp_ids.begin(), vol_comp_ids.end(), comp_id);
        if (it != vol_comp_ids.end())
        {
          const auto& adj_model = GetComponentLiquidModel(comp_id);
          const auto& adj_cell = *adj_model.GetCellPtr();
          const int64_t ir_adj = scint64_t(adj_cell.global_id_);

          // A[ir][ir_adj] += a_j_COE3 * eq_cons_of_mom_j.coeff_sets_[0][i];
          MatSetValue(A_,
                      ir,
                      ir_adj,
                      a_j_COE3 * eq_cons_of_mom_j.coeff_sets_[0][i],
                      ADD_VALUES);
        }
        else
        {
          const auto& bndry_comp_model = component_models_[comp_id];
          const double bndry_p = bndry_comp_model->VarOld("p");
          // b[ir] -= a_j_COE3 * eq_cons_of_mom_j.coeff_sets_[0][i] * bndry_p;
          VecSetValue(b_,
                      ir,
                      -a_j_COE3 * eq_cons_of_mom_j.coeff_sets_[0][i] * bndry_p,
                      ADD_VALUES);
        }
      } // for i
    }   // for j
  }     // for volume components

  VecAssemblyBegin(b_);
  VecAssemblyEnd(b_);
  MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);

  KSPSolve(pressure_solver_setup.ksp, b_, x_);

  pressure_vector_->CopyLocalValues(x_);
  pressure_vector_->CommunicateGhostEntries();

  // Update pressure
  VecDbl p(num_volumetric_components, 0.0);
  for (size_t vol_comp_id : pipe_system_ptr_->VolumeComponentIDs())
  {
    auto& vol_comp_model = GetComponentLiquidModel(vol_comp_id);
    const auto& cell = *vol_comp_model.GetCellPtr();

    const int64_t map_id =
      pressure_vector_->MapGhostToLocal(scint64_t(cell.global_id_));

    vol_comp_model.VarNew("p") = (*pressure_vector_)[map_id];
  }
  t_psystem.TimeSectionEnd();

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

  // Update density and internal energy + other state variables
  // for (const auto& [vol_comp_id, ir] : vol_comp_id_2_row_i_map)
  for (size_t vol_comp_id : pipe_system_ptr_->VolumeComponentIDs())
  {
    auto& vol_model = GetComponentLiquidModel(vol_comp_id);
    const auto& cell = *vol_model.GetCellPtr();
    if (cell.partition_id_ != Chi::mpi.location_id) continue;

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

    const std::vector<std::string> var_names = {
      "rho", "e", "T", "p", "k", "mu"};

    const auto state_map = EvaluateState(var_names, state_specs);

    for (const auto& var_name : var_names)
      vol_model.VarNew(var_name) = state_map.at(var_name);

    // Compute Reynold's number
    const double rho = rho_i_tp1;
    const double u = vol_model.VarNew("u");
    const double Dh = vol_model.HydraulicDiameter();
    const double mu = vol_model.VarNew("mu");

    const double Re = std::fabs(rho * u * Dh / mu);

    vol_model.VarNew("Re") = Re;
  }

  for (size_t vol_comp_id : pipe_system_ptr_->VolumeComponentIDs())
  {
    auto& vol_model = GetComponentLiquidModel(vol_comp_id);
    const auto& cell = *vol_model.GetCellPtr();
    BroadcastStateMap(
      vol_model.VarNames(), vol_model.VarsMapNew(), cell.partition_id_);
  }

  step_time_ = Chi::program_timer.GetTime() - start_time;
  intgl_step_time_ += step_time_;
  step_counter_ += 1.0;

  t_solve.TimeSectionEnd();
}

} // namespace piper