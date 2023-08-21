#include "IncompressibleLiquidPhysics.h"

#include "piper/Piper.h"
#include "piper/utils/CoolPropInterface.h"
#include "piper/components/HardwareComponent.h"

#include "math/chi_math.h"

#include "chi_log.h"

namespace piper
{

void IncompressibleLiquidPhysics::Step()
{
  ChiLogicalErrorIf(Chi::mpi.process_count != 1, "ONLY SERIAL!!!");
  typedef chi_mesh::Vector3 Vec3;
  const double epsilon = 1.0e-8;

  const double dt = 0.001;

  const auto& pipe_system = *pipe_system_ptr_;
  const auto& vol_comp_ids = pipe_system.VolumeComponentIDs();
  const auto& junc_comp_ids = pipe_system.JunctionComponentIDs();

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

  std::map<size_t, EqCoeffs> junction_COM_eqs;
  std::map<size_t, EqCoeffs> volume_COMass_eqs;
  std::map<size_t, EqCoeffs> volume_EOS_eqs;

  for (const size_t junc_comp_id : junc_comp_ids)
  {
    const auto& junction_model = component_models_.at(junc_comp_id);
    const auto& junction_connections = junction_model->ConnectionPoints();
    const size_t I = junction_connections.size();

    VecDbl u_i(I, 0.0);
    VecDbl rho_i(I, 0.0);
    VecDbl Ax_i(I, 0.0);
    VecDbl V_i(I, 0.0);
    VecDbl f_i(I, 0.0);
    VecDbl friction_loss_i(I, 0.0);
    double V_j = 0.0;
    double rho_j = 0.0;

    VecDbl p_i_old(I, 0.0); // for checking

    for (size_t i = 0; i < I; ++i)
    {
      const bool junction_outgoing_at_i =
        junction_model->FlowOrientationRelToConPoint(i) ==
        utils::FlowOrientation::OUTGOING;
      const auto& junction_connection = junction_connections[i];

      const auto& volume_i_model =
        component_models_.at(junction_connection.connected_comp_id_);

      const bool is_boundary =
        volume_i_model->Category() == ComponentCategory::BoundaryLike;

      const double A_i =
        is_boundary ? junction_model->Area() : volume_i_model->Area();
      const double V_half = 0.5 * volume_i_model->Volume();
      const bool volume_i_outgoing =
        volume_i_model->FlowOrientationRelToConPoint(
          junction_connection.connected_comp_connection_point_id_) ==
        utils::FlowOrientation::OUTGOING;
      const chi_mesh::Vector3& vol_orientation =
        volume_i_model->GetOrientation().Vector();

      const double u_i_temp = volume_i_model->VarOld("u");

      u_i[i] = volume_i_outgoing ? -u_i_temp : u_i_temp;
      rho_i[i] = volume_i_model->VarOld("rho");
      Ax_i[i] = junction_outgoing_at_i ? A_i : -A_i;
      V_i[i] = V_half;

      double gravity_force = rho_i[i] * gravity_.Dot(vol_orientation);

      gravity_force *= volume_i_outgoing ? 1.0 : -1.0;
      gravity_force *= junction_outgoing_at_i ? 1.0 : -1.0;

      f_i[i] += gravity_force;

      p_i_old[i] = volume_i_model->VarOld("p");

      V_j += V_half;
      rho_j += V_half * rho_i[i];
    } // for connection i
    rho_j /= V_j;

    const double form_loss = 0.0;

    EqCoeffs eq_cons_of_mom;
    {
      auto& eq = eq_cons_of_mom;

      const double C = dt / (V_j * rho_j);

      eq.coeff_sets_ = {VecDbl(I, 0.0)};
      for (size_t i = 0; i < I; ++i)
        eq.coeff_sets_[0][i] = -C * Ax_i[i];

      eq.id_maps_ = {std::vector<size_t>(I, 0)};
      for (size_t i = 0; i < I; ++i)
        eq.id_maps_[0][i] = junction_connections[i].connected_comp_id_;

      eq.rhs_ = junction_model->VarOld("u");
      for (size_t i = 0; i < I; ++i)
      {
        eq.rhs_ += -C * Ax_i[i] * rho_i[i] * u_i[i] * u_i[i];
        eq.rhs_ += C * V_i[i] * f_i[i];
        eq.rhs_ += C * V_i[i] * friction_loss_i[i];
        eq.rhs_ += form_loss;
      }

      double check = eq.rhs_;
      for (size_t i = 0; i < I; ++i)
        check += eq.coeff_sets_[0][i] * p_i_old[i];

      Chi::log.Log0Verbose1()
        << "Junction \"" << junction_model->Name()
        << "\" delta_u_j=" << check - junction_model->VarOld("u");

      // if (junction_model->Name() == "j1")
      //{
      //   printf("Density %.10e %.10e\n", rho_i[1], rho_i[0]);
      //   Chi::log.Log() << "  rho = " << rho_i[0];
      //   Chi::log.Log() << "  p0 = " << p_i_old[0];
      //   Chi::log.Log() << "  p1 = " << p_i_old[1];
      //   Chi::log.Log() << "  dp = " << p_i_old[1] - p_i_old[0];
      //   Chi::log.Log() << " Adp0 = " << Ax_i[0]*(p_i_old[1] - p_i_old[0]);
      //   Chi::log.Log() << " Adp1 = " << Ax_i[1]*(p_i_old[1] - p_i_old[0]);
      //   Chi::log.Log() << " V f = " << V_i[1] * (f_i[1]);
      //   Chi::log.Log() << " C V f = " << C * V_i[1] * (f_i[1]);
      //   //Chi::log.Log() << "  Axi pi = " << Ax_i[i] * p_i_old[i];
      // }
    }
    junction_COM_eqs[junc_comp_id] = std::move(eq_cons_of_mom);

  } // for junction component

  for (const auto& vol_comp_id : vol_comp_ids)
  {
    auto& vol_comp_model = *component_models_.at(vol_comp_id);
    const auto& vol_connections = vol_comp_model.ConnectionPoints();
    const size_t J = vol_connections.size();

    const double Ai = vol_comp_model.Area();
    const double Vi = vol_comp_model.Volume();
    const double p_t = vol_comp_model.VarOld("p");
    const double rho_t = vol_comp_model.VarOld("rho");
    const double e_t = vol_comp_model.VarOld("e");
    const double q_t = 0.0; // TODO: compute or get somewhere

    VecDbl u_j(J, 0.0);
    VecDbl rho_j(J, 0.0);
    VecDbl e_j(J, 0.0);
    VecDbl p_j(J, 0.0);
    VecDbl gz_j(J, 0.0);
    VecDbl Ax_j(J, 0.0);
    double avg_flowrate = 0.0;

    for (size_t j = 0; j < J; ++j)
    {
      const bool volume_outgoing_at_j =
        vol_comp_model.FlowOrientationRelToConPoint(j) ==
        utils::FlowOrientation::OUTGOING;
      const auto& vol_connection = vol_connections[j];

      const auto& junction_j_model =
        component_models_.at(vol_connection.connected_comp_id_);
      const double A_j = junction_j_model->Area();
      const bool junction_j_outgoing =
        junction_j_model->FlowOrientationRelToConPoint(
          vol_connection.connected_comp_connection_point_id_) ==
        utils::FlowOrientation::OUTGOING;

      const auto& junction_j_connections = junction_j_model->ConnectionPoints();
      const size_t volumeA_id = junction_j_connections.at(0).connected_comp_id_;
      const size_t volumeB_id = junction_j_connections.at(1).connected_comp_id_;

      const auto& volumeA = component_models_.at(volumeA_id);
      const auto& volumeB = component_models_.at(volumeB_id);

      const double u_j_old = junction_j_model->VarOld("u");
      const auto& upwind_volume = (u_j_old >= 0) ? volumeA : volumeB;

      u_j[j] = junction_j_outgoing ? -u_j_old : u_j_old;
      // u_j[j] = u_j_old;
      rho_j[j] = upwind_volume->VarOld("rho");
      e_j[j] = upwind_volume->VarOld("e");
      p_j[j] = upwind_volume->VarOld("p");
      gz_j[j] = junction_j_model->MakeCentroid().Dot(gravity_);
      Ax_j[j] = volume_outgoing_at_j ? A_j : -A_j;

      avg_flowrate += Ax_j[j] * u_j[j];
    } // for connection j
    avg_flowrate /= double(J);

    vol_comp_model.VarOld("u") = avg_flowrate / Ai;

    EqCoeffs eq_cons_of_mass;
    {
      auto& eq = eq_cons_of_mass;

      eq.coeff_sets_ = {VecDbl(J, 0.0)};
      for (size_t j = 0; j < J; ++j)
        eq.coeff_sets_[0][j] = -1.0 * (dt / Vi) * Ax_j[j] * rho_j[j];

      eq.rhs_ = rho_t;
    }

    EqCoeffs eq_cons_of_energy;
    {
      auto& eq = eq_cons_of_energy;

      eq.coeff_sets_ = {VecDbl(J, 0.0)};
      for (size_t j = 0; j < J; ++j)
        eq.coeff_sets_[0][j] =
          -1.0 * (dt / Vi) * Ax_j[j] *
          (rho_j[j] * e_j[j] + rho_j[j] * gz_j[j] + p_j[j]);

      eq.rhs_ = dt * q_t + rho_t * e_t;
    }

    EqCoeffs eq_eos;
    {
      auto& eq = eq_eos;
      const double e_star = (1.0 + epsilon) * e_t;
      const double rho_star = PropSI("rho", "p", p_t, "e", e_star, fluid_name_);

      const double d_rho_e_d_rho =
        (rho_star * e_star - rho_t * e_t) / (rho_star - rho_t);
      eq.coeff_sets_ = {VecDbl(1, d_rho_e_d_rho)};
      eq.rhs_ = rho_t * e_t - d_rho_e_d_rho * rho_t;
    }

    EqCoeffs eq_cons_of_energy2 = eq_cons_of_energy;
    {
      auto& eq = eq_cons_of_energy2;

      const double F = eq_eos.coeff_sets_[0][0];

      for (size_t j = 0; j < J; ++j)
        eq.coeff_sets_[0][j] /= F;

      eq.rhs_ += F * rho_t - rho_t * e_t;
      eq.rhs_ /= F;
    }

    EqCoeffs eq_cons_of_energy3;
    {
      auto& eq = eq_cons_of_energy3;

      eq.coeff_sets_ = eq_cons_of_mass.coeff_sets_;
      for (size_t j = 0; j < J; ++j)
        eq.coeff_sets_[0][j] -= eq_cons_of_energy2.coeff_sets_[0][j];

      eq.rhs_ = eq_cons_of_energy2.rhs_ - eq_cons_of_mass.rhs_;

      Chi::log.Log0Verbose1() << "COE3 rhs=" << eq.rhs_;
      Chi::log.Log0Verbose1() << "COE3 ui=" << avg_flowrate;
    }

    {
      const size_t ir = vol_comp_id_2_row_i_map.at(vol_comp_id);

      b[ir] = eq_cons_of_energy3.rhs_;
      for (size_t j = 0; j < J; ++j)
      {
        const size_t junction_j_id = vol_connections[j].connected_comp_id_;
        const auto& eq_cons_of_mom_j = junction_COM_eqs[junction_j_id];

        const double a_j_COE3 = eq_cons_of_energy3.coeff_sets_[0][j];

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
    }
    volume_COMass_eqs[vol_comp_id] = std::move(eq_cons_of_mass);
    volume_EOS_eqs[vol_comp_id] = std::move(eq_eos);

  } // for volume components

  // chi_math::PrintMatrix(A);
  // chi_math::PrintVector(b);
  p = b;
  chi_math::GaussElimination(A, p, static_cast<int>(num_volumetric_components));
  chi_math::PrintVector(p);

  // Update pressure
  for (const auto& [vol_comp_id, ir] : vol_comp_id_2_row_i_map)
  {
    auto& vol_comp_model = *component_models_.at(vol_comp_id);
    vol_comp_model.VarOld("p") = p[ir];
  }

  // Update junction velocity
  for (const size_t junc_comp_id : junc_comp_ids)
  {
    auto& junction_model = component_models_.at(junc_comp_id);
    const auto& junction_connections = junction_model->ConnectionPoints();
    const size_t I = junction_connections.size();

    const auto& eq_COM = junction_COM_eqs.at(junc_comp_id);

    double u_j_tp1 = eq_COM.rhs_;
    for (size_t i = 0; i < I; ++i)
    {
      const size_t comp_id = eq_COM.id_maps_[0][i];
      const double coeff = eq_COM.coeff_sets_[0][i];

      const auto& vol_comp_model = component_models_.at(comp_id);
      u_j_tp1 += coeff * vol_comp_model->VarOld("p");
    }

    junction_model->VarOld("u") = u_j_tp1;

    Chi::log.Log0Verbose1() << "--------- " << junction_model->Name();
    Chi::log.Log0Verbose1() << "u_j_new = " << u_j_tp1;
  }

  // Update density and internal energy + other state variables
  for (const auto& [vol_comp_id, ir] : vol_comp_id_2_row_i_map)
  {
    auto& vol_comp_model = *component_models_.at(vol_comp_id);
    const auto& vol_connections = vol_comp_model.ConnectionPoints();
    const size_t J = vol_connections.size();

    const auto& eq_COMass = volume_COMass_eqs.at(vol_comp_id);

    Chi::log.Log0Verbose1() << "rho_old = " << vol_comp_model.VarOld("rho");
    double rho_i_tp1 = eq_COMass.rhs_;
    for (size_t j = 0; j < J; ++j)
    {
      const auto& vol_connection = vol_connections[j];
      const auto& junction_j_model =
        component_models_.at(vol_connection.connected_comp_id_);

      rho_i_tp1 += eq_COMass.coeff_sets_[0][j] * junction_j_model->VarOld("u");
    }
    vol_comp_model.VarOld("rho") = rho_i_tp1;

    //const auto& eq_EOS = volume_EOS_eqs.at(vol_comp_id);
    //double rho_e_i_tp1 = eq_EOS.coeff_sets_[0][1] * rho_i_tp1 + eq_EOS.rhs_;

    // const double e_i_tp1 = rho_e_i_tp1 / rho_i_tp1;
    // vol_comp_model.VarNew("e") = e_i_tp1;
    // Chi::log.Log0Verbose1() << "e_new=" << e_i_tp1;

    const std::vector<StateVal> state_specs = {
      {"p", vol_comp_model.VarOld("p")}, {"rho", vol_comp_model.VarOld("rho")}};

    const auto state = EvaluateState(state_specs);

    const std::vector<std::string> var_names = {
      "rho", "e", "T", "p", "h", "s", "k", "Pr", "mu"};

    Chi::log.Log0Verbose1()
      << "************* " << vol_comp_model.Name();
    for (const auto& var_name : var_names)
    {
      vol_comp_model.VarOld(var_name) = state.at(var_name);
      Chi::log.Log0Verbose1() << var_name << " = " << state.at(var_name);
    }
  }
}

} // namespace piper