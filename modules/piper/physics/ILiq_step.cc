#include "IncompressibleLiquidPhysics.h"

#include "piper/Piper.h"
#include "piper/utils/CoolPropInterface.h"
#include "piper/components/HardwareComponent.h"

#include "math/chi_math.h"

#include "chi_log.h"

#include <functional>

namespace piper
{

void IncompressibleLiquidPhysics::Step()
{
  ChiLogicalErrorIf(Chi::mpi.process_count != 1, "ONLY SERIAL!!!");
  typedef chi_mesh::Vector3 Vec3;
  const double epsilon = 1.0e-8;

  const double dt = 0.1;

  const auto& pipe_system = *pipe_system_ptr_;
  const auto& vol_comp_ids = pipe_system.VolumeComponentIDs();
  const auto& jnc_comp_ids = pipe_system.JunctionComponentIDs();

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

  std::map<size_t, EqCoeffs> jnc_COM_eqs;
  std::map<size_t, EqCoeffs> vol_COMass_eqs;
  std::map<size_t, EqCoeffs> vol_EOS_eqs;

  for (const size_t jnc_comp_id : jnc_comp_ids)
  {
    const auto& jnc_model = component_models_.at(jnc_comp_id);
    const size_t I = jnc_model->Connections().size();

    VecDbl u_i(I, 0.0);
    VecDbl rho_i(I, 0.0);
    VecDbl Ax_i(I, 0.0);
    VecDbl V_i(I, 0.0);
    VecDbl f_i(I, 0.0);
    VecDbl friction_loss_i(I, 0.0);
    double V_j = 0.0;
    double rho_j = 0.0;

    double form_loss = 0.0;

    VecDbl p_i_old(I, 0.0); // for checking

    enum class ABDir
    {
      AB,
      BA
    };

    for (auto [i, vol_i_model, connection_i] : jnc_model->Connections())
    {
      const bool is_boundary =
        vol_i_model.Category() == ComponentCategory::BoundaryLike;

      const double A_i = is_boundary ? jnc_model->Area() : vol_i_model.Area();
      const double V_half = 0.5 * vol_i_model.Volume();

      const bool vol_i_outgoing_wrt_i = vol_i_model.IsOutgoingRelToConPoint(
        connection_i.connected_comp_connection_point_id_);

      const ABDir vol_dir = std::invoke(
        [&](size_t node_id)
        {
          ABDir val;
          if (node_id == 0) val = vol_i_outgoing_wrt_i ? ABDir::BA : ABDir::AB;
          if (node_id == 1) val = vol_i_outgoing_wrt_i ? ABDir::AB : ABDir::BA;
          return val;
        },
        /*node_id=*/i);

      const Vec3& vol_orientation = vol_i_model.GetOrientation().Vector();

      const double u_i_temp = vol_i_model.VarOld("u");

      u_i[i] = vol_dir == ABDir::BA ? -u_i_temp : u_i_temp;
      rho_i[i] = vol_i_model.VarOld("rho");
      Ax_i[i] = i == 0 ? -A_i : A_i;
      V_i[i] = V_half;

      const double gravity_force = rho_i[i] * gravity_.Dot(vol_orientation);

      f_i[i] = vol_dir == ABDir::BA ? -gravity_force : gravity_force;

      if (not is_boundary)
      {
        const double f_D = FrictionFactor(vol_i_model);
        const double Dh = vol_i_model.HydraulicDiameter();
        const double rho = rho_i[i];
        const double u = u_i[i];
        const double friction_loss = -f_D * 0.5 * rho * u * u * (1.0 / Dh);

        friction_loss_i[i] =
          vol_dir == ABDir::BA ? -friction_loss : friction_loss;
      }

      p_i_old[i] = vol_i_model.VarOld("p");

      V_j += V_half;
      rho_j += V_half * rho_i[i];
    } // for connection i
    rho_j /= V_j;

    if (jnc_model->Name() == "j1")
    {
      const double u = jnc_model->VarOld("u");
      const double K = 5.0;
      form_loss = -K * 0.5 * rho_j * u * u ;
    }

    EqCoeffs eq_cons_of_mom;
    {
      auto& eq = eq_cons_of_mom;

      const double C = dt / (V_j * rho_j);

      eq.coeff_sets_ = {VecDbl(I, 0.0)};
      for (size_t i = 0; i < I; ++i)
        eq.coeff_sets_[0][i] = -C * Ax_i[i];

      auto& jnc_connections = jnc_model->ConnectionPoints();
      eq.id_maps_ = {std::vector<size_t>(I, 0)};
      for (size_t i = 0; i < I; ++i)
        eq.id_maps_[0][i] = jnc_connections[i].connected_comp_id_;

      eq.rhs_ = jnc_model->VarOld("u");
      for (size_t i = 0; i < I; ++i)
      {
        eq.rhs_ += -C * Ax_i[i] * rho_i[i] * u_i[i] * u_i[i];
        eq.rhs_ += C * V_i[i] * f_i[i];
        eq.rhs_ += C * V_i[i] * friction_loss_i[i];
      }
      eq.rhs_ += C * V_j * form_loss;

      double check = eq.rhs_;
      for (size_t i = 0; i < I; ++i)
        check += eq.coeff_sets_[0][i] * p_i_old[i];

      Chi::log.Log0Verbose1()
        << "Junction \"" << jnc_model->Name()
        << "\" delta_u_j=" << check - jnc_model->VarOld("u");
    }
    jnc_COM_eqs[jnc_comp_id] = std::move(eq_cons_of_mom);

  } // for junction component

  for (const auto& vol_comp_id : vol_comp_ids)
  {
    auto& vol_model = *component_models_.at(vol_comp_id);
    const size_t J = vol_model.Connections().size();

    const double Ai = vol_model.Area();
    const double Vi = vol_model.Volume();
    const double p_t = vol_model.VarOld("p");
    const double rho_t = vol_model.VarOld("rho");
    const double e_t = vol_model.VarOld("e");
    const double q_t = 1.0e6; // TODO: compute or get somewhere

    VecDbl u_j(J, 0.0);
    VecDbl rho_j(J, 0.0);
    VecDbl e_j(J, 0.0);
    VecDbl p_j(J, 0.0);
    VecDbl gz_j(J, 0.0);
    VecDbl form_loss_j(J, 0.0);
    VecDbl Ax_j(J, 0.0);
    double avg_flowrate = 0.0;

    enum class ABDir
    {
      AB,
      BA
    };

    for (auto [j, jnc_j_model, connection_j] : vol_model.Connections())
    {
      const double A_j = jnc_j_model.Area();

      const bool jnc_j_outgoing = jnc_j_model.IsOutgoingRelToConPoint(
        connection_j.connected_comp_connection_point_id_);

      const ABDir jnc_dir = std::invoke(
        [&](size_t node_id)
        {
          ABDir val;
          if (node_id == 0) val = jnc_j_outgoing ? ABDir::BA : ABDir::AB;
          if (node_id == 1) val = jnc_j_outgoing ? ABDir::AB : ABDir::BA;
          return val;
        },
        /*node_id=*/j);

      const auto& junction_j_connections = jnc_j_model.ConnectionPoints();
      const size_t volumeA_id = junction_j_connections.at(0).connected_comp_id_;
      const size_t volumeB_id = junction_j_connections.at(1).connected_comp_id_;

      const auto& volumeA = component_models_.at(volumeA_id);
      const auto& volumeB = component_models_.at(volumeB_id);

      const double u_j_old = jnc_j_model.VarOld("u");
      const auto& upwind_volume = (u_j_old >= 0) ? volumeA : volumeB;

      u_j[j] = jnc_dir == ABDir::AB ? u_j_old : -u_j_old;
      rho_j[j] = upwind_volume->VarOld("rho");
      e_j[j] = upwind_volume->VarOld("e");
      p_j[j] = upwind_volume->VarOld("p");
      gz_j[j] = jnc_j_model.MakeCentroid().Dot(gravity_);
      Ax_j[j] = j == 0 ? -A_j : A_j;

      if (jnc_j_model.Name() == "j1")
        form_loss_j[j] = -5.0 * 0.5 * rho_j[j] * u_j_old*u_j_old * u_j[j];

      avg_flowrate += Ax_j[j] * u_j[j];
    } // for connection j
    avg_flowrate /= double(J);

    vol_model.VarOld("u") = avg_flowrate / Ai;

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
      {
        eq.coeff_sets_[0][j] =
          -1.0 * (dt / Vi) * Ax_j[j] *
          (rho_j[j] * e_j[j] + rho_j[j] * gz_j[j] + p_j[j]);
        eq.rhs_ += (dt/Vi) * form_loss_j[j];
      }

      eq.rhs_ += dt * q_t + rho_t * e_t;
    }

    EqCoeffs eq_eos;
    {
      auto& eq = eq_eos;
      // const double e_star = (1.0 + epsilon) * e_t;
      // const double rho_star = PropSI("rho", "p", p_t, "e", e_star,
      // fluid_name_);

      const double T_star = (1.0 + epsilon) * vol_model.VarOld("T");
      const double e_star = PropSI("e", "p", p_t, "T", T_star, fluid_name_);
      const double rho_star = PropSI("rho", "p", p_t, "T", T_star, fluid_name_);

      const double d_rho_e_d_rho =
        (rho_star * e_star - rho_t * e_t) / (rho_star - rho_t);

      Chi::log.Log0Verbose1() << "d_rho_e_d_rho " << d_rho_e_d_rho << " "
                              << (rho_star - rho_t) << " " << e_star - e_t;

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

    // Assemble terms into pressure system
    {
      const size_t ir = vol_comp_id_2_row_i_map.at(vol_comp_id);

      b[ir] = eq_cons_of_energy3.rhs_;

      for (auto [j, jnc_j_model, connection_j] : vol_model.Connections())
      {
        const auto& eq_cons_of_mom_j = jnc_COM_eqs[jnc_j_model.GetID()];

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
    vol_COMass_eqs[vol_comp_id] = std::move(eq_cons_of_mass);
    vol_EOS_eqs[vol_comp_id] = std::move(eq_eos);

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
    vol_comp_model.VarOld("p") = p[vol_comp_id_2_row_i_map.at(vol_comp_id)];
  }

  // Update junction velocity
  for (const size_t junc_comp_id : jnc_comp_ids)
  {
    auto& junction_model = component_models_.at(junc_comp_id);
    const auto& junction_connections = junction_model->ConnectionPoints();
    const size_t I = junction_connections.size();

    const auto& eq_COM = jnc_COM_eqs.at(junc_comp_id);

    double u_j_tp1 = eq_COM.rhs_;
    for (size_t i = 0; i < I; ++i)
    {
      const size_t comp_id = eq_COM.id_maps_[0][i];
      const double coeff = eq_COM.coeff_sets_[0][i];

      const auto& vol_comp_model = component_models_.at(comp_id);
      u_j_tp1 += coeff * vol_comp_model->VarOld("p");
    }

    junction_model->VarOld("u") = u_j_tp1;
  }

  // Update density and internal energy + other state variables
  for (const auto& [vol_comp_id, ir] : vol_comp_id_2_row_i_map)
  {
    auto& vol_model = *component_models_.at(vol_comp_id);
    const size_t J = vol_model.Connections().size();

    const double Ai = vol_model.Area();

    const auto& eq_COMass = vol_COMass_eqs.at(vol_comp_id);

    double rho_i_tp1 = eq_COMass.rhs_;
    double avg_flowrate = 0.0;
    for (auto [j, jnc_j_model, connection_j] : vol_model.Connections())
    {
      const double A_j = jnc_j_model.Area();

      const bool junction_j_outgoing = jnc_j_model.IsOutgoingRelToConPoint(
        connection_j.connected_comp_connection_point_id_);

      const double u_j_old = jnc_j_model.VarOld("u");

      rho_i_tp1 += eq_COMass.coeff_sets_[0][j] * u_j_old;

      const double ujj = junction_j_outgoing ? u_j_old : -u_j_old;
      const double Axj = j == 0 ? -A_j : A_j;

      avg_flowrate += Axj * ujj;
    }
    vol_model.VarOld("rho") = rho_i_tp1;
    avg_flowrate /= double(J);

    vol_model.VarOld("u") = avg_flowrate / Ai;

    const auto& eq_EOS = vol_EOS_eqs.at(vol_comp_id);
    double rho_e_i_tp1 = eq_EOS.coeff_sets_[0][0] * rho_i_tp1 + eq_EOS.rhs_;

    const double e_i_tp1 = rho_e_i_tp1 / rho_i_tp1;
    vol_model.VarNew("e") = e_i_tp1;

    const std::vector<StateVal> state_specs = {{"p", vol_model.VarOld("p")},
                                               {"e", e_i_tp1}};

    const auto state = EvaluateState(state_specs);

    const std::vector<std::string> var_names = {
      "rho", "e", "T", "p", "h", "s", "k", "Pr", "mu"};

    Chi::log.Log0Verbose1() << "************* " << vol_model.Name();
    for (const auto& var_name : var_names)
    {
      vol_model.VarOld(var_name) = state.at(var_name);
      Chi::log.Log0Verbose1() << var_name << " = " << state.at(var_name);
    }
    Chi::log.Log() << "T=" << vol_model.VarOld("T")
                   << " rho=" << vol_model.VarOld("rho")
                   << " p=" << vol_model.VarOld("p")
                   << " e=" << vol_model.VarOld("e")
                   << " mu=" << vol_model.VarOld("mu")
                   << " u=" << vol_model.VarOld("u");
    // Compute Reynold's number

    const double rho = rho_i_tp1;
    const double u = vol_model.VarOld("u");
    const double Dh = vol_model.HydraulicDiameter();
    const double mu = vol_model.VarOld("mu");

    const double Re = std::fabs(rho * u * Dh / mu);

    vol_model.VarOld("Re") = Re;
  }
}

} // namespace piper