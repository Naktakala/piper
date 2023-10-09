#include "SingleVolumeLiquidModel.h"

#include "piper/physics/LiquidPhysics.h"

#include "physics/TimeSteppers/TimeStepper.h"

#include "mesh/Cell/cell.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

#include <functional>

namespace piper
{

RegisterChiObjectParametersOnly(piper, SingleVolumeLiquidModel);

chi::InputParameters SingleVolumeLiquidModel::GetInputParameters()
{
  chi::InputParameters params = ComponentLiquidModel::GetInputParameters();

  params.SetGeneralDescription(
    "General control volume with liquid physics model.");
  params.SetDocGroup("Piper");
  params.SetObjectType("SingleVolumeLiquidModel");

  params.AddOptionalParameter(
    "volumetric_heat_generation",
    0.0,
    "Externally applied heat generation to the entire volume.");

  return params;
}

SingleVolumeLiquidModel::SingleVolumeLiquidModel(
  const chi::InputParameters& params,
  std::vector<std::unique_ptr<ComponentModel>>& family,
  LiquidPhysics& liquid_physics,
  const HardwareComponent& hardware_component,
  const chi_mesh::Cell* cell,
  const std::vector<std::string>& variable_names)
  : ComponentLiquidModel(
      params, family, liquid_physics, hardware_component, cell, variable_names),
    volumetric_heat_generation_(
      params.GetParamValue<double>("volumetric_heat_generation"))
{
}

void SingleVolumeLiquidModel::AssembleEquations()
{
  typedef std::vector<double> VecDbl;
  typedef chi_mesh::Vector3 Vec3;

  const double dt = physics_.GetTimeStepper().TimeStepSize();
  const Vec3& gravity = physics_.GravityVector();
  const double epsilon = 1.0e-8;

  auto& vol_model = *this;
  const size_t J = vol_model.Connections().size();

  const double Ai = vol_model.Area();
  const double Vi = vol_model.Volume();
  const double p_t = vol_model.VarOld("p");
  const double rho_t = vol_model.VarOld("rho");
  const double e_t = vol_model.VarOld("e");
  const double q_t = volumetric_heat_generation_;

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

    const auto& volumeA = family_.at(volumeA_id);
    const auto& volumeB = family_.at(volumeB_id);

    const double u_j_old = jnc_j_model.VarOld("u");
    const auto& upwind_volume = (u_j_old >= 0) ? volumeA : volumeB;

    u_j[j] = jnc_dir == ABDir::AB ? u_j_old : -u_j_old;
    rho_j[j] = upwind_volume->VarOld("rho");
    e_j[j] = upwind_volume->VarOld("e");
    p_j[j] = upwind_volume->VarOld("p");
    gz_j[j] = jnc_j_model.MakeCentroid().Dot(gravity);
    Ax_j[j] = j == 0 ? -A_j : A_j;

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
      eq.coeff_sets_[0][j] = -1.0 * (dt / Vi) * Ax_j[j] *
                             (rho_j[j] * e_j[j] + rho_j[j] * gz_j[j] + p_j[j]);
      eq.rhs_ += (dt / Vi) * form_loss_j[j];
    }

    eq.rhs_ += dt * q_t + rho_t * e_t;
  }

  EqCoeffs eq_eos;
  {
    auto& eq = eq_eos;

    const double T_star = (1.0 + epsilon) * vol_model.VarOld("T");

    const auto state_map =
      physics_.EvaluateState({"e", "rho"}, {{"p", p_t}, {"T", T_star}});

    const double e_star = state_map.at("e");
    const double rho_star = state_map.at("rho");
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
  }

  equation_coefficients_ = {std::move(eq_cons_of_mass),
                            std::move(eq_eos),
                            std::move(eq_cons_of_energy3)};
}

} // namespace piper