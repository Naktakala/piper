#include "JunctionLiquidModel.h"

#include "piper/physics/LiquidPhysics.h"
#include "piper/utils/Orientation.h"

#include "ChiObjectFactory.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <functional>

namespace piper
{

RegisterChiObjectParametersOnly(piper, JunctionLiquidModel);

chi::InputParameters JunctionLiquidModel::GetInputParameters()
{
  chi::InputParameters params = ComponentLiquidModel::GetInputParameters();

  params.SetGeneralDescription("General junction with liquid physics model.");
  params.SetDocGroup("Piper");
  params.SetObjectType("JunctionLiquidModel");

  params.AddOptionalParameter(
    "form_loss_forward", 0.0, "Forward form loss coefficient.");
  params.AddOptionalParameter(
    "form_loss_reverse", 0.0, "Reverse form loss coefficient.");

  return params;
}

JunctionLiquidModel::JunctionLiquidModel(
  const chi::InputParameters& params,
  std::vector<std::unique_ptr<ComponentModel>>& family,
  LiquidPhysics& liquid_physics,
  const HardwareComponent& hardware_component,
  const chi_mesh::Cell* cell,
  const std::vector<std::string>& variable_names)
  : ComponentLiquidModel(
      params, family, liquid_physics, hardware_component, cell, variable_names),
    forward_loss_coefficient_(
      params.GetParamValue<double>("form_loss_forward")),
    reverse_loss_coefficient_(params.GetParamValue<double>("form_loss_reverse"))
{
}

void JunctionLiquidModel::AssembleEquations()
{
  typedef std::vector<double> VecDbl;
  typedef chi_mesh::Vector3 Vec3;

  const double dt = physics_.DeltaT();
  const Vec3& gravity = physics_.GravityVector();

  auto& jnc_model = *this;
  const size_t I = jnc_model.Connections().size();

  VecDbl u_i(I, 0.0);
  VecDbl rho_i(I, 0.0);
  VecDbl Ax_i(I, 0.0);
  VecDbl V_i(I, 0.0);
  VecDbl f_i(I, 0.0);
  VecDbl friction_loss_i(I, 0.0);
  double V_j = 0.0;
  double rho_j = 0.0;

  VecDbl p_i_old(I, 0.0); // for checking

  enum class ABDir
  {
    AB,
    BA
  };

  for (auto [i, vol_i_model, connection_i] : jnc_model.Connections())
  {
    const bool is_boundary =
      vol_i_model.Category() == ComponentCategory::BoundaryLike;

    const double A_i = is_boundary ? jnc_model.Area() : vol_i_model.Area();
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

    const double gravity_force = rho_i[i] * gravity.Dot(vol_orientation);

    f_i[i] = vol_dir == ABDir::BA ? -gravity_force : gravity_force;

    if (not is_boundary)
    {
      const double f_D = physics_.FrictionFactorFuncion()(vol_i_model);
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

  double form_loss = 0.0;
  {
    const double u = jnc_model.VarOld("u");
    if (u > 0.0) form_loss = -forward_loss_coefficient_ * 0.5 * rho_j * u * u;
    if (u < 0.0) form_loss = reverse_loss_coefficient_ * 0.5 * rho_j * u * u;
  }

  EqCoeffs eq_cons_of_mom;
  {
    auto& eq = eq_cons_of_mom; // shorthand

    const double C = dt / (V_j * rho_j);

    eq.coeff_sets_ = {VecDbl(I, 0.0)};
    for (size_t i = 0; i < I; ++i)
      eq.coeff_sets_[0][i] = -C * Ax_i[i];

    auto& jnc_connections = jnc_model.ConnectionPoints();
    eq.id_maps_ = {std::vector<size_t>(I, 0)};
    for (size_t i = 0; i < I; ++i)
      eq.id_maps_[0][i] = jnc_connections[i].connected_comp_id_;

    eq.rhs_ = jnc_model.VarOld("u");
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

    Chi::log.Log0Verbose1() << "Junction \"" << jnc_model.Name()
                            << "\" delta_u_j=" << check - jnc_model.VarOld("u");
  }

  equation_coefficients_ = {std::move(eq_cons_of_mom)};
}

} // namespace piper