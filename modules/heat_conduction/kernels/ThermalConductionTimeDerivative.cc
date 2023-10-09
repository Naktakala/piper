#include "ThermalConductionTimeDerivative.h"

#include "ChiObjectFactory.h"

namespace hcm
{

RegisterChiObject(hcm, ThermalConductionTimeDerivative);

chi::InputParameters ThermalConductionTimeDerivative::GetInputParameters()
{
  chi::InputParameters params = chi_math::FEMTimeKernel::GetInputParameters();

  params.SetGeneralDescription("A time kernel with material properties density "
                               "and heat capacity as coefficients.");
  params.SetDocGroup("doc_HeatConduction");

  params.AddOptionalParameter("rho_property_name",
                              "density",
                              "Name of the material property for density.");
  params.AddOptionalParameter(
    "cp_property_name",
    "heat_capacity",
    "Name of the material property for heat capacity.");

  return params;
}

ThermalConductionTimeDerivative::ThermalConductionTimeDerivative(
  const chi::InputParameters& params)
  : chi_math::FEMTimeKernel(params),
    rho_property_name_(params.GetParamValue<std::string>("rho_property_name")),
    cp_property_name_(params.GetParamValue<std::string>("cp_property_name")),
    rho_(GetFEMMaterialProperty(rho_property_name_)),
    Cp_(GetFEMMaterialProperty(cp_property_name_))
{
}

double ThermalConductionTimeDerivative::ResidualEntryAtQP()
{
  return rho_[qp_] * Cp_[qp_] * test_values_[i_][qp_] * var_dot_value_[qp_];
}

/**(a/dt) * d/dT (rho Cp T). Let rho Cp T = fgh.
 * d/dx (fgh) = fg dh/dx + h dfg/dx
 *            = fg dh/dx + h (f dg/dx + g df/dx)
 *            = fg dh/dx + f dg/dx h + df/dx g h*/
double ThermalConductionTimeDerivative::JacobianEntryAtQP()
{
  const double a_div_dt = var_dot_dvar_;

  const double test_i = test_values_[i_][qp_];
  const double shape_j = shape_values_[j_][qp_];
  const double h = test_i * var_dot_value_[qp_];

  const double fg_dhdx = rho_[qp_] * Cp_[qp_] * test_i * shape_j * a_div_dt;
  const double f_dgdx_h = rho_[qp_] * Cp_.dvar_[qp_] * h * shape_j;
  const double dfdx_g_h = rho_.dvar_[qp_] * Cp_[qp_] * h * shape_j;

  return fg_dhdx + f_dgdx_h + dfdx_g_h;
}

} // namespace hcm