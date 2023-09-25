#include "CoupledConvectiveHeatFluxBC.h"

#include "ChiObjectFactory.h"

namespace hcm
{

RegisterChiObject(hcm, CoupledConvectiveHeatFluxBC);

chi::InputParameters CoupledConvectiveHeatFluxBC::GetInputParameters()
{
  chi::InputParameters params =
    chi_math::FEMBoundaryCondition::GetInputParameters();

  params.SetGeneralDescription("A boundary condition modeling a convective "
                               "surface but with the parameters coupled to"
                               " other variables.");
  params.SetDocGroup("doc_HeatConduction");

  params.AddOptionalParameter(
    "T_bulk", 1.0, "A value or the name of the variable to couple to.");
  params.SetParameterTypeMismatchAllowed("T_bulk");

  params.AddOptionalParameter(
    "convection_coefficient", 1.0, "The heat transfer coefficient");
  params.SetParameterTypeMismatchAllowed("convection_coefficient");

  return params;
}

CoupledConvectiveHeatFluxBC::CoupledConvectiveHeatFluxBC(
  const chi::InputParameters& params)
  : chi_math::FEMBoundaryCondition(params),
    T_bulk_(GetCoupledField(params.GetParamValue<std::string>("T_bulk"))),
    h_(params.GetParamValue<double>("convection_coefficient"))
{
}

double CoupledConvectiveHeatFluxBC::ResidualEntryAtQP()
{
  return -test_values_[i_][qp_] * h_ * (T_bulk_[qp_] - var_value_[qp_]);
}

double CoupledConvectiveHeatFluxBC::JacobianEntryAtQP()
{
   return -test_values_[i_][qp_] * shape_values_[j_][qp_] * (-h_);
}

} // namespace hcm