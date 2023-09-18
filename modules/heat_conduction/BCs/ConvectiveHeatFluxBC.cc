#include "ConvectiveHeatFluxBC.h"

#include "ChiObjectFactory.h"

namespace hcm
{

RegisterChiObject(hcm, ConvectiveHeatFluxBC);

chi::InputParameters ConvectiveHeatFluxBC::GetInputParameters()
{
  chi::InputParameters params =
    chi_math::FEMBoundaryCondition::GetInputParameters();

  params.SetGeneralDescription(
    "A boundary condition modeling a convective surface.");
  params.SetDocGroup("doc_HeatConduction");

  params.AddRequiredParameter<double>("T_bulk", "Bulk fluid temperature.");
  params.AddRequiredParameter<double>("convection_coefficient",
                                      "The convection coefficient (h).");

  return params;
}

ConvectiveHeatFluxBC::ConvectiveHeatFluxBC(const chi::InputParameters& params)
: chi_math::FEMBoundaryCondition(params),
    T_bulk_(params.GetParamValue<double>("T_bulk")),
    h_(params.GetParamValue<double>("convection_coefficient"))
{
}

double ConvectiveHeatFluxBC::ResidualEntryAtQP()
{
  return -test_values_[i_][qp_] * h_ * (T_bulk_ - var_value_[qp_]);
}

double ConvectiveHeatFluxBC::JacobianEntryAtQP()
{
  return -test_values_[i_][qp_] * shape_values_[j_][qp_] * (-h_);
}

} // namespace hcm