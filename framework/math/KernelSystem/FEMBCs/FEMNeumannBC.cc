#include "FEMNeumannBC.h"

#include "ChiObjectFactory.h"

namespace chi_math
{

RegisterChiObject(chi_math, FEMNeumannBC);

chi::InputParameters FEMNeumannBC::GetInputParameters()
{
  chi::InputParameters params = FEMBoundaryCondition::GetInputParameters();

  params.SetGeneralDescription("Implementation of a general NeumannBC");
  params.SetDocGroup("doc_KernelSystem");

  params.AddOptionalParameter(
    "flux_value", 0.0, "Value of the NeumannBC flux.");

  return params;
}

FEMNeumannBC::FEMNeumannBC(const chi::InputParameters& params)
  : FEMBoundaryCondition(params),
    flux_value_(params.GetParamValue<double>("flux_value"))
{
}

double FEMNeumannBC::ResidualEntryAtQP()
{
  return -test_i_qp_ * flux_value_;
}

} // namespace chi_math