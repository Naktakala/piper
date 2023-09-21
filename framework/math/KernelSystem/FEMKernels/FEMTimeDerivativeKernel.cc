#include "FEMTimeDerivativeKernel.h"

#include "ChiObjectFactory.h"

namespace chi_math
{

RegisterChiObject(chi_math, FEMTimeDerivativeKernel);

chi::InputParameters FEMTimeDerivativeKernel::GetInputParameters()
{
  chi::InputParameters params = FEMTimeKernel::GetInputParameters();

  params.SetGeneralDescription("Time derivative kernel");
  params.SetDocGroup("doc_Kernels");

  return params;
}

FEMTimeDerivativeKernel::FEMTimeDerivativeKernel(
  const chi::InputParameters& params)
  : FEMTimeKernel(params)
{
}

double FEMTimeDerivativeKernel::ResidualEntryAtQP()
{
  return test_values_[i_][qp_] * var_dot_value_[qp_];
}

double FEMTimeDerivativeKernel::JacobianEntryAtQP()
{
  return test_values_[i_][qp_] * shape_values_[j_][qp_] * var_dot_dvar_;
}

} // namespace chi_math