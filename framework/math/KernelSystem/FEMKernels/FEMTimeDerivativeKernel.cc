#include "FEMTimeDerivativeKernel.h"

#include "ChiObjectFactory.h"

namespace chi_math
{

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
  return test_i_qp_ * (var_qp_value_ - var_value_old_qp_)/dt_;
}

double FEMTimeDerivativeKernel::JacobianEntryAtQP()
{
  return test_i_qp_ * shape_j_qp_ / dt_;
}

} // namespace chi_math