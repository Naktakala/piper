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
  //if (time_ > 0.0)
  //  std::cout << old_var_qp_value_[TIME_T] << "\n";
  return test_i_qp_ * (var_qp_value_ - old_var_qp_value_[TIME_T]) / dt_;
}

double FEMTimeDerivativeKernel::JacobianEntryAtQP()
{
  return test_i_qp_ * shape_j_qp_ / dt_;
}

} // namespace chi_math