#include "FEMTimeKernel.h"

#include "math/SpatialDiscretization/FiniteElement/finite_element.h"

namespace chi_math
{

chi::InputParameters FEMTimeKernel::GetInputParameters()
{
  chi::InputParameters params = FEMKernel::GetInputParameters();

  return params;
}

FEMTimeKernel::FEMTimeKernel(const chi::InputParameters& params)
  : FEMKernel(params)
{
}

bool FEMTimeKernel::IsTimeKernel() const { return true; }

} // namespace chi_math