#ifndef CHITECH_FEMTIMEDERIVATIVEKERNEL_H
#define CHITECH_FEMTIMEDERIVATIVEKERNEL_H

#include "FEMTimeKernel.h"

namespace chi_math
{

class FEMTimeDerivativeKernel : public FEMTimeKernel
{
public:
  static chi::InputParameters GetInputParameters();
  explicit FEMTimeDerivativeKernel(const chi::InputParameters& params);

  double ResidualEntryAtQP() override;
  double JacobianEntryAtQP() override;
};

}

#endif // CHITECH_FEMTIMEDERIVATIVEKERNEL_H
