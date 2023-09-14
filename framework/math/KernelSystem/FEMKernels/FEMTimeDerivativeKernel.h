#ifndef CHI_FEMTIMEDERIVATIVEKERNEL_H
#define CHI_FEMTIMEDERIVATIVEKERNEL_H

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

#endif // CHI_FEMTIMEDERIVATIVEKERNEL_H
