#ifndef PIPER_COUPLEDHEATGENERATION_H
#define PIPER_COUPLEDHEATGENERATION_H

#include "math/KernelSystem/FEMKernels/FEMKernel.h"

namespace hcm
{

class CoupledHeatGeneration : public chi_math::FEMKernel
{
public:
  static chi::InputParameters GetInputParameters();
  explicit CoupledHeatGeneration(const chi::InputParameters& params);

  double ResidualEntryAtQP() override;
  double JacobianEntryAtQP() override {return 0.0;}

protected:
  const chi_math::FEMCoupledField& e_gen_;
};

}

#endif // PIPER_COUPLEDHEATGENERATION_H
