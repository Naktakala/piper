#ifndef PIPER_THERMAL_CONDUCTIONKERNEL_H
#define PIPER_THERMAL_CONDUCTIONKERNEL_H

#include "ChiObject.h"
#include "math/KernelSystem/FEMKernels/FEMKernel.h"

namespace hcm
{

class ThermalConductionKernel : public chi_math::FEMKernel
{
public:
  static chi::InputParameters GetInputParameters();
  explicit ThermalConductionKernel(
    const chi::InputParameters& params);

  double ResidualEntryAtQP() override;
  double JacobianEntryAtQP() override;

protected:
  double k_;
};

}

#endif // PIPER_THERMAL_CONDUCTIONKERNEL_H
