#ifndef PIPER_THERMALCONDUCTIONKERNEL2_H
#define PIPER_THERMALCONDUCTIONKERNEL2_H

#include "ChiObject.h"
#include "math/KernelSystem/FEMKernels/FEMKernel.h"
#include "math/KernelSystem/Coupling/FEMMaterialProperty.h"

namespace chi_math
{
class FEMMaterialProperty;
}

namespace hcm
{

class ThermalConductionKernel2 : public chi_math::FEMKernel

{
public:
  static chi::InputParameters GetInputParameters();
  explicit ThermalConductionKernel2(
    const chi::InputParameters& params);

  double ResidualEntryAtQP() override;
  double JacobianEntryAtQP() override;

protected:
  const chi_math::FEMMaterialProperty& k_;
};

}

#endif // PIPER_THERMALCONDUCTIONKERNEL2_H
