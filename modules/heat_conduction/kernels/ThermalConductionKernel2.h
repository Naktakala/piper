#ifndef PIPER_THERMALCONDUCTIONKERNEL2_H
#define PIPER_THERMALCONDUCTIONKERNEL2_H

#include "math/KernelSystem/FEMKernels/FEMKernel.h"

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
  explicit ThermalConductionKernel2(const chi::InputParameters& params);

  double ResidualEntryAtQP() override;
  double JacobianEntryAtQP() override;

protected:
  const std::string k_property_name_;
  const chi_math::FEMMaterialProperty& k_;
};

} // namespace hcm

#endif // PIPER_THERMALCONDUCTIONKERNEL2_H
