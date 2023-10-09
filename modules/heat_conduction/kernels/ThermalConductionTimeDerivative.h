#ifndef PIPER_THERMALCONDUCTIONTIMEDERIVATIVE_H
#define PIPER_THERMALCONDUCTIONTIMEDERIVATIVE_H

#include "math/KernelSystem/FEMKernels/FEMTimeKernel.h"

namespace hcm
{

class ThermalConductionTimeDerivative : public chi_math::FEMTimeKernel
{
public:
  static chi::InputParameters GetInputParameters();
  explicit ThermalConductionTimeDerivative(const chi::InputParameters& params);

  double ResidualEntryAtQP() override;
  double JacobianEntryAtQP() override;

protected:
  const std::string rho_property_name_;
  const std::string cp_property_name_;
  const chi_math::FEMMaterialProperty& rho_;
  const chi_math::FEMMaterialProperty& Cp_;
};

} // namespace hcm

#endif // PIPER_THERMALCONDUCTIONTIMEDERIVATIVE_H
