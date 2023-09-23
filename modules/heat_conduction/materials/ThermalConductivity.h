#ifndef PIPER_THERMALCONDUCTIVITY_H
#define PIPER_THERMALCONDUCTIVITY_H

#include "materials/MaterialProperty.h"

namespace chi_math
{
class FunctionDimAToDimB;
}

namespace hcm
{

class ThermalConductivity : public chi::MaterialProperty
{
public:
  static chi::InputParameters GetInputParameters();
  explicit ThermalConductivity(const chi::InputParameters& params);

  std::vector<std::string> RequiredInputNames() const override;
  double ComputeScalarValue(const std::vector<double>& input_params) const override;

  bool HasDerivative() const override;

protected:
  double constant_value_;
  std::shared_ptr<chi_math::FunctionDimAToDimB> function_;
  std::string temperature_fieldname_;

private:
  std::shared_ptr<chi_math::FunctionDimAToDimB>
  GetFunction(const chi::InputParameters& params);
};

} // namespace hcm

#endif // PIPER_THERMALCONDUCTIVITY_H
