#ifndef CHITECH_CONSTANTMATERIALPROPERTY_H
#define CHITECH_CONSTANTMATERIALPROPERTY_H

#include "MaterialProperty.h"

namespace chi
{

class ConstantMaterialProperty : public MaterialProperty2
{
public:
  static InputParameters GetInputParameters();
  explicit ConstantMaterialProperty(const InputParameters& params);

  std::vector<std::string> RequiredInputNames() const override;
  double
  ComputeScalarValue(const std::vector<double>& input_params) const override;

protected:
  double scalar_value_;
};

} // namespace chi

#endif // CHITECH_CONSTANTMATERIALPROPERTY_H
