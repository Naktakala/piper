#include "ConstantMaterialProperty.h"

#include "ChiObjectFactory.h"

namespace chi
{

RegisterChiObject(chi, ConstantMaterialProperty);

InputParameters ConstantMaterialProperty::GetInputParameters()
{
  InputParameters params = MaterialProperty2::GetInputParameters();

  params.SetGeneralDescription("A scalar material property that is constant.");
  params.SetDocGroup("doc_MaterialProperties");

  params.AddRequiredParameter<double>("scalar_value", "Value of the property.");

  return params;
}

ConstantMaterialProperty::ConstantMaterialProperty(
  const InputParameters& params)
  : MaterialProperty2(params),
    scalar_value_(params.GetParamValue<double>("scalar_value"))
{
}

double ConstantMaterialProperty::ComputeScalarValue(
  const std::vector<double>& input_params) const
{
  return scalar_value_;
}

std::vector<std::string> ConstantMaterialProperty::RequiredInputNames() const
{
  return {};
}

} // namespace chi