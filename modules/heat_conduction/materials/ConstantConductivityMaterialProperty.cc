#include "ConstandConductivityMaterialProperty.h"

#include "ChiObjectFactory.h"

namespace hcm
{

RegisterChiObject(hcm, ConstantConductivityMaterialProperty);

chi::InputParameters ConstantConductivityMaterialProperty::GetInputParameters()
{
  chi::InputParameters params = chi::MaterialProperty::GetInputParameters();

  params.SetGeneralDescription("A material property where the thermal "
                               "conductivity is constant.");
  params.SetDocGroup("doc_HeatConduction");

  params.AddRequiredParameter<double>(
    "k", "Value for the constant thermal conductivity.");

  return params;
}

ConstantConductivityMaterialProperty::ConstantConductivityMaterialProperty(
  const chi::InputParameters& params)
  : chi::MaterialProperty(params), k_(params.GetParamValue<double>("k"))
{
}

} // namespace hcm