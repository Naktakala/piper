#include "MaterialPropertiesData.h"

#include "materials/MaterialProperty.h"

#include "ChiObjectFactory.h"

namespace chi
{

RegisterChiObject(chi, MaterialPropertiesData);

InputParameters MaterialPropertiesData::GetInputParameters()
{
  InputParameters params = ChiObject::GetInputParameters();

  params.SetGeneralDescription(
    "A block where material properties are populated");

  params.AddRequiredParameterArray(
    "properties", "An array of handles to material properties.");

  return params;
}

MaterialPropertiesData::MaterialPropertiesData(const InputParameters& params)
  : ChiObject(params),
    properties_(
      std::move(AssemblePropertiesList(params.GetParam("properties"))))
{
}

std::shared_ptr<MaterialPropertiesData> MaterialPropertiesData::MakeEmptyData()
{
  return std::make_shared<MaterialPropertiesData>();
}

const std::vector<std::shared_ptr<const MaterialProperty2>>&
MaterialPropertiesData::Properties() const
{
  return properties_;
}

std::vector<std::shared_ptr<const MaterialProperty2>>
MaterialPropertiesData::AssemblePropertiesList(
  const ParameterBlock& property_param)
{
  std::vector<std::shared_ptr<const MaterialProperty2>> properties;
  for (const auto& sub_param : property_param)
  {
    const size_t handle = sub_param.GetValue<size_t>();
    auto property_ptr = Chi::GetStackItemPtrAsType<MaterialProperty2>(
      Chi::object_stack, handle, __FUNCTION__);

    properties.push_back(property_ptr);
  }
  return properties;
}

} // namespace chi