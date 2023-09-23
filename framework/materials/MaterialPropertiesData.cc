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
    property_map_(std::move(BuildPropertyMap(params.GetParam("properties"))))
{
}

std::shared_ptr<MaterialPropertiesData> MaterialPropertiesData::MakeEmptyData()
{
  return std::make_shared<MaterialPropertiesData>();
}

const std::map<size_t, std::shared_ptr<const MaterialProperty>>&
MaterialPropertiesData::PropertyMap() const
{
  return property_map_;
}

std::map<size_t, std::shared_ptr<const MaterialProperty>>
MaterialPropertiesData::BuildPropertyMap(const ParameterBlock& property_param)
{
  std::map<size_t, std::shared_ptr<const MaterialProperty>> property_map;
  for (const auto& sub_param : property_param)
  {
    const size_t handle = sub_param.GetValue<size_t>();
    auto property_ptr = Chi::GetStackItemPtrAsType<MaterialProperty>(
      Chi::object_stack, handle, __FUNCTION__);

    property_map[handle] = property_ptr;
  }
  return property_map;
}

} // namespace chi