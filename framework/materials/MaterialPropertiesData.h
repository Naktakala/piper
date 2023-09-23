#ifndef CHITECH_MATERIALPROPERTIESDATA_H
#define CHITECH_MATERIALPROPERTIESDATA_H

#include "ChiObject.h"

namespace chi
{

class MaterialProperty;

class MaterialPropertiesData : public ChiObject
{
public:
  static InputParameters GetInputParameters();
  explicit MaterialPropertiesData(const InputParameters& params);
  MaterialPropertiesData(){}; // Empty constructor

  static std::shared_ptr<MaterialPropertiesData> MakeEmptyData();

  const std::map<size_t, std::shared_ptr<const MaterialProperty>>&
  PropertyMap() const;

protected:
  std::map<size_t, std::shared_ptr<const MaterialProperty>> property_map_;

private:
  /**Builds a map of handles to property_ptrs from a list of handles.*/
  static std::map<size_t, std::shared_ptr<const MaterialProperty>>
  BuildPropertyMap(const ParameterBlock& property_param);
};

} // namespace chi

#endif // CHITECH_MATERIALPROPERTIESDATA_H
