#ifndef CHITECH_MATERIALPROPERTIESDATA_H
#define CHITECH_MATERIALPROPERTIESDATA_H

#include "ChiObject.h"

namespace chi
{

class MaterialProperty2;

class MaterialPropertiesData : public ChiObject
{
public:
  static InputParameters GetInputParameters();
  explicit MaterialPropertiesData(const InputParameters& params);
  MaterialPropertiesData(){}; // Empty constructor

  static std::shared_ptr<MaterialPropertiesData> MakeEmptyData();

  const std::vector<std::shared_ptr<const MaterialProperty2>>&
  Properties() const;

protected:
  const std::vector<std::shared_ptr<const MaterialProperty2>> properties_;

private:
  /**Builds a map of handles to property_ptrs from a list of handles.*/
  static std::vector<std::shared_ptr<const MaterialProperty2>>
  AssemblePropertiesList(const ParameterBlock& property_param);
};

} // namespace chi

#endif // CHITECH_MATERIALPROPERTIESDATA_H
