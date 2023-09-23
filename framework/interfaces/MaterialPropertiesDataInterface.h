#ifndef CHITECH_MATERIALPROPERTIESDATAINTERFACE_H
#define CHITECH_MATERIALPROPERTIESDATAINTERFACE_H

#include "parameters/input_parameters.h"

namespace chi
{

class MaterialPropertiesData;

class MaterialPropertiesDataInterface
{
public:
  static InputParameters GetInputParameters();
  explicit MaterialPropertiesDataInterface(const InputParameters& params);

protected:
  std::shared_ptr<MaterialPropertiesData> material_properties_data_;

private:
  static std::shared_ptr<MaterialPropertiesData>
  GetDataBlock(const InputParameters& params);
};

} // namespace chi

#endif // CHITECH_MATERIALPROPERTIESDATAINTERFACE_H
