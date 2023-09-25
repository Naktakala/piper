#include "KernelSystem.h"

#include "materials/MaterialPropertiesData.h"
#include "materials/MaterialProperty.h"
#include "math/KernelSystem/Coupling/FEMMaterialProperty.h"

namespace chi_math
{

void KernelSystem::PopulateMaterialProperties()
{
  const auto& property_map = material_properties_data_->PropertyMap();
  for (const auto& [handle, property_ptr] : property_map)
    fem_material_properties_.push_back(
      std::make_shared<FEMMaterialProperty>(*property_ptr));
}

} // namespace chi_math