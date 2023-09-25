#ifndef CHITECH_MATERIALPROPERTYINTERFACE_H
#define CHITECH_MATERIALPROPERTYINTERFACE_H

#include "parameters/input_parameters.h"
#include "math/KernelSystem/Coupling/FEMMaterialProperty.h"

namespace chi
{
class MaterialPropertiesData;
}

namespace chi_math
{

/**This interface expects a handle to a block of type MaterialPropertiesData
 * which it uses to establish a reference to that block*/
class FEMMaterialPropertyInterface
{
protected:
  static chi::InputParameters GetInputParameters();
  explicit FEMMaterialPropertyInterface(const chi::InputParameters& params,
                                        const std::vector<int>& mat_id_scope);

  void PreComputeInternalMaterialProperties();
  void PreComputeFaceMaterialProperties();

  const chi_math::FEMMaterialProperty&
  GetFEMMaterialProperty(const std::string& name);

private:
  static const chi::MaterialPropertiesData&
  GetFEMMaterialProperties(const chi::InputParameters& params);

  const size_t fem_data_handle_;

  const chi::MaterialPropertiesData& material_properties_data_;

  const std::vector<int>& mat_id_scope_;

  std::vector<std::unique_ptr<FEMMaterialProperty>> coupled_material_props_;
};

} // namespace chi_math

#endif // CHITECH_MATERIALPROPERTYINTERFACE_H
