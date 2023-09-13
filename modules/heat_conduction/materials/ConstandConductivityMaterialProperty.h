#ifndef PIPER_CONSTANDCONDUCTIVITYMATERIALPROPERTY_H
#define PIPER_CONSTANDCONDUCTIVITYMATERIALPROPERTY_H

#include "materials/chi_material_property.h"

namespace hcm
{

class ConstantConductivityMaterialProperty : public chi::MaterialProperty
{
public:
  static chi::InputParameters GetInputParameters();
  explicit ConstantConductivityMaterialProperty(
    const chi::InputParameters& params);

protected:
  double k_;
};

} // namespace hcm

#endif // PIPER_CONSTANDCONDUCTIVITYMATERIALPROPERTY_H
