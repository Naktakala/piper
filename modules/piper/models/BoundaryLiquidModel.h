#ifndef PIPER_BOUNDARYLIQUIDMODEL_H
#define PIPER_BOUNDARYLIQUIDMODEL_H

#include "ComponentLiquidModel.h"

namespace piper
{

class BoundaryLiquidModel : public ComponentLiquidModel
{
public:
  static chi::InputParameters GetInputParameters();
  BoundaryLiquidModel(const chi::InputParameters& params,
                      std::vector<std::unique_ptr<ComponentModel>>& family,
                      LiquidPhysics& liquid_physics,
                      const HardwareComponent& hardware_component,
                      const chi_mesh::Cell* cell,
                      const std::vector<std::string>& variable_names);
};

} // namespace piper

#endif // PIPER_BOUNDARYLIQUIDMODEL_H
