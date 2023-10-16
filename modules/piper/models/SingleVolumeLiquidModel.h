#ifndef PIPER_SINGLEVOLUMELIQUIDMODEL_H
#define PIPER_SINGLEVOLUMELIQUIDMODEL_H

#include "ComponentLiquidModel.h"
#include "VolumetricComponentModel.h"

namespace piper
{

class SingleVolumeLiquidModel : public ComponentLiquidModel,
                                public VolumetricComponentModel
{
public:
  static chi::InputParameters GetInputParameters();
  SingleVolumeLiquidModel(const chi::InputParameters& params,
                          std::vector<std::unique_ptr<ComponentModel>>& family,
                          LiquidPhysics& liquid_physics,
                          const HardwareComponent& hardware_component,
                          const chi_mesh::Cell* cell,
                          const std::vector<std::string>& variable_names);

  void AssembleEquations() override;
};

} // namespace piper

#endif // PIPER_SINGLEVOLUMELIQUIDMODEL_H
