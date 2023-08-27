#ifndef PIPER_JUNCTIONLIQUIDMODEL_H
#define PIPER_JUNCTIONLIQUIDMODEL_H

#include "ComponentLiquidModel.h"

namespace piper
{

class JunctionLiquidModel : public ComponentLiquidModel
{
public:
  static chi::InputParameters GetInputParameters();
  JunctionLiquidModel(const chi::InputParameters& params,
                      std::vector<std::unique_ptr<ComponentModel>>& family,
                      LiquidPhysics& liquid_physics,
                      const HardwareComponent& hardware_component,
                      const chi_mesh::Cell* cell,
                      const std::vector<std::string>& variable_names);

  void AssembleEquations() override;

protected:
  double forward_loss_coefficient_;
  double reverse_loss_coefficient_;
};

} // namespace piper

#endif // PIPER_JUNCTIONLIQUIDMODEL_H
