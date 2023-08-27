#include "BoundaryLiquidModel.h"

#include "ChiObjectFactory.h"

namespace piper
{

RegisterChiObjectParametersOnly(piper, BoundaryLiquidModel);

chi::InputParameters BoundaryLiquidModel::GetInputParameters()
{
  return ComponentLiquidModel::GetInputParameters();
}

BoundaryLiquidModel::BoundaryLiquidModel(
  const chi::InputParameters& params,
  std::vector<std::unique_ptr<ComponentModel>>& family,
  LiquidPhysics& liquid_physics,
  const HardwareComponent& hardware_component,
  const chi_mesh::Cell* cell,
  const std::vector<std::string>& variable_names)
  : ComponentLiquidModel(
      params, family, liquid_physics, hardware_component, cell, variable_names)
{
}

} // namespace piper