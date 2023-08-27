#include "ComponentLiquidModel.h"

#include "chi_log.h"

namespace piper
{

chi::InputParameters ComponentLiquidModel::GetInputParameters()
{
  return chi::InputParameters{};
}

void ComponentLiquidModel::AssembleEquations()
{
  ChiLogicalError("Method not implemented.");
}

ComponentLiquidModel::ComponentLiquidModel(
  const chi::InputParameters& params,
  std::vector<std::unique_ptr<ComponentModel>>& family,
  LiquidPhysics& liquid_physics,
  const HardwareComponent& hardware_component,
  const chi_mesh::Cell* cell,
  const std::vector<std::string>& variable_names)
  : ComponentModel(family, hardware_component, cell, variable_names),
    physics_(liquid_physics)
{
}

const ComponentLiquidModel::EqCoeffs&
ComponentLiquidModel::GetEquationCoefficients(size_t id) const
{
  ChiLogicalErrorIf(id >= equation_coefficients_.size(), "Invalid equation id");
  return equation_coefficients_[id];
}

ComponentLiquidModel::EqCoeffs&
ComponentLiquidModel::GetEquationCoefficients(size_t id)
{
  ChiLogicalErrorIf(id >= equation_coefficients_.size(), "Invalid equation id");
  return equation_coefficients_[id];
}

} // namespace piper