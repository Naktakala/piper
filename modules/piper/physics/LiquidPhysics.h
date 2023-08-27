#ifndef PIPER_LIQUIDPHYSICS_H
#define PIPER_LIQUIDPHYSICS_H

#include "FluidPhysics.h"
#include "piper/models/ComponentLiquidModel.h"

namespace piper
{

class ComponentModel;

class LiquidPhysics : public FluidPhysics
{
public:
  struct EqCoeffs
  {
    double rhs_;
    std::vector<std::vector<double>> coeff_sets_;
    std::vector<std::vector<size_t>> id_maps_;
  };
  static chi::InputParameters GetInputParameters();
  explicit LiquidPhysics(const chi::InputParameters& params);

  /**Returns the fluid name.*/
  const std::string& FluidName() const;

  void Initialize() override;
  void InitializeUnknowns() override;
  std::vector<std::string>
  MakeVariableNamesList(ComponentCategory hw_comp_category) override;
  void StaticGravityInitializer(const chi::ParameterBlock& params);

  ComponentLiquidModel& GetComponentLiquidModel(size_t component_id);

  void Step() override;
  void Advance() override;

protected:
  typedef std::pair<std::string, double> StateVal;
  typedef std::vector<StateVal> StateValsList;
  typedef std::map<std::string, double> StateValMap;

  StateValMap EvaluateState(const StateValsList& state_vals_list);

  const std::string fluid_name_;
};

} // namespace piper

#endif // PIPER_LIQUIDPHYSICS_H
