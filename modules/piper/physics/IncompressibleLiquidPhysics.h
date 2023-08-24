#ifndef PIPER_INCOMPRESSIBLELIQUIDPHYSICS_H
#define PIPER_INCOMPRESSIBLELIQUIDPHYSICS_H

#include "FluidPhysics.h"

namespace piper
{

class ComponentModel;

class IncompressibleLiquidPhysics : public FluidPhysics
{
public:
  struct EqCoeffs
  {
    double rhs_;
    std::vector<std::vector<double>> coeff_sets_;
    std::vector<std::vector<size_t>> id_maps_;
  };
  static chi::InputParameters GetInputParameters();
  explicit IncompressibleLiquidPhysics(const chi::InputParameters& params);

  void InitializeUnknowns() override;
  std::vector<std::string>
  MakeVariableNamesList(ComponentCategory hw_comp_category) override;
  void StaticGravityInitializer(const chi::ParameterBlock& params);

  void Step() override;

protected:
  typedef std::pair<std::string, double> StateVal;
  typedef std::vector<StateVal> StateValsList;
  typedef std::map<std::string, double> StateValMap;

  StateValMap EvaluateState(const StateValsList& state_vals_list);

  static double FrictionFactor(const ComponentModel& model);

  const std::string fluid_name_;
};

} // namespace piper

#endif // PIPER_INCOMPRESSIBLELIQUIDPHYSICS_H
