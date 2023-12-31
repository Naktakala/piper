#ifndef PIPER_LIQUIDPHYSICS_H
#define PIPER_LIQUIDPHYSICS_H

#include "FluidPhysics.h"
#include "piper/models/ComponentLiquidModel.h"
#include "math/PETScUtils/petsc_utils.h"
#include "math/ParallelVector/GhostedParallelSTLVector.h"

namespace chi_math::PETScUtils
{
struct PETScSolverSetup;
}

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

  typedef std::pair<std::string, double> StateVal;
  typedef std::vector<StateVal> StateValsList;
  typedef std::map<std::string, double> StateValMap;

  StateValMap EvaluateState(const std::vector<std::string>& vals_wanted,
                            const StateValsList& state_vals_list);

protected:
  const std::string fluid_name_;

  Mat A_ = nullptr;
  Vec x_ = nullptr;
  Vec b_ = nullptr;
  chi_math::PETScUtils::PETScSolverSetup pressure_solver_setup;

  std::unique_ptr<chi_math::GhostedParallelSTLVector> pressure_vector_;
};

} // namespace piper

#endif // PIPER_LIQUIDPHYSICS_H
