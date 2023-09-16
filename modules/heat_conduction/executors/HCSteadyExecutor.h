#ifndef PIPER_HCSTEADYEXECUTOR_H
#define PIPER_HCSTEADYEXECUTOR_H

#include "physics/SolverBase/chi_solver.h"
#include "math/UnknownManager/unknown_manager.h"
#include "math/NonLinearSolver/non_linear_solver.h"

#include <petscsnes.h>

namespace chi_mesh
{
class MeshContinuum;
}

namespace chi_math
{
class SpatialDiscretization;
class GhostedParallelVector;
class FEMKernelSystem;
class NonLinearExecutioner;
} // namespace chi_math

namespace hcm
{

class HeatConductionSystem;

class HCSteadyExecutor : public chi_physics::Solver
{
public:
  static chi::InputParameters GetInputParameters();
  explicit HCSteadyExecutor(const chi::InputParameters& params);

  void Initialize() override;
  virtual std::unique_ptr<chi_math::NonLinearExecutioner> SetExecutioner();
  void Execute() override;

protected:
  std::shared_ptr<HeatConductionSystem> hc_system_;

  const double T_initial_value_;

  const chi::ParameterBlock nl_solver_params_;

  chi_math::UnknownManager T_uk_man_;

  std::unique_ptr<chi_math::GhostedParallelVector> T_old_;
  std::shared_ptr<chi_math::FEMKernelSystem> fem_kernel_system_;
  std::unique_ptr<chi_math::NonLinearExecutioner> nl_executioner_;
  std::unique_ptr<chi_math::NonLinearSolver<Mat, Vec, SNES>> nl_solver_;
  std::map<std::string, std::shared_ptr<chi_physics::FieldFunctionGridBased>>
    name_2_ff_map_;
};

} // namespace hcm

#endif // PIPER_HCSTEADYEXECUTOR_H
