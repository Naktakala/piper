#ifndef PIPER_HEATCONDUCTIONSTEADYSTATEEXECUTOR_H
#define PIPER_HEATCONDUCTIONSTEADYSTATEEXECUTOR_H

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
} // namespace chi_math

namespace hcm
{

class HeatConductionSystem;

class HeatConductionSteadyStateExecutor : public chi_physics::Solver
{
public:
  static chi::InputParameters GetInputParameters();
  explicit HeatConductionSteadyStateExecutor(
    const chi::InputParameters& params);

  void Initialize() override;
  void Execute() override;

protected:
  std::shared_ptr<HeatConductionSystem> hc_system_;

  const double T_initial_value_;

  const chi::ParameterBlock solver_params_;

  chi_math::UnknownManager T_uk_man_;

  size_t num_local_dofs_ = 0;
  size_t num_globl_dofs_ = 0;

  std::unique_ptr<chi_math::GhostedParallelVector> T_old_;

  std::unique_ptr<chi_math::FEMKernelSystem> fem_kernel_system_;

  std::unique_ptr<chi_math::NonLinearSolver<Mat, Vec, SNES>> nl_solver_;

  std::map<std::string, std::shared_ptr<chi_physics::FieldFunctionGridBased>>
    name_2_ff_map_;
};

} // namespace hcm

#endif // PIPER_HEATCONDUCTIONSTEADYSTATEEXECUTOR_H
