#ifndef CHITECH_NONLINEAREXECUTIONER_H
#define CHITECH_NONLINEAREXECUTIONER_H

#include "physics/SolverBase/chi_solver.h"
#include "math/NonLinearSolver/non_linear_solver.h"
#include "math/Systems/EquationSystemTimeData.h"

#include <petscsnes.h>

namespace chi_math
{

class ParallelVector;
class ParallelMatrix;
struct ParallelMatrixSparsityPattern;

class EquationSystem;
class EquationSystemTimeData;

class NLSolverFailedException : public std::exception
{
public:
  NLSolverFailedException() = default;
};

/**Abstract Non-Linear executioner class that interfaces with
 * BasicNonLinearSolver.*/
class NonLinearExecutioner : public chi_physics::Solver
{
public:
  /**Returns the number of local DOFs across all unknowns.*/
  virtual int64_t NumLocalDOFs() const;
  /**Returns the number of local DOFs across all unknowns.*/
  virtual int64_t NumGlobalDOFs() const;

  /**Returns a reference to the current solution vector.*/
  virtual ParallelVector& SolutionVector();

  /**Sets the current solution vector.*/
  virtual void SetInitialSolution();

  virtual void ComputeResidual(const ParallelVector& x, ParallelVector& r) = 0;
  virtual void ComputeJacobian(const ParallelVector& x, ParallelMatrix& J) = 0;

  virtual ParallelMatrixSparsityPattern BuildMatrixSparsityPattern() const;

  virtual void
  AddToPreConditionerOptions(std::vector<std::string>& pc_options){};

  void SetTimeData(EquationSystemTimeData time_data);

  void SetModeToTimeOnly();
  void SetModeToNonTimeOnly();
  void SetModeToTimeAndNonTime();

  void Initialize() override;

protected:
  static chi::InputParameters GetInputParameters();
  explicit NonLinearExecutioner(const chi::InputParameters& params);
  static std::shared_ptr<EquationSystem> GetEquationSystem(size_t handle);

  static bool TimeIDListHasID(const std::vector<TimeID>& time_ids, TimeID id);

  std::shared_ptr<EquationSystem> eq_system_;
  const chi::ParameterBlock nl_solver_params_;
  std::unique_ptr<chi_math::NonLinearSolver<Mat, Vec, SNES>> nl_solver_;
};

} // namespace chi_math

#endif // CHITECH_NONLINEAREXECUTIONER_H
