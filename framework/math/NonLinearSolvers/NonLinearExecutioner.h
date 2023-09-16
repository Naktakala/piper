#ifndef CHITECH_NONLINEAREXECUTIONER_H
#define CHITECH_NONLINEAREXECUTIONER_H

#include "ChiObject.h"

namespace chi_math
{

class ParallelVector;
class GhostedParallelVector;
class ParallelMatrix;
struct ParallelMatrixSparsityPattern;

class EquationSystem;
class EquationSystemTimeData;

/**Abstract Non-Linear executioner class that interfaces with
 * BasicNonLinearSolver.*/
class NonLinearExecutioner : public ChiObject
{
public:
  /**Returns the number of local DOFs across all unknowns.*/
  virtual int64_t NumLocalDOFs() const;
  /**Returns the number of local DOFs across all unknowns.*/
  virtual int64_t NumGlobalDOFs() const;

  /**Returns a reference to the current solution vector.*/
  virtual GhostedParallelVector& SolutionVector();
  ///**Returns a reference to the current residual vector.*/
  //virtual ParallelVector& ResidualVector();

  /**Sets the current solution vector.*/
  virtual void SetInitialSolution();

  virtual void ComputeResidual(const GhostedParallelVector& x,
                               ParallelVector& r) = 0;
  virtual void ComputeJacobian(const GhostedParallelVector& x,
                               ParallelMatrix& J) = 0;

  virtual ParallelMatrixSparsityPattern BuildMatrixSparsityPattern() const;

  virtual void
  AddToPreConditionerOptions(std::vector<std::string>& pc_options){};

  void SetTimeData(EquationSystemTimeData time_data);

  virtual void Advance(EquationSystemTimeData time_data);

protected:
  static chi::InputParameters GetInputParameters();
  explicit NonLinearExecutioner(
    const chi::InputParameters& params,
    std::shared_ptr<EquationSystem> equation_system);

  std::shared_ptr<EquationSystem> eq_system_;
};

} // namespace chi_math

#endif // CHITECH_NONLINEAREXECUTIONER_H
