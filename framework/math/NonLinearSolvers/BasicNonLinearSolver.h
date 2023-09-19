#ifndef CHITECH_BASICNONLINEARSOLVER_H
#define CHITECH_BASICNONLINEARSOLVER_H

#include "math/NonLinearSolver/non_linear_solver.h"
#include "math/NonLinearSolver/NonLinearSolverOptions.h"

#include <petscsnes.h>
#include <vector>

namespace chi_math
{

class NonLinearExecutioner;

struct BasicNLSolverContext : public chi_math::NonLinearSolverContext<Vec, SNES>
{
  explicit BasicNLSolverContext(NonLinearExecutioner& executioner);

  virtual ~BasicNLSolverContext() = default;

  NonLinearExecutioner& executioner_;
};

class BasicNonLinearSolver : public chi_math::NonLinearSolver<Mat, Vec, SNES>
{
public:
  explicit BasicNonLinearSolver(NonLinearExecutioner& executioner,
                                const chi::InputParameters& params);

protected:
  void SetMonitor() override;
  void SetPreconditioner() override;

  void SetSystemSize() override;
  void SetSystem() override;
  void SetFunction() override;
  void SetJacobian() override;
  void PostSetupCallback() override;

  void SetInitialGuess() override;

  void PostSolveCallback() override;

  static PetscErrorCode ResidualFunction(SNES snes, Vec x, Vec r, void* ctx);
  static PetscErrorCode
  ComputeJacobian(SNES snes, Vec x, Mat Jmat, Mat Pmat, void* ctx);



private:
  void SetupMatrix(Mat& A);
};

} // namespace chi_math

#endif // CHITECH_BASICNONLINEARSOLVER_H
