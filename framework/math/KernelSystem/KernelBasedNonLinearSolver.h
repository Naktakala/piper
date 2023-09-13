#ifndef PIPER_KERNELBASEDNONLINEARSOLVER_H
#define PIPER_KERNELBASEDNONLINEARSOLVER_H

#include "math/NonLinearSolver/non_linear_solver.h"
#include "math/NonLinearSolver/NonLinearSolverOptions.h"

#include <petscsnes.h>
#include <vector>

namespace chi_math
{

class FEMKernelSystem;

class KernelBasedNonLinearSolver
  : public chi_math::NonLinearSolver<Mat, Vec, SNES>
{
public:
  explicit KernelBasedNonLinearSolver(
    FEMKernelSystem& kernel_system,
    const chi::InputParameters& params);

  void UpdateSolution(const std::vector<double>& stl_vector);

protected:
  struct KernelBasedContext : public chi_math::NonLinearSolverContext<Vec, SNES>
  {
    explicit KernelBasedContext(FEMKernelSystem& kernel_system);

    virtual ~KernelBasedContext() = default;

    FEMKernelSystem& kernel_system_;
  };

  void SetMonitor() override;
  void SetPreconditioner() override;

  void SetSystemSize() override;
  void SetSystem() override;
  void SetFunction() override;
  void SetJacobian() override;
  void PostSetupCallback() override;

  void SetInitialGuess() override;

  void PostSolveCallback() override;

  static PetscErrorCode ResidualFunction(SNES snes, Vec phi, Vec r, void* ctx);
  static PetscErrorCode
  ComputeJacobian(SNES snes, Vec x, Mat Jmat, Mat Pmat, void* ctx);

private:
  void SetupMatrix(Mat& A);
};

} // namespace chi_math

#endif // PIPER_KERNELBASEDNONLINEARSOLVER_H
