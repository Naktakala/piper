#include "BasicNonLinearSolver.h"

#include "math/PETScUtils/petsc_utils.h"
#include "math/NonLinearSolvers/NonLinearExecutioner.h"
#include "math/PETScUtils/petsc_snes_utils.h"
#include "math/ParallelMatrix/ParallelPETScMatrixProxy.h"
#include "math/ParallelVector/ghosted_parallel_vector.h"

#include "chi_runtime.h"
#include "chi_log.h"

#define CheckContext(x)                                                        \
  if (not x)                                                                   \
  throw std::runtime_error(std::string(__PRETTY_FUNCTION__) +                  \
                           ": context casting failure")
#define GetKernelBasedContextPtr(x)                                            \
  std::dynamic_pointer_cast<BasicNLSolverContext>(x);                          \
  CheckContext(x)

namespace chi_math
{

// ##################################################################
BasicNLSolverContext::BasicNLSolverContext(NonLinearExecutioner& executioner)
  : executioner_(executioner)
{
}

// ##################################################################
BasicNonLinearSolver::BasicNonLinearSolver(NonLinearExecutioner& executioner,
                                           const chi::InputParameters& params)
  : chi_math::NonLinearSolver<Mat, Vec, SNES>(
      std::make_shared<BasicNLSolverContext>(executioner), params)
{
}

// ##################################################################
// void BasicNonLinearSolver::UpdateSolution(
//  const std::vector<double>& stl_vector)
//{
//  ChiLogicalErrorIf(not x_,
//                    "Solution vector not initialized. This method can only"
//                    "be made after a call to Setup().");
//  ChiLogicalErrorIf(
//    stl_vector.size() < num_local_dofs_,
//    "The STL-vector's local-size does not match that of the solution
//    vector.");
//
//  PetscScalar* x_raw;
//  VecGetArray(x_, &x_raw);
//  for (size_t i = 0; i < num_local_dofs_; ++i)
//    x_raw[i] = stl_vector[i];
//  VecRestoreArray(x_, &x_raw);
//}

// ##################################################################
void BasicNonLinearSolver::SetMonitor()
{
  SNESMonitorSet(nl_solver_,
                 &chi_math::PETScUtils::BasicSNESMonitor,
                 PETSC_VIEWER_STDOUT_WORLD,
                 nullptr);

  KSP ksp;
  SNESGetKSP(nl_solver_, &ksp);
  KSPMonitorSet(
    ksp, &chi_math::PETScUtils::KSPMonitorStraight, nullptr, nullptr);
}

// ##################################################################
void BasicNonLinearSolver::SetPreconditioner()
{
  if (options_.nl_method_ == "NEWTON" or options_.nl_method_ == "PJFNK")
  {
    auto& pc_params = options_.pc_options_;

    KSP ksp;
    SNESGetKSP(nl_solver_, &ksp);
    PC pc;
    KSPGetPC(ksp, &pc);
    KSPSetPCSide(ksp, PC_RIGHT);

    const auto pc_type_str = pc_params.GetParamValue<std::string>("pc_type");

    PCSetType(pc, pc_type_str.c_str());

    std::vector<std::string> pc_options;
    if (pc_type_str == "hypre")
    {
      const auto pc_hypre_type =
        pc_params.GetParamValue<std::string>("pc_hypre_type");
      PCHYPRESetType(pc, pc_hypre_type.c_str());

      pc_options = {"pc_hypre_boomeramg_strong_threshold 0.7",
                    "pc_hypre_boomeramg_agg_nl 0",
                    "pc_hypre_boomeramg_agg_num_paths 1",
                    "pc_hypre_boomeramg_max_levels 25",
                    //"pc_hypre_boomeramg_coarsen_type HMIS",
                    "pc_hypre_boomeramg_coarsen_type Falgout",
                    //"pc_hypre_boomeramg_interp_type ext+i",
                    "pc_hypre_boomeramg_interp_type classical",
                    "pc_hypre_boomeramg_P_max 4",
                    "pc_hypre_boomeramg_truncfactor 0.0",
                    "pc_mg_galerkin_mat_product_algorithm hypre"};

      auto nl_context_ptr = GetKernelBasedContextPtr(context_ptr_);
      nl_context_ptr->executioner_.AddToPreConditionerOptions(pc_options);
    } // if hypre

    for (const auto& option : pc_options)
      PetscOptionsInsertString(nullptr, ("-" + solver_name_ + option).c_str());

    for (const auto& user_option : pc_params)
    {
      const std::string option_string =
        user_option.Name() + " " + user_option.GetValue<std::string>();
      PetscOptionsInsertString(nullptr,
                               ("-" + solver_name_ + option_string).c_str());
    }

    PCSetFromOptions(pc);
  } // if NEWTON or PJFNK
}

// ##################################################################
void BasicNonLinearSolver::SetSystemSize()
{
  auto kb_context =
    std::static_pointer_cast<BasicNLSolverContext>(context_ptr_);
  num_local_dofs_ = kb_context->executioner_.NumLocalDOFs();
  num_globl_dofs_ = kb_context->executioner_.NumGlobalDOFs();
}

// ##################################################################
void BasicNonLinearSolver::SetSystem()
{
  //============================================= Create the vectors
  x_ = chi_math::PETScUtils::CreateVector(num_local_dofs_, num_globl_dofs_);
  VecDuplicate(x_, &r_);
}

// ##################################################################
void BasicNonLinearSolver::SetFunction()
{
  SNESSetFunction(nl_solver_, r_, ResidualFunction, this);
}

// ##################################################################
void BasicNonLinearSolver::SetJacobian()
{
  if (options_.nl_method_ == "JFNK")
  {
    MatCreateSNESMF(nl_solver_, &J_);
    SNESSetJacobian(nl_solver_, J_, J_, MatMFFDComputeJacobian, this);
    SNESSetUseMatrixFree(nl_solver_, PETSC_FALSE, PETSC_TRUE);
  }
  else if (options_.nl_method_ == "PJFNK")
  {
    SetupMatrix(P_);
    MatCreateSNESMF(nl_solver_, &J_);
    SNESSetJacobian(nl_solver_, J_, P_, ComputeJacobian, this);
    SNESSetUseMatrixFree(nl_solver_, PETSC_TRUE, PETSC_FALSE);
  }
  else if (options_.nl_method_ == "NEWTON")
  {
    SetupMatrix(J_);
    SNESSetJacobian(nl_solver_, J_, J_, ComputeJacobian, this);
  }
  else
    ChiInvalidArgument("Unsupported nl_method \"" + options_.nl_method_ +
                       "\".");

  // MatSNESMFSetReuseBase(J_, PETSC_TRUE);
  SNESSetLagJacobian(nl_solver_, 1);
}

// ##################################################################
void BasicNonLinearSolver::PostSetupCallback()
{
  Chi::log.Log0Verbose1() << "BasicNonLinearSolver setup completed.";
}

// ##################################################################
void BasicNonLinearSolver::SetInitialGuess()
{
  auto context = std::static_pointer_cast<BasicNLSolverContext>(context_ptr_);
  context->executioner_.SetInitialSolution();

  auto& solution_vector = context->executioner_.SolutionVector();

  chi_math::PETScUtils::CopySTLvectorToVec(
    solution_vector.RawValues(), x_, solution_vector.LocalSize());
}

// ##################################################################
PetscErrorCode
BasicNonLinearSolver::ResidualFunction(SNES snes, Vec x, Vec r, void*)
{
  BasicNLSolverContext* nl_context_ptr;
  SNESGetApplicationContext(snes, &nl_context_ptr);

  auto& executioner = nl_context_ptr->executioner_;
  auto solution_vector = executioner.SolutionVector();
  auto residual_vector = solution_vector;

  chi_math::PETScUtils::CopyVecToSTLvector(x,
                                           solution_vector.RawValues(),
                                           solution_vector.LocalSize(),
                                           /*resize_STL=*/false);
  solution_vector.CommunicateGhostEntries();

  residual_vector.Set(0.0);

  executioner.ComputeResidual(solution_vector, residual_vector);

  chi_math::PETScUtils::CopySTLvectorToVec(
    residual_vector.RawValues(), r, residual_vector.LocalSize());

  return 0;
}

// ##################################################################
PetscErrorCode BasicNonLinearSolver::ComputeJacobian(
  SNES snes, Vec x, Mat Jmat, Mat Pmat, void* ctx)
{
  auto solver_ptr = static_cast<BasicNonLinearSolver*>(ctx);
  BasicNLSolverContext* nl_context_ptr;
  SNESGetApplicationContext(snes, &nl_context_ptr);

  auto& executioner = nl_context_ptr->executioner_;
  auto solution_vector = executioner.SolutionVector();

  chi_math::PETScUtils::CopyVecToSTLvector(x,
                                           solution_vector.RawValues(),
                                           solution_vector.LocalSize(),
                                           /*resize_STL=*/false);
  solution_vector.CommunicateGhostEntries();

  auto& options = solver_ptr->options_;

  Mat matrix_to_assemble;

  /*if (options.nl_method_ == "JFNK") Never happens*/

  // For PJFNK The Jacobian is done with Matrix-Free actions using
  // residual evaluations whilst the preconditioner is assembled as the
  // Jacobian
  if (options.nl_method_ == "PJFNK")
  {
    MatMFFDComputeJacobian(solver_ptr->nl_solver_, x, Jmat, nullptr, nullptr);
    matrix_to_assemble = Pmat;
  }
  // For NEWTON the jacobian and the preconditioner is the same therefore
  // we don't need to touch Pmat
  else if (options.nl_method_ == "NEWTON")
    matrix_to_assemble = Jmat;
  else
    ChiInvalidArgument("Unsupported nl_method \"" + options.nl_method_ + "\".");

  MatZeroEntries(matrix_to_assemble);

  ParallelPETScMatrixProxy J_proxy(matrix_to_assemble,
                                   executioner.NumLocalDOFs(),
                                   executioner.NumLocalDOFs(),
                                   executioner.NumGlobalDOFs(),
                                   executioner.NumGlobalDOFs(),
                                   Chi::mpi.comm);

  executioner.ComputeJacobian(solution_vector, J_proxy);

  J_proxy.Assemble();

  return 0;
}

// ##################################################################
void BasicNonLinearSolver::PostSolveCallback()
{
  auto context = std::static_pointer_cast<BasicNLSolverContext>(context_ptr_);

  auto& executioner = context->executioner_;
  auto& solution_vector = executioner.SolutionVector();

  chi_math::PETScUtils::CopyVecToSTLvector(x_,
                                           solution_vector.RawValues(),
                                           solution_vector.LocalSize(),
                                           /*resize_STL=*/false);

  solution_vector.CommunicateGhostEntries();

  // SNESView(nl_solver_, PETSC_VIEWER_STDOUT_WORLD);
}

// ##################################################################
void BasicNonLinearSolver::SetupMatrix(Mat& A)
{
  auto nl_context_ptr = GetKernelBasedContextPtr(context_ptr_);

  auto& executioner = nl_context_ptr->executioner_;

  const auto sparsity_pattern = executioner.BuildMatrixSparsityPattern();

  A = chi_math::PETScUtils::CreateSquareMatrix(executioner.NumLocalDOFs(),
                                               executioner.NumGlobalDOFs());

  chi_math::PETScUtils::InitMatrixSparsity(
    A,
    sparsity_pattern.nodal_nnz_in_diag_,
    sparsity_pattern.nodal_nnz_off_diag_);
}

} // namespace chi_math