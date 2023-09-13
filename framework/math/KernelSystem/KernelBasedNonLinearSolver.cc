#include "KernelBasedNonLinearSolver.h"

#include "math/PETScUtils/petsc_utils.h"
#include "math/KernelSystem/FEMKernelSystem.h"
#include "math/PETScUtils/petsc_snes_utils.h"
#include "math/ParallelMatrix/ParallelPETScMatrixProxy.h"
#include "math/SpatialDiscretization/spatial_discretization.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

#define CheckContext(x)                                                        \
  if (not x)                                                                   \
  throw std::runtime_error(std::string(__PRETTY_FUNCTION__) +                  \
                           ": context casting failure")
#define GetKernelBasedContextPtr(x)                                            \
  std::dynamic_pointer_cast<KernelBasedContext>(x);                            \
  CheckContext(x)

namespace chi_math
{

// ##################################################################
KernelBasedNonLinearSolver::KernelBasedNonLinearSolver(
  FEMKernelSystem& kernel_system, const chi::InputParameters& params)
  : chi_math::NonLinearSolver<Mat, Vec, SNES>(
      std::make_shared<KernelBasedContext>(kernel_system), params)
{
}

// ##################################################################
void KernelBasedNonLinearSolver::UpdateSolution(
  const std::vector<double>& stl_vector)
{
  ChiLogicalErrorIf(not x_,
                    "Solution vector not initialized. This method can only"
                    "be made after a call to Setup().");
  ChiLogicalErrorIf(
    stl_vector.size() < num_local_dofs_,
    "The STL-vector's local-size does not match that of the solution vector.");

  PetscScalar* x_raw;
  VecGetArray(x_, &x_raw);
  for (size_t i = 0; i < num_local_dofs_; ++i)
    x_raw[i] = stl_vector[i];
  VecRestoreArray(x_, &x_raw);
}

// ##################################################################
void KernelBasedNonLinearSolver::SetMonitor()
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
void KernelBasedNonLinearSolver::SetPreconditioner()
{
  if (options_.nl_method_ == "NEWTON" or options_.nl_method_ == "PJFNK")
  {
    auto& pc_params = options_.pc_options_;

    KSP ksp;
    SNESGetKSP(nl_solver_, &ksp);
    PC pc;
    KSPGetPC(ksp, &pc);

    const auto pc_type_str = pc_params.GetParamValue<std::string>("pc_type");

    PCSetType(pc, pc_type_str.c_str());

    std::vector<std::string> pc_options;
    if (pc_type_str == "hypre")
    {
      const auto pc_hypre_type =
        pc_params.GetParamValue<std::string>("pc_hypre_type");
      PCHYPRESetType(pc, pc_hypre_type.c_str());

      pc_options = {
        //"pc_hypre_boomeramg_agg_nl 1",
        "pc_hypre_boomeramg_P_max 4",
        "pc_hypre_boomeramg_grid_sweeps_coarse 1",
        "pc_hypre_boomeramg_max_levels 25",
        "pc_hypre_boomeramg_relax_type_all symmetric-SOR/Jacobi",
        "pc_hypre_boomeramg_coarsen_type HMIS",
        "pc_hypre_boomeramg_interp_type ext+i"};

      auto nl_context_ptr = GetKernelBasedContextPtr(context_ptr_);

      auto& kernel_system = nl_context_ptr->kernel_system_;

      const auto& grid = kernel_system.SDM().Grid();
      if (grid.Attributes() & chi_mesh::DIMENSION_1)
        pc_options.emplace_back("pc_hypre_boomeramg_strong_threshold 0.5");
      if (grid.Attributes() & chi_mesh::DIMENSION_2)
        pc_options.emplace_back("pc_hypre_boomeramg_strong_threshold 0.6");
      if (grid.Attributes() & chi_mesh::DIMENSION_3)
        pc_options.emplace_back("pc_hypre_boomeramg_strong_threshold 0.7");
    }

    for (const auto& option : pc_options)
      PetscOptionsInsertString(nullptr, ("-" + solver_name_ + option).c_str());

    PCSetFromOptions(pc);
  }
}

// ##################################################################
void KernelBasedNonLinearSolver::SetSystemSize()
{
  auto kb_context = std::static_pointer_cast<KernelBasedContext>(context_ptr_);
  num_local_dofs_ = kb_context->kernel_system_.NumLocalDOFs();
  num_globl_dofs_ = kb_context->kernel_system_.NumGlobalDOFs();
}

// ##################################################################
void KernelBasedNonLinearSolver::SetSystem()
{
  //============================================= Create the vectors
  x_ = chi_math::PETScUtils::CreateVector(num_local_dofs_, num_globl_dofs_);
  VecDuplicate(x_, &r_);
}

// ##################################################################
void KernelBasedNonLinearSolver::SetFunction()
{
  SNESSetFunction(nl_solver_, r_, ResidualFunction, this);
}

// ##################################################################
void KernelBasedNonLinearSolver::SetJacobian()
{
  if (options_.nl_method_ == "JFNK")
  {
    MatCreateSNESMF(nl_solver_, &J_);
    SNESSetJacobian(nl_solver_, J_, J_, MatMFFDComputeJacobian, this);
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
}

// ##################################################################
void KernelBasedNonLinearSolver::PostSetupCallback()
{
  Chi::log.Log0Verbose1() << "KernelBasedNonLinearSolver setup completed.";
}

// ##################################################################
void KernelBasedNonLinearSolver::SetInitialGuess()
{
  auto context = std::static_pointer_cast<KernelBasedContext>(context_ptr_);
  context->kernel_system_.SetInitialSolution();

  auto& solution_vector = context->kernel_system_.SolutionVector();

  chi_math::PETScUtils::CopySTLvectorToVec(
    solution_vector.RawValues(), x_, solution_vector.LocalSize());
}

// ##################################################################
PetscErrorCode
KernelBasedNonLinearSolver::ResidualFunction(SNES snes, Vec phi, Vec r, void*)
{
  KernelBasedContext* nl_context_ptr;
  SNESGetApplicationContext(snes, &nl_context_ptr);

  auto& kernel_system = nl_context_ptr->kernel_system_;
  auto& solution_vector = kernel_system.SolutionVector();
  auto& residual_vector = kernel_system.ResidualVector();

  chi_math::PETScUtils::CopyVecToSTLvector(phi,
                                           solution_vector.RawValues(),
                                           solution_vector.LocalSize(),
                                           /*resize_STL=*/false);
  solution_vector.CommunicateGhostEntries();

  residual_vector.Set(0.0);

  kernel_system.ComputeResidual(residual_vector);

  chi_math::PETScUtils::CopySTLvectorToVec(
    residual_vector.RawValues(), r, residual_vector.LocalSize());

  return 0;
}

// ##################################################################
PetscErrorCode KernelBasedNonLinearSolver::ComputeJacobian(
  SNES snes, Vec x, Mat Jmat, Mat Pmat, void* ctx)
{
  auto solver_ptr = static_cast<KernelBasedNonLinearSolver*>(ctx);
  KernelBasedContext* nl_context_ptr;
  SNESGetApplicationContext(snes, &nl_context_ptr);

  auto& kernel_system = nl_context_ptr->kernel_system_;
  auto& solution_vector = kernel_system.SolutionVector();

  chi_math::PETScUtils::CopyVecToSTLvector(x,
                                           solution_vector.RawValues(),
                                           solution_vector.LocalSize(),
                                           /*resize_STL=*/false);
  solution_vector.CommunicateGhostEntries();

  auto& options = solver_ptr->options_;

  Mat ref_mat;
  if (options.nl_method_ == "NEWTON") ref_mat = Jmat;
  else
  {
    ref_mat = Pmat;
    MatMFFDComputeJacobian(solver_ptr->nl_solver_, x, Jmat, nullptr, nullptr);
  }

  MatZeroEntries(ref_mat);

  ParallelPETScMatrixProxy J_proxy(ref_mat,
                                   kernel_system.NumLocalDOFs(),
                                   kernel_system.NumLocalDOFs(),
                                   kernel_system.NumGlobalDOFs(),
                                   kernel_system.NumGlobalDOFs(),
                                   Chi::mpi.comm);

  kernel_system.ComputeJacobian(J_proxy);

  J_proxy.Assemble();

  return 0;
}

// ##################################################################
KernelBasedNonLinearSolver::KernelBasedContext::KernelBasedContext(
  FEMKernelSystem& kernel_system)
  : kernel_system_(kernel_system)
{
}

// ##################################################################
void KernelBasedNonLinearSolver::PostSolveCallback()
{
  auto context = std::static_pointer_cast<KernelBasedContext>(context_ptr_);

  auto& kernel_system = context->kernel_system_;
  auto& solution_vector = kernel_system.SolutionVector();

  chi_math::PETScUtils::CopyVecToSTLvector(
    x_, solution_vector.RawValues(), solution_vector.LocalSize(),
    /*resize_STL=*/false);
}

// ##################################################################
void KernelBasedNonLinearSolver::SetupMatrix(Mat& A)
{
  auto nl_context_ptr = GetKernelBasedContextPtr(context_ptr_);

  auto& kernel_system = nl_context_ptr->kernel_system_;
  auto& sdm = kernel_system.SDM();
  auto& uk_man = kernel_system.UnknownStructure();

  A = chi_math::PETScUtils::CreateSquareMatrix(kernel_system.NumLocalDOFs(),
                                               kernel_system.NumGlobalDOFs());

  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm.BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, uk_man);

  chi_math::PETScUtils::InitMatrixSparsity(
    A, nodal_nnz_in_diag, nodal_nnz_off_diag);
}

} // namespace chi_math