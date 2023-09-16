#include "KernelSystem.h"

#include "chi_log.h"

namespace chi_math
{

KernelSystem::KernelSystem(int64_t num_local_dofs,
                           int64_t num_globl_dofs,
                           const std::vector<int64_t>& ghost_ids,
                           TimeID oldest_time_id /*=TimeID::T_PLUS_1*/)
  : EquationSystem(num_local_dofs, num_globl_dofs, ghost_ids, oldest_time_id),
    active_kernels_(STD_KERNELS | BNDRY_KERNELS)
{
  for (int t = 0; t < num_old_blocks_; ++t)
  {
    auto new_sol_vec = std::make_unique<GhostedParallelVector>(
      num_local_dofs, num_globl_dofs, ghost_ids, Chi::mpi.comm);
    auto new_res_vec = std::make_unique<ParallelVector>(
      num_local_dofs, num_globl_dofs, Chi::mpi.comm);

    old_solution_vectors_.push_back(std::move(new_sol_vec));
    old_residual_vectors_.push_back(std::move(new_res_vec));
  }
}

/**Returns true if time-kernels are active.*/
bool KernelSystem::AreTimeKernelsActive() const
{
  return active_kernels_ & TIME_KERNELS;
}
/**Controls whether time kernels are active or not.*/
void KernelSystem::SetActiveKernels(ActiveKernels active_kernels)
{
  active_kernels_ = active_kernels;
}

void KernelSystem::ComputeTimeResidual(const GhostedParallelVector& x,
                                       ParallelVector& r)
{
  SetActiveKernels(TIME_KERNELS);
  ComputeResidual(x, r);
  SetActiveKernels(STD_KERNELS | BNDRY_KERNELS);
}

void KernelSystem::ComputeTimeJacobian(const GhostedParallelVector& x,
                                       ParallelMatrix& J)
{
  SetActiveKernels(TIME_KERNELS);
  ComputeJacobian(x, J);
  SetActiveKernels(STD_KERNELS | BNDRY_KERNELS);
}

} // namespace chi_math