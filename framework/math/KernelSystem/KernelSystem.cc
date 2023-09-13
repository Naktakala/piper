#include "KernelSystem.h"

namespace chi_math
{

KernelSystem::KernelSystem(int64_t num_local_dofs,
                           int64_t num_globl_dofs,
                           const std::vector<int64_t>& ghost_ids)
  : num_local_dofs_(num_local_dofs),
    num_globl_dofs_(num_globl_dofs),
    solution_vector_(num_local_dofs, num_globl_dofs, ghost_ids, Chi::mpi.comm),
    residual_vector_(num_local_dofs, num_globl_dofs, Chi::mpi.comm)
{
}

int64_t KernelSystem::NumLocalDOFs() const { return num_local_dofs_; }

int64_t KernelSystem::NumGlobalDOFs() const { return num_globl_dofs_; }

} // namespace chi_math