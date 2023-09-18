#include "KernelSystem.h"

#include "chi_log.h"

namespace chi_math
{

KernelSystem::KernelSystem(int64_t num_local_dofs,
                           int64_t num_globl_dofs,
                           const std::vector<int64_t>& ghost_ids,
                           TimeID oldest_time_id)
  : EquationSystem(num_local_dofs, num_globl_dofs, ghost_ids, oldest_time_id)
{
}


} // namespace chi_math