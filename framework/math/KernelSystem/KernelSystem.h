#ifndef CHI_KERNELSYSTEM_H
#define CHI_KERNELSYSTEM_H

#include "math/Systems/EquationSystem.h"
#include "math/ParallelVector/ghosted_parallel_vector.h"
#include "mesh/chi_mesh.h"

namespace chi_math
{

/**Derived from an EquationSystem, this system computes items using Kernels.
* The type of kernels here is still abstract.*/
class KernelSystem : public EquationSystem
{
public:

protected:
  KernelSystem(int64_t num_local_dofs,
               int64_t num_globl_dofs,
               const std::vector<int64_t>& ghost_ids,
               TimeID oldest_time_id);
};

} // namespace chi_math

#endif // CHI_KERNELSYSTEM_H
