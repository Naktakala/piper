#ifndef CHI_KERNELSYSTEM_H
#define CHI_KERNELSYSTEM_H

#include "math/Systems/EquationSystem.h"
#include "math/ParallelVector/ghosted_parallel_vector.h"
#include "mesh/chi_mesh.h"

namespace chi_math
{

enum ActiveKernels : int
{
  NONE = 0,
  TIME_KERNELS = (1 << 0),
  STD_KERNELS = (2 << 1),
  BNDRY_KERNELS = (3 << 2)
};

inline ActiveKernels operator|(const ActiveKernels f1, const ActiveKernels f2)
{
  return static_cast<ActiveKernels>(static_cast<int>(f1) |
                                    static_cast<int>(f2));
}
inline bool operator&(const ActiveKernels f1, const ActiveKernels f2)
{
  return static_cast<int>(f1) & static_cast<int>(f2);
}

/**Abstract KernelSystem.*/
class KernelSystem : public EquationSystem
{
public:
  /**Returns true if time-kernels are active.*/
  bool AreTimeKernelsActive() const;
  /**Controls whether time kernels are active or not.*/
  void SetActiveKernels(ActiveKernels active_kernels);

protected:
  KernelSystem(int64_t num_local_dofs,
               int64_t num_globl_dofs,
               const std::vector<int64_t>& ghost_ids,
               TimeID oldest_time_id);

  void ComputeTimeResidual(const GhostedParallelVector& x,
                           ParallelVector& r) override;
  void ComputeTimeJacobian(const GhostedParallelVector& x,
                           ParallelMatrix& J) override;

  ActiveKernels active_kernels_;
};

} // namespace chi_math

#endif // CHI_KERNELSYSTEM_H
