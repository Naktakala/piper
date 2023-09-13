#ifndef PIPER_KERNELSYSTEM_H
#define PIPER_KERNELSYSTEM_H

#include "math/ParallelVector/ghosted_parallel_vector.h"

namespace chi_math
{

class ParallelMatrix;

/**Abstract KernelSystem.*/
class KernelSystem
{
public:
  int64_t NumLocalDOFs() const;
  int64_t NumGlobalDOFs() const;

  GhostedParallelVector& SolutionVector() { return solution_vector_; }
  ParallelVector& ResidualVector() { return residual_vector_; }

  virtual void SetInitialSolution() {};

  virtual void ComputeResidual(ParallelVector& r) = 0;
  virtual void ComputeJacobian(ParallelMatrix& J) = 0;

protected:
  KernelSystem(int64_t num_local_dofs,
               int64_t num_globl_dofs,
               const std::vector<int64_t>& ghost_ids);

  const int64_t num_local_dofs_;
  const int64_t num_globl_dofs_;
  chi_math::GhostedParallelVector solution_vector_;
  chi_math::ParallelVector residual_vector_;
};

} // namespace chi_math

#endif // PIPER_KERNELSYSTEM_H
