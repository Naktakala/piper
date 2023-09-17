#ifndef CHITECH_EQUATIONSYSTEM_H
#define CHITECH_EQUATIONSYSTEM_H

#include "EquationSystemTimeData.h"
#include "math/ParallelVector/ghosted_parallel_vector.h"
#include "mesh/chi_mesh.h"

namespace chi_math
{

class ParallelMatrix;
struct ParallelMatrixSparsityPattern;

enum class EqTermScope : int
{
  NONE = 0,
  TIME_TERMS = (1 << 0),
  DOMAIN_TERMS = (1 << 2),
  BOUNDARY_TERMS = (1 << 3)
};

inline EqTermScope operator|(const EqTermScope f1, const EqTermScope f2)
{
  return static_cast<EqTermScope>(static_cast<int>(f1) |
                                    static_cast<int>(f2));
}
inline bool operator&(const EqTermScope f1, const EqTermScope f2)
{
  return static_cast<int>(f1) & static_cast<int>(f2);
}

class EquationSystem
{
public:
  typedef chi_mesh::Vector3 Vec3;
  typedef std::vector<Vec3> VecVec3;

  /**Returns the number of local DOFs across all unknowns.*/
  int64_t NumLocalDOFs() const;
  /**Returns the number of local DOFs across all unknowns.*/
  int64_t NumGlobalDOFs() const;

  /**Returns a reference to the current solution vector.*/
  GhostedParallelVector& SolutionVector(TimeID time_id = TimeID::T_PLUS_1);
  /**Returns a reference to the current residual vector.*/
  ParallelVector& ResidualVector(TimeID time_id = TimeID::T_PLUS_1);

  /**Sets the current solution vector.*/
  virtual void SetInitialSolution(){};

  //virtual void ComputeTimeResidual(const GhostedParallelVector& x,
  //                                 ParallelVector& r) = 0;
  //virtual void ComputeTimeJacobian(const GhostedParallelVector& x,
  //                                 ParallelMatrix& J) = 0;
  virtual void ComputeResidual(const GhostedParallelVector& x,
                               ParallelVector& r) = 0;
  virtual void ComputeJacobian(const GhostedParallelVector& x,
                               ParallelMatrix& J) = 0;

  /**Sets the current time data.*/
  void SetTimeData(EquationSystemTimeData time_data);

  /**Returns the current equation terms-scope.*/
  EqTermScope EquationTermsScope() const;
  /**Sets the scope of current equations.*/
  void SetEquationTermsScope(EqTermScope eq_term_scope);

  /**Determines if a particular equation term is active.*/
  bool QueryTermsActive(EqTermScope query_scope) const;

  /**Returns a reference to the system current time data.*/
  const EquationSystemTimeData& GetTimeData() const;

  /**Uses the underlying system to build a sparsity pattern.*/
  virtual ParallelMatrixSparsityPattern BuildMatrixSparsityPattern() const = 0;

  /**Advances the system in time.*/
  void Advance(EquationSystemTimeData time_data,
               const ParallelVector& latest_std_residual);

protected:
  EquationSystem(int64_t num_local_dofs,
                 int64_t num_globl_dofs,
                 const std::vector<int64_t>& ghost_ids,
                 TimeID oldest_time_id);

  const int64_t num_local_dofs_;
  const int64_t num_globl_dofs_;
  const std::vector<int64_t> ghost_ids_;

  GhostedParallelVector main_solution_vector_;

  std::vector<std::unique_ptr<GhostedParallelVector>> old_solution_vectors_;
  std::vector<std::unique_ptr<ParallelVector>> old_residual_vectors_;

  size_t num_old_blocks_;
  EquationSystemTimeData time_data_;

private:
  EqTermScope eq_term_scope_;
};

} // namespace chi_math

#endif // CHITECH_EQUATIONSYSTEM_H
