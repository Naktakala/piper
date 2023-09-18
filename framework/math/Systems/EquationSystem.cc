#include "EquationSystem.h"

#include "chi_log.h"

namespace chi_math
{

EquationSystem::EquationSystem(int64_t num_local_dofs,
                               int64_t num_globl_dofs,
                               const std::vector<int64_t>& ghost_ids,
                               TimeID oldest_time_id)
  : num_local_dofs_(num_local_dofs),
    num_globl_dofs_(num_globl_dofs),
    ghost_ids_(ghost_ids),
    main_solution_vector_(
      num_local_dofs, num_globl_dofs, ghost_ids, Chi::mpi.comm),
    num_old_blocks_(static_cast<int>(oldest_time_id) + 1),
    eq_term_scope_(EqTermScope::DOMAIN_TERMS | EqTermScope::BOUNDARY_TERMS)
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

/**Returns the number of local DOFs across all unknowns.*/
int64_t EquationSystem::NumLocalDOFs() const { return num_local_dofs_; }
/**Returns the number of local DOFs across all unknowns.*/
int64_t EquationSystem::NumGlobalDOFs() const { return num_globl_dofs_; }

/**Returns a reference to the current solution vector.*/
GhostedParallelVector&
EquationSystem::SolutionVector(TimeID time_id /*=TimeID::T_PLUS_1*/)
{
  if (time_id == TimeID::T_PLUS_1) return main_solution_vector_;

  const size_t index = static_cast<int>(time_id);
  ChiLogicalErrorIf(index >= old_solution_vectors_.size(),
                    "Solution vector not available at " + TimeIDName(time_id));

  return *old_solution_vectors_[index];
}

/**Returns a reference to the current residual vector.*/
ParallelVector&
EquationSystem::ResidualVector(TimeID time_id /*=TimeID::T_PLUS_1*/)
{
  ChiInvalidArgumentIf(time_id == TimeID::T_PLUS_1,
                       "Residual at T_PLUS_1 not available.");

  const size_t index = static_cast<int>(time_id);
  ChiLogicalErrorIf(index >= old_residual_vectors_.size(),
                    "Residual vector not available at " + TimeIDName(time_id));

  return *old_residual_vectors_[index];
}

/**Sets the current time data.*/
void EquationSystem::SetTimeData(EquationSystemTimeData time_data)
{
  time_data_ = time_data;
}

/**Returns the current equation terms-scope.*/
EqTermScope EquationSystem::EquationTermsScope() const
{
  return eq_term_scope_;
}

/**Sets the scope of current equations.*/
void EquationSystem::SetEquationTermsScope(EqTermScope eq_term_scope)
{
  eq_term_scope_ = eq_term_scope;
}

bool EquationSystem::QueryTermsActive(EqTermScope query_scope) const
{
  return eq_term_scope_ & query_scope;
}

/**Returns a reference to the system current time data.*/
const EquationSystemTimeData& EquationSystem::GetTimeData() const
{
  return time_data_;
}

/**Advances the system in time.*/
void EquationSystem::Advance(EquationSystemTimeData time_data,
                             const ParallelVector& latest_std_residual)
{
  Chi::log.Log0Error() << "Advanced";
  SetTimeData(time_data);

  const auto& sol_vec = main_solution_vector_;
  const auto& res_vec = latest_std_residual;

  if (num_old_blocks_ > 0)
  {
    {
      auto& vec = old_solution_vectors_;
      auto start = vec.begin();
      vec.insert(start, std::make_unique<GhostedParallelVector>(sol_vec));
      vec.pop_back();
    }
    {
      auto& vec = old_residual_vectors_;
      auto start = vec.begin();
      vec.insert(start, std::make_unique<ParallelVector>(res_vec));
      vec.pop_back();
    }
  }
}

} // namespace chi_math