#include "EquationSystem.h"

#include "math/TimeIntegrators/ImplicitEulerTimeIntegrator.h"

#include "ChiObjectFactory.h"

namespace chi_math
{

chi::InputParameters EquationSystem::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddOptionalParameter("verbosity", 0, "Level of verbosity");
  params.AddOptionalParameter("time_integrator",
                              0,
                              "Handle to a TimeIntegrator object. A default "
                              "ImplicitEuler time integrator will"
                              "be used for transient simulations.");

  return params;
}

EquationSystem::EquationSystem(const chi::InputParameters& params)
  : ChiObject(params),
    verbosity_(params.GetParamValue<int>("verbosity")),
    time_integrator_(InitTimeIntegrator(params)),
    num_solution_histories_(
      time_integrator_->NumberOfSolutionHistoriesRequired()),
    num_residual_histories_(
      time_integrator_->NumberOfResidualHistoriesRequired()),
    time_data_(0.01, 0.0, 1.0),
    eq_term_scope_(EqTermScope::DOMAIN_TERMS | EqTermScope::BOUNDARY_TERMS)
{
}

// ##################################################################
std::shared_ptr<TimeIntegrator>
EquationSystem::InitTimeIntegrator(const chi::InputParameters& params)
{
  const auto& params_at_assignment = params.ParametersAtAssignment();

  size_t handle;
  if (params_at_assignment.Has("time_integrator"))
    handle = params.GetParamValue<size_t>("time_integrator");
  else
  {
    chi::InputParameters valid_params =
      ImplicitEulerTimeIntegrator::GetInputParameters();

    // Has no parameters to set

    auto& factory = ChiObjectFactory::GetInstance();
    handle = factory.MakeRegisteredObjectOfType(
      "chi_math::ImplicitEulerTimeIntegrator", valid_params);
  }

  auto time_integrator = Chi::GetStackItemPtrAsType<TimeIntegrator>(
    Chi::object_stack, handle, "EquationSystem::InitTimeIntegrator");

  return time_integrator;
}

/**Returns a reference to the system's time integrator.*/
TimeIntegrator& EquationSystem::GetTimeIntegrator()
{
  ChiLogicalErrorIf(not time_integrator_, "Time integrator object is nullptr.");
  return *time_integrator_;
}

/**Returns a reference to the system current time data.*/
const EquationSystemTimeData& EquationSystem::GetTimeData() const
{
  return time_data_;
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

/**Returns a reference to the current solution vector.*/
ParallelVector&
EquationSystem::SolutionVector(TimeID time_id /*=TimeID::T_PLUS_1*/)
{
  if (time_id == TimeID::T_PLUS_1) return *main_solution_vector_;

  const size_t index = static_cast<int>(time_id);
  ChiLogicalErrorIf(index >= old_solution_vectors_.size(),
                    "Solution vector not available at " + TimeIDName(time_id));

  return *old_solution_vectors_[index];
}

/**Returns a const reference to the current solution vector.*/
const ParallelVector&
EquationSystem::SolutionVector(TimeID time_id /*=TimeID::T_PLUS_1*/) const
{
  if (time_id == TimeID::T_PLUS_1) return *main_solution_vector_;

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

/**Advances the system in time.*/
void EquationSystem::Advance(
  EquationSystemTimeData time_data,
  std::map<TimeID, const ParallelVector*>& latest_std_residuals)
{
  for (auto& callback : callbacks_on_pre_advance_)
    callback();

  SetTimeData(time_data);

  const auto& sol_vec = main_solution_vector_;

  if (num_solution_histories_ > 0)
  {
    {
      auto& vec = old_solution_vectors_;
      auto start = vec.begin();
      vec.insert(start, sol_vec->MakeCopy());
      vec.pop_back();
    }
  }

  if (num_residual_histories_ > 0)
  {
    for (const auto& [time_id, vec_ptr] : latest_std_residuals)
    {
      const int time_index = static_cast<int>(time_id) + 1;
      if (time_index >= 0 and time_index < num_residual_histories_)
        old_residual_vectors_.at(time_index) = vec_ptr->MakeCopy();
    }
  }
}

chi::ParameterBlock
EquationSystem::GetInfo(const chi::ParameterBlock& params) const
{
  return chi::ParameterBlock{};
}

} // namespace chi_math