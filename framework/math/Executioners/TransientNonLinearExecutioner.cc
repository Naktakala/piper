#include "TransientNonLinearExecutioner.h"

#include "math/Systems/FieldEquationSystem.h"
#include "math/TimeIntegrators/TimeIntegrator.h"

#include "physics/TimeSteppers/TimeStepper.h"
#include "physics/PhysicsEventPublisher.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiObjectFactory.h"

namespace chi_math
{

RegisterChiObject(chi_math, TransientNonLinearExecutioner);

chi::InputParameters TransientNonLinearExecutioner::GetInputParameters()
{
  chi::InputParameters params = NonLinearExecutioner::GetInputParameters();

  return params;
}

TransientNonLinearExecutioner::TransientNonLinearExecutioner(
  const chi::InputParameters& params)
  : NonLinearExecutioner(params)
{
}

void TransientNonLinearExecutioner::ComputeResidual(const ParallelVector& x,
                                                    ParallelVector& r)
{
  t_residual_.TimeSectionBegin();

  auto& time_integrator = eq_system_->GetTimeIntegrator();

  const auto r_time_ids = time_integrator.GetResidualTimeIDsNeeded();

  // Naming convention:
  // Time residual     : r^t, named r_t
  // Non-time residual : r(x), named r_x
  // Time notations    : tp1 = t+1, t=t, tm1 = t-1, tm2 = t-2 etc.

  r_map_.clear();

  r_t_tp1_ = r.MakeClone();
  SetModeToTimeOnly();
  eq_system_->ComputeResidual(x, *r_t_tp1_);

  if (TimeIDListHasID(r_time_ids, TimeID::T_PLUS_1))
  {
    r_x_tp1_ = r.MakeClone();
    SetModeToNonTimeOnly();
    eq_system_->ComputeResidual(x, *r_x_tp1_);
    r_map_[TimeID::T_PLUS_1] = (&(*r_x_tp1_));
  }

  for (const TimeID time_id : r_time_ids)
    if (time_id != TimeID::T_PLUS_1)
      r_map_[time_id] = (&eq_system_->ResidualVector(time_id));

  time_integrator.ComputeResidual(r, *r_t_tp1_, r_map_);

  // No longer needing r_t_tp1
  r_t_tp1_ = nullptr;

  // Delete r_x_tp1 if r_x_t is not required
  if (not TimeIDListHasID(r_time_ids, TimeID::T)) r_x_tp1_ = nullptr;

  t_residual_.TimeSectionEnd();
}

void TransientNonLinearExecutioner::ComputeJacobian(const ParallelVector& x,
                                                    ParallelMatrix& J)
{
  t_jacobian_.TimeSectionBegin();
  SetModeToTimeAndNonTime();
  eq_system_->GetTimeIntegrator().ComputeJacobian(x, J, *eq_system_);
  t_jacobian_.TimeSectionEnd();
}

void TransientNonLinearExecutioner::SetInitialSolution()
{
  if (not initial_solution_set_)
  {
    const double dt = timestepper_->TimeStepSize();
    const double time = timestepper_->Time();

    const double var_dot_dvar =
      eq_system_->GetTimeIntegrator().GetTimeCoefficient(dt);

    eq_system_->SetTimeData({dt, time, var_dot_dvar});

    eq_system_->SetInitialSolution();

    const auto r_time_ids =
      eq_system_->GetTimeIntegrator().GetResidualTimeIDsNeeded();

    if (TimeIDListHasID(r_time_ids, TimeID::T))
    {
      const auto& x = eq_system_->SolutionVector();
      auto r = x.MakeClone();
      SetModeToNonTimeOnly();
      eq_system_->ComputeResidual(x, *r);

      r_map_[TimeID::T_PLUS_1] = &(*r);
      eq_system_->Advance({dt, time, var_dot_dvar}, r_map_);

      eq_system_->UpdateFields();
    }

    initial_solution_set_ = true;
  }
}

void TransientNonLinearExecutioner::Step()
{
  t_solve_.TimeSectionBegin();

  const double dt = timestepper_->TimeStepSize();
  const double time = timestepper_->Time();

  if (print_header_)
    Chi::log.Log() << "Solver \"" + TextName() + "\" " +
                        timestepper_->StringTimeInfo(false);

  const double var_dot_dvar =
    eq_system_->GetTimeIntegrator().GetTimeCoefficient(dt);

  eq_system_->SetTimeData({dt, time + dt, var_dot_dvar});

  nl_solver_->Solve();

  if (print_footer_)
    Chi::log.Log() << nl_solver_->GetConvergedReasonString() << "\n\n";

  t_solve_.TimeSectionEnd();
}

void TransientNonLinearExecutioner::Advance()
{
  const bool last_solve_converged = nl_solver_->IsConverged();
  if (last_solve_converged)
  {
    timestepper_->Advance();

    const double dt = timestepper_->TimeStepSize();
    const double time = timestepper_->Time();
    const size_t t_index = timestepper_->TimeStepIndex();
    const double var_dot_dvar =
      eq_system_->GetTimeIntegrator().GetTimeCoefficient(dt);

    eq_system_->Advance({dt, time + dt, var_dot_dvar}, r_map_);

    eq_system_->UpdateFields();
    eq_system_->OutputFields(static_cast<int>(t_index));
  }
  else
  {
    if (not timestepper_->Adapt(chi_physics::TimeStepStatus::FAILURE))
      throw NLSolverFailedException();
  }
}

void TransientNonLinearExecutioner::Execute()
{
  auto& physics_event_publisher =
    chi_physics::PhysicsEventPublisher::GetInstance();

  if (not initial_solution_set_)
  {
    Chi::log.Log() << timestepper_->StringTimeInfo(true);
    Chi::log.Log() << "Setting initial solution";

    SetInitialSolution();
  }

  try
  {
    while (timestepper_->IsActive())
    {
      physics_event_publisher.SolverStep(*this);
      physics_event_publisher.SolverAdvance(*this);
    }
  }
  catch (const NLSolverFailedException&)
  {
    Chi::log.Log0Error() << "Solver failed: "
                         << timestepper_->StringTimeInfo(false);
  }
}

} // namespace chi_math