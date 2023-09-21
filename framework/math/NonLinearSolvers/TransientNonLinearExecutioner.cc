#include "TransientNonLinearExecutioner.h"

#include "math/Systems/EquationSystem.h"
#include "math/TimeIntegrators/TimeIntegrator.h"

#include "physics/TimeStepControllers/TimeStepController.h"

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
  auto& time_integrator = eq_system_->GetTimeIntegrator();

  const auto r_time_ids = time_integrator.GetResidualTimeIDsNeeded();

  // Naming convention:
  // Time residual     : r^t, named r_t
  // Non-time residual : r(x), named r_x
  // Time notations    : tp1 = t+1, t=t, tm1 = t-1, tm2 = t-2 etc.

  r_map_.clear();

  r_t_tp1_ = r.MakeNewVector();
  SetModeToTimeOnly();
  eq_system_->ComputeResidual(x, *r_t_tp1_);

  if (TimeIDListHasID(r_time_ids, TimeID::T_PLUS_1))
  {
    r_x_tp1_ = r.MakeNewVector();
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
}

void TransientNonLinearExecutioner::ComputeJacobian(const ParallelVector& x,
                                                    ParallelMatrix& J)
{
  SetModeToTimeAndNonTime();
  eq_system_->GetTimeIntegrator().ComputeJacobian(x, J, *eq_system_);
}

void TransientNonLinearExecutioner::SetInitialSolution()
{
  if (time_step_controller_->TimeIndex() == 0)
  {
    Chi::log.Log() << "Setting initial solution";
    const double dt = time_step_controller_->GetTimeStepSize();
    const double time = time_step_controller_->Time();
    const double var_dot_dvar =
      eq_system_->GetTimeIntegrator().GetTimeCoefficient(dt);

    eq_system_->SetTimeData({dt, time, var_dot_dvar});

    eq_system_->SetInitialSolution();

    const auto r_time_ids =
      eq_system_->GetTimeIntegrator().GetResidualTimeIDsNeeded();

    if (TimeIDListHasID(r_time_ids, TimeID::T))
    {
      const auto& x = eq_system_->SolutionVector();
      auto r = x.MakeNewVector();
      SetModeToNonTimeOnly();
      eq_system_->ComputeResidual(x, *r);

      r_map_[TimeID::T_PLUS_1] = &(*r);
      eq_system_->Advance({dt, time, var_dot_dvar},
                          r_map_);

      eq_system_->UpdateFields();
    }
  }
}

void TransientNonLinearExecutioner::Step()
{
  Chi::log.Log() << time_step_controller_->StringTimeInfo();

  const double dt = time_step_controller_->GetTimeStepSize();
  const double time = time_step_controller_->Time();
  const double var_dot_dvar =
    eq_system_->GetTimeIntegrator().GetTimeCoefficient(dt);

  eq_system_->SetTimeData({dt, time, var_dot_dvar});

  nl_solver_->Solve();

  Chi::log.Log() << nl_solver_->GetConvergedReasonString() << "\n\n";
}

void TransientNonLinearExecutioner::Advance()
{
  const bool last_solve_converged = nl_solver_->IsConverged();
  if (last_solve_converged)
  {
    time_step_controller_->Advance();

    const double dt = time_step_controller_->GetTimeStepSize();
    const double time = time_step_controller_->Time();
    const double var_dot_dvar =
      eq_system_->GetTimeIntegrator().GetTimeCoefficient(dt);

    eq_system_->Advance({dt, time, var_dot_dvar}, r_map_);

    eq_system_->UpdateFields();
    time_ = time;
    
  }
  else
  {
    if (not time_step_controller_->Adapt(chi_physics::TimeStepStatus::FAILURE))
    {
      Chi::log.Log0Error() << "Solver failed: "
                           << time_step_controller_->StringTimeInfo();
      Chi::Exit(EXIT_FAILURE);
    }
  }
}

} // namespace chi_math