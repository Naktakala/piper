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
  eq_system_->SetEquationTermsScope(EqTermScope::TIME_TERMS |
                                    EqTermScope::DOMAIN_TERMS |
                                    EqTermScope::BOUNDARY_TERMS);

  auto& time_integrator = eq_system_->GetTimeIntegrator();

  const auto time_ids = time_integrator.GetTimeIDsNeeded();

  auto time_residual = r.MakeNewVector();

  std::vector<const ParallelVector*> residuals;
  for (const TimeID time_id : time_ids)
  {
    if (time_id == TimeID::T_PLUS_1)
    {
      residual_tp1_ = r.MakeNewVector();
      eq_system_->ComputeResidual(x, *residual_tp1_);
      residuals.push_back(&(*residual_tp1_));
    }
    else
      residuals.push_back(&eq_system_->ResidualVector(time_id));
  }

  time_integrator.ComputeResidual(r, *time_residual, residuals);
}

void TransientNonLinearExecutioner::ComputeJacobian(const ParallelVector& x,
                                                    ParallelMatrix& J)
{
  eq_system_->SetEquationTermsScope(EqTermScope::TIME_TERMS |
                                    EqTermScope::DOMAIN_TERMS |
                                    EqTermScope::BOUNDARY_TERMS);
  eq_system_->GetTimeIntegrator().ComputeJacobian(x, J, *eq_system_);
}

void TransientNonLinearExecutioner::Step()
{
  Chi::log.Log() << time_step_controller_->StringTimeInfo();

  SetTimeData(
    {time_step_controller_->GetTimeStepSize(), time_step_controller_->Time()});

  nl_solver_->Solve();

  Chi::log.Log() << nl_solver_->GetConvergedReasonString() << "\n\n";
}

void TransientNonLinearExecutioner::Advance()
{
  const bool last_solve_converged = nl_solver_->IsConverged();
  if (last_solve_converged)
  {
    time_step_controller_->Advance();
    eq_system_->Advance({time_step_controller_->GetTimeStepSize(),
                         time_step_controller_->Time()}, *residual_tp1_);
    eq_system_->UpdateFields();
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