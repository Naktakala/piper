#include "SteadyNonLinearExecutioner.h"

#include "math/Systems/FieldEquationSystem.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

namespace chi_math
{

RegisterChiObject(chi_math, SteadyNonLinearExecutioner);

chi::InputParameters SteadyNonLinearExecutioner::GetInputParameters()
{
  chi::InputParameters params = NonLinearExecutioner::GetInputParameters();

  return params;
}

SteadyNonLinearExecutioner::SteadyNonLinearExecutioner(
  const chi::InputParameters& params)
  : NonLinearExecutioner(params)
{
}

void SteadyNonLinearExecutioner::ComputeResidual(const ParallelVector& x,
                                                 ParallelVector& r)
{
  Chi::log.LogEvent(t_tag_residual_, chi::ChiLog::EventType::EVENT_BEGIN);
  SetModeToNonTimeOnly();
  eq_system_->ComputeResidual(x, r);
  Chi::log.LogEvent(t_tag_residual_, chi::ChiLog::EventType::EVENT_END);
}

void SteadyNonLinearExecutioner::ComputeJacobian(const ParallelVector& x,
                                                 ParallelMatrix& J)
{
  Chi::log.LogEvent(t_tag_jacobian_, chi::ChiLog::EventType::EVENT_BEGIN);
  SetModeToNonTimeOnly();
  eq_system_->ComputeJacobian(x, J);
  Chi::log.LogEvent(t_tag_jacobian_, chi::ChiLog::EventType::EVENT_END);
}

void SteadyNonLinearExecutioner::Execute()
{
  Chi::log.LogEvent(t_tag_solve_, chi::ChiLog::EventType::EVENT_BEGIN);

  if (print_header_)
    Chi::log.Log() << "\nExecuting solver \"" + TextName() + "\"";

  nl_solver_->Solve();
  eq_system_->UpdateFields();

  if (print_footer_)
    Chi::log.Log() << nl_solver_->GetConvergedReasonString() << "\n\n";

  Chi::log.LogEvent(t_tag_solve_, chi::ChiLog::EventType::EVENT_END);

  eq_system_->OutputFields(-1); // -1 = latest

  PrintTimingInfo();
}

} // namespace chi_math