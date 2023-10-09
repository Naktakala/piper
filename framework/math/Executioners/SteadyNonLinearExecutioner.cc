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
  t_residual_.TimeSectionBegin();
  SetModeToNonTimeOnly();
  eq_system_->ComputeResidual(x, r);
  t_residual_.TimeSectionEnd();
}

void SteadyNonLinearExecutioner::ComputeJacobian(const ParallelVector& x,
                                                 ParallelMatrix& J)
{
  t_jacobian_.TimeSectionBegin();
  SetModeToNonTimeOnly();
  eq_system_->ComputeJacobian(x, J);
  t_jacobian_.TimeSectionEnd();
}

void SteadyNonLinearExecutioner::Execute()
{
  t_solve_.TimeSectionBegin();

  if (print_header_)
    Chi::log.Log() << "\nExecuting solver \"" + TextName() + "\"";

  nl_solver_->Solve();
  eq_system_->UpdateFields();

  if (print_footer_)
    Chi::log.Log() << nl_solver_->GetConvergedReasonString() << "\n\n";

  t_solve_.TimeSectionEnd();

  eq_system_->OutputFields(-1); // -1 = latest
}

} // namespace chi_math