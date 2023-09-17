#include "TransientNonLinearExecutioner.h"

#include "math/Systems/EquationSystem.h"
#include "math/TimeIntegrators/TimeIntegrator.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_math
{

chi::InputParameters TransientNonLinearExecutioner::GetInputParameters()
{
  chi::InputParameters params = NonLinearExecutioner::GetInputParameters();

  params.AddRequiredParameter<size_t>("time_integrator",
                                      "Handle to a time integrator.");

  return params;
}

TransientNonLinearExecutioner::TransientNonLinearExecutioner(
  const chi::InputParameters& params,
  std::shared_ptr<EquationSystem> equation_system)
  : NonLinearExecutioner(params, std::move(equation_system)),
    time_integrator_(Chi::GetStackItemPtrAsType<TimeIntegrator>(
      Chi::object_stack,
      params.GetParamValue<size_t>("time_integrator"),
      __FUNCTION__))
{
}

void TransientNonLinearExecutioner::ComputeResidual(
  const GhostedParallelVector& x, ParallelVector& r)
{
  eq_system_->SetEquationTermsScope(EqTermScope::TIME_TERMS |
                                    EqTermScope::DOMAIN_TERMS |
                                    EqTermScope::BOUNDARY_TERMS);
  ParallelVector time_residual = r;

  const auto time_ids = time_integrator_->GetTimeIDsNeeded();

  std::vector<const ParallelVector*> residuals;
  for (const TimeID time_id : time_ids)
  {
    if (time_id == TimeID::T_PLUS_1)
    {
      residual_tp1_ = std::make_unique<ParallelVector>(r);
      eq_system_->ComputeResidual(x, *residual_tp1_);
      residuals.push_back(&(*residual_tp1_));
    }
    else
      residuals.push_back(&eq_system_->ResidualVector(time_id));
  }

  time_integrator_->ComputeResidual(r, time_residual, residuals);
}

void TransientNonLinearExecutioner::ComputeJacobian(
  const GhostedParallelVector& x, ParallelMatrix& J)
{
  eq_system_->SetEquationTermsScope(EqTermScope::TIME_TERMS |
                                    EqTermScope::DOMAIN_TERMS |
                                    EqTermScope::BOUNDARY_TERMS);
  time_integrator_->ComputeJacobian(x, J, *eq_system_);
}

void TransientNonLinearExecutioner::Advance(EquationSystemTimeData time_data)
{
  eq_system_->Advance(time_data, *residual_tp1_);
}

} // namespace chi_math