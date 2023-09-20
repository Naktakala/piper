#include "SteadyNonLinearExecutioner.h"

#include "math/Systems/EquationSystem.h"

#include "ChiObjectFactory.h"

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
  eq_system_->SetEquationTermsScope(EqTermScope::DOMAIN_TERMS |
                                    EqTermScope::BOUNDARY_TERMS);
  eq_system_->ComputeResidual(x, r);
}

void SteadyNonLinearExecutioner::ComputeJacobian(const ParallelVector& x,
                                                 ParallelMatrix& J)
{
  eq_system_->SetEquationTermsScope(EqTermScope::DOMAIN_TERMS |
                                    EqTermScope::BOUNDARY_TERMS);
  eq_system_->ComputeJacobian(x, J);
}

void SteadyNonLinearExecutioner::Execute()
{
  nl_solver_->Solve();
  eq_system_->UpdateFields();
}

} // namespace chi_math