#include "SteadyNonLinearExecutioner.h"

#include "math/Systems/EquationSystem.h"

namespace chi_math
{

chi::InputParameters SteadyNonLinearExecutioner::GetInputParameters()
{
  chi::InputParameters params = NonLinearExecutioner::GetInputParameters();

  return params;
}

SteadyNonLinearExecutioner::SteadyNonLinearExecutioner(
  const chi::InputParameters& params,
  std::shared_ptr<EquationSystem> equation_system)
  : NonLinearExecutioner(params, std::move(equation_system))
{
}

void SteadyNonLinearExecutioner::ComputeResidual(const GhostedParallelVector& x,
                                                 ParallelVector& r)
{
  eq_system_->ComputeResidual(x, r);
}

void SteadyNonLinearExecutioner::ComputeJacobian(const GhostedParallelVector& x,
                                                 ParallelMatrix& J)
{
  eq_system_->ComputeJacobian(x, J);
}

} // namespace chi_math