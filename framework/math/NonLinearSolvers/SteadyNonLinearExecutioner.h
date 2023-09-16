#ifndef CHITECH_STEADYNONLINEAREXECUTIONER_H
#define CHITECH_STEADYNONLINEAREXECUTIONER_H

#include "NonLinearExecutioner.h"

namespace chi_math
{

class SteadyNonLinearExecutioner : public NonLinearExecutioner
{
public:
  static chi::InputParameters GetInputParameters();
  explicit SteadyNonLinearExecutioner(
    const chi::InputParameters& params,
    std::shared_ptr<EquationSystem> equation_system);

  void ComputeResidual(const GhostedParallelVector& x,
                       ParallelVector& r) override;
  void ComputeJacobian(const GhostedParallelVector& x,
                       ParallelMatrix& J) override;
};

} // namespace chi_math

#endif // CHITECH_STEADYNONLINEAREXECUTIONER_H
