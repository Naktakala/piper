#ifndef CHITECH_STEADYNONLINEAREXECUTIONER_H
#define CHITECH_STEADYNONLINEAREXECUTIONER_H

#include "NonLinearExecutioner.h"

namespace chi_math
{

class SteadyNonLinearExecutioner : public NonLinearExecutioner
{
public:
  static chi::InputParameters GetInputParameters();
  explicit SteadyNonLinearExecutioner(const chi::InputParameters& params);

  void ComputeResidual(const ParallelVector& x, ParallelVector& r) override;
  void ComputeJacobian(const ParallelVector& x, ParallelMatrix& J) override;

  void Execute() override;
};

} // namespace chi_math

#endif // CHITECH_STEADYNONLINEAREXECUTIONER_H
