#ifndef CHITECH_TRANSIENTNONLINEAREXECUTIONER_H
#define CHITECH_TRANSIENTNONLINEAREXECUTIONER_H

#include "NonLinearExecutioner.h"

namespace chi_math
{

class TimeIntegrator;

class TransientNonLinearExecutioner : public NonLinearExecutioner
{
public:
  static chi::InputParameters GetInputParameters();
  explicit TransientNonLinearExecutioner(
    const chi::InputParameters& params,
    std::shared_ptr<EquationSystem> equation_system);

  void ComputeResidual(const ParallelVector& x,
                       ParallelVector& r) override;
  void ComputeJacobian(const ParallelVector& x,
                       ParallelMatrix& J) override;

  void Advance(EquationSystemTimeData time_data) override;

protected:
  std::shared_ptr<TimeIntegrator> time_integrator_;

  std::unique_ptr<ParallelVector> residual_tp1_;
};

} // namespace chi_math

#endif // CHITECH_TRANSIENTNONLINEAREXECUTIONER_H
