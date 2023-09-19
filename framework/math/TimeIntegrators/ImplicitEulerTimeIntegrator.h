#ifndef CHITECH_IMPLICITEULERTIMEINTEGRATOR_H
#define CHITECH_IMPLICITEULERTIMEINTEGRATOR_H

#include "TimeIntegrator.h"

namespace chi_math
{

class ImplicitEulerTimeIntegrator : public TimeIntegrator
{
public:
  static chi::InputParameters GetInputParameters();
  explicit ImplicitEulerTimeIntegrator(const chi::InputParameters& params);

  std::vector<TimeID> GetTimeIDsNeeded() const override;

  void ComputeResidual(
    ParallelVector& r,
    const ParallelVector& time_residual,
    const std::vector<const ParallelVector*>& std_residuals) override;

  void ComputeJacobian(const ParallelVector& x,
                       ParallelMatrix& J,
                       EquationSystem& equation_system) override;
};

} // namespace chi_math

#endif // CHITECH_IMPLICITEULERTIMEINTEGRATOR_H
