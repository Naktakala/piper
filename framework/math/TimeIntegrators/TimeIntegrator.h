#ifndef CHITECH_TIMEINTEGRATOR_H
#define CHITECH_TIMEINTEGRATOR_H

#include "ChiObject.h"
#include "math/Systems/EquationSystemTimeData.h"

namespace chi_math
{

class ParallelVector;
class ParallelMatrix;
class EquationSystem;

/**Base class for time integrators*/
class TimeIntegrator : public ChiObject
{
public:
  virtual std::vector<TimeID> GetTimeIDsNeeded() const = 0;
  virtual size_t NumberOfSolutionHistoriesRequired() const = 0;
  virtual size_t NumberOfResidualHistoriesRequired() const = 0;

  virtual void
  ComputeResidual(ParallelVector& r,
                  const ParallelVector& time_residual,
                  const std::vector<const ParallelVector*>& std_residuals) = 0;
  virtual void ComputeJacobian(const ParallelVector& x,
                               ParallelMatrix& J,
                               EquationSystem& equation_system) = 0;

protected:
  static chi::InputParameters GetInputParameters();
  explicit TimeIntegrator(const chi::InputParameters& params);
};

} // namespace chi_math

#endif // CHITECH_TIMEINTEGRATOR_H
