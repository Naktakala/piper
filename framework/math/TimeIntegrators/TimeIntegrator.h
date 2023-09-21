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
  /**For some time integrations we can compute the time and non-time
  * residual, at t+1, at the same time. This function allows us to control
  * whether or not that is the case.*/
  virtual std::vector<TimeID> GetResidualTimeIDsNeeded() const = 0;
  virtual size_t NumberOfSolutionHistoriesRequired() const = 0;
  virtual size_t NumberOfResidualHistoriesRequired() const = 0;
  virtual double GetTimeCoefficient(double dt) const = 0;

  virtual void
  ComputeResidual(ParallelVector& r,
                  const ParallelVector& time_residual,
                  const std::map<TimeID, const ParallelVector*>& std_residuals) = 0;
  virtual void ComputeJacobian(const ParallelVector& x,
                               ParallelMatrix& J,
                               EquationSystem& equation_system) = 0;

protected:
  static chi::InputParameters GetInputParameters();
  explicit TimeIntegrator(const chi::InputParameters& params);
};

} // namespace chi_math

#endif // CHITECH_TIMEINTEGRATOR_H
