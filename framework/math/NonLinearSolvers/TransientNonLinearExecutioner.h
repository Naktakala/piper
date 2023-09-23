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
  explicit TransientNonLinearExecutioner(const chi::InputParameters& params);

  void ComputeResidual(const ParallelVector& x, ParallelVector& r) override;
  void ComputeJacobian(const ParallelVector& x, ParallelMatrix& J) override;

  /**Sets the current solution vector.*/
  virtual void SetInitialSolution() override;

  void Step() override;
  void Advance() override;
  void Execute() override;

protected:
  std::unique_ptr<ParallelVector> r_t_tp1_;
  std::unique_ptr<ParallelVector> r_x_tp1_;
  std::map<TimeID, const ParallelVector*> r_map_;
  bool initial_solution_set_ = false;
};

} // namespace chi_math

#endif // CHITECH_TRANSIENTNONLINEAREXECUTIONER_H
