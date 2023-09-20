#ifndef CHITECH_FEMDIRICHLETBC_H
#define CHITECH_FEMDIRICHLETBC_H

#include "FEMBoundaryCondition.h"

namespace chi_math
{

class FEMDirichletBC : public FEMBoundaryCondition
{
public:
  static chi::InputParameters GetInputParameters();
  explicit FEMDirichletBC(const chi::InputParameters& params);

  bool IsDirichlet() const override;

  double ComputeLocalResidual(size_t f, uint32_t i) override;
  double ComputeLocalJacobian(size_t f, uint32_t i, uint32_t j) override;

  double ResidualEntryAtQP() override;

  double BCValue() const;
  bool AllowApplyBeforeSolve() const;

protected:
  double bc_value_ = 0.0;

private:
  bool apply_before_solve_ = true;
};

} // namespace chi_math

#endif // CHITECH_FEMDIRICHLETBC_H
