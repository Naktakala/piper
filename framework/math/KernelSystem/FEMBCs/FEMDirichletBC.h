#ifndef PIPER_FEMDIRICHLETBC_H
#define PIPER_FEMDIRICHLETBC_H

#include "FEMBoundaryCondition.h"

namespace chi_math
{

class FEMDirichletBC : public FEMBoundaryCondition
{
public:
  static chi::InputParameters GetInputParameters();
  explicit FEMDirichletBC(const chi::InputParameters& params);

  bool IsDirichlet() const override;

  double ComputeLocalResidual(uint32_t i) override;

  double ResidualEntryAtQP() override;

  double BCValue() const;
  bool AllowApplyBeforeSolve() const;

protected:
  double bc_value_ = 0.0;

  double var_value_ = 0.0;
  chi_mesh::Vector3 node_xyz_;

private:
  bool apply_before_solve_ = true;
};

} // namespace chi_math

#endif // PIPER_FEMDIRICHLETBC_H
