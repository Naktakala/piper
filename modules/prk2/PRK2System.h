#ifndef PIPER_PRK2SYSTEM_H
#define PIPER_PRK2SYSTEM_H

#include "math/Systems/EquationSystem.h"
#include "math/dynamic_matrix.h"
#include "math/dynamic_vector.h"

namespace prk2
{

class PRK2System : public chi_math::EquationSystem
{
public:
  static chi::InputParameters GetInputParameters();
  explicit PRK2System(const chi::InputParameters& params);

  int64_t NumLocalDOFs() const override { return num_local_dofs_; }
  int64_t NumGlobalDOFs() const override { return num_global_dofs_; }

  void ComputeResidual(const chi_math::ParallelVector& x,
                       chi_math::ParallelVector& r) override;
  void ComputeJacobian(const chi_math::ParallelVector& x,
                       chi_math::ParallelMatrix& J) override;

  chi_math::ParallelMatrixSparsityPattern
  BuildMatrixSparsityPattern() const override;

  void SetInitialSolution() override;

  void PreAdvanceCallback();

  void SetProperties(const chi::ParameterBlock& params);

  chi::ParameterBlock GetInfo(const chi::ParameterBlock& params) const override;

private:
  void InitializeSystem();
  std::vector<double> lambdas_;
  std::vector<double> betas_;
  const double beta_;
  double gen_time_;
  double rho_;
  double source_strength_;

  size_t num_precursors_;

  const int64_t num_local_dofs_;
  const int64_t num_global_dofs_;

  chi_math::DynamicMatrix<double> A_;
  chi_math::DynamicVector<double> x_initial_, q_;

  double population_ = 0.0;
  double period_tph_ = 1.0e6;
};

} // namespace prk2

#endif // PIPER_PRK2SYSTEM_H
