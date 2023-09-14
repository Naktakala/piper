#ifndef CHI_FEMTIMEKERNEL_H
#define CHI_FEMTIMEKERNEL_H

#include "FEMKernel.h"

namespace chi_math
{

struct FEMTimeKernelRefData : public FEMKernelRefData
{
  FEMTimeKernelRefData(
    const std::shared_ptr<const finite_element::InternalQuadraturePointData>&
      qp_data,
    const VecDbl& var_values,
    const VecVec3& var_grad_values,
    const VecDbl& var_values_old,
    double dt,
    double time,
    double time_old);

  const std::vector<double>& var_values_old_;
  double time_ = 0.0;
  double time_old_ = 0.0;
  double dt_ = 1.0;
};

class FEMTimeKernel : public FEMKernel
{
public:
  static chi::InputParameters GetInputParameters();
  explicit FEMTimeKernel(const chi::InputParameters& params);

  void SetTimeReferenceData(FEMTimeKernelRefData* ref_data_ptr);

  double ComputeLocalResidual(uint32_t i) override;
  double ComputeLocalJacobian(uint32_t i, uint32_t j) override;

protected:
  double var_value_old_qp_ = 0.0;
  double time_ = 0.0;
  double time_old_ = 0.0;
  double dt_ = 1.0;

  FEMTimeKernelRefData* time_ref_data_ptr_ = nullptr;
};

} // namespace chi_math

#endif // CHI_FEMTIMEKERNEL_H
