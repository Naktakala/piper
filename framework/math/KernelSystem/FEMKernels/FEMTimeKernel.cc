#include "FEMTimeKernel.h"

#include "math/SpatialDiscretization/FiniteElement/finite_element.h"

namespace chi_math
{

FEMTimeKernelRefData::FEMTimeKernelRefData(
  const std::shared_ptr<const finite_element::InternalQuadraturePointData>&
    qp_data,
  const VecDbl& var_values,
  const VecVec3& var_grad_values,
  const VecDbl& var_values_old,
  double dt,
  double time,
  double time_old)
  : FEMKernelRefData(qp_data, var_values, var_grad_values),
    var_values_old_(var_values_old),
    dt_(dt),
    time_(time),
    time_old_(time_old)
{
}

chi::InputParameters FEMTimeKernel::GetInputParameters()
{
  chi::InputParameters params = FEMKernel::GetInputParameters();

  return params;
}

FEMTimeKernel::FEMTimeKernel(const chi::InputParameters& params)
  : FEMKernel(params)
{
}

void FEMTimeKernel::SetTimeReferenceData(FEMTimeKernelRefData* ref_data_ptr)
{
  FEMKernel::SetReferenceData(ref_data_ptr_);
  time_ref_data_ptr_ = ref_data_ptr;

  time_ = time_ref_data_ptr_->time_;
  time_old_ = time_ref_data_ptr_->time_old_;
  dt_ = time_ref_data_ptr_->dt_;
}

double FEMTimeKernel::ComputeLocalResidual(uint32_t i)
{
  double local_r = 0.0;
  for (uint32_t qp : ref_data_ptr_->qp_indices_)
  {
    test_i_qp_ = ref_data_ptr_->shape_values_[i][qp];
    test_grad_i_qp_ = ref_data_ptr_->shape_grad_values_[i][qp];
    var_qp_value_ = ref_data_ptr_->var_qp_values_[qp];
    var_grad_qp_value_ = ref_data_ptr_->var_grad_qp_values_[qp];
    qp_xyz_ = ref_data_ptr_->qpoints_xyz_[qp];

    var_value_old_qp_ = time_ref_data_ptr_->var_values_old_[qp];

    local_r += ref_data_ptr_->JxW_values_[qp] * ResidualEntryAtQP();
  }

  return local_r;
}

double FEMTimeKernel::ComputeLocalJacobian(uint32_t i, uint32_t j)
{
  double local_j = 0.0;
  for (uint32_t qp : ref_data_ptr_->qp_indices_)
  {
    test_i_qp_ = ref_data_ptr_->shape_values_[i][qp];
    test_grad_i_qp_ = ref_data_ptr_->shape_grad_values_[i][qp];
    var_qp_value_ = ref_data_ptr_->var_qp_values_[qp];
    var_grad_qp_value_ = ref_data_ptr_->var_grad_qp_values_[qp];
    qp_xyz_ = ref_data_ptr_->qpoints_xyz_[qp];

    shape_j_qp_ = ref_data_ptr_->shape_values_[j][qp];
    shape_grad_j_qp_ = ref_data_ptr_->shape_grad_values_[j][qp];

    local_j += ref_data_ptr_->JxW_values_[qp] * JacobianEntryAtQP();
  }

  return local_j;
}

} // namespace chi_math