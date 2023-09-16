#include "FEMKernel.h"

#include "math/Systems/EquationSystemTimeData.h"
#include "math/SpatialDiscretization/FiniteElement/finite_element.h"

#include "chi_log.h"

namespace chi_math
{

/**Constructor for reference data used by a kernel.*/
FEMKernelRefData::FEMKernelRefData(
  const EquationSystemTimeData& time_data,
  const std::shared_ptr<const finite_element::InternalQuadraturePointData>&
    qp_data_ptr_,
  VecDbl var_qp_values,
  VecVec3 var_grad_qp_values,
  MatDbl old_var_qp_values)
  : time_data_(time_data),
    qp_data_ptr_(qp_data_ptr_),
    qp_indices_(qp_data_ptr_->QuadraturePointIndices()),
    qpoints_xyz_(qp_data_ptr_->QPointsXYZ()),
    shape_values_(qp_data_ptr_->ShapeValues()),
    shape_grad_values_(qp_data_ptr_->ShapeGradValues()),
    JxW_values_(qp_data_ptr_->JxW_Values()),
    var_qp_values_(std::move(var_qp_values)),
    var_grad_qp_values_(std::move(var_grad_qp_values)),
    old_var_qp_values_(std::move(old_var_qp_values))
{
}

chi::InputParameters FEMKernel::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddOptionalParameterArray(
    "mat_ids",
    std::vector<chi::ParameterBlock>{},
    "A list of material-ids to which this block will be"
    "restricted");

  return params;
}

FEMKernel::FEMKernel(const chi::InputParameters& params) : ChiObject(params)
{
  const auto& user_params = params.ParametersAtAssignment();

  if (user_params.Has("mat_ids"))
    mat_ids_ = params.GetParamVectorValue<int>("mat_ids");
}

void FEMKernel::SetReferenceData(
  std::shared_ptr<FEMKernelRefData>& ref_data_ptr)
{
  ref_data_ptr_ = ref_data_ptr;
}

const std::vector<int>& FEMKernel::GetMaterialIDScope() const
{
  return mat_ids_;
}

bool FEMKernel::IsTimeKernel() const
{
  return false;
}

double FEMKernel::ComputeLocalResidual(uint32_t i)
{
  dt_ = ref_data_ptr_->time_data_.dt_;
  time_ = ref_data_ptr_->time_data_.time_;

  const size_t max_t = ref_data_ptr_->old_var_qp_values_.size();
  const auto& old_var_qp_values = ref_data_ptr_->old_var_qp_values_;
  old_var_qp_value_.assign(max_t, 0.0);

  double local_r = 0.0;
  for (uint32_t qp : ref_data_ptr_->qp_indices_)
  {
    test_i_qp_ = ref_data_ptr_->shape_values_[i][qp];
    test_grad_i_qp_ = ref_data_ptr_->shape_grad_values_[i][qp];
    var_qp_value_ = ref_data_ptr_->var_qp_values_[qp];
    var_grad_qp_value_ = ref_data_ptr_->var_grad_qp_values_[qp];
    qp_xyz_ = ref_data_ptr_->qpoints_xyz_[qp];

    for (size_t t=0; t<max_t; ++t)
      old_var_qp_value_[t] = old_var_qp_values[t][qp];

    local_r += ref_data_ptr_->JxW_values_[qp] * ResidualEntryAtQP();
  }

  return local_r;
}

double FEMKernel::ComputeLocalJacobian(uint32_t i, uint32_t j)
{
  dt_ = ref_data_ptr_->time_data_.dt_;
  time_ = ref_data_ptr_->time_data_.time_;

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

double FEMKernel::ResidualEntryAtQP() { return 0.0; }
double FEMKernel::JacobianEntryAtQP() { return 0.0; }

} // namespace chi_math