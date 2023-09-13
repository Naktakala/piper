#include "FEMKernel.h"

#include "math/SpatialDiscretization/FiniteElement/finite_element.h"

#include "chi_log.h"

namespace chi_math
{

FEMKernelRefData::FEMKernelRefData(
  const finite_element::InternalQuadraturePointData& qp_data,
  const VecDbl& var_values,
  const VecVec3& var_grad_values)
  : qp_indices_(qp_data.QuadraturePointIndices()),
    qpoints_xyz_(qp_data.QPointsXYZ()),
    shape_values_(qp_data.ShapeValues()),
    shape_grad_values_(qp_data.ShapeGradValues()),
    JxW_values_(qp_data.JxW_Values()),
    var_values_(var_values),
    var_grad_values_(var_grad_values)
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

double FEMKernel::ComputeLocalResidual(uint32_t i)
{
  double local_r = 0.0;
  for (uint32_t qp : ref_data_ptr_->qp_indices_)
  {
    test_i_qp_ = ref_data_ptr_->shape_values_[i][qp];
    test_grad_i_qp_ = ref_data_ptr_->shape_grad_values_[i][qp];
    var_value_qp_ = ref_data_ptr_->var_values_[qp];
    var_grad_value_qp_ = ref_data_ptr_->var_grad_values_[qp];
    qp_xyz_ = ref_data_ptr_->qpoints_xyz_[qp];

    local_r += ref_data_ptr_->JxW_values_[qp] * ResidualEntryAtQP();
  }

  return local_r;
}

double FEMKernel::ComputeLocalJacobian(uint32_t i, uint32_t j)
{
  double local_j = 0.0;
  for (uint32_t qp : ref_data_ptr_->qp_indices_)
  {
    test_i_qp_ = ref_data_ptr_->shape_values_[i][qp];
    test_grad_i_qp_ = ref_data_ptr_->shape_grad_values_[i][qp];
    var_value_qp_ = ref_data_ptr_->var_values_[qp];
    var_grad_value_qp_ = ref_data_ptr_->var_grad_values_[qp];
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