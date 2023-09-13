#include "FEMBoundaryCondition.h"

#include "math/SpatialDiscretization/FiniteElement/finite_element.h"

#include "chi_log.h"

namespace chi_math
{

FEMBCRefData::FEMBCRefData(
  const finite_element::FaceQuadraturePointData& qp_data,
  const VecDbl& var_values,
  const VecVec3& var_grad_values,
  const VecDbl& var_nodal_values,
  const VecVec3& node_locations)
  : qp_indices_(qp_data.QuadraturePointIndices()),
    qpoints_xyz_(qp_data.QPointsXYZ()),
    shape_values_(qp_data.ShapeValues()),
    shape_grad_values_(qp_data.ShapeGradValues()),
    JxW_values_(qp_data.JxW_Values()),
    var_values_(var_values),
    var_grad_values_(var_grad_values),
    var_nodal_values_(var_nodal_values),
    node_locations_(node_locations)
{
}

chi::InputParameters FEMBoundaryCondition::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddOptionalParameterArray("boundaries",
                                   std::vector<std::string>{},
                                   "A list of boundary names on which this "
                                   "boundary condition mush be supplied.");

  return params;
}

FEMBoundaryCondition::FEMBoundaryCondition(const chi::InputParameters& params)
  : ChiObject(params),
    boundary_scope_(params.GetParamVectorValue<std::string>("boundaries"))
{
}

void FEMBoundaryCondition::SetReferenceData(
  std::shared_ptr<FEMBCRefData>& ref_data_ptr)
{
  ref_data_ptr_ = ref_data_ptr;
}

const std::vector<std::string>& FEMBoundaryCondition::GetBoundaryScope() const
{
  return boundary_scope_;
}

double FEMBoundaryCondition::ComputeLocalResidual(uint32_t i)
{
  i_ = i;
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

double FEMBoundaryCondition::ComputeLocalJacobian(uint32_t i, uint32_t j)
{
  i_ = i;
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

double FEMBoundaryCondition::ResidualEntryAtQP() { return 0.0; }
double FEMBoundaryCondition::JacobianEntryAtQP() { return 0.0; }

} // namespace chi_math