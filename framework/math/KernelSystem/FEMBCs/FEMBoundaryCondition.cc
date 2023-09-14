#include "FEMBoundaryCondition.h"

#include "math/SpatialDiscretization/FiniteElement/finite_element.h"

#include "chi_log.h"

namespace chi_math
{

FEMBCRefData::FEMBCRefData(
  const std::shared_ptr<const finite_element::FaceQuadraturePointData>& qp_data,
  VecDbl var_qp_values,
  VecVec3 var_grad_qp_values,
  const VecDbl& var_nodal_values,
  const VecVec3& node_locations)
  : qp_data_(qp_data),
    qp_indices_(qp_data->QuadraturePointIndices()),
    qpoints_xyz_(qp_data->QPointsXYZ()),
    shape_values_(qp_data->ShapeValues()),
    shape_grad_values_(qp_data->ShapeGradValues()),
    JxW_values_(qp_data->JxW_Values()),
    normals_(qp_data->Normals()),
    var_qp_values_(std::move(var_qp_values)),
    var_grad_qp_values_(std::move(var_grad_qp_values)),
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

bool FEMBoundaryCondition::IsDirichlet() const { return false; }

void FEMBoundaryCondition::SetFaceReferenceData(
  FaceID face_index, std::shared_ptr<FEMBCRefData>& ref_data_ptr)
{
  face_id_2_ref_data_ptr_map_[face_index] = ref_data_ptr;
}

const std::vector<std::string>& FEMBoundaryCondition::GetBoundaryScope() const
{
  return boundary_scope_;
}

double FEMBoundaryCondition::ComputeLocalResidual(size_t f, uint32_t i)
{
  auto& ref_data = *face_id_2_ref_data_ptr_map_.at(f);
  i_ = i;
  double local_r = 0.0;
  for (uint32_t qp : ref_data.qp_indices_)
  {
    test_i_qp_ = ref_data.shape_values_[i][qp];
    test_grad_i_qp_ = ref_data.shape_grad_values_[i][qp];
    var_value_qp_ = ref_data.var_qp_values_[qp];
    var_grad_value_qp_ = ref_data.var_grad_qp_values_[qp];
    qp_xyz_ = ref_data.qpoints_xyz_[qp];
    normal_qp_ = ref_data.normals_[qp];

    local_r += ref_data.JxW_values_[qp] * ResidualEntryAtQP();
  }

  return local_r;
}

double FEMBoundaryCondition::ComputeLocalJacobian(size_t f, uint32_t i, uint32_t j)
{
  auto& ref_data = *face_id_2_ref_data_ptr_map_.at(f);
  i_ = i;
  double local_j = 0.0;
  for (uint32_t qp : ref_data.qp_indices_)
  {
    test_i_qp_ = ref_data.shape_values_[i][qp];
    test_grad_i_qp_ = ref_data.shape_grad_values_[i][qp];
    var_value_qp_ = ref_data.var_qp_values_[qp];
    var_grad_value_qp_ = ref_data.var_grad_qp_values_[qp];
    qp_xyz_ = ref_data.qpoints_xyz_[qp];
    normal_qp_ = ref_data.normals_[qp];

    shape_j_qp_ = ref_data.shape_values_[j][qp];
    shape_grad_j_qp_ = ref_data.shape_grad_values_[j][qp];

    local_j += ref_data.JxW_values_[qp] * JacobianEntryAtQP();
  }

  return local_j;
}

double FEMBoundaryCondition::ResidualEntryAtQP() { return 0.0; }
double FEMBoundaryCondition::JacobianEntryAtQP() { return 0.0; }

} // namespace chi_math