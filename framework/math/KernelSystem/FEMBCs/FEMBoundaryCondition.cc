#include "FEMBoundaryCondition.h"

#include "math/Systems/EquationSystemTimeData.h"
#include "math/KernelSystem/FEMKernelSystem.h"

#include "chi_log.h"

namespace chi_math
{

chi::InputParameters FEMBoundaryCondition::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<std::string>("type", "The bc-kernel type.");

  params.AddOptionalParameterArray("boundaries",
                                   std::vector<std::string>{},
                                   "A list of boundary names on which this "
                                   "boundary condition mush be supplied.");

  params.AddRequiredParameter<std::string>(
    "var", "Name of the variable this BC acts on.");

  params.AddOptionalParameter(
    "var_component", 0, "On which component this BC is active.");

  params.AddRequiredParameter<size_t>("fem_data_handle",
                                      "Handle to a FEMKernelSystemData block.");

  return params;
}

FEMBoundaryCondition::FEMBoundaryCondition(const chi::InputParameters& params)
  : ChiObject(params),
    var_name_component_(params.GetParamValue<std::string>("var"),
                        params.GetParamValue<uint32_t>("var_component")),
    boundary_scope_(params.GetParamVectorValue<std::string>("boundaries")),
    fem_data_(Chi::GetStackItem<FEMKernelSystemData>(
      Chi::object_stack,
      params.GetParamValue<size_t>("fem_data_handle"),
      __FUNCTION__)),
    dt_(fem_data_.time_data_.dt_),
    time_(fem_data_.time_data_.time_),
    JxW_values_(fem_data_.face_qp_data_.JxW_Values()),
    test_values_(fem_data_.face_qp_data_.ShapeValues()),
    test_grad_values_(fem_data_.face_qp_data_.ShapeGradValues()),
    shape_values_(fem_data_.face_qp_data_.ShapeValues()),
    shape_grad_values_(fem_data_.face_qp_data_.ShapeGradValues()),
    var_value_(fem_data_.face_var_qp_values_),
    var_grad_value_(fem_data_.face_var_grad_qp_values_),
    qp_xyz_(fem_data_.face_qp_data_.QPointsXYZ()),
    normal_(fem_data_.face_qp_data_.Normals()),
    nodal_var_values_(fem_data_.nodal_var_values_),
    node_locations_(fem_data_.node_locations_)
{
}

bool FEMBoundaryCondition::IsDirichlet() const { return false; }

const std::pair<std::string, uint32_t>&
FEMBoundaryCondition::ActiveVariableAndComponent() const
{
  return var_name_component_;
}

const std::vector<std::string>& FEMBoundaryCondition::GetBoundaryScope() const
{
  return boundary_scope_;
}

double FEMBoundaryCondition::ComputeLocalResidual(size_t f, uint32_t i)
{
  i_ = i;
  const size_t num_qpoints = var_value_.size();

  double local_r = 0.0;

  for (qp_ = 0; qp_ < num_qpoints; ++qp_)
    local_r += JxW_values_[qp_] * ResidualEntryAtQP();

  return local_r;
}

double
FEMBoundaryCondition::ComputeLocalJacobian(size_t f, uint32_t i, uint32_t j)
{
  i_ = i;
  j_ = j;
  const size_t num_qpoints = var_value_.size();

  double local_j = 0.0;

  for (qp_ = 0; qp_ < num_qpoints; ++qp_)
    local_j += JxW_values_[qp_] * JacobianEntryAtQP();

  return local_j;
}

double FEMBoundaryCondition::ResidualEntryAtQP() { return 0.0; }
double FEMBoundaryCondition::JacobianEntryAtQP() { return 0.0; }

} // namespace chi_math