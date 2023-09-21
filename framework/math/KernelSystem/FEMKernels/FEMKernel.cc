#include "FEMKernel.h"

#include "math/Systems/EquationSystemTimeData.h"
#include "math/KernelSystem/FEMKernelSystem.h"

#include "chi_log.h"

namespace chi_math
{

chi::InputParameters FEMKernel::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<std::string>("type", "The kernel type.");

  params.AddOptionalParameterArray(
    "mat_ids",
    std::vector<chi::ParameterBlock>{},
    "A list of material-ids to which this block will be"
    "restricted");

  params.AddRequiredParameter<std::string>(
    "var", "Name of the variable this kernel acts on.");
  params.AddOptionalParameter(
    "var_component", 0, "On which component this kernel is active.");

  params.AddRequiredParameter<size_t>("fem_data_handle",
                                      "Handle to a FEMKernelSystemData block.");

  return params;
}

FEMKernel::FEMKernel(const chi::InputParameters& params)
  : ChiObject(params),
    var_name_component_(params.GetParamValue<std::string>("var"),
      params.GetParamValue<uint32_t>("var_component")),
    fem_data_(Chi::GetStackItem<FEMKernelSystemData>(
      Chi::object_stack,
      params.GetParamValue<size_t>("fem_data_handle"),
      __FUNCTION__)),
    dt_(fem_data_.time_data_.dt_),
    time_(fem_data_.time_data_.time_),
    var_dot_dvar_(fem_data_.time_data_.var_dot_dvar_),
    JxW_values_(fem_data_.qp_data_.JxW_Values()),
    test_values_(fem_data_.qp_data_.ShapeValues()),
    test_grad_values_(fem_data_.qp_data_.ShapeGradValues()),
    shape_values_(fem_data_.qp_data_.ShapeValues()),
    shape_grad_values_(fem_data_.qp_data_.ShapeGradValues()),
    var_value_(fem_data_.var_qp_values_),
    var_grad_value_(fem_data_.var_grad_qp_values_),
    var_dot_value_(fem_data_.var_dot_qp_values_),
    qp_xyz_(fem_data_.qp_data_.QPointsXYZ())
{
  const auto& user_params = params.ParametersAtAssignment();

  if (user_params.Has("mat_ids"))
    mat_ids_ = params.GetParamVectorValue<int>("mat_ids");
}

const std::vector<int>& FEMKernel::GetMaterialIDScope() const
{
  return mat_ids_;
}

bool FEMKernel::IsTimeKernel() const { return false; }

const std::pair<std::string, uint32_t>&
FEMKernel::ActiveVariableAndComponent() const
{
  return var_name_component_;
}

double FEMKernel::ComputeLocalResidual(uint32_t i)
{
  i_ = i;
  const size_t num_qpoints = var_value_.size();

  double local_r = 0.0;

  for (qp_ = 0; qp_ < num_qpoints; ++qp_)
    local_r += JxW_values_[qp_] * ResidualEntryAtQP();

  return local_r;
}

double FEMKernel::ComputeLocalJacobian(uint32_t i, uint32_t j)
{
  i_ = i;
  j_ = j;
  const size_t num_qpoints = var_value_.size();

  double local_j = 0.0;

  for (qp_ = 0; qp_ < num_qpoints; ++qp_)
    local_j += JxW_values_[qp_] * JacobianEntryAtQP();

  return local_j;
}

double FEMKernel::ResidualEntryAtQP() { return 0.0; }
double FEMKernel::JacobianEntryAtQP() { return 0.0; }

} // namespace chi_math