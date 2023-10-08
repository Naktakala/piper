#include "FEMKernel.h"

#include "math/Systems/EquationSystemTimeData.h"
#include "math/KernelSystem/KernelSystem.h"
#include "math/KernelSystem/Coupling/FEMMaterialProperty.h"
#include "materials/MaterialPropertiesData.h"

#include "chi_log.h"

namespace chi_math
{

chi::InputParameters FEMKernel::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();
  params += chi::MaterialIDScopeInterface::GetInputParameters();
  params += chi_math::CoupledFieldInterface::GetInputParameters();
  params += chi_math::FEMMaterialPropertyInterface::GetInputParameters();

  params.AddRequiredParameter<std::string>("type", "The kernel type.");

  params.AddRequiredParameter<std::string>(
    "var", "Name of the variable this kernel acts on.");
  params.AddOptionalParameter(
    "var_component", 0, "On which component this kernel is active.");

  params.AddRequiredParameter<size_t>("fem_data_handle",
                                      "Handle to a FEMKernelSystemData block.");

  params.AddOptionalParameter(
    "properties_at_qps",
    true,
    "Flag, if false, will evaluate material properties at the centroid instead "
    "of at quadrature points.");

  return params;
}

FEMKernel::FEMKernel(const chi::InputParameters& params)
  : ChiObject(params),
    chi::MaterialIDScopeInterface(params),
    chi_math::CoupledFieldInterface(params),
    chi_math::FEMMaterialPropertyInterface(params, GetMaterialIDScope()),
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
    coord_(fem_data_.coord_qp_values_),
    qp_xyz_(fem_data_.qp_data_.QPointsXYZ())
{
}

void FEMKernel::PreComputeValues()
{
  PreComputeInternalCoupledFields();
  PreComputeInternalMaterialProperties();
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
    local_r += coord_[qp_] *JxW_values_[qp_] * ResidualEntryAtQP();

  return local_r;
}

double FEMKernel::ComputeLocalJacobian(uint32_t i, uint32_t j)
{
  i_ = i;
  j_ = j;
  const size_t num_qpoints = var_value_.size();

  double local_j = 0.0;

  for (qp_ = 0; qp_ < num_qpoints; ++qp_)
    local_j += coord_[qp_] * JxW_values_[qp_] * JacobianEntryAtQP();

  return local_j;
}

} // namespace chi_math