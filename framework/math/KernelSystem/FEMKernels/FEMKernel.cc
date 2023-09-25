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

const chi_math::FEMMaterialProperty&
FEMKernel::GetFEMMaterialProperty(const std::string& name)
{
  const auto& kernel_mat_scope = this->GetMaterialIDScope();

  for (const auto& property_ptr : fem_data_.fem_material_properties_)
  {
    if (property_ptr->TextName() == name)
    {
      const auto& property_mat_scope = property_ptr->GetMaterialIDScope();
      if (property_mat_scope.empty()) return *property_ptr;
      else if (kernel_mat_scope.empty())
      {
        ChiLogicalError("Kernel defined on all materials but material property "
                        "\"" +
                        name + "\" is not.");
      }
      else
      {
        for (int mi : kernel_mat_scope)
        {
          bool found = false;
          for (int pmi : property_mat_scope)
          {
            if (mi == pmi)
            {
              found = true;
              break;
            }
          }//for pmi

          ChiLogicalErrorIf(not found,
                            "Kernel is defined on material id " +
                              std::to_string(mi) + " but material property \"" +
                              name + "\" is not.");
        }// for mi
      }
      return *property_ptr;
    } // if property name found
  } // for property in map

  ChiLogicalError("Kernel required parameter \"" + name + "\" not found.");
}

} // namespace chi_math