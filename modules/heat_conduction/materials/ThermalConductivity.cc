#include "ThermalConductivity.h"

#include "math/Functions/function_dimA_to_dimB.h"

#include "ChiObjectFactory.h"
#include "chi_log.h"

namespace hcm
{

RegisterChiObject(hcm, ThermalConductivity);

chi::InputParameters ThermalConductivity::GetInputParameters()
{
  chi::InputParameters params = chi::MaterialProperty::GetInputParameters();

  params.SetGeneralDescription("Material property for thermal conductivity.");

  params.AddOptionalParameter(
    "constant_value", 1.0, "Constant value for the thermal conductivity.");

  params.AddOptionalParameter("value_function", 0, "Handle to a function.");

  params.AddOptionalParameter(
    "temperature_fieldname", "T", "Default fieldname for temperature");

  return params;
}

ThermalConductivity::ThermalConductivity(const chi::InputParameters& params)
  : chi::MaterialProperty(params),
    constant_value_(params.GetParamValue<double>("constant_value")),
    function_(GetFunction(params))
{
  const auto& user_params = params.ParametersAtAssignment();
  if (user_params.Has("constant_value") and user_params.Has("function"))
    Chi::log.Log0Warning()
      << "ThermalConductivity: parameters \"constant_value\""
         "and \"function\" have both been specified. Defaulting to the "
         "function.";
}

std::shared_ptr<chi_math::FunctionDimAToDimB>
ThermalConductivity::GetFunction(const chi::InputParameters& params)
{
  if (params.ParametersAtAssignment().Has("value_function"))
  {
    const size_t handle = params.GetParamValue<size_t>("value_function");
    auto function_ptr =
      Chi::GetStackItemPtrAsType<chi_math::FunctionDimAToDimB>(
        Chi::object_stack, handle, __FUNCTION__);

    const size_t in_dim = function_ptr->InputDimension();
    const size_t out_dim = function_ptr->OutputDimension();

    size_t required_dim = 1;
    if (IsPositionDependent()) required_dim += 3;
    if (IsTimeDependent()) required_dim += 1;

    ChiInvalidArgumentIf(
      in_dim != required_dim,
      "The supplied function has input dimension " + std::to_string(in_dim) +
        " but this material property requires and input dimension of " +
        std::to_string(required_dim));

    ChiInvalidArgumentIf(
      out_dim != 1,
      "The supplied function has output dimension " + std::to_string(out_dim) +
        " but this material property requires and output dimension of 1.");

    return function_ptr;
  }
  else { return nullptr; }
}

std::vector<std::string> ThermalConductivity::RequiredInputNames() const
{
  if (not function_) return {};
  else
    return {temperature_fieldname_};
}

double ThermalConductivity::ComputeScalarValue(
  const std::vector<double>& input_params) const
{
  if (not function_) return constant_value_;
  else
  {
    auto output = function_->Evaluate(input_params);
    return output.front();
  }
}

bool ThermalConductivity::HasDerivative() const
{
  return function_ != nullptr;
}

} // namespace hcm