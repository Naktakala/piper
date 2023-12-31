#include "ThermalConductionKernel2.h"

#include "ChiObjectFactory.h"

namespace hcm
{

RegisterChiObject(hcm, ThermalConductionKernel2);

chi::InputParameters ThermalConductionKernel2::GetInputParameters()
{
  chi::InputParameters params = FEMKernel::GetInputParameters();

  params.SetGeneralDescription("Kernel for classical heat conduction"
                               "driven by a conduction coefficient.");
  params.SetDocGroup("doc_HeatConduction");

  params.AddOptionalParameter(
    "k_property_name", "k", "Property name of the thermal conductivity");

  return params;
}

ThermalConductionKernel2::ThermalConductionKernel2(
  const chi::InputParameters& params)
  : chi_math::FEMKernel(params),
    k_property_name_(params.GetParamValue<std::string>("k_property_name")),
    k_(GetFEMMaterialProperty(k_property_name_))
{
}

double ThermalConductionKernel2::ResidualEntryAtQP()
{
  return k_[qp_] * test_grad_values_[i_][qp_].Dot(var_grad_value_[qp_]);
}

double ThermalConductionKernel2::JacobianEntryAtQP()
{
  double entry =
    k_[qp_] * test_grad_values_[i_][qp_].Dot(shape_grad_values_[j_][qp_]);
  if (k_.HasDerivative())
    entry += k_.dvar_[qp_] * shape_values_[j_][qp_] *
             test_grad_values_[i_][qp_].Dot(var_grad_value_[qp_]);
  return entry;
}

} // namespace hcm