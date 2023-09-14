#include "ThermalConductionKernel.h"

#include "ChiObjectFactory.h"

namespace hcm
{

RegisterChiObject(hcm, ThermalConductionKernel);

chi::InputParameters ThermalConductionKernel::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.SetGeneralDescription("A material property where the thermal "
                               "conductivity is constant.");
  params.SetDocGroup("doc_HeatConduction");

  params.AddRequiredParameter<double>(
    "k", "Value for the constant thermal conductivity.");

  return params;
}

ThermalConductionKernel::ThermalConductionKernel(
  const chi::InputParameters& params)
  : chi_math::FEMKernel(params),
    k_(params.GetParamValue<double>("k"))
{
}

double ThermalConductionKernel::ResidualEntryAtQP()
{
  return k_ * test_grad_i_qp_.Dot(var_grad_qp_value_);
}

double ThermalConductionKernel::JacobianEntryAtQP()
{
  return k_ * test_grad_i_qp_.Dot(shape_grad_j_qp_);
}

} // namespace hcm