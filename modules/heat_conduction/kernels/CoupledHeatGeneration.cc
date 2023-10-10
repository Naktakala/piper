#include "CoupledHeatGeneration.h"

#include "ChiObjectFactory.h"

namespace hcm
{

RegisterChiObject(hcm, CoupledHeatGeneration);

chi::InputParameters CoupledHeatGeneration::GetInputParameters()
{
  chi::InputParameters params = chi_math::FEMKernel::GetInputParameters();

  params.SetGeneralDescription(
    "Heat generation kernel with strength from another field.");
  params.SetDocGroup("doc_HeatConduction");

  params.AddOptionalParameter(
    "e_gen", 1.0, "A value or the name of the variable to couple to.");
  params.SetParameterTypeMismatchAllowed("e_gen");

  return params;
}

CoupledHeatGeneration::CoupledHeatGeneration(const chi::InputParameters& params)
: chi_math::FEMKernel(params),
    e_gen_(GetCoupledField(params.GetParamValue<std::string>("e_gen")))
{}

double CoupledHeatGeneration::ResidualEntryAtQP()
{
  return -e_gen_[qp_] * test_values_[i_][qp_];
}

} // namespace hcm