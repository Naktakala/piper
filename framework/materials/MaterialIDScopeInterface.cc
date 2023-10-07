#include "MaterialIDScopeInterface.h"

namespace chi
{

InputParameters MaterialIDScopeInterface::GetInputParameters()
{
  InputParameters params;

  params.AddOptionalParameterArray(
    "mat_ids",
    std::vector<chi::ParameterBlock>{},
    "A list of material-ids to which this block will be"
    "restricted");

  return params;
}

MaterialIDScopeInterface::MaterialIDScopeInterface(
  const InputParameters& params)
  : mat_ids_((params.ParametersAtAssignment().Has("mat_ids"))
               ? params.GetParamVectorValue<int>("mat_ids")
               : std::vector<int>{})
{
}

const std::vector<int>& MaterialIDScopeInterface::GetMaterialIDScope() const
{
  return mat_ids_;
}

} // namespace chi