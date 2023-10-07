#include "MaterialProperty.h"

namespace chi
{

InputParameters MaterialProperty2::GetInputParameters()
{
  InputParameters params = ChiObject::GetInputParameters();
  params += MaterialIDScopeInterface::GetInputParameters();

  params.AddRequiredParameter<std::string>(
    "name", "Text name associated with this property");

  params.AddOptionalParameter(
    "position_dependent",
    false,
    "Flag controlling whether this property is position dependent.");

  params.AddOptionalParameter(
    "time_dependent",
    false,
    "Flag controlling whether this property is time dependent.");

  return params;
}

MaterialProperty2::MaterialProperty2(const chi::InputParameters& params)
  : ChiObject(params),
    MaterialIDScopeInterface(params),
    name_(params.GetParamValue<std::string>("name")),
    position_dependent_(params.GetParamValue<bool>("position_dependent")),
    time_dependent_(params.GetParamValue<bool>("time_dependent"))
{
}

const std::string& MaterialProperty2::TextName() const { return name_; }

bool MaterialProperty2::IsPositionDependent() const
{
  return position_dependent_;
}

bool MaterialProperty2::IsTimeDependent() const { return time_dependent_; }
} // namespace chi
