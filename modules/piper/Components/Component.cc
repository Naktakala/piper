#include "Component.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace piper
{

std::string ComponentTypeName(ComponentType type)
{
  // clang-format off
  switch (type)
  {
    case (ComponentType::BoundaryLike): return "BoundaryLike";
    case (ComponentType::Volumetric):   return "Volumetric";
    case (ComponentType::JunctionLike): return "JunctionLike";
    default: throw std::logic_error("in ComponentTypeName");
  }
  // clang-format on
}

chi::InputParameters Component::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.SetGeneralDescription("Base class for all components");
  params.SetDocGroup("Piper");

  params.AddRequiredParameter<std::string>(
    "name", "Name associated with this component.");

  params.AddOptionalParameterBlock(
    "orientation",
    chi::ParameterBlock(),
    "The orientation of the component. See `piper::Orientation`");

  params.LinkParameterToBlock("orientation", "piper::Orientation");

  return params;
}

Component::Component(const chi::InputParameters& params,
                     ComponentType type,
                     std::vector<utils::Connection> connection_points)
  : ChiObject(params),
    orientation_(
      params.ParametersAtAssignment().Has("orientation")
        ? Orientation::MakeOrientation(params.GetParam("orientation"))
        : Orientation::MakeOrientation(chi::ParameterBlock())),
    connection_points_(std::move(connection_points)),
    type_(type),
    name_(params.GetParamValue<std::string>("name"))
{
  Chi::log.Log0Verbose1() << "Component \"" << Name() << "\" created. "
                          << "orientation={azimuthal_angle="
                          << orientation_.Varphi()
                          << ", polar_anlge=" << orientation_.Theta() << "}";
}

const std::string& Component::Name() const { return name_; }

const std::vector<utils::Connection>& Component::ConnectionPoints() const
{
  return connection_points_;
}
std::vector<utils::Connection>& Component::ConnectionPoints()
{
  return connection_points_;
}

const Orientation& Component::GetOrientation() const { return orientation_; }

size_t Component::GetID() const { return id_; }
void Component::SetID(size_t id) { id_ = id; }
ComponentType Component::Type() const { return type_; }

} // namespace piper