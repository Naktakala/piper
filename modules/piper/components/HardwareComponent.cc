#include "HardwareComponent.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace piper
{

std::string ComponentCategoryName(ComponentCategory category)
{
  // clang-format off
  switch (category)
  {
    case (ComponentCategory::BoundaryLike): return "BoundaryLike";
    case (ComponentCategory::Volumetric):   return "Volumetric";
    case (ComponentCategory::JunctionLike): return "JunctionLike";
    default: throw std::logic_error("in ComponentCategoryName");
  }
  // clang-format on
}

// ##################################################################
chi::InputParameters HardwareComponent::GetInputParameters()
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

// ##################################################################
HardwareComponent::HardwareComponent(
  const chi::InputParameters& params,
  ComponentCategory type,
  std::vector<utils::Connection> connection_points)
  : ChiObject(params),
    orientation_(
      params.ParametersAtAssignment().Has("orientation")
        ? Orientation::MakeOrientation(params.GetParam("orientation"))
        : Orientation::MakeOrientation(chi::ParameterBlock())),
    connection_points_(std::move(connection_points)),
    category_(type),
    name_(params.GetParamValue<std::string>("name"))
{
  Chi::log.Log0Verbose1() << "Component \"" << Name() << "\" created. "
                          << "orientation={azimuthal_angle="
                          << orientation_.Varphi()
                          << ", polar_anlge=" << orientation_.Theta() << "}";
}

// ##################################################################
const std::string& HardwareComponent::Name() const { return name_; }

// ##################################################################
const std::vector<utils::Connection>&
HardwareComponent::ConnectionPoints() const
{
  return connection_points_;
}
std::vector<utils::Connection>& HardwareComponent::ConnectionPoints()
{
  return connection_points_;
}

// ##################################################################
const Orientation& HardwareComponent::GetOrientation() const
{
  return orientation_;
}

size_t HardwareComponent::GetID() const { return id_; }
void HardwareComponent::SetID(size_t id) { id_ = id; }
ComponentCategory HardwareComponent::Category() const { return category_; }

// ##################################################################
const std::string&
HardwareComponent::ConnectionPointName(size_t con_point_id) const
{
  ChiLogicalErrorIf(con_point_id >= connection_points_.size(),
                    "Invalid con_point_id=" + std::to_string(con_point_id) +
                      " for list of size " +
                      std::to_string(connection_points_.size()) + ".");
  return connection_points_[con_point_id].name_;
}

// ##################################################################
/**Returns the connection point id given the name.*/
size_t
HardwareComponent::ConnectionPointID(const std::string& con_point_name) const
{
  auto functor = [&con_point_name](const utils::Connection& conn)
  { return conn.name_ == con_point_name; };
  auto it =
    std::find_if(connection_points_.begin(), connection_points_.end(), functor);
  ChiInvalidArgumentIf(it == connection_points_.end(),
                       Name() + ": No connection point with name \"" +
                         con_point_name + "\"");

  return std::distance(connection_points_.begin(), it);
}

// ##################################################################
/**Returns true if the component's flow orientation is outgoing relative
 * to the connection point.*/
bool HardwareComponent::IsOutgoingRelToConPoint(size_t con_point_id) const
{
  return FlowOrientationRelToConPoint(con_point_id) ==
         utils::FlowOrientation::OUTGOING;
}

// ##################################################################
const chi_mesh::Vector3& HardwareComponent::GetRootNodePosition() const
{
  ChiLogicalErrorIf(connection_points_.empty(), "Empty nodes");
  return connection_points_.front().position_;
}

// ##################################################################
chi_mesh::Vector3 HardwareComponent::MakeCentroid() const
{
  ChiLogicalErrorIf(connection_points_.empty(), "Empty nodes");

  chi_mesh::Vector3 centroid;
  for (const auto& node : connection_points_)
    centroid += node.position_;
  centroid /= double(connection_points_.size());

  return centroid;
}

} // namespace piper