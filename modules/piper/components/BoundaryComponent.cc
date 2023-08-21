#include "BoundaryComponent.h"

#include "ChiObjectFactory.h"

namespace piper
{

RegisterChiObject(piper, BoundaryComponent);

chi::InputParameters BoundaryComponent::GetInputParameters()
{
  chi::InputParameters params = HardwareComponent::GetInputParameters();

  params.SetGeneralDescription(
    "A boundary component is a proxy component for eventual boundary "
    "conditions. "
    "It can be used as a placeholder for boundaries where the boundary "
    "condition "
    "will be different depending on the physics modelled in the system.");
  params.SetDocGroup("Piper");

  return params;
}

BoundaryComponent::BoundaryComponent(const chi::InputParameters& params)
  : HardwareComponent(params,
                      ComponentCategory::BoundaryLike,
                      {utils::Connection{"to_or_from"}})
{
}

/**Returns the component's flow orientation relative to the connection
 * point.*/
utils::FlowOrientation
BoundaryComponent::FlowOrientationRelToConPoint(size_t con_point_id) const
{
  if (con_point_id == 0) return utils::FlowOrientation::OUTGOING;
  else if (con_point_id == 1)
    return utils::FlowOrientation::INCOMING;
  else
    ChiInvalidArgument("Invalid con_point_id " + std::to_string(con_point_id));
}

/**We do the same here that we do for single junctions.*/
void BoundaryComponent::Nodalize(size_t connection_point_id,
                                 const chi_mesh::Vector3& datum)
{
  for (auto& connection : connection_points_)
    connection.position_ = datum;
}

} // namespace piper