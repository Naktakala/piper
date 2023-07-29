#include "BoundaryComponent.h"

#include "ChiObjectFactory.h"

namespace piper
{

RegisterChiObject(piper, BoundaryComponent);

chi::InputParameters BoundaryComponent::GetInputParameters()
{
  chi::InputParameters params = Component::GetInputParameters();

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
  : Component(
      params, ComponentType::BoundaryLike, {utils::Connection{"to_or_from"}})
{
}

/**We do the same here that we do for single junctions.*/
void BoundaryComponent::Nodalize(std::string connection_point_name,
                                 const chi_mesh::Vector3& datum)
{
  for (auto& connection : connection_points_)
    connection.position_ = datum;
}

} // namespace piper