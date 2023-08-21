#ifndef PIPER_BOUNDARYCOMPONENT_H
#define PIPER_BOUNDARYCOMPONENT_H

#include "HardwareComponent.h"

namespace piper
{

/**A boundary component is a proxy component for eventual boundary conditions.
* It can be used as a placeholder for boundaries where the boundary condition
* will be different depending on the physics modelled in the system.*/
class BoundaryComponent : public HardwareComponent
{
public:
  static chi::InputParameters GetInputParameters();
  explicit BoundaryComponent(const chi::InputParameters& params);

  utils::FlowOrientation
  FlowOrientationRelToConPoint(size_t con_point_id) const override;

  void Nodalize(size_t connection_point_id,
                const chi_mesh::Vector3& datum) override;
};

}

#endif // PIPER_BOUNDARYCOMPONENT_H
