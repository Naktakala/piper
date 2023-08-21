#ifndef PIPER_SINGLEJUNCTION_H
#define PIPER_SINGLEJUNCTION_H

#include "HardwareComponent.h"

namespace piper
{
class HardwareComponent;

class SingleJunction : public HardwareComponent
{
public:
  static chi::InputParameters GetInputParameters();
  explicit SingleJunction(const chi::InputParameters& params);

  utils::FlowOrientation
  FlowOrientationRelToConPoint(size_t con_point_id) const override;

  void Nodalize(size_t connection_point_id,
                const chi_mesh::Vector3& datum) override;

  void EstablishHardwareReferenceData(
    const std::vector<std::shared_ptr<HardwareComponent>>& hardware_components)
    override;

  double Area() const override;

protected:
  double A_; /// Flow area
  chi::ParameterBlock area_param_;
};

} // namespace piper

#endif // PIPER_SINGLEJUNCTION_H
