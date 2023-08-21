#ifndef PIPER_SINGLEVOLUME_H
#define PIPER_SINGLEVOLUME_H

#include "HardwareComponent.h"

namespace piper
{

class SingleVolume : public HardwareComponent
{
public:
  static chi::InputParameters GetInputParameters();
  explicit SingleVolume(const chi::InputParameters& params);

  utils::FlowOrientation
  FlowOrientationRelToConPoint(size_t con_point_id) const override;

  void Nodalize(size_t connection_point_id,
                const chi_mesh::Vector3& datum) override;

  double Area() const override;
  double Volume() const override;

protected:
  double Dh_; /// Hydraulic diameter
  double A_; /// Flow area
  double length_;
  double roughness_;
};

}

#endif // PIPER_SINGLEVOLUME_H
