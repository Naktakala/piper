#ifndef PIPER_SINGLEVOLUME_H
#define PIPER_SINGLEVOLUME_H

#include "Component.h"

namespace piper
{

class SingleVolume : public Component
{
public:
  static chi::InputParameters GetInputParameters();
  explicit SingleVolume(const chi::InputParameters& params);

  void Nodalize(std::string connection_point_name,
                const chi_mesh::Vector3& datum) override;

protected:
  double Dh_; /// Hydraulic diameter
  double A_; /// Flow area
  double length_;
  double roughness_;
};

}

#endif // PIPER_SINGLEVOLUME_H
