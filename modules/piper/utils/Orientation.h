#ifndef PIPER_ORIENTATION_H
#define PIPER_ORIENTATION_H

#include "ChiObject.h"
#include "mesh/chi_mesh.h"

namespace piper
{

class Orientation : public ChiObject
{
public:
  static chi::InputParameters GetInputParameters();
  explicit Orientation(const chi::InputParameters& params);

  static Orientation MakeOrientation(const chi::ParameterBlock& param_block);

  double Varphi() const;
  double Theta() const;
  const chi_mesh::Vector3& Vector() const;

protected:
  double varphi_ = 0.0;
  double theta_ = 0.0;
  chi_mesh::Vector3 vector_;
};

}

#endif // PIPER_ORIENTATION_H
