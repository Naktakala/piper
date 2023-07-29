#ifndef PIPER_ORIENTATION_H
#define PIPER_ORIENTATION_H

#include "ChiObject.h"

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

protected:
  double varphi_ = 0.0;
  double theta_ = 0.0;
};

}

#endif // PIPER_ORIENTATION_H
