#ifndef PIPER_SINGLEJUNCTION_H
#define PIPER_SINGLEJUNCTION_H

#include "Component.h"

namespace piper
{
class Component;

class SingleJunction : public Component
{
public:
  static chi::InputParameters GetInputParameters();
  explicit SingleJunction(const chi::InputParameters& params);

  void Nodalize(std::string connection_point_name,
                const chi_mesh::Vector3& datum) override;

protected:
};

}

#endif // PIPER_SINGLEJUNCTION_H
