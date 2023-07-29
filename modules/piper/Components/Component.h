#ifndef PIPER_COMPONENT_H
#define PIPER_COMPONENT_H

#include "ChiObject.h"
#include "piper/utils/Orientation.h"
#include "piper/utils/Connections.h"

namespace piper
{
enum class ComponentType : int
{
  BoundaryLike = 0,
  Volumetric = 1,
  JunctionLike = 2
};

std::string ComponentTypeName(ComponentType type);

/**Base class for all components*/
class Component : public ChiObject
{
public:
  const std::string& Name() const;
  const std::vector<utils::Connection>& ConnectionPoints() const;
  std::vector<utils::Connection>& ConnectionPoints();
  const Orientation& GetOrientation() const;
  size_t GetID() const;
  void SetID(size_t id);
  ComponentType Type() const;

  virtual void Nodalize(std::string connection_point_name,
                        const chi_mesh::Vector3& datum) = 0;

protected:
  // input items
  Orientation orientation_;

  // constructed or initialization items

  /**These will be junctions. */
  std::vector<utils::Connection> connection_points_;

  static chi::InputParameters GetInputParameters();
  explicit Component(const chi::InputParameters& params,
                     ComponentType type,
                     std::vector<utils::Connection> connection_points);

private:
  const ComponentType type_;
  const std::string name_;
  size_t id_ = 0;
};

} // namespace piper

#endif // PIPER_COMPONENT_H
