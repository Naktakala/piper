#ifndef PIPER_HARDWARECOMPONENT_H
#define PIPER_HARDWARECOMPONENT_H

#include "ChiObject.h"
#include "piper/utils/utils.h"
#include "piper/utils/Orientation.h"
#include "piper/utils/Connections.h"

namespace piper
{

std::string ComponentCategoryName(ComponentCategory category);

/**Base class for all components*/
class HardwareComponent : public ChiObject
{
public:
  const std::string& Name() const;

  const std::vector<utils::Connection>& ConnectionPoints() const;
  std::vector<utils::Connection>& ConnectionPoints();

  const Orientation& GetOrientation() const;

  size_t GetID() const;
  void SetID(size_t id);

  ComponentCategory Category() const;

  const std::string& ConnectionPointName(size_t con_point_id) const;
  size_t ConnectionPointID(const std::string& con_point_name) const;

  /**Returns the component's flow orientation relative to the connection
   * point.*/
  virtual utils::FlowOrientation
  FlowOrientationRelToConPoint(size_t con_point_id) const = 0;

  virtual void Nodalize(size_t connection_point_id,
                        const chi_mesh::Vector3& datum) = 0;

  const chi_mesh::Vector3& GetRootNodePosition() const;
  virtual chi_mesh::Vector3 MakeCentroid() const;

  virtual double Area() const { return 0.0; }
  virtual double Volume() const { return 0.0; }

  /**This function gets called in Piper right after the hardware components
   * are connected to each other. It allows for setting things like "the min of
   * adjoining areas".*/
  virtual void EstablishHardwareReferenceData(
    const std::vector<std::shared_ptr<HardwareComponent>>&
      hardware_components){};

protected:
  // input items
  Orientation orientation_;

  // constructed or initialization items

  /**These will be junctions. */
  std::vector<utils::Connection> connection_points_;

  static chi::InputParameters GetInputParameters();
  explicit HardwareComponent(const chi::InputParameters& params,
                             ComponentCategory type,
                             std::vector<utils::Connection> connection_points);

private:
  const ComponentCategory category_;
  const std::string name_;
  size_t id_ = 0;
};

} // namespace piper

#endif // PIPER_HARDWARECOMPONENT_H
