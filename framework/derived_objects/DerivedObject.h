#ifndef CHITECH_DERIVEDOBJECT_H
#define CHITECH_DERIVEDOBJECT_H

#include "ChiObject.h"
#include "event_system/EventSubscriber.h"

namespace chi_mesh
{
struct Vector3;
}

namespace chi
{

class DerivedObject : public ChiObject,
                      public EventSubscriber
{
public:
  /**Response function to events.*/
  void ReceiveEventUpdate(const Event& event) override;

  /**Calls the base ChiObject's method and adds a subscription to
   * `chi_physics::PhysicsEventPublisher` singleton.*/
  void PushOntoStack(std::shared_ptr<ChiObject>& new_object) override;

  /**Update function specific to each derived object.*/
  virtual void Update(const Event& event) = 0;

  /**Returns the name of the derived object.*/
  const std::string& Name() const;

  /**Overrideable method to evaluate the object's value at the given position.*/
  virtual double SpatialValue(const chi_mesh::Vector3& position) const;
protected:
  static InputParameters GetInputParameters();
  explicit DerivedObject(const InputParameters& params);

  const std::string name_;
  std::vector<std::string> subscribed_events_for_execution_;
  std::string solvername_filter_;
};

}

#endif // CHITECH_DERIVEDOBJECT_H
