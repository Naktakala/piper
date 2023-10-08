#ifndef CHITECH_EVENTHOOK_H
#define CHITECH_EVENTHOOK_H

#include "ChiObject.h"
#include "event_system/EventSubscriber.h"

namespace chi
{

class EventHook : public ChiObject, public EventSubscriber
{
public:
  static InputParameters GetInputParameters();
  explicit EventHook(const InputParameters& params);

  void PushOntoStack(std::shared_ptr<ChiObject>& new_object) override;

  void ReceiveEventUpdate(const Event& event) override;

  const std::string& Name() const {return name_;}

  virtual void Execute(const Event& event) = 0;

private:
  const std::string name_;
  const std::vector<std::string> subscribed_events_for_execution_;
};

}

#endif // CHITECH_EVENTHOOK_H
