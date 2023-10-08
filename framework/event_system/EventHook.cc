#include "EventHook.h"

#include "event_system/SystemWideEventPublisher.h"
#include "event_system/Event.h"

#include "chi_log.h"

namespace chi
{

InputParameters EventHook::GetInputParameters()
{
  InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<std::string>("name", "Event hook name");

  params.AddOptionalParameterArray(
    "execute_on",
    std::vector<std::string>{"SolverInitialized",
                             "SolverPreAdvance"
                             "SolverAdvanced",
                             "SolverExecuted",
                             "ProgramExecuted"},
    "List of events at which the event-hook will execute.");

  return params;
}

EventHook::EventHook(const InputParameters& params)
  : ChiObject(params),
    name_(params.GetParamValue<std::string>("name")),
    subscribed_events_for_execution_(
      params.GetParamVectorValue<std::string>("execute_on"))
{
}

void EventHook::PushOntoStack(std::shared_ptr<ChiObject>& new_object)
{
  auto evhook_ptr = std::dynamic_pointer_cast<EventHook>(new_object);
  ChiLogicalErrorIf(not evhook_ptr,
                    "Failed to cast new object to chi::EventHook");

  Chi::object_stack.push_back(new_object);
  new_object->SetStackID(Chi::postprocessor_stack.size() - 1);

  auto new_subscriber =
    std::dynamic_pointer_cast<chi::EventSubscriber>(evhook_ptr);

  ChiLogicalErrorIf(
    not new_subscriber,
    "Failure to cast chi::PostProcessor to chi::EventSubscriber");

  auto& publisher = chi::SystemWideEventPublisher::GetInstance();
  publisher.AddSubscriber(new_subscriber);
}

void EventHook::ReceiveEventUpdate(const Event& event)
{
  auto it = std::find(subscribed_events_for_execution_.begin(),
                      subscribed_events_for_execution_.end(),
                      event.Name());

  if (it != subscribed_events_for_execution_.end())
  {

    Execute(event);
    if (Chi::log.GetVerbosity() >= 1)
      Chi::log.Log0Verbose1() << "Event hook \"" << Name()
                              << "\" executed on "
                                 "event \""
                              << event.Name() << "\".";
  }
}

} // namespace chi