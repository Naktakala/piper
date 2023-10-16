#include "DerivedObject.h"

#include "physics/PhysicsEventPublisher.h"
#include "event_system/Event.h"

#include "chi_log.h"

namespace chi
{

InputParameters DerivedObject::GetInputParameters()
{
  InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<std::string>("name",
                                           "Name of the derived object");

  params.AddOptionalParameterArray(
    "execute_on",
    std::vector<std::string>{
      "SolverInitialized", "SolverAdvanced", "SolverExecuted"},
    "List of events at which the derived object will update its state.");

  params.AddOptionalParameter(
    "solvername_filter",
    "",
    "Controls update events to only execute on the relevant solver's"
    "event calls.");

  return params;
}

DerivedObject::DerivedObject(const InputParameters& params)
  : ChiObject(params),
    name_(params.GetParamValue<std::string>("name")),
    subscribed_events_for_execution_(
      params.GetParamVectorValue<std::string>("execute_on"))
{
}

void DerivedObject::PushOntoStack(std::shared_ptr<ChiObject>& new_object)
{
  ChiObject::PushOntoStack(new_object);

  auto new_subscriber =
    std::dynamic_pointer_cast<chi::EventSubscriber>(new_object);

  ChiLogicalErrorIf(
    not new_subscriber,
    "Failure to cast chi::DerivedObject to chi::EventSubscriber");

  auto& publisher = chi_physics::PhysicsEventPublisher::GetInstance();
  publisher.AddSubscriber(new_subscriber);
}

void DerivedObject::ReceiveEventUpdate(const Event& event)
{
  auto it = std::find(subscribed_events_for_execution_.begin(),
                      subscribed_events_for_execution_.end(),
                      event.Name());

  if (it != subscribed_events_for_execution_.end())
  {
    if (event.Code() >= 31 and event.Code() <= 38 and
        not solvername_filter_.empty())
    {
      if (event.Parameters().GetParamValue<std::string>("solver_name") !=
          solvername_filter_)
        return;
    }

    Update(event);
    if (Chi::log.GetVerbosity() >= 1)
      Chi::log.Log0Verbose1() << "DerivedObject \"" << Name()
                              << "\" executed on "
                                 "event \""
                              << event.Name() << "\".";
  }
}

const std::string& DerivedObject::Name() const { return name_; }

double DerivedObject::SpatialValue(const chi_mesh::Vector3& position) const
{
  ChiLogicalError("Method not available");
}

} // namespace chi