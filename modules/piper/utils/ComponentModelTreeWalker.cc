#include "ComponentModelTreeWalker.h"

#include "piper/models/ComponentModel.h"
#include "piper/components/HardwareComponent.h"

#include "chi_log.h"

namespace piper
{

ComponentModelTreeWalker::ComponentModelTreeWalker(
  ModelMap& model_map,
  size_t starting_component_id,
  const OperatorFunction& operator_function)
  : model_map_(model_map),
    starting_component_id_(starting_component_id),
    op_(operator_function)
{
}

void ComponentModelTreeWalker::Execute()
{
  std::set<size_t> comp_ids_visited;

  ChiLogicalErrorIf(starting_component_id_ >= model_map_.size(),
                    "Bad starting component id");
  auto& starting_model_ptr = model_map_[starting_component_id_];

  Walk(*starting_model_ptr, comp_ids_visited);
}

void ComponentModelTreeWalker::Execute(
  ModelMap& model_map,
  size_t starting_component_id,
  const OperatorFunction& operator_function)
{
  ComponentModelTreeWalker walker(
    model_map, starting_component_id, operator_function);

  walker.Execute();
}

// NOLINTBEGIN(misc-no-recursion)
void ComponentModelTreeWalker::Walk(ComponentModel& component,
                                    std::set<size_t>& comp_ids_visited)
{
  const auto& hw_comp = component.GetHardwareComponent();
  if (comp_ids_visited.find(hw_comp.GetID()) == comp_ids_visited.end())
  {
    op_(component);

    comp_ids_visited.insert(hw_comp.GetID());

    const auto& connections = hw_comp.ConnectionPoints();
    for (const auto& connection : connections)
    {
      ChiLogicalErrorIf(connection.connected_comp_id_ >= model_map_.size(),
                        "Bad starting component id");
      auto& connected_model_ptr = model_map_[connection.connected_comp_id_];

      Walk(*connected_model_ptr, comp_ids_visited);
    }
  }
}
// NOLINTEND(misc-no-recursion)

} // namespace piper