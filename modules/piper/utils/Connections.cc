#include "Connections.h"

#include "piper/components/HardwareComponent.h"

namespace piper::utils
{

// NOLINTBEGIN(misc-no-recursion)
void Nodalize(size_t component_id,
              size_t connection_id,
              const chi_mesh::Vector3& datum,
              std::vector<std::shared_ptr<HardwareComponent>>& components,
              std::set<size_t>& components_visited)
{
  if (components_visited.find(component_id) != components_visited.end()) return;

  auto& component = components.at(component_id);
  component->Nodalize(connection_id, datum);

  components_visited.insert(component_id);

  for (const auto& connection : component->ConnectionPoints())
    Nodalize(connection.connected_comp_id_,
             connection.connected_comp_connection_point_id_,
             connection.position_,
             components,
             components_visited);
}
// NOLINTEND(misc-no-recursion)

} // namespace piper::utils