#include "Connections.h"

#include "piper/Components/Component.h"

namespace piper::utils
{

// NOLINTBEGIN(misc-no-recursion)
void Nodalize(const std::string& component_name,
              const std::string& connection_name,
              const chi_mesh::Vector3& datum,
              std::map<std::string, std::shared_ptr<Component>>& components,
              std::set<std::string>& components_visited)
{
  if (components_visited.find(component_name) != components_visited.end())
    return;

  auto& component = components.at(component_name);
  component->Nodalize(connection_name, datum);

  components_visited.insert(component_name);

  for (const auto& connection : component->ConnectionPoints())
    Nodalize(connection.connected_comp_name_,
             connection.connected_comp_connection_point_name_,
             connection.position_,
             components,
             components_visited);
}
// NOLINTEND(misc-no-recursion)

} // namespace piper