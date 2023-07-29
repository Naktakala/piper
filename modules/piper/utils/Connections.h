#ifndef PIPER_CONNECTIONS_H
#define PIPER_CONNECTIONS_H

#include "mesh/chi_mesh.h"
#include <string>
#include <map>

namespace piper
{
class Component;
}

namespace piper::utils
{

struct Connection
{
  std::string name_;
  std::string connected_comp_name_;
  std::string connected_comp_connection_point_name_;
  chi_mesh::Vector3 position_;
};

/**Recursively nodalizes connections*/
void Nodalize(const std::string& component_name,
              const std::string& connection_name,
              const chi_mesh::Vector3& datum,
              std::map<std::string, std::shared_ptr<Component>>& components,
              std::set<std::string>& components_visited);

} // namespace piper

#endif // PIPER_CONNECTIONS_H
