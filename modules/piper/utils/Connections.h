#ifndef PIPER_CONNECTIONS_H
#define PIPER_CONNECTIONS_H

#include "mesh/chi_mesh.h"
#include <string>
#include <map>

namespace piper
{
class HardwareComponent;
}

namespace piper::utils
{

enum class FlowOrientation
{
  INCOMING, OUTGOING
};

struct Connection
{
  std::string name_;
  std::string connected_comp_name_;
  std::string connected_comp_connection_point_name_;
  size_t connected_comp_id_ = -1;
  size_t connected_comp_connection_point_id_ = -1;
  chi_mesh::Vector3 position_;
};

/**Recursively nodalizes connections*/
void Nodalize(size_t component_id,
              size_t connection_id,
              const chi_mesh::Vector3& datum,
              std::vector<std::shared_ptr<HardwareComponent>>& components,
              std::set<size_t>& components_visited);

} // namespace piper::utils

#endif // PIPER_CONNECTIONS_H
