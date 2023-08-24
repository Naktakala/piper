#ifndef PIPER_COMPONENT_MODEL_H
#define PIPER_COMPONENT_MODEL_H

#include <vector>
#include <string>
#include <map>

#include "piper/utils/utils.h"
#include "piper/utils/Connections.h"
#include "piper/models/ConnectionInterface.h"

#include "mesh/chi_mesh.h"

namespace chi_mesh
{
class Cell;
}

namespace piper::utils
{
struct Connection;
}

namespace piper
{

class HardwareComponent;
class Orientation;

class ComponentModel
{
public:
  ComponentModel(std::vector<std::unique_ptr<ComponentModel>>& family,
                 const HardwareComponent& hardware_component,
                 const chi_mesh::Cell* cell,
                 const std::vector<std::string>& variable_names);

  const HardwareComponent& GetHardwareComponent() const;

  const std::string& Name() const;
  const std::vector<utils::Connection>& ConnectionPoints() const;
  const Orientation& GetOrientation() const;
  size_t GetID() const;

  ComponentCategory Category() const;
  double Area() const;
  double Volume() const;
  double Length() const;
  double HydraulicDiameter() const;

  /**Returns true if the component's flow orientation is outgoing relative
  * to the connection point.*/
  bool IsOutgoingRelToConPoint(size_t con_point_id) const;

  /**Returns the component's flow orientation relative to the connection
   * point.*/
  utils::FlowOrientation
  FlowOrientationRelToConPoint(size_t con_point_id) const;

  const chi_mesh::Vector3& GetRootNodePosition() const;
  chi_mesh::Vector3 MakeCentroid() const;

  std::vector<std::unique_ptr<ComponentModel>>& Family();
  ConnectionInterface Connections();

  const std::vector<std::string>& VarNames() const;

  const double& VarOld(const std::string& name) const;
  double& VarOld(const std::string& name);

  const double& VarNew(const std::string& name) const;
  double& VarNew(const std::string& name);

protected:
  static std::map<std::string, double>
  MakeVariablesMap(const std::vector<std::string>& variable_names);

  std::vector<std::unique_ptr<ComponentModel>>& family_;
  const HardwareComponent* hardware_component_;
  const chi_mesh::Cell* cell_ptr_;
  const std::vector<std::string> variable_names_;
  std::map<std::string, double> vars_old_;
  std::map<std::string, double> vars_new_;
};

} // namespace piper

#endif // PIPER_COMPONENT_MODEL_H
