#include "ComponentModel.h"

#include "piper/components/HardwareComponent.h"

#include "chi_log.h"

namespace piper
{

ComponentModel::ComponentModel(const HardwareComponent& hardware_component,
                               const chi_mesh::Cell* cell_ptr,
                               const std::vector<std::string>& variable_names)
  : hardware_component_(&hardware_component),
    cell_ptr_(cell_ptr),
    variable_names_(variable_names),
    vars_old_(MakeVariablesMap(variable_names)),
    vars_new_(MakeVariablesMap(variable_names))
{
}

std::map<std::string, double>
ComponentModel::MakeVariablesMap(const std::vector<std::string>& variable_names)
{
  std::map<std::string, double> vars_map;

  for (const auto& var_name : variable_names)
    vars_map.insert(std::make_pair(var_name, 0.0));

  return vars_map;
}

const HardwareComponent& ComponentModel::GetHardwareComponent() const
{
  return *hardware_component_;
}

const std::string& ComponentModel::Name() const
{
  return hardware_component_->Name();
}

const std::vector<utils::Connection>& ComponentModel::ConnectionPoints() const
{
  return hardware_component_->ConnectionPoints();
}

const Orientation& ComponentModel::GetOrientation() const
{
  return hardware_component_->GetOrientation();
}

size_t ComponentModel::GetID() const { return hardware_component_->GetID(); }

ComponentCategory ComponentModel::Category() const
{
  return hardware_component_->Category();
}

double ComponentModel::Area() const { return hardware_component_->Area(); }

double ComponentModel::Volume() const { return hardware_component_->Volume(); }

utils::FlowOrientation
ComponentModel::FlowOrientationRelToConPoint(size_t con_point_id) const
{
  return hardware_component_->FlowOrientationRelToConPoint(con_point_id);
}

const chi_mesh::Vector3& ComponentModel::GetRootNodePosition() const
{
  return hardware_component_->GetRootNodePosition();
}

chi_mesh::Vector3 ComponentModel::MakeCentroid() const
{
  return hardware_component_->MakeCentroid();
}

// ##################################################################
const std::vector<std::string>& ComponentModel::VarNames() const
{
  return variable_names_;
}

const double& ComponentModel::VarOld(const std::string& name) const
{
  auto& ref_map = vars_old_;
  auto it = ref_map.find(name);

  if (it == ref_map.end())
  {
    const auto& comp_name = hardware_component_->Name();
    ChiInvalidArgument("Invalid variable name \"" + name +
                       "\" for component model of "
                       "component \"" +
                       comp_name + "\"");
  }

  return (*it).second;
}

double& ComponentModel::VarOld(const std::string& name)
{
  auto& ref_map = vars_old_;
  auto it = ref_map.find(name);

  if (it == ref_map.end())
  {
    const auto& comp_name = hardware_component_->Name();
    ChiInvalidArgument("Invalid variable name \"" + name +
                       "\" for component model of "
                       "component \"" +
                       comp_name + "\"");
  }

  return (*it).second;
}

const double& ComponentModel::VarNew(const std::string& name) const
{
  auto& ref_map = vars_new_;
  auto it = ref_map.find(name);

  if (it == ref_map.end())
  {
    const auto& comp_name = hardware_component_->Name();
    ChiInvalidArgument("Invalid variable name \"" + name +
                       "\" for component model of "
                       "component \"" +
                       comp_name + "\"");
  }

  return (*it).second;
}

double& ComponentModel::VarNew(const std::string& name)
{
  auto& ref_map = vars_new_;
  auto it = ref_map.find(name);

  if (it == ref_map.end())
  {
    const auto& comp_name = hardware_component_->Name();
    ChiInvalidArgument("Invalid variable name \"" + name +
                       "\" for component model of "
                       "component \"" +
                       comp_name + "\"");
  }

  return (*it).second;
}
} // namespace piper