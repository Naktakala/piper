#include "Piper.h"

#include "ChiObjectFactory.h"
#include "piper/Components/HardwareComponent.h"


#include "chi_runtime.h"
#include "chi_log.h"

#define MsgJuncCompNonExistent(jname, comp_name)                               \
  ("Junction \"" + jname + "\" is specified to connect to component \"" +      \
   comp_name + "\". But it does not exist.")

#define MsgJuncCompConnNonExistent(jname, comp_name, comp_conn_name)           \
  ("Junction \"" + jname + "\" is specified to connect to component \"" +      \
   comp_name + "\" at connection point \"" + comp_conn_name +                  \
   "\". This connection point was not found on the component.")

// NOLINTBEGIN(performance-inefficient-string-concatenation)
namespace piper
{

RegisterChiObject(piper, Piper);

chi::InputParameters Piper::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.SetGeneralDescription("A module for simulating pipes or channels.");
  params.SetDocGroup("Piper");

  params.AddOptionalParameter("name", "PiperSystem", "System name.");

  params.AddRequiredParameterArray("components",
                                   "List of volumetric components.");
  params.AddRequiredParameterArray("connections",
                                   "List of junction components.");

  params.AddOptionalParameter(
    "root_component",
    std::string(),
    "The name of the component to bind to the system's datum. Each component "
    "has a local datum of {0,0,0}. If this parameter is supplied this "
    "component's global datum will be"
    " assigned to the system's datum (i.e., the parameter \"datum\". If this "
    "parameters is not supplied, then the first component's datum will be at "
    "\"datum\" and all other components will be assembled relative to this.");
  params.AddOptionalParameterArray(
    "datum",
    std::vector<double>{0, 0, 0},
    "A 3 component array of doubles specifying the system origin. The "
    "root_component will have it's first connection point attached to this "
    "datum.");

  params.AddOptionalParameter("print_nodalization",
                              false,
                              "If set to true, will print the nodalization of "
                              "the system in the initialization step.");

  return params;
}

Piper::Piper(const chi::InputParameters& params)
  : ChiObject(params),
    system_name_(params.GetParamValue<std::string>("name")),
    root_component_name_(params.GetParamValue<std::string>("root_component")),
    root_component_id_(0),
    datum_(chi_mesh::Vector3(params.GetParamVectorValue<double>("datum"))),
    print_nodalization_(params.GetParamValue<bool>("print_nodalization"))
{
  auto component_handles = params.GetParamVectorValue<size_t>("components");
  auto junction_handles = params.GetParamVectorValue<size_t>("connections");

  ChiInvalidArgumentIf(component_handles.empty(),
                       "No components added to system");
  ChiInvalidArgumentIf(junction_handles.empty(),
                       "No junctions added to system");

  //=================================================================
  /**Lambda for adding components to a map from chi_object handles.*/
  auto AddComponents =
    [this](const std::vector<size_t>& handles,
           const std::vector<ComponentCategory>& allowed_types,
           const std::string& param_name)
  {
    for (size_t handle : handles)
    {
      auto component_ptr_ = Chi::GetStackItemPtrAsType<HardwareComponent>(
        Chi::object_stack, handle, "Piper");

      const auto& name = component_ptr_->Name();
      const auto& type = component_ptr_->Category();

      if (hw_comp_name_2_id_map_.count(name) != 0)
        ChiInvalidArgument("Duplicate component name \"" + name + "\"");

      const size_t component_id = hardware_components_.size();
      component_ptr_->SetID(component_id);

      if (std::find(allowed_types.begin(), allowed_types.end(), type) ==
          allowed_types.end())
      {
        std::string allowed_type_names;
        for (const auto& allowed_type : allowed_types)
        {
          allowed_type_names.append("\"");
          allowed_type_names.append(ComponentCategoryName(allowed_type));
          allowed_type_names.append("\"");
          if (allowed_type != allowed_types.back())
            allowed_type_names.append(" ");
        }
        ChiLogicalError("Component \"" + name + "\" is of the wrong type \"" +
                        ComponentCategoryName(type) +
                        "\" when used for parameter \"" + param_name +
                        "\". Allowed types are " + allowed_type_names);
      }

      hw_comp_name_2_id_map_.insert(std::make_pair(name, component_id));
      hardware_components_.push_back(component_ptr_);

      // We don't need to check for duplicates here since it is already
      // asserted above.
      // clang-format off
      switch (type)
      {
        case ComponentCategory::BoundaryLike:
          boundary_component_ids_.push_back(component_id);break;
        case ComponentCategory::Volumetric:
          volume_component_ids_.push_back(component_id); break;
        case ComponentCategory::JunctionLike:
          junction_component_ids_.push_back(component_id);break;
      }
      // clang-format on
    }
  };
  //=================================================================
  hw_comp_name_2_id_map_.clear();
  hardware_components_.clear();
  boundary_component_ids_.clear();
  volume_component_ids_.clear();
  junction_component_ids_.clear();

  AddComponents(
    component_handles,
    {ComponentCategory::BoundaryLike, ComponentCategory::Volumetric},
    "components");
  AddComponents(
    junction_handles, {ComponentCategory::JunctionLike}, "connections");

  ConnectComponents();

  Chi::log.Log0Verbose1() << "Piper system \"" << SystemName() << "\" created.";
}

const std::string& Piper::SystemName() const { return system_name_; }
const std::vector<std::shared_ptr<HardwareComponent>>&
Piper::HardwareComponents() const
{
  return hardware_components_;
}

const std::vector<size_t>& Piper::BoundaryComponentIDs() const
{
  return boundary_component_ids_;
}
const std::vector<size_t>& Piper::VolumeComponentIDs() const
{
  return volume_component_ids_;
}
const std::vector<size_t>& Piper::JunctionComponentIDs() const
{
  return junction_component_ids_;
}

size_t Piper::MapHWCompName2ID(const std::string& name) const
{
  auto it = hw_comp_name_2_id_map_.find(name);
  ChiLogicalErrorIf(it == hw_comp_name_2_id_map_.end(),
                    "ID not found for component name \"" + name + "\".");

  return it->second;
}

size_t Piper::RootComponentID() const { return root_component_id_; }

} // namespace piper

// NOLINTEND(performance-inefficient-string-concatenation)