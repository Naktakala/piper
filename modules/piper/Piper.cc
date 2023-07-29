#include "Piper.h"

#include "ChiObjectFactory.h"
#include "piper/Components/Component.h"

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

  params.AddOptionalParameter(
    "name", "PiperSystem", "A name to associated with this system.");

  params.AddRequiredParameterArray("components",
                                   "List of volumetric components");
  params.AddRequiredParameterArray("connections",
                                   "List of junction components");

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
    "A 3 component array of doubles specifying the system origin.");

  params.AddOptionalParameter("print_nodalization",
                              false,
                              "If set to true, will print the nodalization of "
                              "the system in the initialization step.");

  return params;
}

Piper::Piper(const chi::InputParameters& params)
  : chi_physics::Solver(params),
    system_name_(params.GetParamValue<std::string>("name")),
    root_component_(params.GetParamValue<std::string>("root_component")),
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
  size_t component_id = 0;
  /**Lambda for adding components to a map from chi_object handles.*/
  auto AddComponents =
    [this, &component_id](const std::vector<size_t>& handles,
                          const std::vector<ComponentType>& allowed_types,
                          const std::string& param_name)
  {
    for (size_t handle : handles)
    {
      auto component_ptr_ = Chi::GetStackItemPtrAsType<Component>(
        Chi::object_stack, handle, "Piper");

      const auto& name = component_ptr_->Name();
      const auto& type = component_ptr_->Type();

      if (components_.count(name) != 0)
        ChiInvalidArgument("Duplicate component name \"" + name + "\"");

      component_ptr_->SetID(component_id++);

      if (std::find(allowed_types.begin(), allowed_types.end(), type) ==
          allowed_types.end())
      {
        std::string allowed_type_names;
        for (const auto& allowed_type : allowed_types)
        {
          allowed_type_names.append("\"");
          allowed_type_names.append(ComponentTypeName(allowed_type));
          allowed_type_names.append("\"");
          if (allowed_type != allowed_types.back())
            allowed_type_names.append(" ");
        }
        ChiLogicalError("Component \"" + name + "\" is of the wrong type \"" +
                        ComponentTypeName(type) +
                        "\" when used for parameter \"" + param_name +
                        "\". Allowed types are " + allowed_type_names);
      }

      components_[name] = component_ptr_;

      // We don't need to check for duplicates here since it is already
      // asserted above.
      // clang-format off
      switch (type)
      {
        case ComponentType::BoundaryLike:
          boundary_component_names_.push_back(name);break;
        case ComponentType::Volumetric:
          volume_component_names_.push_back(name); break;
        case ComponentType::JunctionLike:
          junction_component_names_.push_back(name);break;
      }
      // clang-format on
    }
  };
  //=================================================================
  components_.clear();
  AddComponents(component_handles,
                {ComponentType::BoundaryLike, ComponentType::Volumetric},
                "components");
  AddComponents(junction_handles, {ComponentType::JunctionLike}, "connections");

  Chi::log.Log0Verbose1() << "Piper system \"" << SystemName() << "\" created.";
}

const std::string& Piper::SystemName() const { return system_name_; }



} // namespace piper

// NOLINTEND(performance-inefficient-string-concatenation)