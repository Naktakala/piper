#include "Piper.h"

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

void Piper::Initialize()
{
  Chi::mpi.Barrier();

  //============================================= Build connectivity
  Chi::log.Log0Verbose1() << "Build connectivity";
  // Junctions are the only components that require fully resolved
  // connections at input-time. Therefore, we use them to complete
  // the connectivity of the volume components
  for (const auto& jname : junction_component_names_)
  {
    const auto& junction = components_[jname];
    const auto& conpoints = junction->ConnectionPoints();

    for (auto& connection : conpoints)
    {
      const std::string& comp_name = connection.connected_comp_name_;
      const std::string& comp_conn_name =
        connection.connected_comp_connection_point_name_;

      // Check validity of connected component names
      ChiLogicalErrorIf(components_.count(comp_name) == 0,
                        MsgJuncCompNonExistent(jname, comp_name));

      // Check validity of the connected component connection point
      auto& volume_component = components_[comp_name];
      bool vol_conpoint_found = false;
      for (auto& vol_conpoint : volume_component->ConnectionPoints())
        if (vol_conpoint.name_ == comp_conn_name)
        {
          vol_conpoint.connected_comp_name_ = jname;
          vol_conpoint.connected_comp_connection_point_name_ = connection.name_;
          vol_conpoint_found = true;
        }

      ChiLogicalErrorIf(
        not vol_conpoint_found,
        MsgJuncCompConnNonExistent(jname, comp_name, comp_conn_name));
    }
  } // for each junction

  // Now run through each volume component and make sure its connected up
  for (const auto& [name, component] : components_)
  {
    const auto& conpoints = component->ConnectionPoints();
    for (auto& connection : conpoints)
    {
      ChiLogicalErrorIf(
        connection.connected_comp_name_.empty() or
          connection.connected_comp_connection_point_name_.empty(),
        "Component \"" + name +
          "\" has no connection at required connection point \"" +
          connection.name_ + "\".");
    }
  }

  //============================================= Fill in nodalization

  // We need a root component, if it was not specified then we will
  // the volume component with id 0, and use that one
  if (root_component_.empty())
  {
    for (const auto& [name, component] : components_)
      if (component->GetID() == 0) root_component_ = name;
  }
  else
    ChiInvalidArgumentIf(
      components_.count(root_component_) == 0,
      "The specified root component name \"" + root_component_ +
        "\" does not exist in the list of volume components."
        " This parameter must be a valid volume component name");

  std::string root_component_datum_name =
    components_[root_component_]->ConnectionPoints().front().name_;

  // Nodalize
  std::set<std::string> components_visited;
  utils::Nodalize(root_component_,
                  root_component_datum_name,
                  datum_,
                  components_,
                  components_visited);

  // Print nodalization
  if (print_nodalization_)
  {
    std::stringstream outstr;
    for (const auto& [comp_name, component] : components_)
    {
      outstr << comp_name << "\n";

      for (const auto& connection : component->ConnectionPoints())
        outstr << connection.name_ << " pos=" << connection.position_.PrintStr()
               << " connected to " << connection.connected_comp_name_ << ":"
               << connection.connected_comp_connection_point_name_ << "\n";
    }
    Chi::log.Log() << outstr.str();
  }
}

}

// NOLINTEND(performance-inefficient-string-concatenation)