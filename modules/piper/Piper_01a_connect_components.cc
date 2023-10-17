#include "Piper.h"

#include "piper/components/HardwareComponent.h"

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

void Piper::ConnectComponents()
{
  //============================================= Build connectivity
  Chi::log.Log0Verbose1() << "Build connectivity";
  // Junctions are the only components that require fully resolved
  // connections at input-time. Therefore, we use them to complete
  // the connectivity of the volume components
  for (const auto& jid : junction_component_ids_)
  {
    auto& junction = *hardware_components_[jid];
    const auto& jname = junction.Name();
    auto& conpoints = junction.ConnectionPoints();

    size_t j = 0;
    for (auto& connection : conpoints)
    {
      const std::string& comp_name = connection.connected_comp_name_;
      const std::string& comp_conn_name =
        connection.connected_comp_connection_point_name_;

      // Check validity of connected component names
      auto it = hw_comp_name_2_id_map_.find(comp_name);
      ChiLogicalErrorIf(it == hw_comp_name_2_id_map_.end(),
                        MsgJuncCompNonExistent(jname, comp_name));
      const size_t comp_id = it->second;

      // Check validity of the connected component connection point
      auto& vol_component = *hardware_components_[comp_id];
      auto& vol_comp_conn_points = vol_component.ConnectionPoints();

      auto functor = [&comp_conn_name](const utils::Connection& vol_conpoint)
      { return vol_conpoint.name_ == comp_conn_name; };

      auto con_it = std::find_if(
        vol_comp_conn_points.begin(), vol_comp_conn_points.end(), functor);

      ChiLogicalErrorIf(
        con_it == vol_comp_conn_points.end(),
        MsgJuncCompConnNonExistent(jname, comp_name, comp_conn_name));
      const size_t comp_conn_id =
        std::distance(vol_comp_conn_points.begin(), con_it);

      // Complete junction's connection point info
      connection.connected_comp_id_ = comp_id;
      connection.connected_comp_connection_point_id_ = comp_conn_id;

      // Populate connected components connection point info
      con_it->connected_comp_name_ = jname;
      con_it->connected_comp_connection_point_name_ = comp_conn_name;
      con_it->connected_comp_id_ = jid;
      con_it->connected_comp_connection_point_id_ = j;

      ++j;
    }
  } // for each junction

  // Now run through each volume component and make sure its connected up
  for (const auto& component_ptr : hardware_components_)
  {
    const auto& conpoints = component_ptr->ConnectionPoints();
    for (auto& connection : conpoints)
    {
      ChiLogicalErrorIf(
        connection.connected_comp_name_.empty() or
          connection.connected_comp_connection_point_name_.empty(),
        "Component \"" + component_ptr->Name() +
          "\" has no connection at required connection point \"" +
          connection.name_ + "\".");
    }
  }

  //============================================= Fill in nodalization

  //=================================== Set root component
  // We need a root component, if it was not specified then we will use
  // the first volume/boundary component with id 0, and use that one
  if (root_component_name_.empty())
    root_component_id_ = hardware_components_.front()->GetID();
  else
  {
    auto it = hw_comp_name_2_id_map_.find(root_component_name_);

    ChiInvalidArgumentIf(it == hw_comp_name_2_id_map_.end(),
                         "The specified root component name \"" +
                           root_component_name_ +
                           "\" does not exist in the list of components.");
    auto candidate_ptr = hardware_components_[it->second];
    ChiInvalidArgumentIf(
      candidate_ptr->Category() != ComponentCategory::Volumetric,
      "The specified root component name \"" + root_component_name_ +
        "\" must be a volume component, i.e., not a boundary or a junction.");

    root_component_id_ = candidate_ptr->GetID();
  }

  //=================================== Recursively Nodalize
  std::set<size_t> components_visited;
  utils::Nodalize(
    root_component_id_, 0, datum_, hardware_components_, components_visited);

  //=================================== Perform post connection operations
  for (auto& component_ptr : hardware_components_)
    component_ptr->EstablishHardwareReferenceData(hardware_components_);

  //=================================== Print nodalization if needed
  if (print_nodalization_)
  {
    std::stringstream outstr;
    for (const auto& component_ptr : hardware_components_)
    {
      outstr << component_ptr->Name() << "\n";

      for (const auto& connection : component_ptr->ConnectionPoints())
        outstr << connection.name_ << " pos=" << connection.position_.PrintStr()
               << " connected to " << connection.connected_comp_name_ << ":"
               << connection.connected_comp_connection_point_name_ << "\n";
    }
    Chi::log.Log() << outstr.str();
  }
}

} // namespace piper

// NOLINTEND(performance-inefficient-string-concatenation)