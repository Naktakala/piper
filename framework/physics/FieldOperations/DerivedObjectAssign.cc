#include "DerivedObjectAssign.h"

#include "derived_objects/DerivedObject.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"

#include "ChiObjectFactory.h"

#include "chi_runtime.h"

namespace chi_physics::field_operations
{

RegisterChiObject(chi_physics::field_operations, DerivedObjectAssign);

chi::InputParameters DerivedObjectAssign::GetInputParameters()
{
  chi::InputParameters params = FieldOperation::GetInputParameters();
  params += GridBasedFieldFunctionInterface::GetInputParameters();

  params.SetGeneralDescription("Assigns the field-value to the spatial values"
                               " associated with a derived object.");

  params.AddRequiredParameter<size_t>("derived_object",
                                      "Handle to a derived object");

  return params;
}

DerivedObjectAssign::DerivedObjectAssign(const chi::InputParameters& params)
  : FieldOperation(params),
    GridBasedFieldFunctionInterface(params),
    derived_object_(Chi::GetStackItem<chi::DerivedObject>(
      Chi::object_stack,
      params.GetParamValue<size_t>("derived_object"),
      __FUNCTION__))
{
}

void DerivedObjectAssign::Execute()
{
  auto gridbased_ff = GetGridBasedFieldFunction();

  ChiLogicalErrorIf(not gridbased_ff, "Invalid field function");

  const auto& sdm = gridbased_ff->GetSpatialDiscretization();
  const auto& grid = sdm.Grid();
  auto& local_data = gridbased_ff->FieldVector();

  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto& node_locations = cell_mapping.GetNodeLocations();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t dof_map =
        sdm.MapDOFLocal(cell, 0, gridbased_ff->GetUnknownManager(), 0, 0);

      local_data[dof_map] = derived_object_.SpatialValue(node_locations[i]);
    }
  } // for cell
}

} // namespace chi_physics::field_operations