#include "HeatGeneration.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"

#include "ChiObjectFactory.h"

namespace piper
{

// ##################################################################
chi::InputParameters HeatGeneration::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<std::string>("type", "Private type name");

  return params;
}

HeatGeneration::HeatGeneration(const chi::InputParameters& params)
  : ChiObject(params)
{
}

// ##################################################################
RegisterChiObject(piper, ConstantHeatGeneration);

chi::InputParameters ConstantHeatGeneration::GetInputParameters()
{
  chi::InputParameters params = HeatGeneration::GetInputParameters();

  params.AddOptionalParameter(
    "value", 0.0, "Strength of the volumetric heat source");

  return params;
}

ConstantHeatGeneration::ConstantHeatGeneration(
  const chi::InputParameters& params)
  : HeatGeneration(params), value_(params.GetParamValue<double>("value"))
{
}

double ConstantHeatGeneration::GetValue(const chi_mesh::Cell& cell) const
{
  return value_;
}

// ##################################################################
RegisterChiObject(piper, CoupledFieldHeatGeneration);

chi::InputParameters CoupledFieldHeatGeneration::GetInputParameters()
{
  chi::InputParameters params = HeatGeneration::GetInputParameters();

  params.AddRequiredParameter<std::string>("field_name",
                                           "Name of the field to couple to");

  return params;
}

CoupledFieldHeatGeneration::CoupledFieldHeatGeneration(
  const chi::InputParameters& params)
  : HeatGeneration(params),
    coupled_field_(
      GetCoupledField(params.GetParamValue<std::string>("field_name")))
{
}

const chi_physics::FieldFunction&
CoupledFieldHeatGeneration::GetCoupledField(const std::string& field_name)
{
  for (auto& ff : Chi::field_function_stack)
    if (ff->TextName() == field_name) return *ff;

  ChiInvalidArgument("Field function \"" + field_name + "\" not found");
}

double CoupledFieldHeatGeneration::GetValue(const chi_mesh::Cell& cell) const
{
  auto grid_ff =
    dynamic_cast<const chi_physics::FieldFunctionGridBased*>(&coupled_field_);

  if (grid_ff)
  {
    const auto& sdm = grid_ff->GetSpatialDiscretization();
    const auto& grid = sdm.Grid();

    const size_t num_local_cells = grid.local_cells.size();

    ChiLogicalErrorIf(cell.local_id_ >= num_local_cells,
                      "Cannot couple to a field based on a different mesh."
                      "Grid does not have the required cell.");

    const auto& ref_cell = grid.local_cells[cell.local_id_];

    ChiLogicalErrorIf(&ref_cell != &cell,
                      "Cannot couple to a field based on a different mesh.");

    if (sdm.Type() == chi_math::SpatialDiscretizationType::FINITE_VOLUME)
    {
      const int64_t dof_map = sdm.MapDOFLocal(cell, 0);
      return grid_ff->FieldVectorRead()[dof_map];
    }
    else { return grid_ff->Evaluate(cell, cell.centroid_, 0); }
  }

  return coupled_field_.Evaluate(cell, cell.centroid_, 0);
}

} // namespace piper