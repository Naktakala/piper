#include "HeatFlux.h"

#include "piper/models/ComponentModel.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

namespace piper
{

// ##################################################################
chi::InputParameters HeatFlux::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<std::string>("type", "Private type name");

  return params;
}

HeatFlux::HeatFlux(const chi::InputParameters& params) : ChiObject(params) {}

// ##################################################################
RegisterChiObject(piper, ConstantHeatFlux);

chi::InputParameters ConstantHeatFlux::GetInputParameters()
{
  chi::InputParameters params = HeatFlux::GetInputParameters();

  params.AddOptionalParameter(
    "value", 0.0, "Strength of the volumetric heat source");

  return params;
}

ConstantHeatFlux::ConstantHeatFlux(const chi::InputParameters& params)
  : HeatFlux(params), value_(params.GetParamValue<double>("value"))
{
}

double ConstantHeatFlux::GetValue(
  const chi_mesh::Cell& cell,
  const std::map<std::string, double>& variable_reference_map_,
  double length) const
{
  return value_;
}

double ConstantHeatFlux::ComputeHCoeff(const ComponentModel& model,
                                       const chi_mesh::Vector3& gravity)
{
  return 0.0;
}

// ##################################################################
RegisterChiObject(piper, ConvectionHeatFlux);

chi::InputParameters ConvectionHeatFlux::GetInputParameters()
{
  chi::InputParameters params = HeatFlux::GetInputParameters();

  params.AddOptionalParameter("surface_temperature_value",
                              0.0,
                              "Constant value for the surface temperature.");

  params.AddOptionalParameter("surface_temperature_field",
                              "T_surf",
                              "Name of the surface_temperature_field");

  params.AddOptionalParameter(
    "wetted_perimeter",
    0.0,
    "Wetted perimetere, used to compute surface area.");

  return params;
}

ConvectionHeatFlux::ConvectionHeatFlux(const chi::InputParameters& params)
  : HeatFlux(params),
    surface_temperature_value_(
      params.GetParamValue<double>("surface_temperature_value")),
    surface_temperature_field_(GetSurfaceTemperatureField(params)),
    wetted_perimeter_(params.GetParamValue<double>("wetted_perimeter"))
{
}

std::shared_ptr<const chi_physics::FieldFunction>
ConvectionHeatFlux::GetSurfaceTemperatureField(
  const chi::InputParameters& params)
{
  if (not params.ParametersAtAssignment().Has("surface_temperature_field"))
    return nullptr;

  const std::string field_name =
    params.GetParamValue<std::string>("surface_temperature_field");
  for (auto& ff : Chi::field_function_stack)
    if (ff->TextName() == field_name) return ff;

  ChiInvalidArgument("Field function \"" + field_name + "\" not found");
}

double ConvectionHeatFlux::GetValue(
  const chi_mesh::Cell& cell,
  const std::map<std::string, double>& variable_reference_map_,
  double length) const
{
  auto hcoeff = variable_reference_map_.find("hcoeff");
  ChiInvalidArgumentIf(hcoeff == variable_reference_map_.end(),
                       "Component variable \"hcoeff\" not found.");

  auto T_bulk = variable_reference_map_.find("T");
  ChiInvalidArgumentIf(T_bulk == variable_reference_map_.end(),
                       "Component variable \"T\" not found.");

  const double area = wetted_perimeter_ * length;

  const double T_surf = EvaluateSurfaceTemperature(cell);

  return hcoeff->second * area * (T_surf - T_bulk->second);
}

double
ConvectionHeatFlux::EvaluateSurfaceTemperature(const chi_mesh::Cell& cell) const
{
  double T_surf = surface_temperature_value_;
  if (surface_temperature_field_ == nullptr) return T_surf;

  auto grid_ff = dynamic_cast<const chi_physics::FieldFunctionGridBased*>(
    &(*surface_temperature_field_));

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
      T_surf = grid_ff->FieldVectorRead()[dof_map];
    }
    else { T_surf = grid_ff->Evaluate(cell, cell.centroid_, 0); }
  }
  else
    T_surf = surface_temperature_field_->Evaluate(cell, cell.centroid_, 0);

  return T_surf;
}

double ConvectionHeatFlux::ComputeHCoeff(const ComponentModel& model,
                                         const chi_mesh::Vector3& gravity)
{
  const double g = gravity.Norm();
  const double beta = model.VarOld("beta");

  const double T_surf = EvaluateSurfaceTemperature(*model.GetCellPtr());
  const double T_bulk = model.VarOld("T");
  const double mu = model.VarOld("mu");
  const double rho = model.VarOld("rho");
  const double nu = mu / rho;
  const double Pr = model.VarOld("Pr");
  const double k = model.VarOld("k");

  const double Dh = model.HydraulicDiameter();

  const double Gr =
    std::fabs(
    g * std::fabs(beta) * (T_surf - T_bulk) * pow(Dh, 3.0) / pow(nu, 2.0));

  const double Ra = Gr * Pr;

  double h_natural;
  if (Ra < 1.0e9 and Ra > 0.0)
  {
    h_natural =
      (k / Dh) * (0.68 + 0.67 * pow(Ra, 0.25) /
                           pow(1.0 + pow(0.492 / Pr, 9.0 / 16.0), 4.0 / 9.0));
  }
  else if (Ra < 1.0e12)
  {
    h_natural = (k / Dh) * pow(0.825 + 0.387 * pow(Ra, 1.0 / 6.0) /
                                         pow(1.0 + pow(0.492 / Pr, 9.0 / 16.0),
                                             8.0 / 27.0),
                               2.0);
  }
  else
    h_natural = 0.0;

  const double Re = model.VarOld("Re");

  const double h_forced = (k / Dh) * 0.023 * pow(Re, 0.8) * pow(Pr, 0.4);

  return std::max(h_natural, h_forced);
}

} // namespace piper