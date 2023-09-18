#include "HeatConductionSystem.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

#include "ChiObjectFactory.h"

namespace hcm
{

RegisterChiObject(hcm, HeatConductionSystem);

chi::InputParameters HeatConductionSystem::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.SetGeneralDescription(
    "A system for the simulation of heat conduction");
  params.SetDocGroup("doc_HeatConduction");

  params.AddOptionalParameter(
    "name", "HeatConductionSystem", "Name of the system.");

  params.AddOptionalParameter(
    "temperature_ff_name",
    "T",
    "Name to give to the field function holding the temperature.");

  params.AddRequiredParameterArray("kernels",
                                   "A list of kernel in parameter form.");
  params.AddRequiredParameterArray("bcs",
                                   "A list of boundary condition handles.");

  return params;
}

HeatConductionSystem::HeatConductionSystem(const chi::InputParameters& params)
  : ChiObject(params),
    temperature_ff_name_(
      params.GetParamValue<std::string>("temperature_ff_name"))
{
  try
  {
    grid_ptr_ = chi_mesh::GetCurrentHandler().GetGrid();
  }
  catch (const std::exception& e)
  {
    ChiLogicalError("Failed to get grid. Make sure a mesh has been created.");
  }

  sdm_ptr_ = chi_math::SpatialDiscretization_PWLC::New(Grid());

  const auto volume_kernels_input = params.GetParam("kernels");

  for (const auto& kernel_input : volume_kernels_input)
    volume_kernel_inputs_.push_back(kernel_input);

  const auto bc_kernels_input = params.GetParam("bcs");

  for (const auto& bc_input : bc_kernels_input)
    boundary_condition_inputs_.push_back(bc_input);
}

const std::string& HeatConductionSystem::GetTemperatureFFName() const
{
  return temperature_ff_name_;
}

const chi_mesh::MeshContinuum& HeatConductionSystem::Grid() const
{
  ChiLogicalErrorIf(not grid_ptr_,
                    "No grid active on solver. Called Initialize?");

  return *grid_ptr_;
}

const chi_math::SpatialDiscretization& HeatConductionSystem::SDM() const
{
  ChiLogicalErrorIf(
    not sdm_ptr_,
    "No spatial discretization active on solver. Called Initialize?");

  return *sdm_ptr_;
}

chi_math::SpatialDiscretizationPtr& HeatConductionSystem::SDMPtr()
{
  ChiLogicalErrorIf(
    not sdm_ptr_,
    "No spatial discretization active on solver. Called Initialize?");

  return sdm_ptr_;
}

const std::vector<chi::ParameterBlock>&
HeatConductionSystem::VolumeKernelInputs() const
{
  return volume_kernel_inputs_;
}

std::vector<chi::ParameterBlock>&
HeatConductionSystem::BoundaryConditionInputs()
{
  return boundary_condition_inputs_;
}

} // namespace hcm