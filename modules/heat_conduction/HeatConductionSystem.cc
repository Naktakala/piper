#include "HeatConductionSystem.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"
#include "math/KernelSystem/FEMKernels/FEMKernel.h"
#include "math/KernelSystem/FEMBCs/FEMBoundaryCondition.h"

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

  params.AddRequiredParameterArray("kernels", "A list of kernel handles.");
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

  const auto kernel_handles = params.GetParamVectorValue<size_t>("kernels");
  const auto bc_handles = params.GetParamVectorValue<size_t>("bcs");

  for (size_t kernel_handle : kernel_handles)
  {
    auto kernel_ptr = Chi::GetStackItemPtrAsType<chi_math::FEMKernel>(
      Chi::object_stack, kernel_handle, __FUNCTION__);

    volume_kernels_.push_back(kernel_ptr);
  }

  for (size_t bc_handle : bc_handles)
  {
    auto bc_ptr = Chi::GetStackItemPtrAsType<chi_math::FEMBoundaryCondition>(
      Chi::object_stack, bc_handle, __FUNCTION__);

    boundary_conditions_.push_back(bc_ptr);
  }
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

std::vector<chi_math::FEMKernelPtr>&
HeatConductionSystem::VolumeKernels()
{
  return volume_kernels_;
}

std::vector<chi_math::FEMBoundaryConditionPtr>&
HeatConductionSystem::BoundaryConditions()
{
  return boundary_conditions_;
}

} // namespace hcm