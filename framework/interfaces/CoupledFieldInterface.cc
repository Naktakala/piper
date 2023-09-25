#include "CoupledFieldInterface.h"

#include "math/KernelSystem/KernelSystem.h"
#include "math/KernelSystem/Coupling/FEMCoupledField.h"
#include "math/Containers/MultiFieldContainer.h"
#include "math/KernelSystem/Coupling/FFContainerFEMCoupledVariable.h"
#include "math/KernelSystem/Coupling/FieldFunctionFEMCoupledField.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_math
{

chi::InputParameters CoupledFieldInterface::GetInputParameters()
{
  chi::InputParameters params;

  params.AddOptionalParameterArray(
    "multifieldcontainer_handles",
    std::vector<size_t>{},
    "A list of handles to MultiFieldContainer that should be searched for "
    "coupling before searching the regular field functions.");

  // implicit parameter fem_data_handle

  return params;
}

CoupledFieldInterface::CoupledFieldInterface(const chi::InputParameters& params)
  : multifieldcontainer_handles_(
      params.GetParamVectorValue<size_t>("multifieldcontainer_handles")),
    fem_data_handle_(params.GetParamValue<size_t>("fem_data_handle"))
{
}

const FEMCoupledField&
CoupledFieldInterface::GetCoupledField(const std::string& field_name)
{
  // First check if we already have it
  for (const auto& ptr : coupled_field_ptrs_)
    if (ptr->FieldName() == field_name) return *ptr;

  // Since we are going to have to create a new one we can pre-emptively
  // grab the fem-data object
  auto fem_data_ptr = Chi::GetStackItemPtrAsType<FEMKernelSystemData>(
    Chi::object_stack, fem_data_handle_, __PRETTY_FUNCTION__);

  // Now check any MultiFieldContainer objects
  for (size_t handle : multifieldcontainer_handles_)
  {
    auto container_ptr =
      Chi::GetStackItemPtrAsType<chi_math::MultiFieldContainer>(
        Chi::object_stack, handle, __PRETTY_FUNCTION__);

    size_t field_index = 0;
    for (const auto& field_info : *container_ptr)
    {
      if (field_info.field_->TextName() == field_name)
      {
        auto new_coupling = std::make_unique<FFContainerFEMCoupledVariable>(
          field_name, *fem_data_ptr, *container_ptr, field_index);

        coupled_field_ptrs_.push_back(std::move(new_coupling));

        return *coupled_field_ptrs_.back();
      }
      ++field_index;
    } // for field in container
  }   // for container

  // Now check regular field functions
  for (const auto& field_ptr : Chi::field_function_stack)
  {
    if (field_ptr->TextName() == field_name)
    {
      auto field_grid_ptr =
        std::dynamic_pointer_cast<chi_physics::FieldFunctionGridBased>(
          field_ptr);

      ChiLogicalErrorIf(not field_grid_ptr,
                        "Failed to cast field \"" + field_name +
                          "\" to FieldFunctionGridBased.");

      auto new_coupling = std::make_unique<FieldFunctionFEMCoupledField>(
        field_name, *fem_data_ptr, *field_grid_ptr);



      coupled_field_ptrs_.push_back(std::move(new_coupling));

      return *coupled_field_ptrs_.back();
    }
  } // for field in stack

  ChiLogicalError("Coupled field with name \"" + field_name +
                  "\" could not be found.");
}

void CoupledFieldInterface::PreComputeInternalCoupledFields()
{
  for (auto& coupled_field_ptr : coupled_field_ptrs_)
    coupled_field_ptr->ComputeFieldInternalQPValues();
}

void CoupledFieldInterface::PreComputeFaceCoupledFields()
{
  for (auto& coupled_field_ptr : coupled_field_ptrs_)
    coupled_field_ptr->ComputeFieldFaceQPValues();
}

} // namespace chi_math