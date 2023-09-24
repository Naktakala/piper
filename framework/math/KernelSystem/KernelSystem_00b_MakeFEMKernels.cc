#include "KernelSystem.h"

#include "FEMKernels/FEMKernel.h"
#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

namespace chi_math
{

std::vector<FEMKernelPtr>
KernelSystem::MakeFEMKernels(const chi::ParameterBlock& volume_kernels_inputs,
                             size_t fem_data_handle)
{
  auto& object_factory = ChiObjectFactory::GetInstance();

  std::vector<FEMKernelPtr> kernels;
  for (auto kernel_input : volume_kernels_inputs) // making a copy
  {
    ChiInvalidArgumentIf(not kernel_input.Has("type"),
                         "A kernel input is missing the \"type\" parameter");

    kernel_input.AddParameter("fem_data_handle", fem_data_handle);

    const auto obj_type = kernel_input.GetParamValue<std::string>("type");
    chi::InputParameters in_params =
      object_factory.GetRegisteredObjectParameters(obj_type);
    in_params.AssignParameters(kernel_input);

    const size_t kernel_handle =
      object_factory.MakeRegisteredObjectOfType(obj_type, in_params);

    auto kernel = Chi::GetStackItemPtrAsType<FEMKernel>(
      Chi::object_stack, kernel_handle, __FUNCTION__);

    kernels.push_back(kernel);
  }

  return kernels;
}

/**Obtains the kernels associated with a material.*/
std::vector<FEMKernelPtr> KernelSystem::GetMaterialKernels(int mat_id)
{
  auto iter = matid_2_volume_kernels_map_.find(mat_id);

  ChiLogicalErrorIf(iter == matid_2_volume_kernels_map_.end(),
                    "No kernel for material id " + std::to_string(mat_id));

  const auto& kernels = iter->second;

  ChiLogicalErrorIf(kernels.empty(),
                    "No kernel for material id " + std::to_string(mat_id));

  std::vector<std::shared_ptr<FEMKernel>> filtered_kernels;
  const bool time_terms_active = QueryTermsActive(EqTermScope::TIME_TERMS);
  const bool domain_terms_active = QueryTermsActive(EqTermScope::DOMAIN_TERMS);
  // const auto& current_field =
  // field_block_info_.at(current_field_index_).field_;
  const auto& current_field =
    primary_fields_container_->GetFieldBlockInfo(current_field_index_).field_;

  for (const auto& kernel : kernels)
  {
    bool kernel_active = false;
    if (kernel->IsTimeKernel() and time_terms_active) kernel_active = true;
    if (not kernel->IsTimeKernel() and domain_terms_active)
      kernel_active = true;

    if (kernel->ActiveVariableAndComponent().first != current_field->TextName())
      kernel_active = false;
    if (kernel->ActiveVariableAndComponent().second != current_field_component_)
      kernel_active = false;

    if (kernel_active) filtered_kernels.push_back(kernel);
  }

  if (time_terms_active and filtered_kernels.empty())
    ChiLogicalError("No time kernels in system");

  if (filtered_kernels.empty())
    ChiLogicalError("No kernels for material id " + std::to_string(mat_id) +
                    " on field \"" + current_field->TextName() + "\"");

  return filtered_kernels;
}

} // namespace chi_math