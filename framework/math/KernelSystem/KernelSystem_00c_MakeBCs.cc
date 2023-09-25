#include "KernelSystem.h"

#include "FEMBCs/FEMBoundaryCondition.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

namespace chi_math
{

std::vector<FEMBoundaryConditionPtr>
KernelSystem::MakeBCs(const chi::ParameterBlock& boundary_condition_inputs,
                      size_t fem_data_handle)
{
  auto& object_factory = ChiObjectFactory::GetInstance();
  std::vector<FEMBoundaryConditionPtr> boundary_conditions;
  for (auto bc_input : boundary_condition_inputs) // making a copy
  {
    ChiInvalidArgumentIf(
      not bc_input.Has("type"),
      "A BoundaryCondition input is missing the \"type\" parameter");

    bc_input.AddParameter("fem_data_handle", fem_data_handle);
    bc_input.AddParameter(
      "multifieldcontainer_handles",
      std::vector<size_t>{primary_fields_container_->StackID()});

    const auto obj_type = bc_input.GetParamValue<std::string>("type");
    chi::InputParameters in_params =
      object_factory.GetRegisteredObjectParameters(obj_type);
    in_params.AssignParameters(bc_input);

    const size_t bc_handle =
      object_factory.MakeRegisteredObjectOfType(obj_type, in_params);

    auto boundary_condition = Chi::GetStackItemPtrAsType<FEMBoundaryCondition>(
      Chi::object_stack, bc_handle, __FUNCTION__);

    boundary_conditions.push_back(boundary_condition);
  }

  return boundary_conditions;
}

/**Obtains the boundary kernel associated with a boundary id.*/
FEMBoundaryConditionPtr KernelSystem::GetBoundaryCondition(uint64_t boundary_id)
{
  if (not(EquationTermsScope() & EqTermScope::BOUNDARY_TERMS)) return nullptr;

  const auto& current_field =
    primary_fields_container_->GetFieldBlockInfo(current_field_index_).field_;
  const auto& field_name = current_field->TextName();

  auto& bid2bc_map =
    varname_comp_2_bid2bc_map_.at({field_name, current_field_component_});

  auto bc_it = bid2bc_map.find(boundary_id);

  if (bc_it == bid2bc_map.end()) return nullptr;

  const auto& bndry_condition = bc_it->second;

  const auto& bc_var_comp = bndry_condition->ActiveVariableAndComponent();

  if (bc_var_comp.first != current_field->TextName() or
      bc_var_comp.second != current_field_component_)
    return nullptr;

  return bc_it->second;
}

} // namespace chi_math