#include "FEMMaterialPropertyInterface.h"

#include "math/KernelSystem/KernelSystem.h"
#include "materials/MaterialPropertiesData.h"
#include "materials/MaterialProperty.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_math
{

chi::InputParameters FEMMaterialPropertyInterface::GetInputParameters()
{
  chi::InputParameters params;

  // implicit parameter fem_data_handle

  return params;
}

FEMMaterialPropertyInterface::FEMMaterialPropertyInterface(
  const chi::InputParameters& params, const std::vector<int>& mat_id_scope)
  : fem_data_handle_(params.GetParamValue<size_t>("fem_data_handle")),
    material_properties_data_(GetFEMMaterialProperties(params)),
    mat_id_scope_(mat_id_scope)
{
}

const chi::MaterialPropertiesData&
FEMMaterialPropertyInterface::GetFEMMaterialProperties(
  const chi::InputParameters& params)
{
  auto& fem_data = Chi::GetStackItem<FEMKernelSystemData>(
    Chi::object_stack,
    params.GetParamValue<size_t>("fem_data_handle"),
    __FUNCTION__);

  return fem_data.mat_props_data_;
}

const chi_math::FEMMaterialProperty&
FEMMaterialPropertyInterface::GetFEMMaterialProperty(const std::string& name)
{
  // Since we are going to have to create a new one we can pre-emptively
  // grab the fem-data object
  const auto& fem_data = Chi::GetStackItem<FEMKernelSystemData>(
    Chi::object_stack, fem_data_handle_, __PRETTY_FUNCTION__);

  std::vector<const chi::MaterialProperty*> candidates;
  // First we just check which properties match the name
  for (const auto& prop : material_properties_data_.Properties())
    if (prop->TextName() == name) candidates.push_back(&(*prop));

  ChiLogicalErrorIf(candidates.size() > 1,
                    "More than property with name \"" + name + "\"");
  ChiLogicalErrorIf(candidates.empty(),
                    "No property with name \"" + name + "\" found.");

  // Making copy
  auto prop_mat_scope = candidates.front()->GetMaterialIDScope();
  const auto& this_mat_scope = mat_id_scope_;

  // If the prop_mat_scope is empty it means it applies to all materials,
  // therefore regardless of the object's scope this will be a match, and
  // vice versa
  if (prop_mat_scope.empty() or this_mat_scope.empty())
  {
    auto new_prop =
      std::make_unique<FEMMaterialProperty>(*candidates.front(), fem_data);
    coupled_material_props_.push_back(std::move(new_prop));
    return *coupled_material_props_.back();
  }

  for (int mat_id : this_mat_scope)
    for (int prop_mat_id : prop_mat_scope)
      if (mat_id == prop_mat_id)
      {
        auto new_prop =
          std::make_unique<FEMMaterialProperty>(*candidates.front(), fem_data);
        coupled_material_props_.push_back(std::move(new_prop));
        return *coupled_material_props_.back();
      }

  std::stringstream outstr;
  for (int mat_id : this_mat_scope)
    outstr << mat_id << " ";
  ChiLogicalError("Could not find material property defined on any of the "
                  "materials : " +
                  outstr.str());
}

void FEMMaterialPropertyInterface::PreComputeInternalMaterialProperties()
{
  for (const auto& prop : coupled_material_props_)
    prop->PreComputeInternalQPValues();
}

void FEMMaterialPropertyInterface::PreComputeFaceMaterialProperties()
{
  for (const auto& prop : coupled_material_props_)
    prop->PreComputeFaceQPValues();
}

} // namespace chi_math