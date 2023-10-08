#include "FieldEquationSystem.h"

#include "math/ParallelVector/GhostedParallelSTLVector.h"
#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/TimeIntegrators/ImplicitEulerTimeIntegrator.h"
#include "math/ParallelMatrix/ParallelMatrix.h"

#include "materials/MaterialPropertiesData.h"

#include "chi_log.h"

namespace chi_math
{

chi::InputParameters FieldEquationSystem::GetInputParameters()
{
  chi::InputParameters params = EquationSystem::GetInputParameters();

  params.AddRequiredParameterArray(
    "fields", "An array of FieldFunctionGridBased handles or names.");

  params.AddOptionalParameter("material_properties",
                              0,
                              "Optional handle to a parameter block of type"
                              "chi::MaterialPropertiesData.");

  params.AddOptionalParameter(
    "output_filename_base",
    "",
    "If present, VTK output will be exported to files with this base name.");

  return params;
}

FieldEquationSystem::FieldEquationSystem(const chi::InputParameters& params)
  : EquationSystem(params),
    material_properties_data_(GetOrMakeMaterialPropertiesData(params)),
    output_file_base_(
      params.GetParamValue<std::string>("output_filename_base")),
    primary_fields_container_(
      std::move(MakeMultifieldContainer(params.GetParam("fields")))),

    num_local_dofs_(primary_fields_container_->TotalNumLocalDOFs()),
    num_globl_dofs_(primary_fields_container_->TotalNumGlobalDOFs())
{
  main_solution_vector_ = std::move(MakeSolutionVector());
  old_solution_vectors_ = std::move(InitSolutionHistory());
  old_residual_vectors_ = std::move(InitResidualHistory());
}

// ##################################################################
const chi::MaterialPropertiesData&
FieldEquationSystem::GetOrMakeMaterialPropertiesData(
  const chi::InputParameters& params)
{
  if (params.ParametersAtAssignment().Has("material_properties"))
  {
    return Chi::GetStackItem<chi::MaterialPropertiesData>(
      /*stack=*/Chi::object_stack,
      /*handle=*/params.GetParamValue<size_t>("material_properties"),
      /*calling_function_name=*/__FUNCTION__);
  }
  else
  {
    auto empty_data = chi::MaterialPropertiesData::MakeEmptyData();
    Chi::object_stack.push_back(empty_data);
    const size_t handle = Chi::object_stack.size() - 1;
    empty_data->SetStackID(handle);

    return *empty_data;
  }
}

// ##################################################################
std::unique_ptr<ParallelVector> FieldEquationSystem::MakeSolutionVector()
{
  auto ghost_ids = primary_fields_container_->GetSystemGhostIDs();

  auto new_vec = std::make_unique<GhostedParallelSTLVector>(
    num_local_dofs_, num_globl_dofs_, ghost_ids, Chi::mpi.comm);

  return new_vec;
}

// ##################################################################
std::vector<std::unique_ptr<ParallelVector>>
FieldEquationSystem::InitSolutionHistory()
{
  std::vector<std::unique_ptr<ParallelVector>> old_solution_vectors;
  for (size_t t = 0; t < num_solution_histories_; ++t)
  {
    auto vec = main_solution_vector_->MakeClone();
    old_solution_vectors.push_back(std::move(vec));
  }

  return old_solution_vectors;
}

// ##################################################################
std::vector<std::unique_ptr<ParallelVector>>
FieldEquationSystem::InitResidualHistory()
{
  std::vector<std::unique_ptr<ParallelVector>> old_residual_vectors;
  for (size_t t = 0; t < num_residual_histories_; ++t)
  {
    auto vec = std::make_unique<ParallelSTLVector>(
      num_local_dofs_, num_globl_dofs_, Chi::mpi.comm);
    old_residual_vectors.push_back(std::move(vec));
  }

  return old_residual_vectors;
}

/**Returns the number of local DOFs across all unknowns.*/
int64_t FieldEquationSystem::NumLocalDOFs() const { return num_local_dofs_; }
/**Returns the number of local DOFs across all unknowns.*/
int64_t FieldEquationSystem::NumGlobalDOFs() const { return num_globl_dofs_; }

std::shared_ptr<MultiFieldContainer>
FieldEquationSystem::MakeMultifieldContainer(const chi::ParameterBlock& params)
{
  auto valid_params = MultiFieldContainer::GetInputParameters();
  chi::ParameterBlock param_block;
  param_block.AddParameter(params);

  valid_params.SetErrorOriginScope("chi_math::MultifieldContainer");
  valid_params.AssignParameters(param_block);

  auto container_ptr = std::make_shared<MultiFieldContainer>(valid_params);

  Chi::object_stack.push_back(container_ptr);

  const size_t handle = Chi::object_stack.size() - 1;

  container_ptr->SetStackID(handle);

  return container_ptr;
}

/**Uses the underlying system to build a sparsity pattern.*/
ParallelMatrixSparsityPattern
FieldEquationSystem::BuildMatrixSparsityPattern() const
{
  return primary_fields_container_->BuildMatrixSparsityPattern();
}

/**Updates the fields.*/
void FieldEquationSystem::UpdateFields()
{
  const auto& x = *main_solution_vector_;
  for (auto& field_info : *primary_fields_container_)
  {
    std::vector<double> local_solution(field_info.num_local_dofs_, 0.0);
    for (int64_t i = 0; i < field_info.num_local_dofs_; ++i)
      local_solution[i] = x[field_info.local_offset_ + i];

    field_info.field_->UpdateFieldVector(local_solution);
  }
}

/**Output fields to VTK. The filename passed via the options will be used
 * plus a time index (if transient). If the file basename is empty then
 * this method will not do anything.*/
void FieldEquationSystem::OutputFields(int time_index)
{
  if (output_file_base_.empty()) return;

  std::vector<std::shared_ptr<const chi_physics::FieldFunctionGridBased>>
    ff_list;
  for (const auto& field_info : *primary_fields_container_)
    ff_list.push_back(field_info.field_);

  const std::string file_base =
    (time_index < 0) ? output_file_base_
                     : output_file_base_ + "_" + std::to_string(time_index);

  chi_physics::FieldFunctionGridBased::ExportMultipleToVTK(file_base, ff_list);
}

} // namespace chi_math