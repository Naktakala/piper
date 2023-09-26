#include "EquationSystem.h"

#include "math/ParallelVector/GhostedParallelSTLVector.h"
#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/spatial_discretization.h"
#include "math/TimeIntegrators/ImplicitEulerTimeIntegrator.h"
#include "math/ParallelMatrix/ParallelMatrix.h"

#include "materials/MaterialPropertiesData.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

namespace chi_math
{

chi::InputParameters EquationSystem::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameterArray(
    "fields", "An array of FieldFunctionGridBased handles or names.");

  params.AddOptionalParameter("time_integrator",
                              0,
                              "Handle to a TimeIntegrator object. A default "
                              "ImplicitEuler time integrator will"
                              "be used for transient simulations.");
  params.AddOptionalParameter("verbosity", 0, "Level of verbosity");
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

EquationSystem::EquationSystem(const chi::InputParameters& params)
  : ChiObject(params),
    material_properties_data_(GetOrMakeMaterialPropertiesData(params)),
    verbosity_(params.GetParamValue<int>("verbosity")),
    output_file_base_(
      params.GetParamValue<std::string>("output_filename_base")),

    primary_fields_container_(
      std::move(MakeMultifieldContainer(params.GetParam("fields")))),

    num_local_dofs_(primary_fields_container_->TotalNumLocalDOFs()),
    num_globl_dofs_(primary_fields_container_->TotalNumGlobalDOFs()),
    main_solution_vector_(std::move(MakeSolutionVector())),
    time_integrator_(InitTimeIntegrator(params)),
    num_solution_histories_(
      time_integrator_->NumberOfSolutionHistoriesRequired()),
    num_residual_histories_(
      time_integrator_->NumberOfResidualHistoriesRequired()),
    old_solution_vectors_(std::move(InitSolutionHistory())),
    old_residual_vectors_(std::move(InitResidualHistory())),
    eq_term_scope_(EqTermScope::DOMAIN_TERMS | EqTermScope::BOUNDARY_TERMS),
    time_data_(0.01, 0.0, 1.0)
{
}

// ##################################################################
const chi::MaterialPropertiesData&
EquationSystem::GetOrMakeMaterialPropertiesData(
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
std::unique_ptr<ParallelVector> EquationSystem::MakeSolutionVector()
{
  auto ghost_ids = primary_fields_container_->GetSystemGhostIDs();

  auto new_vec = std::make_unique<GhostedParallelSTLVector>(
    num_local_dofs_, num_globl_dofs_, ghost_ids, Chi::mpi.comm);

  return new_vec;
}

// ##################################################################
std::shared_ptr<TimeIntegrator>
EquationSystem::InitTimeIntegrator(const chi::InputParameters& params)
{
  const auto& params_at_assignment = params.ParametersAtAssignment();

  size_t handle;
  if (params_at_assignment.Has("time_integrator"))
    handle = params.GetParamValue<size_t>("time_integrator");
  else
  {
    chi::InputParameters valid_params =
      ImplicitEulerTimeIntegrator::GetInputParameters();

    // Has no parameters to set

    auto& factory = ChiObjectFactory::GetInstance();
    handle = factory.MakeRegisteredObjectOfType(
      "chi_math::ImplicitEulerTimeIntegrator", valid_params);
  }

  auto time_integrator = Chi::GetStackItemPtrAsType<TimeIntegrator>(
    Chi::object_stack, handle, "EquationSystem::InitTimeIntegrator");

  return time_integrator;
}

// ##################################################################
std::vector<std::unique_ptr<ParallelVector>>
EquationSystem::InitSolutionHistory()
{
  std::vector<std::unique_ptr<ParallelVector>> old_solution_vectors;
  for (size_t t = 0; t < num_solution_histories_; ++t)
  {
    auto vec = main_solution_vector_->MakeNewVector();
    old_solution_vectors.push_back(std::move(vec));
  }

  return old_solution_vectors;
}

// ##################################################################
std::vector<std::unique_ptr<ParallelVector>>
EquationSystem::InitResidualHistory()
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
int64_t EquationSystem::NumLocalDOFs() const { return num_local_dofs_; }
/**Returns the number of local DOFs across all unknowns.*/
int64_t EquationSystem::NumGlobalDOFs() const { return num_globl_dofs_; }

/**Returns a reference to the current solution vector.*/
ParallelVector&
EquationSystem::SolutionVector(TimeID time_id /*=TimeID::T_PLUS_1*/)
{
  if (time_id == TimeID::T_PLUS_1) return *main_solution_vector_;

  const size_t index = static_cast<int>(time_id);
  ChiLogicalErrorIf(index >= old_solution_vectors_.size(),
                    "Solution vector not available at " + TimeIDName(time_id));

  return *old_solution_vectors_[index];
}

/**Returns a reference to the current residual vector.*/
ParallelVector&
EquationSystem::ResidualVector(TimeID time_id /*=TimeID::T_PLUS_1*/)
{
  ChiInvalidArgumentIf(time_id == TimeID::T_PLUS_1,
                       "Residual at T_PLUS_1 not available.");

  const size_t index = static_cast<int>(time_id);
  ChiLogicalErrorIf(index >= old_residual_vectors_.size(),
                    "Residual vector not available at " + TimeIDName(time_id));

  return *old_residual_vectors_[index];
}

/**Returns a reference to the system's time integrator.*/
TimeIntegrator& EquationSystem::GetTimeIntegrator()
{
  ChiLogicalErrorIf(not time_integrator_, "Time integrator object is nullptr.");
  return *time_integrator_;
}

/**Sets the current time data.*/
void EquationSystem::SetTimeData(EquationSystemTimeData time_data)
{
  time_data_ = time_data;
}

/**Returns the current equation terms-scope.*/
EqTermScope EquationSystem::EquationTermsScope() const
{
  return eq_term_scope_;
}

/**Sets the scope of current equations.*/
void EquationSystem::SetEquationTermsScope(EqTermScope eq_term_scope)
{
  eq_term_scope_ = eq_term_scope;
}

bool EquationSystem::QueryTermsActive(EqTermScope query_scope) const
{
  return eq_term_scope_ & query_scope;
}

/**Returns a reference to the system current time data.*/
const EquationSystemTimeData& EquationSystem::GetTimeData() const
{
  return time_data_;
}

std::shared_ptr<MultiFieldContainer>
EquationSystem::MakeMultifieldContainer(const chi::ParameterBlock& params)
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
ParallelMatrixSparsityPattern EquationSystem::BuildMatrixSparsityPattern() const
{
  return primary_fields_container_->BuildMatrixSparsityPattern();
}

/**Updates the fields.*/
void EquationSystem::UpdateFields()
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
void EquationSystem::OutputFields(int time_index)
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

/**Advances the system in time.*/
void EquationSystem::Advance(
  EquationSystemTimeData time_data,
  std::map<TimeID, const ParallelVector*>& latest_std_residuals)
{
  SetTimeData(time_data);

  const auto& sol_vec = main_solution_vector_;

  if (num_solution_histories_ > 0)
  {
    {
      auto& vec = old_solution_vectors_;
      auto start = vec.begin();
      vec.insert(start, sol_vec->MakeCopy());
      vec.pop_back();
    }
  }

  if (num_residual_histories_ > 0)
  {
    for (const auto& [time_id, vec_ptr] : latest_std_residuals)
    {
      const int time_index = static_cast<int>(time_id) + 1;
      if (time_index >= 0 and time_index < num_residual_histories_)
        old_residual_vectors_.at(static_cast<int>(time_id) + 1) =
          vec_ptr->MakeCopy();
    }
  }
}

} // namespace chi_math