#include "EquationSystem.h"

#include "math/ParallelVector/GhostedParallelSTLVector.h"
#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/spatial_discretization.h"
#include "math/TimeIntegrators/ImplicitEulerTimeIntegrator.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

#include <numeric>



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

  return params;
}

EquationSystem::EquationSystem(const chi::InputParameters& params)
  : ChiObject(params),
    verbosity_(params.GetParamValue<int>("verbosity")),
    field_block_info_(
      std::move(MakeFieldBlockInfo(BuildFieldList(params.GetParam("fields"))))),
    num_local_dofs_(DetermineNumLocalDofs(field_block_info_)),
    num_globl_dofs_(DetermineNumGlobalDofs(field_block_info_)),
    main_solution_vector_(std::move(MakeSolutionVector())),
    time_integrator_(InitTimeIntegrator(params)),
    num_solution_histories_(
      time_integrator_->NumberOfSolutionHistoriesRequired()),
    num_residual_histories_(
      time_integrator_->NumberOfResidualHistoriesRequired()),
    old_solution_vectors_(std::move(InitSolutionHistory())),
    old_residual_vectors_(std::move(InitResidualHistory())),
    eq_term_scope_(EqTermScope::DOMAIN_TERMS | EqTermScope::BOUNDARY_TERMS)
{
}

// ##################################################################
/**Makes a list of grid-base field functions given either an array of
 * warehouse handles or field-function names.*/
FieldList EquationSystem::BuildFieldList(const chi::ParameterBlock& param_array)
{
  ChiInvalidArgumentIf(
    param_array.Type() != chi::ParameterBlockType::ARRAY,
    "\"variables\" parameter for EquationSystem must be of type ARRAY.");

  ChiInvalidArgumentIf(
    param_array.NumParameters() == 0,
    "\"variables\" parameter for EquationSystem must not be empty.");

  std::vector<std::string> field_names_processed;
  auto CheckDuplicateName = [&field_names_processed](const std::string& ff_name)
  {
    return std::any_of(field_names_processed.begin(),
                       field_names_processed.end(),
                       [&ff_name](const std::string& name)
                       { return name == ff_name; });
  };

  const auto values_type = param_array.GetParam(0).Type();
  if (values_type == chi::ParameterBlockType::INTEGER)
  {
    FieldList field_list;
    for (const auto& index_param : param_array)
    {
      const size_t handle = index_param.GetValue<size_t>();
      auto field_base_ptr =
        Chi::GetStackItemPtr(Chi::field_function_stack, handle, __FUNCTION__);

      auto field_ptr =
        std::dynamic_pointer_cast<chi_physics::FieldFunctionGridBased>(
          field_base_ptr);

      ChiLogicalErrorIf(not field_ptr,
                        "Field function \"" + field_base_ptr->TextName() +
                          "\" is not of required type FieldFunctionGridBased.");

      ChiInvalidArgumentIf(CheckDuplicateName(field_ptr->TextName()),
                           "Duplicate field \"" + field_ptr->TextName() +
                             "\" encountered.");

      field_list.push_back(field_ptr);
      field_names_processed.push_back(field_ptr->TextName());
    }

    return field_list;
  }
  else if (values_type == chi::ParameterBlockType::STRING)
  {
    FieldList field_list;
    for (const auto& index_param : param_array)
    {
      const std::string ff_name = index_param.GetValue<std::string>();

      std::shared_ptr<chi_physics::FieldFunctionGridBased> field_ptr = nullptr;

      for (auto& ff : Chi::field_function_stack)
        if (ff->TextName() == ff_name)
        {
          field_ptr =
            std::dynamic_pointer_cast<chi_physics::FieldFunctionGridBased>(ff);
          ChiLogicalErrorIf(
            not field_ptr,
            "Field function \"" + ff_name +
              "\" is not of required type FieldFunctionGridBased.");
        }

      ChiInvalidArgumentIf(not field_ptr,
                           "Field function \"" + ff_name +
                             "\" was not found in field-function warehouse.");

      ChiInvalidArgumentIf(CheckDuplicateName(field_ptr->TextName()),
                           "Duplicate field \"" + field_ptr->TextName() +
                             "\" encountered.");

      field_list.push_back(field_ptr);
      field_names_processed.push_back(field_ptr->TextName());
    }

    return field_list;
  }
  else
    ChiInvalidArgument(
      "The elements of the \"variables\" parameter for "
      "EquationSystem can only be INTEGER or STRING. This is either a handle "
      "to a FieldFunctionGridBased or a name to an existing one.");
}

// ##################################################################
int64_t EquationSystem::DetermineNumLocalDofs(
  const std::vector<FieldBlockInfo>& field_block_info)
{
  size_t num_dofs = 0;
  for (const auto& field_block : field_block_info)
    num_dofs += field_block.num_local_dofs_;

  return static_cast<int64_t>(num_dofs);
}
// ##################################################################
int64_t EquationSystem::DetermineNumGlobalDofs(
  const std::vector<FieldBlockInfo>& field_block_info)
{
  size_t num_dofs = 0;
  for (const auto& field_block : field_block_info)
    num_dofs += field_block.num_global_dofs_;

  return static_cast<int64_t>(num_dofs);
}

// ##################################################################
std::unique_ptr<ParallelVector> EquationSystem::MakeSolutionVector()
{
  std::vector<int64_t> ghost_ids;
  for (const auto& info : field_block_info_)
    for (const int64_t block_ghost_id : info.ghost_ids_)
      ghost_ids.push_back(MapBlockGlobalIDToSystem(info, block_ghost_id));

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

/**Updates the fields.*/
void EquationSystem::UpdateFields()
{
  const auto& x = *main_solution_vector_;
  for (auto& field_info : field_block_info_)
  {
    std::vector<double> local_solution(field_info.num_local_dofs_, 0.0);
    for (int64_t i = 0; i < field_info.num_local_dofs_; ++i)
      local_solution[i] = x[field_info.local_offset_ + i];

    field_info.field_->UpdateFieldVector(local_solution);
  }
}

/**Advances the system in time.*/
void EquationSystem::Advance(EquationSystemTimeData time_data,
                             const ParallelVector& latest_std_residual)
{
  SetTimeData(time_data);

  const auto& sol_vec = main_solution_vector_;
  const auto& res_vec = latest_std_residual;

  if (num_solution_histories_ > 0)
  {
    {
      auto& vec = old_solution_vectors_;
      auto start = vec.begin();
      vec.insert(start, sol_vec->MakeCopy());
      vec.pop_back();
    }
    {
      auto& vec = old_residual_vectors_;
      auto start = vec.begin();
      vec.insert(start, res_vec.MakeCopy());
      vec.pop_back();
    }
  }
}

} // namespace chi_math