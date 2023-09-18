#include "HCSteadyExecutor.h"

#include "heat_conduction/HeatConductionSystem.h"

#include "math/ParallelVector/ghosted_parallel_vector.h"
#include "math/SpatialDiscretization/spatial_discretization.h"

#include "math/Systems/EquationSystemTimeData.h"
#include "math/KernelSystem/FEMKernelSystem.h"
#include "math/NonLinearSolvers/SteadyNonLinearExecutioner.h"
#include "math/NonLinearSolvers/BasicNonLinearSolver.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "ChiObjectFactory.h"

namespace hcm
{

RegisterChiObject(hcm, HCSteadyExecutor);

chi::InputParameters HCSteadyExecutor::GetInputParameters()
{
  chi::InputParameters params = chi_physics::Solver::GetInputParameters();

  params.SetGeneralDescription("Steady state executor for a conduction system");
  params.SetDocGroup("doc_HeatConduction");

  params.ChangeExistingParamToOptional("name", "HCSteadyExecutor");

  params.AddRequiredParameter<size_t>("conduction_system",
                                      "Handle to a conduction system.");

  params.AddOptionalParameter("temperature_initial_value",
                              0.0,
                              "Initial value used for the temperature.");

  params.AddOptionalParameterBlock(
    "solver_params", chi::ParameterBlock{}, "Parameters to pass to solver.");
  params.LinkParameterToBlock("solver_params",
                              "chi_math::NonLinearSolverOptions");

  return params;
}

HCSteadyExecutor::HCSteadyExecutor(const chi::InputParameters& params)
  : chi_physics::Solver(params),
    hc_system_(Chi::GetStackItemPtrAsType<HeatConductionSystem>(
      Chi::object_stack,
      params.GetParamValue<size_t>("conduction_system"),
      __FUNCTION__)),
    T_initial_value_(params.GetParamValue<double>("temperature_initial_value")),
    nl_solver_params_(params.GetParam("solver_params")),
    T_uk_man_({chi_math::Unknown{chi_math::UnknownType::SCALAR}}),
    T_old_(nullptr),
    nl_solver_(nullptr),
    fem_kernel_system_(nullptr)
{
}

void HCSteadyExecutor::Initialize()
{
  const auto& sdm = hc_system_->SDM();
  const size_t num_local_dofs = sdm.GetNumGlobalDOFs(T_uk_man_);
  const size_t num_globl_dofs = sdm.GetNumGlobalDOFs(T_uk_man_);
  const auto ghost_ids = sdm.GetGhostDOFIndices(T_uk_man_);

  T_old_ = std::make_unique<chi_math::GhostedParallelVector>(
    num_local_dofs, num_globl_dofs, ghost_ids, Chi::mpi.comm);
  T_old_->Set(T_initial_value_);
  T_old_->CommunicateGhostEntries();

  fem_kernel_system_ = std::make_shared<chi_math::FEMKernelSystem>(
    sdm,
    T_uk_man_,
    hc_system_->VolumeKernelInputs(),
    hc_system_->BoundaryConditionInputs(),
    chi_math::TimeID::T_PLUS_1);

  nl_executioner_ = SetExecutioner();

  chi::InputParameters nl_params =
    chi_math::NonLinearSolverOptions::GetInputParameters();
  nl_params.AssignParameters(nl_solver_params_);
  nl_solver_ = std::make_unique<chi_math::BasicNonLinearSolver>(
    *nl_executioner_, nl_params);

  auto T_field_function = std::make_shared<chi_physics::FieldFunctionGridBased>(
    hc_system_->GetTemperatureFFName(),
    hc_system_->SDMPtr(),
    T_uk_man_.GetUnknown(0));

  Chi::field_function_stack.push_back(T_field_function);
  this->field_functions_.push_back(T_field_function);

  name_2_ff_map_[hc_system_->GetTemperatureFFName()] = T_field_function;

  nl_solver_->Setup();
}

std::unique_ptr<chi_math::NonLinearExecutioner>
HCSteadyExecutor::SetExecutioner()
{
  auto executioner_in_params =
    chi_math::SteadyNonLinearExecutioner::GetInputParameters();
  auto executioner = std::make_unique<chi_math::SteadyNonLinearExecutioner>(
    executioner_in_params, fem_kernel_system_);

  return executioner;
}

void HCSteadyExecutor::Execute()
{
  nl_solver_->Solve();
  auto& solution = fem_kernel_system_->SolutionVector();

  auto T_ff = name_2_ff_map_.at(hc_system_->GetTemperatureFFName());

  T_ff->UpdateFieldVector(solution.RawValues());
}

} // namespace hcm