#include "HCTransientExecutor.h"

#include "heat_conduction/HeatConductionSystem.h"

#include "math/ParallelVector/ghosted_parallel_vector.h"
#include "math/SpatialDiscretization/spatial_discretization.h"

#include "math/KernelSystem/FEMKernelSystem.h"
#include "math/NonLinearSolvers/TransientNonLinearExecutioner.h"
#include "math/NonLinearSolvers/BasicNonLinearSolver.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "physics/TimeStepControllers/TimeStepController.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"
#include "logging/stringstream_color.h"

namespace hcm
{

RegisterChiObject(hcm, HCTransientExecutor);

chi::InputParameters HCTransientExecutor::GetInputParameters()
{
  chi::InputParameters params = HCSteadyExecutor::GetInputParameters();
  params += chi_math::TransientNonLinearExecutioner::GetInputParameters();

  params.SetGeneralDescription("Transient executor for a conduction system");
  params.SetDocGroup("doc_HeatConduction");

  params.ChangeExistingParamToOptional("name", "HCTransientExecutor");

  return params;
}

HCTransientExecutor::HCTransientExecutor(const chi::InputParameters& params)
  : HCSteadyExecutor(params),
    time_integrator_params_(params.GetParam("time_integrator"))
{
}

void HCTransientExecutor::Initialize()
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
    hc_system_->VolumeKernels(),
    hc_system_->BoundaryConditions(),
    chi_math::TimeID::T);

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
HCTransientExecutor::SetExecutioner()
{
  auto executioner_in_params =
    chi_math::TransientNonLinearExecutioner::GetInputParameters();
  executioner_in_params.SetErrorOriginScope("TransientNonLinearExecutioner");
  chi::ParameterBlock params;
  params.AddParameter(time_integrator_params_);

  executioner_in_params.AssignParameters(params);

  auto executioner = std::make_unique<chi_math::TransientNonLinearExecutioner>(
    executioner_in_params, fem_kernel_system_);

  return executioner;
}

void HCTransientExecutor::Step()
{
  Chi::log.Log() << time_step_controller_->StringTimeInfo();

  nl_executioner_->SetTimeData({time_step_controller_->GetTimeStepSize(),
                                time_step_controller_->Time()});

  nl_solver_->Solve();

  Chi::log.Log() << nl_solver_->GetConvergedReasonString();
}

void HCTransientExecutor::Advance()
{
  const bool last_solve_converged = nl_solver_->IsConverged();
  if (last_solve_converged)
  {
    time_step_controller_->Advance();
    nl_executioner_->Advance({time_step_controller_->GetTimeStepSize(),
                              time_step_controller_->Time()});

    auto& solution = fem_kernel_system_->SolutionVector();
    auto T_ff = name_2_ff_map_.at(hc_system_->GetTemperatureFFName());
    T_ff->UpdateFieldVector(solution.RawValues());
  }
  else
  {
    if (not time_step_controller_->Adapt(chi_physics::TimeStepStatus::FAILURE))
    {
      Chi::log.Log0Error() << "Solver failed: "
                           << time_step_controller_->StringTimeInfo();
      Chi::Exit(EXIT_FAILURE);
    }


  }
}

} // namespace hcm