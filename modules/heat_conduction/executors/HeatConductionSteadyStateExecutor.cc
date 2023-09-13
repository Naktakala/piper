#include "HeatConductionSteadyStateExecutor.h"

#include "heat_conduction/HeatConductionSystem.h"

#include "math/ParallelVector/ghosted_parallel_vector.h"
#include "math/SpatialDiscretization/spatial_discretization.h"
#include "math/KernelSystem/FEMKernelSystem.h"
#include "math/KernelSystem/KernelBasedNonLinearSolver.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "ChiObjectFactory.h"

namespace hcm
{

RegisterChiObject(hcm, HeatConductionSteadyStateExecutor);

chi::InputParameters HeatConductionSteadyStateExecutor::GetInputParameters()
{
  chi::InputParameters params = chi_physics::Solver::GetInputParameters();

  params.SetGeneralDescription("Steady state executor for a conduction system");
  params.SetDocGroup("doc_HeatConduction");

  params.ChangeExistingParamToOptional("name",
                                       "HeatConductionSteadyStateExecutor");

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

HeatConductionSteadyStateExecutor::HeatConductionSteadyStateExecutor(
  const chi::InputParameters& params)
  : chi_physics::Solver(params),
    hc_system_(Chi::GetStackItemPtrAsType<HeatConductionSystem>(
      Chi::object_stack,
      params.GetParamValue<size_t>("conduction_system"),
      __FUNCTION__)),
    T_initial_value_(params.GetParamValue<double>("temperature_initial_value")),
    solver_params_(params.GetParam("solver_params")),
    T_uk_man_({chi_math::Unknown{chi_math::UnknownType::SCALAR}}),
    T_old_(nullptr),
    nl_solver_(nullptr),
    fem_kernel_system_(nullptr)
{
}

void HeatConductionSteadyStateExecutor::Initialize()
{
  const auto& sdm = hc_system_->SDM();
  num_local_dofs_ = sdm.GetNumGlobalDOFs(T_uk_man_);
  num_globl_dofs_ = sdm.GetNumGlobalDOFs(T_uk_man_);
  const auto ghost_ids = sdm.GetGhostDOFIndices(T_uk_man_);

  T_old_ = std::make_unique<chi_math::GhostedParallelVector>(
    num_local_dofs_, num_globl_dofs_, ghost_ids, Chi::mpi.comm);
  T_old_->Set(T_initial_value_);
  T_old_->CommunicateGhostEntries();

  fem_kernel_system_ = std::make_unique<chi_math::FEMKernelSystem>(
    sdm,
    T_uk_man_,
    hc_system_->VolumeKernels(),
    hc_system_->BoundaryConditions());

  chi::InputParameters nl_params =
    chi_math::NonLinearSolverOptions::GetInputParameters();
  nl_params.AssignParameters(solver_params_);
  nl_solver_ = std::make_unique<chi_math::KernelBasedNonLinearSolver>(
    *fem_kernel_system_, nl_params);

  auto T_field_function = std::make_shared<chi_physics::FieldFunctionGridBased>(
    hc_system_->GetTemperatureFFName(),
    hc_system_->SDMPtr(),
    T_uk_man_.GetUnknown(0));

  Chi::field_function_stack.push_back(T_field_function);
  this->field_functions_.push_back(T_field_function);

  name_2_ff_map_[hc_system_->GetTemperatureFFName()] = T_field_function;

  nl_solver_->Setup();
}

void HeatConductionSteadyStateExecutor::Execute()
{
  nl_solver_->Solve();
  auto& solution = fem_kernel_system_->SolutionVector();

  auto T_ff = name_2_ff_map_.at(hc_system_->GetTemperatureFFName());

  T_ff->UpdateFieldVector(solution.RawValues());
}

} // namespace hcm