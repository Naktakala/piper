#include "NonLinearExecutioner.h"

#include "math/Systems/EquationSystem.h"
#include "math/ParallelMatrix/ParallelMatrix.h"
#include "math/NonLinearSolvers/BasicNonLinearSolver.h"

#include "chi_log.h"

namespace chi_math
{

chi::InputParameters NonLinearExecutioner::GetInputParameters()
{
  chi::InputParameters params = chi_physics::Solver::GetInputParameters();

  params.AddRequiredParameter<size_t>("system", "Handle to an EquationSystem");

  params.AddOptionalParameterBlock(
    "solver_params", chi::ParameterBlock{}, "Parameters to pass to solver.");
  params.LinkParameterToBlock("solver_params",
                              "chi_math::NonLinearSolverOptions");

  params.AddOptionalParameter(
    "print_timing_info", false, "If true, will print timing info.");

  return params;
}

NonLinearExecutioner::NonLinearExecutioner(const chi::InputParameters& params)
  : chi_physics::Solver(params),
    print_timing_info_(params.GetParamValue<bool>("print_timing_info")),
    eq_system_(GetEquationSystem(params.GetParamValue<size_t>("system"))),
    nl_solver_params_(params.GetParam("solver_params"))
{
}

std::shared_ptr<EquationSystem>
NonLinearExecutioner::GetEquationSystem(size_t handle)
{
  return Chi::GetStackItemPtrAsType<EquationSystem>(
    Chi::object_stack, handle, "NonLinearExecutioner::GetEquationSystem");
}

int64_t NonLinearExecutioner::NumLocalDOFs() const
{
  return eq_system_->NumLocalDOFs();
}

int64_t NonLinearExecutioner::NumGlobalDOFs() const
{
  return eq_system_->NumGlobalDOFs();
}

ParallelVector& NonLinearExecutioner::SolutionVector()
{
  return eq_system_->SolutionVector();
}

void NonLinearExecutioner::SetInitialSolution()
{
  eq_system_->SetInitialSolution();
}

ParallelMatrixSparsityPattern
NonLinearExecutioner::BuildMatrixSparsityPattern() const
{
  return eq_system_->BuildMatrixSparsityPattern();
}

void NonLinearExecutioner::SetTimeData(EquationSystemTimeData time_data)
{
  eq_system_->SetTimeData(time_data);
}

void NonLinearExecutioner::SetModeToTimeOnly()
{
  eq_system_->SetEquationTermsScope(EqTermScope::TIME_TERMS);
}
void NonLinearExecutioner::SetModeToNonTimeOnly()
{
  eq_system_->SetEquationTermsScope(EqTermScope::DOMAIN_TERMS |
                                    EqTermScope::BOUNDARY_TERMS);
}
void NonLinearExecutioner::SetModeToTimeAndNonTime()
{
  eq_system_->SetEquationTermsScope(EqTermScope::TIME_TERMS |
                                    EqTermScope::DOMAIN_TERMS |
                                    EqTermScope::BOUNDARY_TERMS);
}

bool NonLinearExecutioner::TimeIDListHasID(const std::vector<TimeID>& time_ids,
                                           TimeID id)
{
  return std::find(time_ids.begin(), time_ids.end(), id) != time_ids.end();
}

void NonLinearExecutioner::Initialize()
{
  auto nl_params = chi_math::NonLinearSolverOptions::GetInputParameters();
  nl_params.AssignParameters(nl_solver_params_);
  nl_solver_ =
    std::make_unique<chi_math::BasicNonLinearSolver>(*this, nl_params);

  nl_solver_->Setup();
}

void NonLinearExecutioner::PrintTimingInfo() const
{
  if (not print_timing_info_) return;

  const size_t residual_ev_tag =
    Chi::log.GetExistingRepeatingEventTag("ComputeResidual");
  const size_t jacobian_ev_tag =
    Chi::log.GetExistingRepeatingEventTag("ComputeJacobian");

  std::stringstream outstr;
  outstr << "ComputeResidual: ";
  outstr << Chi::log.ProcessEvent(residual_ev_tag,
                                  chi::ChiLog::EventOperation::TOTAL_DURATION) /
              1.0e6;
  outstr << " ("
         << Chi::log.ProcessEvent(
              residual_ev_tag,
              chi::ChiLog::EventOperation::NUMBER_OF_OCCURRENCES)
         << ")";

  outstr << " avg="
         << Chi::log.ProcessEvent(
              residual_ev_tag, chi::ChiLog::EventOperation::AVERAGE_DURATION);
  outstr << "\n";

  outstr << "ComputeJacobian: ";
  outstr << Chi::log.ProcessEvent(jacobian_ev_tag,
                                  chi::ChiLog::EventOperation::TOTAL_DURATION) /
              1.0e6;
  outstr << " ("
         << Chi::log.ProcessEvent(
              jacobian_ev_tag,
              chi::ChiLog::EventOperation::NUMBER_OF_OCCURRENCES)
         << ")";
  outstr << " avg="
         << Chi::log.ProcessEvent(
              jacobian_ev_tag, chi::ChiLog::EventOperation::AVERAGE_DURATION);
  outstr << "\n";

  for (size_t k = 0; k < 10; ++k)
  {
    const size_t t_tag =
      Chi::log.GetExistingRepeatingEventTag("Timing_" + std::to_string(k));

    outstr << "Timing_" << k << " " << t_tag << " average "
           << Chi::log.ProcessEvent(
                t_tag, chi::ChiLog::EventOperation::TOTAL_DURATION) /
                1.0e3;
    outstr << "\n";
  }

  Chi::log.Log() << outstr.str();
}

} // namespace chi_math