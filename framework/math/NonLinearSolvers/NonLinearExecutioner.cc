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

  params.AddOptionalParameter(
    "print_nl_residual",
    true,
    "If false, no non-linear residual will be printed.");
  params.AddOptionalParameter(
    "print_l_residual", true, "If false, no linear residual will be printed.");
  params.AddOptionalParameter(
    "print_header",
    true,
    "If false, steady executioner will not print header messages like "
    "\"Executing...\", whilst"
    "transient executioners will not print time-info at each timestep.");

  params.AddOptionalParameter(
    "print_footer",
    true,
    "Disables convergence reason printing and stuff below it.");

  return params;
}

NonLinearExecutioner::NonLinearExecutioner(const chi::InputParameters& params)
  : chi_physics::Solver(params),
    t_tag_residual_(
      Chi::log.GetRepeatingEventTag(TextName() + "ComputeResidual")),
    t_tag_jacobian_(
      Chi::log.GetRepeatingEventTag(TextName() + "ComputeJacobian")),
    t_tag_solve_(Chi::log.GetRepeatingEventTag(TextName() + "Solve")),
    print_timing_info_(params.GetParamValue<bool>("print_timing_info")),
    print_nl_residual_(params.GetParamValue<bool>("print_nl_residual")),
    print_l_residual_(params.GetParamValue<bool>("print_l_residual")),
    print_headers_(params.GetParamValue<bool>("print_header")),
    print_footers_(params.GetParamValue<bool>("print_footer")),
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

  constexpr auto DUR = chi::ChiLog::EventOperation::TOTAL_DURATION;
  constexpr auto NUM = chi::ChiLog::EventOperation::NUMBER_OF_OCCURRENCES;
  constexpr auto AVG = chi::ChiLog::EventOperation::AVERAGE_DURATION;

  std::stringstream outstr;
  outstr << "ComputeResidual: ";
  outstr << Chi::log.ProcessEvent(t_tag_residual_, DUR) / 1.0e6;
  outstr << " (" << Chi::log.ProcessEvent(t_tag_residual_, NUM) << ")";
  outstr << " avg=" << Chi::log.ProcessEvent(t_tag_residual_, AVG);
  outstr << "\n";

  outstr << "ComputeJacobian: ";
  outstr << Chi::log.ProcessEvent(t_tag_jacobian_, DUR) / 1.0e6;
  outstr << " (" << Chi::log.ProcessEvent(t_tag_jacobian_, NUM) << ")";
  outstr << " avg=" << Chi::log.ProcessEvent(t_tag_jacobian_, AVG);
  outstr << "\n";

  outstr << "Solve: ";
  outstr << Chi::log.ProcessEvent(t_tag_solve_, DUR) / 1.0e6;
  outstr << " (" << Chi::log.ProcessEvent(t_tag_solve_, NUM) << ")";
  outstr << " avg=" << Chi::log.ProcessEvent(t_tag_solve_, AVG);
  outstr << "\n";

  Chi::log.Log() << "\n" + TextName() + " timing info:\n" + outstr.str();
}

} // namespace chi_math