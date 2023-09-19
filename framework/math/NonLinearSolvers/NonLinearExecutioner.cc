#include "NonLinearExecutioner.h"

#include "math/Systems/EquationSystem.h"
#include "math/ParallelMatrix/ParallelMatrix.h"

namespace chi_math
{

chi::InputParameters NonLinearExecutioner::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  return params;
}

NonLinearExecutioner::NonLinearExecutioner(
  const chi::InputParameters& params,
  std::shared_ptr<EquationSystem> equation_system)
  : ChiObject(params), eq_system_(std::move(equation_system))
{
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

void NonLinearExecutioner::Advance(EquationSystemTimeData time_data) {}

} // namespace chi_math