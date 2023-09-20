#include "ImplicitEulerTimeIntegrator.h"

#include "math/ParallelVector/ParallelVector.h"
#include "math/Systems/EquationSystem.h"

#include "ChiObjectFactory.h"

namespace chi_math
{

RegisterChiObject(chi_math, ImplicitEulerTimeIntegrator);

chi::InputParameters ImplicitEulerTimeIntegrator::GetInputParameters()
{
  chi::InputParameters params = TimeIntegrator::GetInputParameters();

  params.SetGeneralDescription("Basic Implicit Euler time integration");
  params.SetDocGroup("doc_TimeIntegrators");

  return params;
}

ImplicitEulerTimeIntegrator::ImplicitEulerTimeIntegrator(
  const chi::InputParameters& params)
  : TimeIntegrator(params)
{
}

std::vector<TimeID> ImplicitEulerTimeIntegrator::GetTimeIDsNeeded() const
{
  return {TimeID::T_PLUS_1};
}

size_t ImplicitEulerTimeIntegrator::NumberOfSolutionHistoriesRequired() const
{
  return 1;
}

size_t ImplicitEulerTimeIntegrator::NumberOfResidualHistoriesRequired() const
{
  return 0;
}

void ImplicitEulerTimeIntegrator::ComputeResidual(
  ParallelVector& r,
  const ParallelVector& time_residual,
  const std::vector<const ParallelVector*>& std_residuals)
{
  ChiLogicalErrorIf(std_residuals.size() != 1, "std_residuals too short.");

  const auto& r_t_tp1 = time_residual;
  const auto& r_x_tp1 = *std_residuals[static_cast<int>(TimeID::T_PLUS_1) + 1];

  r += r_t_tp1;
  r += r_x_tp1;
}

void ImplicitEulerTimeIntegrator::ComputeJacobian(
  const ParallelVector& x,
  ParallelMatrix& J,
  EquationSystem& equation_system)
{
  equation_system.ComputeJacobian(x, J);
}


} // namespace chi_math