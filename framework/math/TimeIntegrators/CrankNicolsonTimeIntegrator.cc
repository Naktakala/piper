#include "CrankNicolsonTimeIntegrator.h"

#include "math/ParallelVector/ParallelVector.h"
#include "math/Systems/EquationSystem.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

namespace chi_math
{

RegisterChiObject(chi_math, CrankNicolsonTimeIntegrator);

chi::InputParameters CrankNicolsonTimeIntegrator::GetInputParameters()
{
  chi::InputParameters params = TimeIntegrator::GetInputParameters();

  params.SetGeneralDescription("Basic Crank-Nicolson time integration");
  params.SetDocGroup("doc_TimeIntegrators");

  return params;
}

CrankNicolsonTimeIntegrator::CrankNicolsonTimeIntegrator(
  const chi::InputParameters& params)
  : TimeIntegrator(params)
{
}

std::vector<TimeID> CrankNicolsonTimeIntegrator::GetResidualTimeIDsNeeded() const
{
  return {TimeID::T_PLUS_1, TimeID::T};
}

size_t CrankNicolsonTimeIntegrator::NumberOfSolutionHistoriesRequired() const
{
  return 1;
}

size_t CrankNicolsonTimeIntegrator::NumberOfResidualHistoriesRequired() const
{
  return 1;
}

double CrankNicolsonTimeIntegrator::GetTimeCoefficient(double dt) const
{
  return 2.0/dt;
}

void CrankNicolsonTimeIntegrator::ComputeResidual(
  ParallelVector& r,
  const ParallelVector& time_residual,
  const std::map<TimeID, const ParallelVector*>& std_residuals)
{
  ChiLogicalErrorIf(std_residuals.size() != 2, "std_residuals too short.");

  const auto& r_t_tp1 = time_residual;
  const auto& r_x_tp1 = *std_residuals.at(TimeID::T_PLUS_1);
  const auto& r_x_t = *std_residuals.at(TimeID::T);

  r += r_t_tp1;
  r.PlusAY(r_x_tp1, 0.5);
  r.PlusAY(r_x_t, 0.5);
}

void CrankNicolsonTimeIntegrator::ComputeJacobian(
  const ParallelVector& x,
  ParallelMatrix& J, EquationSystem& equation_system)
{
  equation_system.ComputeJacobian(x, J);
}

} // namespace chi_math