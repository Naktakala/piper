#include "FEMRobinBoundaryCondition.h"

namespace chi_math
{

chi::InputParameters FEMRobinBoundaryCondition::GetInputParameters()
{
  chi::InputParameters params = FEMBoundaryCondition::GetInputParameters();

  return params;
}

FEMRobinBoundaryCondition::FEMRobinBoundaryCondition(
  const chi::InputParameters& params)
  : FEMBoundaryCondition(params)
{
}

bool FEMRobinBoundaryCondition::IsDirichlet() const { return false; }

} // namespace chi_math