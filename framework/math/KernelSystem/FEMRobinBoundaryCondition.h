#ifndef PIPER_FEMROBINBOUNDARYCONDITION_H
#define PIPER_FEMROBINBOUNDARYCONDITION_H

#include "math/KernelSystem/FEMBCs/FEMBoundaryCondition.h"

namespace chi_math
{

class FEMRobinBoundaryCondition : public FEMBoundaryCondition
{
public:
  static chi::InputParameters GetInputParameters();
  explicit FEMRobinBoundaryCondition(const chi::InputParameters& params);

  bool IsDirichlet() const override;
};

} // namespace chi_math

#endif // PIPER_FEMROBINBOUNDARYCONDITION_H
