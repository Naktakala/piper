#ifndef CHITECH_FEMNEUMANNBC_H
#define CHITECH_FEMNEUMANNBC_H

#include "FEMBoundaryCondition.h"

namespace chi_math
{

class FEMNeumannBC : public FEMBoundaryCondition
{
public:
  static chi::InputParameters GetInputParameters();
  explicit FEMNeumannBC(const chi::InputParameters& params);

  double ResidualEntryAtQP() override;

public:
  double flux_value_;
};

}

#endif // CHITECH_FEMNEUMANNBC_H
