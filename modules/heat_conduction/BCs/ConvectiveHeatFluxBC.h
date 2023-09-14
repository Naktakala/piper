#ifndef PIPER_CONVECTIVEHEATFLUXBC_H
#define PIPER_CONVECTIVEHEATFLUXBC_H

#include "math/KernelSystem/FEMBCs/FEMBoundaryCondition.h"

namespace hcm
{

class ConvectiveHeatFluxBC : public chi_math::FEMBoundaryCondition
{
public:
  static chi::InputParameters GetInputParameters();
  explicit ConvectiveHeatFluxBC(const chi::InputParameters& params);

  double ResidualEntryAtQP() override;
  double JacobianEntryAtQP() override;

protected:
  double T_bulk_;
  double h_;
};

}

#endif // PIPER_CONVECTIVEHEATFLUXBC_H
