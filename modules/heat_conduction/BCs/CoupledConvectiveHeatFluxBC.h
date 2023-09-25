#ifndef PIPER_COUPLEDCONVECTIVEHEATFLUXBC_H
#define PIPER_COUPLEDCONVECTIVEHEATFLUXBC_H

#include "math/KernelSystem/FEMBCs/FEMBoundaryCondition.h"

namespace hcm
{

class CoupledConvectiveHeatFluxBC : public chi_math::FEMBoundaryCondition
{
public:
  static chi::InputParameters GetInputParameters();
  explicit CoupledConvectiveHeatFluxBC(const chi::InputParameters& params);

  double ResidualEntryAtQP() override;
  double JacobianEntryAtQP() override;

protected:
  const chi_math::FEMCoupledField& T_bulk_;
  const double h_;
};

}

#endif // PIPER_COUPLEDCONVECTIVEHEATFLUXBC_H
