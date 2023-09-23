#ifndef CHITECH_SINKSOURCEFEMKERNEL_H
#define CHITECH_SINKSOURCEFEMKERNEL_H

#include "FEMKernel.h"

namespace chi_math
{

/**Kernel for implementing either a sink or a source.*/
class SinkSourceFEMKernel : public FEMKernel
{
public:
  static chi::InputParameters GetInputParameters();
  explicit SinkSourceFEMKernel(const chi::InputParameters& params);

  double ResidualEntryAtQP() override;
  double JacobianEntryAtQP() override {return 0.0;}

protected:
  double value_;
};

}

#endif // CHITECH_SINKSOURCEFEMKERNEL_H
