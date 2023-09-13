#ifndef PIPER_SINKSOURCEFEMKERNEL_H
#define PIPER_SINKSOURCEFEMKERNEL_H

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

protected:
  double value_;
};

}

#endif // PIPER_SINKSOURCEFEMKERNEL_H
