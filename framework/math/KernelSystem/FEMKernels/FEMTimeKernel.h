#ifndef CHITECH_FEMTIMEKERNEL_H
#define CHITECH_FEMTIMEKERNEL_H

#include "FEMKernel.h"

namespace chi_math
{

class FEMTimeKernel : public FEMKernel
{
public:
  static chi::InputParameters GetInputParameters();
  explicit FEMTimeKernel(const chi::InputParameters& params);

  bool IsTimeKernel() const override;

protected:
  const size_t TIME_T = 0;
  const size_t TIME_T_MINUS_1 = 1;
  const size_t TIME_T_MINUS_2 = 2;
  const size_t TIME_T_MINUS_3 = 3;
};

} // namespace chi_math

#endif // CHITECH_FEMTIMEKERNEL_H
