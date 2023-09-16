#ifndef PIPER_HCTRANSIENTEXECUTOR_H
#define PIPER_HCTRANSIENTEXECUTOR_H

#include "HCSteadyExecutor.h"

namespace hcm
{

class HCTransientExecutor : public HCSteadyExecutor
{
public:
  static chi::InputParameters GetInputParameters();
  explicit HCTransientExecutor(const chi::InputParameters& params);

  void Initialize() override;
  std::unique_ptr<chi_math::NonLinearExecutioner> SetExecutioner() override;
  void Step() override;
  void Advance() override;

protected:
  const chi::ParameterBlock time_integrator_params_;
};

} // namespace hcm

#endif // PIPER_HCTRANSIENTEXECUTOR_H
