#include "TimeIntegrator.h"

namespace chi_math
{

chi::InputParameters TimeIntegrator::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  return params;
}

TimeIntegrator::TimeIntegrator(const chi::InputParameters& params)
  : ChiObject(params)
{
}

} // namespace chi_math