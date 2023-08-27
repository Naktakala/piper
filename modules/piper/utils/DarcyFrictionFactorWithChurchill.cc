#include "piper/models/ComponentModel.h"

namespace piper
{

double DarcyFrictionFactorWithChurchill(const ComponentModel& model)
{
  const double Re = std::max(model.VarOld("Re"), 10.0);
  const double e = model.Roughness();
  const double D = model.HydraulicDiameter();

  const double a = pow(2.457 * log(pow(7.0 / Re, -0.9) + 0.27 * e / D), 16.0);
  const double b = pow(37530 / Re, 16.0);

  const double fT =
    8.0 * pow(pow(8.0 / Re, 12.0) + pow(a + b, -1.5), 1.0 / 12.0);

  if (Re < 2000.0)
  {
    const double f = 64.0 / Re;

    return std::min(1.0e6, f);
  }
  else if (Re < 4000.0)
  {
    const double c = (Re - 2000.0) / (4000.0 - 2000.0);
    return (1.0 - c) * (64.0 / Re) + c * fT;
  }
  else
    return fT;
}

} // namespace piper