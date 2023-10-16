#include "CoolPropInterface.h"

#include "CoolProp.h"

#include "chi_log.h"

namespace piper
{

std::string TranslateChiPropertyName2CoolProp(const std::string& chi_name)
{
  // clang-format off
  if (chi_name == "rho") return "Dmass";
  else if (chi_name == "e"  ) return "Umass";
  else if (chi_name == "T"  ) return "T";
  else if (chi_name == "p"  ) return "P";
  else if (chi_name == "h"  ) return "Hmass";
  else if (chi_name == "s"  ) return "Smass";
  else if (chi_name == "k"  ) return "conductivity";
  else if (chi_name == "Pr" ) return "Prandtl";
  else if (chi_name == "mu" ) return "viscosity";
  else if (chi_name == "Cp" ) return "Cp0mass";
  else if (chi_name == "beta") return "isobaric_expansion_coefficient";
  else
    ChiLogicalError("Translation of \"" + chi_name +"\" not defined");
  // clang-format on
}

double PropSI(const std::string& output_name,
              const std::string& v0n,
              double v0v,
              const std::string& v1n,
              double v1v,
              const std::string& fluid_name)
{
  return CoolProp::PropsSI(TranslateChiPropertyName2CoolProp(output_name),
                           TranslateChiPropertyName2CoolProp(v0n), v0v,
                           TranslateChiPropertyName2CoolProp(v1n), v1v,
                           fluid_name);
}

} // namespace piper