#include "IncompressibleLiquidPhysics.h"

#include "piper/utils/CoolPropInterface.h"

#include "chi_log.h"

#include "CoolProp.h"

namespace piper
{
// ###################################################################
IncompressibleLiquidPhysics::StateValMap
IncompressibleLiquidPhysics::EvaluateState(const StateValsList& state_vals_list)
{
  ChiInvalidArgumentIf(state_vals_list.size() < 2,
                       "At least 2 state variables are required.");
  const auto& [v0n, v0v] = state_vals_list[0];
  const auto& [v1n, v1v] = state_vals_list[1];

  std::vector<std::string> vals_wanted = {
    "rho", "e", "T", "p", "h", "s", "k", "Pr", "mu", "Cp"};

  StateValMap outvals;
  for (const auto& chi_name : vals_wanted)
  {
    outvals.insert(std::make_pair(
      chi_name,
      CoolProp::PropsSI(TranslateChiPropertyName2CoolProp(chi_name),
                        TranslateChiPropertyName2CoolProp(v0n),
                        v0v,
                        TranslateChiPropertyName2CoolProp(v1n),
                        v1v,
                        fluid_name_)));
  }

  return outvals;
}
} // namespace piper