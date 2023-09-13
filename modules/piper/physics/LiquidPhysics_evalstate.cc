#include "LiquidPhysics.h"

#include "piper/utils/CoolPropInterface.h"

#include "chi_log.h"

#include "CoolProp.h"

namespace piper
{

// ###################################################################
LiquidPhysics::StateValMap
LiquidPhysics::EvaluateState(const std::vector<std::string>& vals_wanted,
                             const StateValsList& state_vals_list)
{
  ChiInvalidArgumentIf(state_vals_list.size() < 2,
                       "At least 2 state variables are required.");
  const auto& [v0n, v0v] = state_vals_list[0];
  const auto& [v1n, v1v] = state_vals_list[1];

  std::vector<std::string> propnames;
  for (const auto& chiname : vals_wanted)
    propnames.push_back(TranslateChiPropertyName2CoolProp(chiname));

  const auto v0nn = TranslateChiPropertyName2CoolProp(v0n);
  const auto v1nn = TranslateChiPropertyName2CoolProp(v1n);

  StateValMap outvals;
  const auto vals = CoolProp::PropsSImulti(
    propnames, v0nn, {v0v}, v1nn, {v1v}, "", {fluid_name_}, {1.0});

  size_t k = 0;
  for (const auto& chiname : vals_wanted)
    outvals.insert(std::make_pair(chiname, vals[0][k++]));

  return outvals;
}
} // namespace piper