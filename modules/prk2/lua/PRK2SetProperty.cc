#include "../PRK2System.h"

#include "console/chi_console.h"
#include "lua/chi_lua.h"

namespace prk2::lua_utils
{

int chiPRK2SetProperty(lua_State* L);

RegisterLuaFunctionAsIs(chiPRK2SetProperty);

int chiPRK2SetProperty(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckIntegerValue(fname, L, 1);

  const size_t handle = lua_tointeger(L, 1);

  auto& prk =
    Chi::GetStackItem<prk2::PRK2System>(Chi::object_stack, handle, fname);

  LuaCheckTableValue(fname, L, 2);

  auto params = chi_lua::TableParserAsParameterBlock::ParseTable(L, 2);

  prk.SetProperties(params);

  return 0;
}

} // namespace prk2::lua_utils