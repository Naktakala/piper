#include "../NetworkServer.h"

#include "console/chi_console.h"

#include "chi_runtime.h"
#include "chi_lua.h"

namespace chi::network_server_utils
{

int chiNetworkServerShutdown(lua_State* L);
int chiNetworkServerSynchronize(lua_State* L);

RegisterLuaFunctionAsIs(chiNetworkServerShutdown);
RegisterLuaFunctionAsIs(chiNetworkServerSynchronize);

/**Shuts down the network server pointed to by the handle.
* \param handle int Handle to the network server*/
int chiNetworkServerShutdown(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckIntegerValue(fname, L, 1);

  const size_t handle = lua_tointeger(L, 1);
  auto& server =
    Chi::GetStackItem<chi::NetworkServer>(Chi::object_stack, handle, fname);

  server.ShutdownServer();

  return 0;
}

/***/
int chiNetworkServerSynchronize(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckIntegerValue(fname, L, 1);

  const size_t handle = lua_tointeger(L, 1);
  auto& server =
    Chi::GetStackItem<chi::NetworkServer>(Chi::object_stack, handle, fname);

  server.Synchronize();

  return 0;
}


} // namespace chi::network_server_utils