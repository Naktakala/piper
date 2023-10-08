#include "LuaEventHook.h"

#include "console/chi_console.h"
#include "chi_lua.h"

#include "event_system/Event.h"

#include "ChiObjectFactory.h"

namespace chi
{

RegisterChiObject(chi, LuaEventHook);

InputParameters LuaEventHook::GetInputParameters()
{
  InputParameters params = EventHook::GetInputParameters();

  params.SetGeneralDescription(
    "An event hook that will call a lua function when it"
    "responds to its subscribed events. The function takes 2 arguments, the "
    "event name, and a table of parameters which holds a multitude of "
    "parameters depending on the event");

  params.AddRequiredParameter<std::string>(
    "lua_function_name", "Text name of the lua function to call");

  return params;
}

LuaEventHook::LuaEventHook(const InputParameters& params)
  : EventHook(params),
    lua_function_name_(params.GetParamValue<std::string>("lua_function_name"))
{
}

void LuaEventHook::Execute(const Event& event)
{
  const std::string fname = __PRETTY_FUNCTION__;
  lua_State* L = Chi::console.GetConsoleState();
  lua_getglobal(L, lua_function_name_.c_str());

  ChiLogicalErrorIf(not lua_isfunction(L, -1),
                    std::string("Attempted to access lua-function, ") +
                      lua_function_name_ +
                      ", but it seems the function could "
                      "not be retrieved.");

  lua_pushstring(L, event.Name().c_str());
  chi_lua::PushParameterBlock(L, event.Parameters(), /*level=*/1);

  // 2 arguments, 0 results, 0=original error object
  if (lua_pcall(L, 2, 0, 0) != 0)
    throw std::logic_error(fname + " attempted to call lua-function \"" +
                           lua_function_name_ + "\" but the call failed. " +
                           lua_tostring(L, -1));
}


} // namespace chi