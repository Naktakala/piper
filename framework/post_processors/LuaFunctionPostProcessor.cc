#include "LuaFunctionPostProcessor.h"

#include "chi_lua.h"
#include "console/chi_console.h"

#include "event_system/Event.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

namespace chi
{

RegisterChiObject(chi, LuaFunctionPostProcessor);

InputParameters LuaFunctionPostProcessor::GetInputParameters()
{
  InputParameters params = PostProcessor::GetInputParameters();

  params.SetGeneralDescription(
    "PostProcessor that executes a lua function at events.");
  params.SetDocGroup("doc_PostProcessors");

  params.AddRequiredParameter<std::string>("lua_function_name",
                                           "Name of the Lua-function");

  return params;
}

LuaFunctionPostProcessor::LuaFunctionPostProcessor(const InputParameters& params)
: chi::PostProcessor(params, PPType::SCALAR),
    lua_function_name_(params.GetParamValue<std::string>("lua_function_name"))
{}

void LuaFunctionPostProcessor::Execute(const Event& event_context)
{
  const std::string fname = __PRETTY_FUNCTION__;
  lua_State* L = Chi::console.GetConsoleState();
  lua_getglobal(L, lua_function_name_.c_str());

  ChiLogicalErrorIf(not lua_isfunction(L, -1),
                    std::string("Attempted to access lua-function, ") +
                      lua_function_name_ +
                      ", but it seems the function could "
                      "not be retrieved.");

  chi_lua::PushParameterBlock(L, event_context.Parameters(), /*level=*/1);

  std::vector<double> result;
  // 1 arguments, 1 result (table), 0=original error object
  if (lua_pcall(L, 1, 1, 0) == 0)
  {
    //LuaCheckTableValue(fname, L, -1);
    value_ = chi_lua::StackItemToParameterBlock(L, -1);
  }
  else
    throw std::logic_error(fname + " attempted to call lua-function, " +
                           lua_function_name_ + ", but the call failed. " +
                           lua_tostring(L, -1));
  lua_pop(L, lua_gettop(L));
}

} // namespace chi