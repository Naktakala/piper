#ifndef CHITECH_LUAEVENTHOOK_H
#define CHITECH_LUAEVENTHOOK_H

#include "EventHook.h"

namespace chi
{

class LuaEventHook : public EventHook
{
public:
  static InputParameters GetInputParameters();
  explicit LuaEventHook(const InputParameters& params);

  void Execute(const Event& event) override;

private:
  const std::string lua_function_name_;
};

} // namespace chi

#endif // CHITECH_LUAEVENTHOOK_H
