#ifndef CHITECH_LUAFUNCTIONPOSTPROCESSOR_H
#define CHITECH_LUAFUNCTIONPOSTPROCESSOR_H

#include "post_processors/PostProcessor.h"

namespace chi
{

class LuaFunctionPostProcessor : public PostProcessor
{
public:
  static InputParameters GetInputParameters();
  explicit LuaFunctionPostProcessor(const InputParameters& params);

  void Execute(const Event& event_context) override;

protected:
  const std::string lua_function_name_;
};

}

#endif // CHITECH_LUAFUNCTIONPOSTPROCESSOR_H
