#include "LiquidPhysics.h"

namespace piper
{

void LiquidPhysics::Advance()
{
  SetTime(Time() + DeltaT());

  for (auto& model : component_models_)
    model->AdvanceNewToOld();
}

}