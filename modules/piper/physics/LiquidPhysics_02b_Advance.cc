#include "LiquidPhysics.h"

namespace piper
{

void LiquidPhysics::Advance()
{
  time_ = Time() + DeltaT();

  for (auto& model : component_models_)
    model->AdvanceNewToOld();
}

}