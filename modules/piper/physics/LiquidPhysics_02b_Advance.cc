#include "LiquidPhysics.h"

namespace piper
{

void LiquidPhysics::Advance()
{
  //time_ = Time() + TimeStepSize();

  for (auto& model : component_models_)
    model->AdvanceNewToOld();
}

}