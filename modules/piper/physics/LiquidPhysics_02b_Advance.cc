#include "LiquidPhysics.h"

#include "physics/TimeSteppers/TimeStepper.h"

namespace piper
{

void LiquidPhysics::Advance()
{
  timestepper_->Advance();

  for (auto& model : component_models_)
    model->AdvanceNewToOld();

  UpdateFieldFunctions();
}

}