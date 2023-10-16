#include "LiquidPhysics.h"

#include "piper/Piper.h"
#include "piper/components/HardwareComponent.h"
#include "piper/MeshGenerators/PiperMeshGenerator.h"
#include "piper/models/BoundaryLiquidModel.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace piper
{

// ###################################################################
std::vector<std::string>
LiquidPhysics::MakeVariableNamesList(ComponentCategory hw_comp_category)
{
  std::vector<std::string> variable_names;
// These switches will let the compiler tell us when we miss cases in
// the enum
#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wswitch-enum"
  switch (hw_comp_category)
  {
    case ComponentCategory::BoundaryLike:
    case ComponentCategory::Volumetric:
    {
      variable_names = {
        "rho", "e", "T", "p", "k", "mu", "u", "Pr", "Re", "hcoeff", "beta"};
      break;
    }
    case ComponentCategory::JunctionLike:
    {
      variable_names = {"u"};
      break;
    }
  } // switch on hw_comp_category
#pragma GCC diagnostic pop

  return variable_names;
}

// ###################################################################
ComponentLiquidModel&
LiquidPhysics::GetComponentLiquidModel(size_t component_id)
{
  ComponentModel* model_ptr = component_models_.at(component_id).get();
  auto liquid_model_ptr = dynamic_cast<ComponentLiquidModel*>(model_ptr);

  ChiLogicalErrorIf(not liquid_model_ptr, "Casting error");
  return *liquid_model_ptr;
}

} // namespace piper
