#include "LiquidPhysics.h"

#include "CoolProp.h"
#include "utils/chi_utils.h"
#include "piper/Piper.h"
#include "piper/components/HardwareComponent.h"
#include "piper/MeshGenerators/PiperMeshGenerator.h"


#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiObjectFactory.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace piper
{

RegisterChiObject(piper, LiquidPhysics);

// ###################################################################
chi::InputParameters LiquidPhysics::GetInputParameters()
{
  chi::InputParameters params = FluidPhysics::GetInputParameters();

  params.SetGeneralDescription("Physics for single incompressible liquids");
  params.SetDocGroup("Piper");

  //============================================= Fluid name
  // clang-format off
  params.AddRequiredParameter<std::string>("fluid_name",
    "Name of the fluid. All of the incompressible fluids in CoolProp are "
    "supported");
  // clang-format on

  const auto incompressible_list_pure = chi::StringSplit(
    CoolProp::get_global_param_string("incompressible_list_pure"), ",");

  using namespace chi_data_types;
  params.ConstrainParameterRange(
    "fluid_name", AllowableRangeList::New(incompressible_list_pure));

  return params;
}

// ###################################################################
LiquidPhysics::LiquidPhysics(const chi::InputParameters& params)
  : FluidPhysics(params),
    fluid_name_(params.GetParamValue<std::string>("fluid_name"))
{
  Chi::log.Log() << "IncompressibleLiquidPhysics created";
}


} // namespace piper