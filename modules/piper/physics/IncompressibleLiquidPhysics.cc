#include "IncompressibleLiquidPhysics.h"

#include "CoolProp.h"
#include "utils/chi_utils.h"
#include "piper/Piper.h"
#include "piper/components/HardwareComponent.h"
#include "piper/MeshGenerators/PiperMeshGenerator.h"
#include "piper/models/ComponentModel.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiObjectFactory.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace piper
{

RegisterChiObject(piper, IncompressibleLiquidPhysics);

// ###################################################################
chi::InputParameters IncompressibleLiquidPhysics::GetInputParameters()
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

  const auto& incompressible_list_pure = chi::StringSplit(
    CoolProp::get_global_param_string("incompressible_list_pure"), ",");

  using namespace chi_data_types;
  params.ConstrainParameterRange(
    "fluid_name", AllowableRangeList::New(incompressible_list_pure));

  return params;
}

// ###################################################################
IncompressibleLiquidPhysics::IncompressibleLiquidPhysics(
  const chi::InputParameters& params)
  : FluidPhysics(params),
    fluid_name_(params.GetParamValue<std::string>("fluid_name"))
{
  Chi::log.Log() << "IncompressibleLiquidPhysics created";
}

// ###################################################################
void IncompressibleLiquidPhysics::InitializeUnknowns()
{
  const auto& grid = Grid();
  const auto& pipe_system = PipeSystem();

  auto& hw_components = pipe_system.HardwareComponents();
  const auto& cell_map = mesh_generator_->GetVolumeComponent2CellGIDMap();

  //================================================== Instantiate models
  for (auto& hw_component_ptr : hw_components)
  {
    const auto& hw_component = *hw_component_ptr;
    const auto hw_comp_category = hw_component.Category();

    const chi_mesh::Cell* cell_ptr = nullptr;
    auto cell_iter = cell_map.find(hw_component.Name());
    if (cell_iter != cell_map.end())
    {
      const uint64_t cell_gid = (*cell_iter).second;
      cell_ptr = &(grid.cells[cell_gid]);
    }

    const auto variable_names = MakeVariableNamesList(hw_comp_category);

    auto component_model = std::make_unique<ComponentModel>(
      component_models_, hw_component, cell_ptr, variable_names);

    component_models_.push_back(std::move(component_model));
  }

  //================================================== Execute the initializer
  initializer_param_block_.RequireParameter("type");
  const std::string initializer_type =
    initializer_param_block_.GetParamValue<std::string>("type");

  Chi::mpi.Barrier();
  Chi::log.Log() << "Executing initializer";
  Chi::mpi.Barrier();
  if (initializer_type == "StaticGravity")
  {
    StaticGravityInitializer(initializer_param_block_);
  }
  else
    ChiInvalidArgument("Unsupported initializer type \"" + initializer_type +
                       "\".");
}

// ###################################################################
std::vector<std::string> IncompressibleLiquidPhysics::MakeVariableNamesList(
  ComponentCategory hw_comp_category)
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
        "rho", "e", "T", "p", "h", "s", "k", "Pr", "mu", "u", "Re"};
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

double IncompressibleLiquidPhysics::FrictionFactor(const ComponentModel& model)
{
  const double Re = std::max(model.VarOld("Re"), 10.0);
  const double e = model.GetHardwareComponent().Roughness();
  const double D = model.HydraulicDiameter();

  const double a = pow(2.457 * log(pow(7.0 / Re, -0.9) + 0.27 * e / D), 16.0);
  const double b = pow(37530 / Re, 16.0);

  const double fT =
    8.0 * pow(pow(8.0 / Re, 12.0) + pow(a + b, -1.5), 1.0 / 12.0);

  if (Re < 2000.0)
  {
    const double f = 64.0 / Re;

    return std::min(1.0e6, f);
  }
  else if (Re < 4000.0)
  {
    const double c = (Re - 2000.0) / (4000.0 - 2000.0);
    return (1.0 - c) * (64.0 / Re) + c * fT;
  }
  else
    return fT;
}

} // namespace piper