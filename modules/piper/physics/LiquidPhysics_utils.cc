#include "LiquidPhysics.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "piper/Piper.h"
#include "piper/components/HardwareComponent.h"
#include "piper/MeshGenerators/PiperMeshGenerator.h"
#include "piper/models/BoundaryLiquidModel.h"
#include "piper/models/JunctionLiquidModel.h"
#include "piper/models/SingleVolumeLiquidModel.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace piper
{

// ###################################################################
const std::string& LiquidPhysics::FluidName() const { return fluid_name_; }

// ###################################################################
void LiquidPhysics::InitializeUnknowns()
{
  const auto& grid = Grid();
  const auto& pipe_system = PipeSystem();

  auto& hw_components = pipe_system.HardwareComponents();
  const auto& cell_map = mesh_generator_->GetVolumeComponent2CellGIDMap();

  //================================================== Instantiate models
  for (auto& hw_component_ptr : hw_components)
  {
    const auto& hw_component = *hw_component_ptr;
    const auto& hw_comp_typename = hw_component.TypeName();
    const auto hw_comp_category = hw_component.Category();
    const auto& comp_name = hw_component.Name();

    const chi_mesh::Cell* cell_ptr = nullptr;
    auto cell_iter = cell_map.find(hw_component.Name());
    if (cell_iter != cell_map.end())
      cell_ptr = &(grid.cells[(*cell_iter).second]);

    const auto variable_names = MakeVariableNamesList(hw_comp_category);

    const std::string supported_typenames =
      "\"BoundaryComponent\", \"SingleVolume\", \"SingleJunction\"";
    if (hw_comp_typename == "BoundaryComponent")
    {
      auto valid_params = BoundaryLiquidModel::GetInputParameters();
      valid_params.AssignParameters(compononent_model_parameters_[comp_name]);

      component_models_.push_back(
        std::make_unique<BoundaryLiquidModel>(valid_params,
                                              component_models_,
                                              *this,
                                              hw_component,
                                              cell_ptr,
                                              variable_names));
    }
    else if (hw_comp_typename == "SingleVolume")
    {
      auto valid_params = SingleVolumeLiquidModel::GetInputParameters();
      valid_params.AssignParameters(compononent_model_parameters_[comp_name]);

      component_models_.push_back(
        std::make_unique<SingleVolumeLiquidModel>(valid_params,
                                                  component_models_,
                                                  *this,
                                                  hw_component,
                                                  cell_ptr,
                                                  variable_names));
    }
    else if (hw_comp_typename == "SingleJunction")
    {
      auto valid_params = JunctionLiquidModel::GetInputParameters();
      valid_params.AssignParameters(compononent_model_parameters_[comp_name]);

      component_models_.push_back(
        std::make_unique<JunctionLiquidModel>(valid_params,
                                              component_models_,
                                              *this,
                                              hw_component,
                                              cell_ptr,
                                              variable_names));
    }
    else
      ChiLogicalError("LiquidPhysics has no available model for hardware "
                      "component type \"" +
                      hw_comp_typename +
                      "\". Supported "
                      "types are " +
                      supported_typenames);
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
