#include "LiquidPhysics.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "piper/MeshGenerators/PiperMeshGenerator.h"
#include "piper/Piper.h"
#include "piper/components/HardwareComponent.h"
#include "piper/models/JunctionLiquidModel.h"
#include "piper/models/SingleVolumeLiquidModel.h"
#include "piper/models/BoundaryLiquidModel.h"

#include "math/ParallelVector/GhostedParallelSTLVector.h"

#include "mesh/MeshHandler/chi_meshhandler.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiObjectFactory.h"

#define scint64_t static_cast<int64_t>

namespace piper
{

void LiquidPhysics::Initialize()
{
  Chi::log.Log() << "Making mesh";
  Chi::mpi.Barrier();
  //this->MakeMesh();
  grid_ptr_ = chi_mesh::GetCurrentHandler().GetGrid();
  Chi::log.Log() << "Initializing Unknowns";
  Chi::mpi.Barrier();
  this->InitializeUnknowns();
}

// ###################################################################
void LiquidPhysics::InitializeUnknowns()
{
  const auto& grid = Grid();
  const auto& pipe_system = PipeSystem();

  auto& hw_components = pipe_system.HardwareComponents();
  const auto& cell_map = pipe_system.GetVolumeComponent2CellGIDMap();

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

  //================================================== Create field functions
  const auto variable_names =
    MakeVariableNamesList(ComponentCategory::Volumetric);
  auto& factory = ChiObjectFactory::GetInstance();
  for (const auto& var_name : variable_names)
  {
    chi::ParameterBlock params;
    params.AddParameter("name", var_name);
    params.AddParameter("sdm_type", "FV");

    const size_t ff_handle = factory.MakeRegisteredObjectOfType(
      "chi_physics::FieldFunctionGridBased", params);

    auto ff = Chi::GetStackItemPtrAsType<chi_physics::FieldFunctionGridBased>(
      Chi::field_function_stack, ff_handle, __FUNCTION__);

    this->field_functions_.push_back(ff);

    varname_2_ff_map_[var_name] = ff;
  } // for each variable name


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

  UpdateFieldFunctions();

  //================================================== Develop linear system
  const int64_t num_local_dofs = scint64_t(grid_ptr_->local_cells.size());
  const int64_t num_global_dofs =
    scint64_t(grid_ptr_->GetGlobalNumberOfCells());

  // const int64_t num_global_dofs =
  //   scint64_t(grid_ptr_->GetGlobalNumberOfCells());
  // const int64_t num_local_dofs =
  //   Chi::mpi.location_id == 0 ? num_global_dofs : 0;

  A_ =
    chi_math::PETScUtils::CreateSquareMatrix(num_local_dofs, num_global_dofs);

  chi_math::PETScUtils::InitMatrixSparsity(A_, 3, 3);

  x_ = chi_math::PETScUtils::CreateVector(num_local_dofs, num_global_dofs);
  b_ = chi_math::PETScUtils::CreateVector(num_local_dofs, num_global_dofs);

  pressure_solver_setup = chi_math::PETScUtils::CreateCommonKrylovSolverSetup(
    A_, "PressureSystem", KSPPREONLY, PCLU);

  KSPSetInitialGuessNonzero(pressure_solver_setup.ksp, PETSC_FALSE);
  KSPSetConvergenceTest(
    pressure_solver_setup.ksp, &KSPConvergedDefault, nullptr, nullptr);
  KSPMonitorCancel(pressure_solver_setup.ksp);

  std::vector<int64_t> ghost_ids;
  for (size_t vol_comp_id : pipe_system_ptr_->VolumeComponentIDs())
  {
    auto& vol_comp_model = GetComponentLiquidModel(vol_comp_id);
    const auto& cell = *vol_comp_model.GetCellPtr();
    ghost_ids.push_back(scint64_t(cell.global_id_));
  }
  pressure_vector_ = std::make_unique<chi_math::GhostedParallelSTLVector>(
    num_local_dofs, num_global_dofs, ghost_ids, Chi::mpi.comm);
}

} // namespace piper