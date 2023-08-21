#include "FluidPhysics.h"

#include "piper/MeshGenerators/PiperMeshGenerator.h"

#include "mesh/MeshHandler/chi_meshhandler.h"

#include "chi_log.h"

namespace piper
{

chi::InputParameters FluidPhysics::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddOptionalParameter(
    "min_cells_per_processor",
    10,
    "The minimum amount of cells per processor. If there are fewer than this "
    "many cells, all the cells will be loaded onto a single process.");

  params.AddRequiredParameterBlock(
    "initializer",
    "A parameter block for an initializer."); // TODO: Improve

  params.AddOptionalParameterArray("gravity_vector",
                                   std::vector<double>{0.0,0.0,-9.81},
                                   "The direction vector for gravity.");

  return params;
}

FluidPhysics::FluidPhysics(const chi::InputParameters& params)
  : ChiObject(params),
    min_cells_per_processor_(
      params.GetParamValue<size_t>("min_cells_per_processor")),
    initializer_param_block_(params.GetParam("initializer")),
    pipe_system_ptr_(nullptr),
    mesh_generator_(nullptr),
    grid_ptr_(nullptr),
    gravity_(params.GetParamVectorValue<double>("gravity_vector"))
{
}

// ###################################################################
void FluidPhysics::MakeMesh()
{
  chi::ParameterBlock params;
  auto valid_params = PiperMeshGenerator::GetInputParameters();
  valid_params.AssignParameters(params);
  mesh_generator_ = std::make_unique<PiperMeshGenerator>(valid_params);

  mesh_generator_->SetPipeSystem(PipeSystem());
  mesh_generator_->Execute();

  grid_ptr_ = chi_mesh::GetCurrentHandler().GetGrid();
}

void FluidPhysics::SetPipeSystem(const Piper& pipe_system)
{
  pipe_system_ptr_ = &pipe_system;
}

const Piper& FluidPhysics::PipeSystem() const
{
  ChiLogicalErrorIf(not pipe_system_ptr_, "Uninitialized pipe system");
  return *pipe_system_ptr_;
}

const chi_mesh::Vector3& FluidPhysics::GravityVector() const
{
  return gravity_;
}

const chi_mesh::MeshContinuum& FluidPhysics::Grid() const
{
  ChiLogicalErrorIf(not grid_ptr_, "Uninitialized mesh");
  return *grid_ptr_;
}

} // namespace piper