#include "FluidPhysics.h"

#include "piper/MeshGenerators/PiperMeshGenerator.h"
#include "piper/utils/utils.h"

#include "mesh/MeshHandler/chi_meshhandler.h"

#include "chi_log.h"

#include <iomanip>

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
                                   std::vector<double>{0.0, 0.0, -9.81},
                                   "The direction vector for gravity.");

  params.AddOptionalParameter("dt", 0.01, "Time step size.");
  params.AddOptionalParameter("time", 0.0, "Current time state.");
  params.AddOptionalParameter("end_time", 1.0, "End time.");
  params.AddOptionalParameter("max_time_steps",
                              -1,
                              "Maximum number of time steps to perform. Use a "
                              "negative number to disable.");
  params.AddOptionalParameterArray(
    "model_parameters",
    std::vector<chi::ParameterBlock>{},
    "List of parameters to set for individual components");

  using namespace chi_data_types;
  params.ConstrainParameterRange(
    "dt", AllowableRangeLowLimit::New(0.0, /*low_closed=*/false));

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
    gravity_(params.GetParamVectorValue<double>("gravity_vector")),
    friction_factor_function_(&DarcyFrictionFactorWithChurchill),
    dt_(params.GetParamValue<double>("dt")),
    time_(params.GetParamValue<double>("time")),
    end_time_(params.GetParamValue<double>("end_time")),
    max_time_steps_(params.GetParamValue<int>("max_time_steps")),
    compononent_model_parameters_(
      ConstructComponentModelParametersMap(params.GetParam("model_parameters")))
{
}

// ###################################################################
std::map<std::string, chi::ParameterBlock>
FluidPhysics::ConstructComponentModelParametersMap(
  const chi::ParameterBlock& params)
{
  std::map<std::string, chi::ParameterBlock> param_map;
  for (const auto& param : params)
  {
    const auto& component_name = param.GetParamValue<std::string>("comp_name");
    auto& block = param_map[component_name];
    for (const auto& sub_param : param)
      if (sub_param.Name() != "comp_name") block.AddParameter(sub_param);
  }
  return param_map;
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

double FluidPhysics::DeltaT() const { return dt_; }
void FluidPhysics::SetTimeStep(double dt) { dt_ = dt; }
double FluidPhysics::Time() const { return time_; }
void FluidPhysics::SetTime(double time) { time_ = time; }
double FluidPhysics::EndTime() const { return end_time_; }
void FluidPhysics::SetEndTime(double end_time) { end_time_ = end_time; }

std::function<double(const ComponentModel&)>
FluidPhysics::FrictionFactorFuncion() const
{
  return friction_factor_function_;
}

void FluidPhysics::Execute()
{
  t_index_ = 0;
  while (Time() < EndTime())
  {
    Step();
    Advance();
    ++t_index_;

    std::stringstream outstr;
    outstr << "Timestep " << t_index_ << " dt=" << std::scientific
           << std::setprecision(2) << DeltaT() << " time=" << std::scientific
           << std::setprecision(2) << Time();

    Chi::log.Log() << outstr.str();

    if (max_time_steps_ >= 0 and t_index_ >= max_time_steps_) break;
  }
}

} // namespace piper