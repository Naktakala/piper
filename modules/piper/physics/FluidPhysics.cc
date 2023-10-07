#include "FluidPhysics.h"

#include "piper/Piper.h"
#include "piper/MeshGenerators/PiperMeshGenerator.h"
#include "piper/utils/utils.h"

#include "physics/PhysicsEventPublisher.h"
#include "physics/TimeSteppers/TimeStepper.h"

#include "mesh/MeshHandler/chi_meshhandler.h"

#include "chi_log.h"

#include <iomanip>

namespace piper
{

chi::InputParameters FluidPhysics::GetInputParameters()
{
  chi::InputParameters params = chi_physics::Solver::GetInputParameters();

  params.AddOptionalParameter(
    "min_cells_per_processor",
    10,
    "The minimum amount of cells per processor. If there are fewer than this "
    "many cells, all the cells will be loaded onto a single process.");

  params.AddRequiredParameter<size_t>("pipe_system",
                                      "Handle to a Piper object.");

  params.AddRequiredParameterBlock(
    "initializer",
    "A parameter block for an initializer."); // TODO: Improve

  params.AddOptionalParameterArray("gravity_vector",
                                   std::vector<double>{0.0, 0.0, -9.81},
                                   "The direction vector for gravity.");

  params.AddOptionalParameterArray(
    "model_parameters",
    std::vector<chi::ParameterBlock>{},
    "List of parameters to set for individual components");

  return params;
}

FluidPhysics::FluidPhysics(const chi::InputParameters& params)
  : chi_physics::Solver(params),
    min_cells_per_processor_(
      params.GetParamValue<size_t>("min_cells_per_processor")),
    initializer_param_block_(params.GetParam("initializer")),
    pipe_system_ptr_(
      &Chi::GetStackItem<Piper>(Chi::object_stack,
                                params.GetParamValue<size_t>("pipe_system"),
                                __FUNCTION__)),
    mesh_generator_(nullptr),
    grid_ptr_(nullptr),
    gravity_(params.GetParamVectorValue<double>("gravity_vector")),
    friction_factor_function_(&DarcyFrictionFactorWithChurchill),
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

std::function<double(const ComponentModel&)>
FluidPhysics::FrictionFactorFuncion() const
{
  return friction_factor_function_;
}

void FluidPhysics::Execute()
{
  const size_t tag_jnc_time = Chi::log.GetRepeatingEventTag("JNC_TIME");
  const size_t tag_vol_time = Chi::log.GetRepeatingEventTag("VOL_TIME");
  const size_t tag_slv_time = Chi::log.GetRepeatingEventTag("SLV_TIME");
  const std::vector<size_t> tags = {Chi::log.GetRepeatingEventTag("TIMING0"),
                                    Chi::log.GetRepeatingEventTag("TIMING1"),
                                    Chi::log.GetRepeatingEventTag("TIMING2"),
                                    Chi::log.GetRepeatingEventTag("TIMING3"),
                                    Chi::log.GetRepeatingEventTag("TIMING4")};

  size_t t_index = 0;
  while (timestepper_->Time() < timestepper_->EndTime())
  {
    auto& physics_event_publisher =
      chi_physics::PhysicsEventPublisher::GetInstance();
    physics_event_publisher.SolverStep(*this);
    physics_event_publisher.SolverAdvance(*this);
    ++t_index;

    const double dt = timestepper_->TimeStepSize();
    const double time = timestepper_->Time();
    std::stringstream outstr;
    outstr << "Timestep " << t_index << " dt=" << std::scientific
           << std::setprecision(2) << dt << " time=" << std::scientific
           << std::setprecision(2) << time << std::fixed
           << " step_time=" << step_time_ << "ms";

    Chi::log.Log() << outstr.str();

    const auto max_time_steps = timestepper_->MaxTimeSteps();

    if (max_time_steps >= 0 and t_index >= max_time_steps) break;
  }

  const double jnc_time = Chi::log.ProcessEvent(
    tag_jnc_time, chi::ChiLog::EventOperation::AVERAGE_DURATION);
  const double vol_time = Chi::log.ProcessEvent(
    tag_vol_time, chi::ChiLog::EventOperation::AVERAGE_DURATION);
  const double slv_time = Chi::log.ProcessEvent(
    tag_slv_time, chi::ChiLog::EventOperation::AVERAGE_DURATION);
  Chi::log.Log() << "Timing junction assembly: " << jnc_time * 1000 << "ms";
  Chi::log.Log() << "Timing volume assembly  : " << vol_time * 1000 << "ms";
  Chi::log.Log() << "Timing solve            : " << slv_time * 1000 << "ms";
  size_t k = 0;
  for (size_t tag : tags)
    Chi::log.Log() << "Timing " << k++ << " : "
                   << Chi::log.ProcessEvent(
                        tag, chi::ChiLog::EventOperation::AVERAGE_DURATION) *
                        1000
                   << "ms";
}

chi::ParameterBlock
FluidPhysics::GetInfo(const chi::ParameterBlock& params) const
{
  const auto info_name = params.GetParamValue<std::string>("name");
  if (info_name == "ComponentVariable")
  {
    params.RequireParameter("component_name");
    params.RequireParameter("var_name");

    const auto component_name =
      params.GetParamValue<std::string>("component_name");
    const auto var_name = params.GetParamValue<std::string>("var_name");

    const size_t comp_id = pipe_system_ptr_->MapHWCompName2ID(component_name);
    const auto& comp_model = component_models_.at(comp_id);

    return chi::ParameterBlock("", comp_model->VarOld(var_name));
    // auto item = com
  }
  else
    ChiInvalidArgument("Unsupported information label \"" + info_name + "\".");

  return chi::ParameterBlock("", 0.0);
}

} // namespace piper