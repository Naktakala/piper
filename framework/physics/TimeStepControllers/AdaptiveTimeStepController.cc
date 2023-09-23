#include "AdaptiveTimeStepController.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

namespace chi_physics
{

RegisterChiObject(chi_physics, AdaptiveTimeStepController);

chi::InputParameters AdaptiveTimeStepController::GetInputParameters()
{
  chi::InputParameters params = TimeStepController::GetInputParameters();

  params.SetGeneralDescription("General adaptive timestep controller");
  params.SetDocGroup("doc_TimeStepControllers");

  params.AddOptionalParameter("decrease_ratio",
                              0.5,
                              "Ratio by which the timestep will be multiplied "
                              "when adapting to a failure.");

  params.AddOptionalParameter(
    "increase_ratio",
    2.0,
    "Ratio by which dt will be multiplied when attempting to increase the "
    "timestep beyond the iteration window.");

  params.AddOptionalParameter(
    "iteration_window_size",
    5,
    "Number of iterations before trying to increase the timestep.");

  params.AddOptionalParameter("dt_min", 1.0e-4, "Minimum timestep size.");
  params.AddOptionalParameter("dt_max",
                              0.0,
                              "Optional. Maximum time-step size. If not "
                              "supplied, will default to initial dt");

  return params;
}

AdaptiveTimeStepController::AdaptiveTimeStepController(
  const chi::InputParameters& params)
  : TimeStepController(params),
    decrease_ratio_(params.GetParamValue<double>("decrease_ratio")),
    increase_ratio_(params.GetParamValue<double>("increase_ratio")),
    iteration_window_size_(
      params.GetParamValue<size_t>("iteration_window_size")),
    dt_min_(params.GetParamValue<double>("dt_min")),
    dt_max_(params.ParametersAtAssignment().Has("dt_max")
              ? params.GetParamValue<double>("dt_max")
              : dt_),

    last_change_index_(0)
{
}

double AdaptiveTimeStepController::GetTimeStepSize()
{
  if (t_index_ > (last_change_index_ + iteration_window_size_))
  {
    dt_ *= increase_ratio_;
    dt_ = std::min(dt_, dt_max_);
  }

  if (not AtEndTime(time_ + dt_)) dt_ = std::max(dt_, dt_min_);
  else
    dt_ = StepSizeLimitedToEndTime();

  return dt_;
}

bool AdaptiveTimeStepController::Adapt(TimeStepStatus time_step_status)
{
  if (time_step_status == TimeStepStatus::SUCCESS or
      time_step_status == TimeStepStatus::NEUTRAL)
    return true;
  if (time_step_status == TimeStepStatus::FAILURE)
  {
    if (dt_ > dt_min_)
    {
      dt_ *= decrease_ratio_;
      last_change_index_ = t_index_;
      return true;
    }
    else
      return false;
  }

  return false;
}

} // namespace chi_physics