#ifndef CHITECH_ADAPTIVETIMESTEPCONTROLLER_H
#define CHITECH_ADAPTIVETIMESTEPCONTROLLER_H

#include "physics/TimeStepControllers/TimeStepController.h"

namespace chi_physics
{

class AdaptiveTimeStepController : public TimeStepController
{
public:
  static chi::InputParameters GetInputParameters();
  explicit AdaptiveTimeStepController(const chi::InputParameters& params);

  double GetTimeStepSize() override;
  bool Adapt(TimeStepStatus time_step_status) override;

protected:
  // parameters
  double decrease_ratio_;
  double increase_ratio_;
  size_t iteration_window_size_;
  double dt_min_;
  double dt_max_;

  // runtime
  size_t last_change_index_;
};

}

#endif // CHITECH_ADAPTIVETIMESTEPCONTROLLER_H
