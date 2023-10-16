#ifndef PIPER_VOLUMETRICCOMPONENTMODEL_H
#define PIPER_VOLUMETRICCOMPONENTMODEL_H

#include "ComponentModel.h"
#include "parameters/input_parameters.h"

namespace piper
{

class HeatGeneration;
class HeatFlux;

class VolumetricComponentModel
{
protected:
  static chi::InputParameters GetInputParameters();
  explicit VolumetricComponentModel(const chi::ParameterBlock& params);

  std::vector<std::shared_ptr<HeatGeneration>> heat_generation_sub_models_;
  std::vector<std::shared_ptr<HeatFlux>> heat_flux_sub_models_;
};

} // namespace piper

#endif // PIPER_VOLUMETRICCOMPONENTMODEL_H
