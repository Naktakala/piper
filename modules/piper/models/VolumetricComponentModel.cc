#include "VolumetricComponentModel.h"

#include "Submodels/HeatGeneration.h"
#include "Submodels/HeatFlux.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

namespace piper
{

chi::InputParameters VolumetricComponentModel::GetInputParameters()
{
  chi::InputParameters params;

  params.AddOptionalParameterBlock(
    "heat_generation_models", chi::ParameterBlock{}, "Blocks of input");
  params.SetParameterTypeMismatchAllowed("heat_generation_models");

  params.AddOptionalParameterBlock(
    "heat_flux_models", chi::ParameterBlock{}, "Blocks of input");
  params.SetParameterTypeMismatchAllowed("heat_flux_models");

  return params;
}

VolumetricComponentModel::VolumetricComponentModel(
  const chi::ParameterBlock& params)
{
  auto& object_factory = ChiObjectFactory::GetInstance();

  if (params.Has("heat_generation_models"))
  {
    const auto& heat_models_param = params.GetParam("heat_generation_models");
    heat_models_param.RequireBlockTypeIs(chi::ParameterBlockType::ARRAY);

    for (const auto& heat_params : heat_models_param)
    {
      heat_params.RequireParameter("type");

      const auto obj_type = heat_params.GetParamValue<std::string>("type");
      chi::InputParameters in_params =
        object_factory.GetRegisteredObjectParameters(obj_type);
      in_params.AssignParameters(heat_params);

      const size_t handle =
        object_factory.MakeRegisteredObjectOfType(obj_type, in_params);

      auto heat_gen_obj = Chi::GetStackItemPtrAsType<HeatGeneration>(
        Chi::object_stack, handle, __FUNCTION__);

      heat_generation_sub_models_.push_back(heat_gen_obj);
    } // for heat param
  }   // heat generation models

  if (params.Has("heat_flux_models"))
  {
    const auto& flux_models_param = params.GetParam("heat_flux_models");
    flux_models_param.RequireBlockTypeIs(chi::ParameterBlockType::ARRAY);

    for (const auto& flux_params : flux_models_param)
    {
      flux_params.RequireParameter("type");

      const auto& obj_type = flux_params.GetParamValue<std::string>("type");

      const size_t handle =
        object_factory.MakeRegisteredObjectOfType(obj_type, flux_params);

      auto flux_obj = Chi::GetStackItemPtrAsType<HeatFlux>(
        Chi::object_stack, handle, __FUNCTION__);

      heat_flux_sub_models_.push_back(flux_obj);
    }// for flux param
  }// heat flux models
}

} // namespace piper