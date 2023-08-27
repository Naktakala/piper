#include "SingleJunction.h"

#include "ChiObjectFactory.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace piper
{

RegisterChiObject(piper, SingleJunction);

chi::InputParameters SingleJunction::GetInputParameters()
{
  chi::InputParameters params = HardwareComponent::GetInputParameters();

  params.SetGeneralDescription(
    "A junction that connects one component to another");
  params.SetDocGroup("Piper");

  params.AddRequiredParameterArray("from",
                                   "An table containing two strings, the "
                                   "component name, and the connection name");
  params.AddRequiredParameterArray("to",
                                   "An table containing two strings, the "
                                   "component name, and the connection name");

  chi::ParameterBlock default_area_param;
  default_area_param.SetBlockName("A");
  default_area_param.AddParameter("type", "use_from_area");

  params.AddOptionalParameterBlock("A", default_area_param, "Flow area.");
  params.SetParameterTypeMismatchAllowed("A");

  return params;
}

SingleJunction::SingleJunction(const chi::InputParameters& params)
  : HardwareComponent(params,
                      "SingleJunction",
                      ComponentCategory::JunctionLike,
                      /*connection_points=*/
                      {utils::Connection{"from"}, utils::Connection{"to"}}),
    A_(0.0),
    area_param_(params.GetParam("A"))
{
  auto from_params = params.GetParamVectorValue<std::string>("from");
  auto to_params = params.GetParamVectorValue<std::string>("to");

  ChiLogicalErrorIf(from_params.size() != 2,
                    "\"from\" parameter requires 2 values.");
  ChiLogicalErrorIf(to_params.size() != 2,
                    "\"to\" parameter requires 2 values.");

  connection_points_ = {
    utils::Connection{"from", from_params[0], from_params[1]},
    utils::Connection{"to", to_params[0], to_params[1]}};

  Chi::log.Log0Verbose1() << Name() << ": connections: "
                          << "from=" << from_params[0] << "," << from_params[1]
                          << " to=" << to_params[0] << "," << to_params[1];
}

/**Returns the component's flow orientation relative to the connection
 * point.*/
utils::FlowOrientation
SingleJunction::FlowOrientationRelToConPoint(size_t con_point_id) const
{
  if (con_point_id == 0) return utils::FlowOrientation::OUTGOING;
  else if (con_point_id == 1)
    return utils::FlowOrientation::INCOMING;
  else
    ChiInvalidArgument("Invalid con_point_id " + std::to_string(con_point_id));
}

/**Since junctions don't really have a position we simply set them equal to the
 * datum. This needs to be revisited.*/
void SingleJunction::Nodalize(size_t connection_point_id,
                              const chi_mesh::Vector3& datum)
{
  for (auto& connection : connection_points_)
    connection.position_ = datum;
}

void SingleJunction::EstablishHardwareReferenceData(
  const std::vector<std::shared_ptr<HardwareComponent>>& hardware_components)
{
  if (area_param_.Type() == chi::ParameterBlockType::FLOAT)
  {
    A_ = area_param_.GetValue<double>();

    Chi::log.Log0Verbose1()
      << "Junction \"" + Name() + "\" area specified as A=" << A_;
  }
  else if (area_param_.Type() == chi::ParameterBlockType::BLOCK)
  {
    area_param_.SetErrorOriginScope("Junction \"" + Name() + "\"");
    area_param_.RequireParameter("type");
    const std::string type = area_param_.GetParamValue<std::string>("type");

    const size_t from_comp_id = connection_points_.front().connected_comp_id_;
    const size_t to_comp_id = connection_points_.back().connected_comp_id_;

    const auto& from_comp = hardware_components.at(from_comp_id);
    const auto& to_comp = hardware_components.at(to_comp_id);

    if (type == "use_from_area")
    {
      if (from_comp->Category() == ComponentCategory::BoundaryLike)
        A_ = to_comp->Area();
      else
        A_ = from_comp->Area();
    }
    else if (type == "use_to_area")
    {
      if (to_comp->Category() == ComponentCategory::BoundaryLike)
        A_ = from_comp->Area();
      else
        A_ = to_comp->Area();
    }
    else
      ChiInvalidArgument("Invalid area specification type \"" + type + "\".");

    Chi::log.Log0Verbose1() << "Junction \"" + Name() + "\" area spec: " << type
                            << " set to A=" << A_;
  }
  else
    ChiInvalidArgument("Invalid parameter type \"" +
                       chi::ParameterBlockTypeName(area_param_.Type()) +
                       "\" supplied for area parameter of junction \"" +
                       Name() + "\".");
}

double SingleJunction::Area() const { return A_; }

} // namespace piper