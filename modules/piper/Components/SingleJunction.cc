#include "SingleJunction.h"

#include "ChiObjectFactory.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace piper
{

RegisterChiObject(piper, SingleJunction);

chi::InputParameters SingleJunction::GetInputParameters()
{
  chi::InputParameters params = Component::GetInputParameters();

  params.SetGeneralDescription(
    "A junction that connects one component to another");
  params.SetDocGroup("Piper");

  params.AddRequiredParameterArray("from",
                                   "An table containing two strings, the "
                                   "component name, and the connection name");
  params.AddRequiredParameterArray("to",
                                   "An table containing two strings, the "
                                   "component name, and the connection name");

  return params;
}

SingleJunction::SingleJunction(const chi::InputParameters& params)
  : Component(params,
              ComponentType::JunctionLike,
              /*connection_points=*/
              {utils::Connection{"from"}, utils::Connection{"to"}})
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

/**Since junctions don't really have a position we simply set them equal to the
 * datum. This needs to be revisited.*/
void SingleJunction::Nodalize(std::string connection_point_name,
                              const chi_mesh::Vector3& datum)
{
  for (auto& connection : connection_points_)
    connection.position_ = datum;
}

} // namespace piper