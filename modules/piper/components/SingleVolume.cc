#include "SingleVolume.h"

#include "ChiObjectFactory.h"

namespace piper
{

RegisterChiObject(piper, SingleVolume);

chi::InputParameters SingleVolume::GetInputParameters()
{
  chi::InputParameters params = HardwareComponent::GetInputParameters();

  params.SetGeneralDescription(
    "Single control volume component with an inlet and an outlet");
  params.SetDocGroup("Piper");

  params.AddRequiredParameter<double>("Dh", "Hydraulic diameter.");
  params.AddRequiredParameter<double>("A", "Flow area.");
  params.AddRequiredParameter<double>("length", "Length.");

  params.AddOptionalParameter("roughness", 0.0, "Wall roughness");

  return params;
}

SingleVolume::SingleVolume(const chi::InputParameters& params)
  : HardwareComponent(
      params,
      "SingleVolume",
      ComponentCategory::Volumetric,
      /*connection_points=*/
      {utils::Connection{"inlet"}, utils::Connection{"outlet"}}),
    Dh_(params.GetParamValue<double>("Dh")),
    A_(params.GetParamValue<double>("A")),
    length_(params.GetParamValue<double>("length")),
    roughness_(params.GetParamValue<double>("roughness"))
{
}

/**Returns the component's flow orientation relative to the connection
 * point.*/
utils::FlowOrientation
SingleVolume::FlowOrientationRelToConPoint(size_t con_point_id) const
{
  if (con_point_id == 0) return utils::FlowOrientation::OUTGOING;
  else if (con_point_id == 1)
    return utils::FlowOrientation::INCOMING;
  else
    ChiInvalidArgument("Invalid con_point_id " + std::to_string(con_point_id));
}

/**The nodalization for a single volume is quite simple. Given the datum point,
 * the other point is just along the length vector.*/
void SingleVolume::Nodalize(size_t connection_point_id,
                            const chi_mesh::Vector3& datum)
{
  const size_t ci = connection_point_id;
  const size_t cj = ci == 0 ? 1 : 0;

  connection_points_[ci].position_ = datum;

  const double varphi = orientation_.Varphi();
  const double theta = orientation_.Theta();

  const auto length_vector =
    length_ * chi_mesh::Vector3(
                cos(varphi) * sin(theta), sin(varphi) * sin(theta), cos(theta));

  if (ci == 0) connection_points_[cj].position_ = datum + length_vector;
  else
    connection_points_[cj].position_ = datum - length_vector;
}

double SingleVolume::Area() const { return A_; }

double SingleVolume::Volume() const { return A_ * length_; }

double SingleVolume::Length() const { return length_; }

double SingleVolume::HydraulicDiameter() const { return Dh_; }

} // namespace piper