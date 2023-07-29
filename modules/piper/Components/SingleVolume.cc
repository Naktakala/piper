#include "SingleVolume.h"

#include "ChiObjectFactory.h"

namespace piper
{

RegisterChiObject(piper, SingleVolume);

chi::InputParameters SingleVolume::GetInputParameters()
{
  chi::InputParameters params = Component::GetInputParameters();

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
  : Component(params,
              ComponentType::Volumetric,
              /*connection_points=*/
              {utils::Connection{"inlet"}, utils::Connection{"outlet"}}),
    Dh_(params.GetParamValue<double>("Dh")),
    A_(params.GetParamValue<double>("A")),
    length_(params.GetParamValue<double>("length")),
    roughness_(params.GetParamValue<double>("roughness"))
{
}

/**The nodalization for a single volume is quite simple. Given the datum point,
 * the other point is just along the length vector.*/
void SingleVolume::Nodalize(std::string connection_point_name,
                            const chi_mesh::Vector3& datum)
{
  size_t ci = 0; // cp = connection point
  for (const auto& connection_point : connection_points_)
  {
    if (connection_point.name_ == connection_point_name) break;
    ++ci;
  }
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

} // namespace piper