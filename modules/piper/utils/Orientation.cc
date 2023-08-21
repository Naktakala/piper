#include "Orientation.h"

#include "ChiObjectFactory.h"
#include <math.h>

namespace piper
{

RegisterChiObjectParametersOnly(piper, Orientation);

chi::InputParameters Orientation::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.SetGeneralDescription(
    "A multi-functional object to establish an orientation");
  params.SetDocGroup("Piper");

  params.AddOptionalParameter("azimuthal",
                              0.0,
                              "Azimuthal angle in degrees. Measured with "
                              "respect to the x-axis and the right-hand rule");
  params.AddOptionalParameter("polar",
                              0.0,
                              "Polar angle in degrees. Measured with respect "
                              "to the z-axis and the right-hand rule.");

  return params;
}

Orientation::Orientation(const chi::InputParameters& params)
  : ChiObject(params),
    varphi_(params.GetParamValue<double>("azimuthal") * M_PI / 180.0),
    theta_(params.GetParamValue<double>("polar") * M_PI / 180.0),
    vector_(sin(theta_) * cos(varphi_), sin(theta_) * sin(varphi_), cos(theta_))
{
}

Orientation Orientation::MakeOrientation(const chi::ParameterBlock& param_block)
{
  chi::InputParameters params = Orientation::GetInputParameters();

  params.AssignParameters(param_block);

  return Orientation(params);
}

double Orientation::Varphi() const { return varphi_; }

double Orientation::Theta() const { return theta_; }

const chi_mesh::Vector3& Orientation::Vector() const
{
  return vector_;
}

} // namespace piper