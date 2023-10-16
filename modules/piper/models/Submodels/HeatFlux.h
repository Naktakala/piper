#ifndef PIPER_HEATFLUX_H
#define PIPER_HEATFLUX_H

#include "ChiObject.h"

namespace chi_mesh
{
class Cell;
struct Vector3;
} // namespace chi_mesh

namespace chi_physics
{
class FieldFunction;
}

namespace piper
{

class ComponentModel;

// ##################################################################
class HeatFlux : public ChiObject
{
public:
  virtual double
  GetValue(const chi_mesh::Cell& cell,
           const std::map<std::string, double>& variable_reference_map_,
           double length) const = 0;

  virtual double ComputeHCoeff(const ComponentModel& model,
                               const chi_mesh::Vector3& gravity) = 0;

protected:
  static chi::InputParameters GetInputParameters();
  explicit HeatFlux(const chi::InputParameters& params);

  double convection_coefficient_ = 0.0;
};

// ##################################################################
class ConstantHeatFlux : public HeatFlux
{
public:
  static chi::InputParameters GetInputParameters();
  explicit ConstantHeatFlux(const chi::InputParameters& params);

  double GetValue(const chi_mesh::Cell& cell,
                  const std::map<std::string, double>& variable_reference_map_,
                  double length) const override;

  double ComputeHCoeff(const ComponentModel& model,
                       const chi_mesh::Vector3& gravity) override;

protected:
  double value_;
};

// ##################################################################
class ConvectionHeatFlux : public HeatFlux
{
public:
  static chi::InputParameters GetInputParameters();
  explicit ConvectionHeatFlux(const chi::InputParameters& params);

  double GetValue(const chi_mesh::Cell& cell,
                  const std::map<std::string, double>& variable_reference_map_,
                  double length) const override;


  double ComputeHCoeff(const ComponentModel& model,
                       const chi_mesh::Vector3& gravity) override;

protected:
  static std::shared_ptr<const chi_physics::FieldFunction>
  GetSurfaceTemperatureField(const chi::InputParameters& params);
  double EvaluateSurfaceTemperature(const chi_mesh::Cell& cell) const;


  double surface_temperature_value_;
  std::shared_ptr<const chi_physics::FieldFunction> surface_temperature_field_;
  double wetted_perimeter_;
};
} // namespace piper

#endif // PIPER_HEATFLUX_H
