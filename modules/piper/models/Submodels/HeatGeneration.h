#ifndef PIPER_HEATGENERATION_H
#define PIPER_HEATGENERATION_H

#include "ChiObject.h"

namespace chi_mesh
{
class Cell;
}

namespace chi_physics
{
class FieldFunction;
}

namespace piper
{
// ##################################################################
class HeatGeneration : public ChiObject
{
public:
  virtual double GetValue(const chi_mesh::Cell& cell) const = 0;

  virtual ~HeatGeneration() = default;

protected:
  static chi::InputParameters GetInputParameters();
  explicit HeatGeneration(const chi::InputParameters& params);
};

// ##################################################################
class ConstantHeatGeneration : public HeatGeneration
{
public:
  static chi::InputParameters GetInputParameters();
  explicit ConstantHeatGeneration(const chi::InputParameters& params);

  double GetValue(const chi_mesh::Cell& cell) const override;

protected:
  double value_;
};

// ##################################################################
class CoupledFieldHeatGeneration : public HeatGeneration
{
public:
  static chi::InputParameters GetInputParameters();
  explicit CoupledFieldHeatGeneration(const chi::InputParameters& params);

  double GetValue(const chi_mesh::Cell& cell) const override;

protected:
  static const chi_physics::FieldFunction&
  GetCoupledField(const std::string& field_name);

  const chi_physics::FieldFunction& coupled_field_;
};

} // namespace piper

#endif // PIPER_HEATGENERATION_H
