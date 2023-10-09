#ifndef PIPER_HEATCAPACITY_H
#define PIPER_HEATCAPACITY_H

#include "materials/MaterialProperty.h"

namespace chi_math
{
class FunctionDimAToDimB;
}

namespace hcm
{

class HeatCapacity : public chi::MaterialProperty2
{
public:
  static chi::InputParameters GetInputParameters();
  explicit HeatCapacity(const chi::InputParameters& params);

  std::vector<std::string> RequiredInputNames() const override;
  double
  ComputeScalarValue(const chi_mesh::Vector3& position,
                     double time,
                     double param) const override;
  double ComputeScalarValueSlope(const chi_mesh::Vector3& position,
                                 double time,
                                 double param) const override;

  bool HasDerivative() const override;

protected:
  double constant_value_;
  std::shared_ptr<chi_math::FunctionDimAToDimB> function_;
  std::string temperature_fieldname_;

private:
  std::shared_ptr<chi_math::FunctionDimAToDimB>
  GetFunction(const chi::InputParameters& params);
};

} // namespace hcm

#endif // PIPER_HEATCAPACITY_H
