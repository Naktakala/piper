#ifndef CHITECH_MATERIALPROPERTY_H
#define CHITECH_MATERIALPROPERTY_H

#include "ChiObject.h"
#include "MaterialIDScopeInterface.h"

namespace chi
{

/**Base class for a material property.*/
class MaterialProperty2 : public ChiObject, public MaterialIDScopeInterface
{
public:
  static chi::InputParameters GetInputParameters();
  explicit MaterialProperty2(const chi::InputParameters& params);

  const std::string& TextName() const;

  virtual std::vector<std::string> RequiredInputNames() const = 0;

  virtual double
  ComputeScalarValue(const std::vector<double>& input_params) const = 0;
  virtual double ComputeScalarValueSlope(
    const std::vector<double>& input_params) const {return 0.0;}

  bool IsPositionDependent() const;
  bool IsTimeDependent() const;
  virtual bool HasDerivative() const { return false; }

  virtual ~MaterialProperty2() = default;

private:
  const std::string name_;
  bool position_dependent_;
  bool time_dependent_;
};

} // namespace chi

#endif // CHITECH_MATERIALPROPERTY_H
