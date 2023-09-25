#ifndef CHITECH_FEMMATERIALPROPERTY_H
#define CHITECH_FEMMATERIALPROPERTY_H

#include <vector>
#include <cstddef>
#include <string>

namespace chi_mesh
{
class Vector3;
}

namespace chi
{
class MaterialProperty;
}

namespace chi_math
{

class FEMKernelSystemData;

/**An interface class that allows simple material properties to interface
 * with FEMKernels.*/
class FEMMaterialProperty
{
public:
  explicit FEMMaterialProperty(const chi::MaterialProperty& property,
                               const FEMKernelSystemData& fem_data);

  double operator[](size_t qp_index) const;

  struct DerivativeHelper
  {
    const std::vector<double>& derivative_qp_values_;
    double operator[](size_t qp_index) const;
  } dvar_{derivative_qp_values_};

  bool HasDerivative() const;

  void PreComputeInternalQPValues();
  void PreComputeFaceQPValues();

protected:
  double Evaluate(const chi_mesh::Vector3& xyz,
                  double t,
                  const std::vector<double>& v);
  double EvaluateSlope(const chi_mesh::Vector3& xyz,
                       double t,
                       const std::vector<double>& v);

  const chi::MaterialProperty& property_;
  const FEMKernelSystemData& fem_data_;

  std::vector<double> qp_values_;
  std::vector<double> derivative_qp_values_;
};

} // namespace chi_math

#endif // CHITECH_FEMMATERIALPROPERTY_H
