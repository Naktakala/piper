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

/**An interface class that allows simple material properties to interface
 * with FEMKernels.*/
class FEMMaterialProperty
{
public:
  explicit FEMMaterialProperty(const chi::MaterialProperty& property);

  const std::string& TextName() const;
  const std::vector<int>& GetMaterialIDScope() const;

  double operator[](size_t qp_index) const;

  struct DerivativeHelper
  {
    const std::vector<double>& derivative_qp_values_;
    double operator[](size_t qp_index) const;
  } dvar_{derivative_qp_values_};

  bool HasDerivative() const;

  bool IsActive(int material_id) const;

  void SetQPValues(std::vector<double> qp_values,
                   std::vector<double> derivative_qp_values)
  {
    qp_values_ = std::move(qp_values);
    derivative_qp_values_ = std::move(derivative_qp_values);
  }

  double Evaluate(const chi_mesh::Vector3& xyz,
                  double t,
                  const std::vector<double>& v);

protected:
  const chi::MaterialProperty& property_;

  std::vector<double> qp_values_;
  std::vector<double> derivative_qp_values_;
};

} // namespace chi_math

#endif // CHITECH_FEMMATERIALPROPERTY_H
