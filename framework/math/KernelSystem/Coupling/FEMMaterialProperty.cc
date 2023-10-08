#include "FEMMaterialProperty.h"

#include "math/KernelSystem/KernelSystem.h"
#include "materials/MaterialProperty.h"
#include "mesh/chi_mesh.h"

#include "chi_log.h"

namespace chi_math
{

FEMMaterialProperty::FEMMaterialProperty(const chi::MaterialProperty2& property,
                                         const FEMKernelSystemData& fem_data)
  : property_(property), fem_data_(fem_data)
{
}

double FEMMaterialProperty::operator[](size_t qp_index) const
{
  return qp_values_[qp_index];
}

double FEMMaterialProperty::DerivativeHelper::operator[](size_t qp_index) const
{
  return derivative_qp_values_[qp_index];
}

bool FEMMaterialProperty::HasDerivative() const
{
  return property_.HasDerivative();
}

double FEMMaterialProperty::Evaluate(const chi_mesh::Vector3& xyz,
                                     double t,
                                     const std::vector<double>& v)
{
  std::vector<double> inputs;
  if (property_.IsPositionDependent()) inputs = {xyz.x, xyz.y, xyz.z};
  if (property_.IsTimeDependent()) inputs.push_back(t);

  for (double value : v)
    inputs.push_back(value);
  return property_.ComputeScalarValue(inputs);
}

double FEMMaterialProperty::EvaluateSlope(const chi_mesh::Vector3& xyz,
                                          double t,
                                          const std::vector<double>& v)
{
  std::vector<double> inputs;
  if (property_.IsPositionDependent()) inputs = {xyz.x, xyz.y, xyz.z};
  if (property_.IsTimeDependent()) inputs.push_back(t);

  for (double value : v)
    inputs.push_back(value);
  return property_.ComputeScalarValueSlope(inputs);
}

void FEMMaterialProperty::PreComputeInternalQPValues()
{
  const double time = fem_data_.time_data_.time_;

  const auto& qp_data = fem_data_.cell_data_.qp_data_;
  const auto& var_qp_values = fem_data_.cell_data_.var_qp_values_;

  const auto& qp_indices = qp_data.QuadraturePointIndices();
  const auto& qp_xyz = qp_data.QPointsXYZ();
  qp_values_.assign(qp_indices.size(), 0.0);
  derivative_qp_values_.assign(qp_indices.size(), 0.0);

  for (size_t qp : qp_indices)
  {
    const double var_value_qp = var_qp_values[qp];
    qp_values_[qp] = this->Evaluate(qp_xyz[qp], time, {var_value_qp});

    if (property_.HasDerivative())
      derivative_qp_values_[qp] =
        this->EvaluateSlope(qp_xyz[qp], time, {var_value_qp});
  }
}

void FEMMaterialProperty::PreComputeFaceQPValues()
{
  const double time = fem_data_.time_data_.time_;

  const auto& qp_data = fem_data_.face_data_.qp_data_;
  const auto& var_qp_values = fem_data_.face_data_.var_qp_values_;

  const auto& qp_indices = qp_data.QuadraturePointIndices();
  const auto& qp_xyz = qp_data.QPointsXYZ();
  qp_values_.assign(qp_indices.size(), 0.0);
  derivative_qp_values_.assign(qp_indices.size(), 0.0);

  for (size_t qp : qp_indices)
  {
    const double var_value_qp = var_qp_values[qp];
    qp_values_[qp] = this->Evaluate(qp_xyz[qp], time, {var_value_qp});

    if (property_.HasDerivative())
    {
      derivative_qp_values_[qp] =
        this->EvaluateSlope(qp_xyz[qp], time, {var_value_qp});
    }
  }
}

} // namespace chi_math