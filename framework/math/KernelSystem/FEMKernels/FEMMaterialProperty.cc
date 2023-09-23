#include "FEMMaterialProperty.h"

#include "materials/MaterialProperty.h"
#include "mesh/chi_mesh.h"

#include "chi_log.h"

namespace chi_math
{

FEMMaterialProperty::FEMMaterialProperty(const chi::MaterialProperty& property)
  : property_(property)
{
}

const std::string& FEMMaterialProperty::TextName() const
{
  return property_.TextName();
}

const std::vector<int>& FEMMaterialProperty::GetMaterialIDScope() const
{
  return property_.GetMaterialIDScope();
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

bool FEMMaterialProperty::IsActive(int material_id) const
{
  const auto& mat_scope = property_.GetMaterialIDScope();
  if (mat_scope.empty()) return true;
  return std::find(mat_scope.begin(), mat_scope.end(), material_id) !=
         mat_scope.end();
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

} // namespace chi_math