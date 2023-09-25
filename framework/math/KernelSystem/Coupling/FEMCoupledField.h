#ifndef CHITECH_FEMCOUPLEDVARIABLE_H
#define CHITECH_FEMCOUPLEDVARIABLE_H

#include <vector>
#include <string>

namespace chi_math
{

class FEMKernelSystemData;

class FEMCoupledField
{
public:
  double operator[](size_t qp_index) const {return qp_values_.at(qp_index);}

  virtual void ComputeFieldInternalQPValues() = 0;
  virtual void ComputeFieldFaceQPValues() = 0;
  virtual ~FEMCoupledField() = default;

  const std::string& FieldName() const;

protected:
  explicit FEMCoupledField(const std::string& field_name,
                              const FEMKernelSystemData& fem_data);

  const FEMKernelSystemData& fem_data_;

  std::vector<double> qp_values_;

private:
  const std::string& field_name_;
};

} // namespace chi_math

#endif // CHITECH_FEMCOUPLEDVARIABLE_H
