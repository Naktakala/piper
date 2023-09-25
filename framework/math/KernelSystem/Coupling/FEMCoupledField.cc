#include "FEMCoupledField.h"

namespace chi_math
{

FEMCoupledField::FEMCoupledField(const std::string& field_name,
                                       const FEMKernelSystemData& fem_data)
  : fem_data_(fem_data), field_name_(field_name)
{
}

const std::string& FEMCoupledField::FieldName() const { return field_name_; }

} // namespace chi_math