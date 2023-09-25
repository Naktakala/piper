#ifndef CHITECH_NONSYSTEMFEMCOUPLEDVARIABLE_H
#define CHITECH_NONSYSTEMFEMCOUPLEDVARIABLE_H

#include "math/KernelSystem/Coupling/FEMCoupledField.h"

#include <memory>

namespace chi_physics
{
class FieldFunctionGridBased;
}

namespace chi_mesh
{
class Cell;
}

namespace chi_math
{

class SpatialDiscretization;
class UnknownManager;
class CellMapping;

class FieldFunctionFEMCoupledField : public FEMCoupledField
{
public:
  FieldFunctionFEMCoupledField(
    const std::string& field_name,
    const FEMKernelSystemData& fem_data,
    const chi_physics::FieldFunctionGridBased& field);

  void ComputeFieldInternalQPValues() override;
  void ComputeFieldFaceQPValues() override;


protected:
  void SharedSDMComputeFieldInternalQPValues(const chi_mesh::Cell& cell,
                                     const chi_math::CellMapping& cell_mapping);
  void NonSharedSDMComputeFieldInternalQPValues(const chi_mesh::Cell& cell);

  void SharedSDMComputeFieldFaceQPValues(const chi_mesh::Cell& cell,
                                             const chi_math::CellMapping& cell_mapping);
  void NonSharedSDMComputeFieldFaceQPValues(const chi_mesh::Cell& cell);

  const chi_physics::FieldFunctionGridBased& field_;
  const chi_math::SpatialDiscretization& sdm_;
  const chi_math::UnknownManager& uk_man_;
};

} // namespace chi_math

#endif // CHITECH_NONSYSTEMFEMCOUPLEDVARIABLE_H
