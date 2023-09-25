#ifndef CHITECH_SYSTEMFEMCOUPLEDVARIABLE_H
#define CHITECH_SYSTEMFEMCOUPLEDVARIABLE_H

#include "FEMCoupledField.h"

namespace chi_mesh
{
class Cell;
}

namespace chi_math
{

class MultiFieldContainer;
class SpatialDiscretization;
class CellMapping;

class FFContainerFEMCoupledVariable : public FEMCoupledField
{
public:
  FFContainerFEMCoupledVariable(const std::string& field_name,
                           const FEMKernelSystemData& fem_data,
                           const MultiFieldContainer& multifield_container,
                           size_t field_index);

  void ComputeFieldInternalQPValues() override;
  void ComputeFieldFaceQPValues() override;

protected:
  void SharedSDMComputeFieldInternalQPValues(const chi_mesh::Cell& cell,
                                     const chi_math::CellMapping& cell_mapping);
  void NonSharedSDMComputeFieldInternalQPValues(const chi_mesh::Cell& cell);

  void SharedSDMComputeFieldFaceQPValues(const chi_mesh::Cell& cell,
                                         const chi_math::CellMapping& cell_mapping);
  void NonSharedSDMComputeFieldFaceQPValues(const chi_mesh::Cell& cell);

  const MultiFieldContainer& container_;
  const size_t field_index_;
  const chi_math::SpatialDiscretization& sdm_;
};

}

#endif // CHITECH_SYSTEMFEMCOUPLEDVARIABLE_H
