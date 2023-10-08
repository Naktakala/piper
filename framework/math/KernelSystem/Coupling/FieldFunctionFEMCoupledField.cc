#include "FieldFunctionFEMCoupledField.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/KernelSystem/KernelSystem.h"

#include "chi_log.h"

#define cint64_t const int64_t

namespace chi_math
{

FieldFunctionFEMCoupledField::FieldFunctionFEMCoupledField(
  const std::string& field_name,
  const FEMKernelSystemData& fem_data,
  const chi_physics::FieldFunctionGridBased& field)
  : FEMCoupledField(field_name, fem_data),
    field_(field),
    sdm_(field_.GetSpatialDiscretization()),
    uk_man_(field_.GetUnknownManager())
{
}

// ##################################################################
void FieldFunctionFEMCoupledField::ComputeFieldInternalQPValues()
{
  const auto& cell = *fem_data_.cell_data_.cell_ptr_;
  const auto& cell_mapping = sdm_.GetCellMapping(cell);

  if (&cell_mapping == fem_data_.cell_data_.cell_mapping_ptr_)
    SharedSDMComputeFieldInternalQPValues(cell, cell_mapping);
  else
    NonSharedSDMComputeFieldInternalQPValues(cell);
}

void FieldFunctionFEMCoupledField::SharedSDMComputeFieldInternalQPValues(
  const chi_mesh::Cell& cell, const chi_math::CellMapping& cell_mapping)
{
  const size_t num_nodes = cell_mapping.NumNodes();
  const auto& field_data = field_.FieldVectorRead();
  const auto& shape_values = fem_data_.cell_data_.qp_data_.ShapeValues();

  const auto& qp_indices =
    fem_data_.cell_data_.qp_data_.QuadraturePointIndices();
  qp_values_.assign(qp_indices.size(), 0.0);
  for (uint32_t qp : qp_indices)
  {
    for (size_t i = 0; i < num_nodes; ++i)
    {
      cint64_t dof_id = sdm_.MapDOFLocal(cell, i, uk_man_, 0, 0);
      qp_values_[qp] += shape_values[i][qp] * field_data[dof_id];
    } // for i
  }   // for qp
}

void FieldFunctionFEMCoupledField::NonSharedSDMComputeFieldInternalQPValues(
  const chi_mesh::Cell& cell)
{
  const auto& cell_mapping = sdm_.GetCellMapping(cell);
  const size_t num_nodes = cell_mapping.NumNodes();
  const auto& field_data = field_.FieldVectorRead();
  std::vector<double> shape_values(num_nodes, 0.0);

  const auto& qp_indices =
    fem_data_.cell_data_.qp_data_.QuadraturePointIndices();
  const auto& qp_xyz = fem_data_.cell_data_.qp_data_.QPointsXYZ();
  qp_values_.assign(qp_indices.size(), 0.0);
  for (uint32_t qp : qp_indices)
  {
    cell_mapping.ShapeValues(qp_xyz[qp], shape_values);
    for (size_t i = 0; i < num_nodes; ++i)
    {
      cint64_t dof_id = sdm_.MapDOFLocal(cell, i, uk_man_, 0, 0);
      qp_values_[qp] += shape_values[i] * field_data[dof_id];
    } // for i
  }   // for qp
}

// ##################################################################
void FieldFunctionFEMCoupledField::ComputeFieldFaceQPValues()
{
  const auto& cell = *fem_data_.cell_data_.cell_ptr_;
  const auto& cell_mapping = sdm_.GetCellMapping(cell);

  if (&cell_mapping == &(*fem_data_.cell_data_.cell_mapping_ptr_))
    SharedSDMComputeFieldFaceQPValues(cell, cell_mapping);
  else
    NonSharedSDMComputeFieldFaceQPValues(cell);
}

void FieldFunctionFEMCoupledField::SharedSDMComputeFieldFaceQPValues(
  const chi_mesh::Cell& cell, const chi_math::CellMapping& cell_mapping)
{
  const size_t num_nodes = cell_mapping.NumNodes();
  const auto& field_data = field_.FieldVectorRead();
  const auto& shape_values = fem_data_.face_data_.qp_data_.ShapeValues();

  const auto& qp_indices =
    fem_data_.face_data_.qp_data_.QuadraturePointIndices();

  if (num_nodes == 1)
  {
    cint64_t dof_id = sdm_.MapDOFLocal(cell, 0, uk_man_, 0, 0);
    qp_values_.assign(qp_indices.size(), field_data[dof_id]);
    return;
  }

  qp_values_.assign(qp_indices.size(), 0.0);
  for (uint32_t qp : qp_indices)
  {
    for (size_t i = 0; i < num_nodes; ++i)
    {
      cint64_t dof_id = sdm_.MapDOFLocal(cell, i, uk_man_, 0, 0);
      qp_values_[qp] += shape_values[i][qp] * field_data[dof_id];
    } // for i
  }   // for qp
}

void FieldFunctionFEMCoupledField::NonSharedSDMComputeFieldFaceQPValues(
  const chi_mesh::Cell& cell)
{
  const auto& cell_mapping = sdm_.GetCellMapping(cell);
  const size_t num_nodes = cell_mapping.NumNodes();
  const auto& field_data = field_.FieldVectorRead();
  std::vector<double> shape_values(num_nodes, 0.0);

  const auto& qp_indices =
    fem_data_.face_data_.qp_data_.QuadraturePointIndices();

  if (num_nodes == 1)
  {
    cint64_t dof_id = sdm_.MapDOFLocal(cell, 0, uk_man_, 0, 0);
    qp_values_.assign(qp_indices.size(), field_data[dof_id]);
    return;
  }

  const auto& qp_xyz = fem_data_.face_data_.qp_data_.QPointsXYZ();
  qp_values_.assign(qp_indices.size(), 0.0);
  for (uint32_t qp : qp_indices)
  {
    cell_mapping.ShapeValues(qp_xyz[qp], shape_values);
    for (size_t i = 0; i < num_nodes; ++i)
    {
      cint64_t dof_id = sdm_.MapDOFLocal(cell, i, uk_man_, 0, 0);
      qp_values_[qp] += shape_values[i] * field_data[dof_id];
    } // for i
  }   // for qp
}

} // namespace chi_math