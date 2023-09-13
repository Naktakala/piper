#include "FEMKernelSystem.h"

#include "math/SpatialDiscretization/spatial_discretization.h"
#include "math/KernelSystem/FEMKernels/FEMKernel.h"
#include "math/KernelSystem/FEMBCs/FEMDirichletBC.h"
#include "math/ParallelMatrix/ParallelMatrix.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"

#define scint64_t(x) static_cast<int64_t>(x)
#define cint64_t const int64_t

namespace chi_math
{

// ##################################################################
void FEMKernelSystem::InitCellKernelData(const chi_mesh::Cell& cell)
{
  cur_cell_data.cur_cell_ptr_ = &cell;

  const auto& cell_mapping = sdm_.GetCellMapping(cell);
  const size_t num_nodes = cell_mapping.NumNodes();

  cur_cell_data.cell_mapping_ptr_ = &cell_mapping;
  cur_cell_data.node_id_sets_ = sdm_.MakeCellInternalAndBndryNodeIDs(cell);
  cur_cell_data.node_locations_ = cell_mapping.GetNodeLocations();

  auto& dof_map = cur_cell_data.dof_map_;
  dof_map.assign(num_nodes, 0);
  for (size_t i = 0; i < num_nodes; ++i)
    dof_map[i] = sdm_.MapDOF(cell, i, uk_man_, 0, 0);

  std::vector<double> local_x(num_nodes, 0.0);
  for (size_t i = 0; i < num_nodes; ++i)
  {
    cint64_t dof_id = sdm_.MapDOFLocal(cell, i, uk_man_, 0, 0);
    local_x[i] = solution_vector_[dof_id];
  }
  cur_cell_data.local_x_ = std::move(local_x);
}

// ##################################################################
std::vector<std::shared_ptr<FEMKernel>>
FEMKernelSystem::SetupCellInternalKernels(const chi_mesh::Cell& cell)
{
  const auto& cell_mapping = *cur_cell_data.cell_mapping_ptr_;

  cur_cell_data.cell_qp_data_ =
    std::make_unique<CellQPData>(cell_mapping.MakeVolumeQuadraturePointData());

  const auto& kernels = GetMaterialKernels(cell.material_id_);
  SetupInternalKernelsRefData(kernels);

  return kernels;
}

// ##################################################################
std::vector<std::pair<size_t, FEMBoundaryConditionPtr>>
FEMKernelSystem::SetupCellBCKernels(const chi_mesh::Cell& cell)
{
  std::vector<std::pair<size_t, FEMBoundaryConditionPtr>> bndry_conditions;
  const size_t num_faces = cell.faces_.size();
  for (size_t f = 0; f < num_faces; ++f)
  {
    const auto& face = cell.faces_[f];
    if (face.has_neighbor_) continue;

    auto bc_it = bid_2_boundary_conditions_map_.find(face.neighbor_id_);
    if (bc_it == bid_2_boundary_conditions_map_.end())
      continue; // Default natural BC

    cur_cell_data.current_face_index_ = f;
    auto bndry_condition = bc_it->second;
    SetupBoundaryConditionRefData(*bndry_condition);
    bndry_conditions.emplace_back(f, bndry_condition);
  } // for face

  return bndry_conditions;
}

// ##################################################################
void FEMKernelSystem::SetupInternalKernelsRefData(
  const std::vector<FEMKernelPtr>& kernels)
{
  const auto& cell_mapping = *cur_cell_data.cell_mapping_ptr_;
  const auto& local_x = cur_cell_data.local_x_;
  const size_t num_nodes = cell_mapping.NumNodes();

  //====================================== Compute variable values on each qp
  const auto& qp_data = *cur_cell_data.cell_qp_data_;
  const auto& qp_indices = qp_data.QuadraturePointIndices();
  const size_t num_qps = qp_indices.size();

  const auto& shape_values = qp_data.ShapeValues();
  const auto& shape_grad_values = qp_data.ShapeGradValues();

  std::vector<double> var_values(num_qps, 0.0);
  std::vector<chi_mesh::Vector3> var_grad_values(num_qps);

  for (size_t j = 0; j < num_nodes; ++j)
    for (uint32_t qp : qp_indices)
    {
      var_values[qp] += shape_values[j][qp] * local_x[j];
      var_grad_values[qp] += shape_grad_values[j][qp] * local_x[j];
    }

  cur_cell_data.cell_var_values_ = std::move(var_values);
  cur_cell_data.cell_var_grad_values_ = std::move(var_grad_values);

  //====================================== Set reference data
  auto ref_data_ptr =
    std::make_shared<FEMKernelRefData>(*cur_cell_data.cell_qp_data_,
                                       cur_cell_data.cell_var_values_,
                                       cur_cell_data.cell_var_grad_values_);

  for (const auto& kernel : kernels)
    kernel->SetReferenceData(ref_data_ptr);
}

// ##################################################################
void FEMKernelSystem::SetupBoundaryConditionRefData(
  chi_math::FEMBoundaryCondition& bndry_condition)
{
  const auto& cell_mapping = *cur_cell_data.cell_mapping_ptr_;
  const size_t f = cur_cell_data.current_face_index_;
  const auto& local_x = cur_cell_data.local_x_;

  const size_t num_nodes = cell_mapping.NumNodes();

  //====================================== For non-dirichlet BCs we
  //                                       compute variable values on each qp
  if (not bndry_condition.IsDirichlet())
  {
    cur_cell_data.face_qp_data_ =
      std::make_unique<FaceQPData>(cell_mapping.MakeFaceQuadraturePointData(f));
    const auto& face_qp_data = *cur_cell_data.face_qp_data_;
    const auto& face_qp_indices = face_qp_data.QuadraturePointIndices();
    const size_t face_num_qps = face_qp_indices.size();

    const auto& face_shape_values = face_qp_data.ShapeValues();
    const auto& face_shape_grad_values = face_qp_data.ShapeGradValues();

    std::vector<double> bc_var_values(face_num_qps, 0.0);
    std::vector<chi_mesh::Vector3> bc_var_grad_values(face_num_qps);

    for (size_t j = 0; j < num_nodes; ++j)
      for (uint32_t qp : face_qp_indices)
      {
        bc_var_values[qp] += face_shape_values[j][qp] * local_x[j];
        bc_var_grad_values[qp] += face_shape_grad_values[j][qp] * local_x[j];
      }

    cur_cell_data.face_var_values_ = std::move(bc_var_values);
    cur_cell_data.face_var_grad_values_ = std::move(bc_var_grad_values);
  }
  else
  {
    cur_cell_data.face_qp_data_ = std::make_unique<FaceQPData>();
    cur_cell_data.face_qp_data_->InitializeEmpty();
  }

  //====================================== Set reference data
  auto ref_bc_data_ptr =
    std::make_shared<FEMBCRefData>(*cur_cell_data.face_qp_data_,
                                   cur_cell_data.face_var_values_,
                                   cur_cell_data.face_var_grad_values_,
                                   cur_cell_data.local_x_,
                                   cur_cell_data.node_locations_);

  bndry_condition.SetReferenceData(ref_bc_data_ptr);
}

} // namespace chi_math