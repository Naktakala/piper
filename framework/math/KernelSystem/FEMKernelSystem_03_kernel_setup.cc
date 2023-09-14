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
/**Initializes cell data prior kernel and BC setup.*/
void FEMKernelSystem::InitCellData(const chi_mesh::Cell& cell)
{
  const auto& cell_mapping = sdm_.GetCellMapping(cell);
  const size_t num_nodes = cell_mapping.NumNodes();

  cur_cell_data.cell_mapping_ptr_ = &cell_mapping;
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
/**Prepares all the necessary data for internal kernels.*/
std::vector<std::shared_ptr<FEMKernel>>
FEMKernelSystem::SetupCellInternalKernels(const chi_mesh::Cell& cell)
{
  const auto& cell_mapping = *cur_cell_data.cell_mapping_ptr_;

  auto cell_qp_data_ptr = cell_mapping.NewVolumeQuadraturePointData();

  const auto& kernels = GetMaterialKernels(cell.material_id_);

  const auto& local_x = cur_cell_data.local_x_;
  const size_t num_nodes = cell_mapping.NumNodes();

  //====================================== Compute variable values on each qp
  const auto& qp_data = *cell_qp_data_ptr;
  const auto& qp_indices = qp_data.QuadraturePointIndices();
  const size_t num_qps = qp_indices.size();

  const auto& shape_values = qp_data.ShapeValues();
  const auto& shape_grad_values = qp_data.ShapeGradValues();

  std::vector<double> cell_var_qp_values(num_qps, 0.0);
  std::vector<chi_mesh::Vector3> cell_var_grad_qp_values(num_qps);

  for (size_t j = 0; j < num_nodes; ++j)
    for (uint32_t qp : qp_indices)
    {
      cell_var_qp_values[qp] += shape_values[j][qp] * local_x[j];
      cell_var_grad_qp_values[qp] += shape_grad_values[j][qp] * local_x[j];
    }

  //====================================== Set reference data
  auto ref_data_ptr = std::make_shared<FEMKernelRefData>(
    cell_qp_data_ptr, cell_var_qp_values, cell_var_grad_qp_values);

  for (const auto& kernel : kernels)
    kernel->SetReferenceData(ref_data_ptr);

  return kernels;
}

// ##################################################################
/**Prepares all the necessary data for boundary kernels.*/
std::vector<std::pair<size_t, FEMBoundaryConditionPtr>>
FEMKernelSystem::SetupCellBCKernels(const chi_mesh::Cell& cell)
{
  std::vector<std::pair<size_t, FEMBoundaryConditionPtr>> bndry_conditions;

  const auto& cell_mapping = *cur_cell_data.cell_mapping_ptr_;

  const auto& local_x = cur_cell_data.local_x_;

  const size_t num_nodes = cell_mapping.NumNodes();

  const size_t num_faces = cell.faces_.size();
  for (size_t f = 0; f < num_faces; ++f)
  {
    const auto& face = cell.faces_[f];
    if (face.has_neighbor_) continue;

    auto bc_it = bid_2_BCKernel_map_.find(face.neighbor_id_);

    if (bc_it == bid_2_BCKernel_map_.end()) continue; // Natural BC

    auto bndry_condition = bc_it->second;
    bndry_conditions.emplace_back(f, bndry_condition);

    std::shared_ptr<FaceQPData> face_qp_data_ptr = nullptr;

    std::vector<double> bc_var_qp_values;
    std::vector<chi_mesh::Vector3> bc_var_grad_qp_values;

    //====================================== For non-dirichlet BCs we compute
    //                                       variable values and gradients on
    //                                       each qp
    if (not bndry_condition->IsDirichlet())
    {
      face_qp_data_ptr = cell_mapping.NewFaceQuadraturePointData(f);
      const auto& face_qp_data = *face_qp_data_ptr;
      const auto& face_qp_indices = face_qp_data.QuadraturePointIndices();
      const size_t face_num_qps = face_qp_indices.size();

      const auto& face_shape_values = face_qp_data.ShapeValues();
      const auto& face_shape_grad_values = face_qp_data.ShapeGradValues();

      bc_var_qp_values.assign(face_num_qps, 0.0);
      bc_var_grad_qp_values.assign(face_num_qps, {});

      for (size_t j = 0; j < num_nodes; ++j)
        for (uint32_t qp : face_qp_indices)
        {
          bc_var_qp_values[qp] += face_shape_values[j][qp] * local_x[j];
          bc_var_grad_qp_values[qp] += face_shape_grad_values[j][qp] * local_x[j];
        }
    }
    //====================================== Dirichlet BCs do not require
    //                                       quadrature point info, so we skip
    //                                       it
    else
    {
      face_qp_data_ptr = std::make_unique<FaceQPData>();
      face_qp_data_ptr->InitializeEmpty();
    }

    //====================================== Set reference data
    auto ref_bc_data_ptr =
      std::make_shared<FEMBCRefData>(face_qp_data_ptr,
                                     bc_var_qp_values,
                                     bc_var_grad_qp_values,
                                     cur_cell_data.local_x_,
                                     cur_cell_data.node_locations_);

    bndry_condition->SetFaceReferenceData(f, ref_bc_data_ptr);
  } // for face

  return bndry_conditions;
}

} // namespace chi_math