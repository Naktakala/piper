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
void FEMKernelSystem::InitCellData(const GhostedParallelVector& x,
                                   const chi_mesh::Cell& cell)
{
  const auto& cell_mapping = sdm_.GetCellMapping(cell);
  const size_t num_nodes = cell_mapping.NumNodes();

  cur_cell_data.cell_mapping_ptr_ = &cell_mapping;
  cur_cell_data.node_locations_ = cell_mapping.GetNodeLocations();

  auto& dof_map = cur_cell_data.dof_map_;
  dof_map.assign(num_nodes, 0);
  for (size_t i = 0; i < num_nodes; ++i)
    dof_map[i] = sdm_.MapDOF(cell, i, uk_man_, 0, 0);

  const size_t max_t = AreTimeKernelsActive() ? num_old_blocks_ : 0;

  VecDbl local_x(num_nodes, 0.0);
  MatDbl old_local_x(max_t, VecDbl(num_nodes, 0.0));

  for (size_t i = 0; i < num_nodes; ++i)
  {
    cint64_t dof_id = sdm_.MapDOFLocal(cell, i, uk_man_, 0, 0);
    local_x[i] = x[dof_id];
    for (size_t t = 0; t < max_t; ++t)
      old_local_x[t][i] = old_solution_vectors_[t]->operator[](dof_id);
  }

  cur_cell_data.local_x_ = std::move(local_x);
  cur_cell_data.old_local_x_ = std::move(old_local_x);
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
  const auto& old_local_x = cur_cell_data.old_local_x_;
  const size_t num_nodes = cell_mapping.NumNodes();

  //====================================== Compute variable values on each qp
  const auto& qp_data = *cell_qp_data_ptr;
  const auto& qp_indices = qp_data.QuadraturePointIndices();
  const size_t num_qps = qp_indices.size();

  const auto& shape_values = qp_data.ShapeValues();
  const auto& shape_grad_values = qp_data.ShapeGradValues();

  VecDbl var_qp_values(num_qps, 0.0);
  VecVec3 var_grad_qp_values(num_qps);

  for (size_t j = 0; j < num_nodes; ++j)
    for (uint32_t qp : qp_indices)
    {
      var_qp_values[qp] += shape_values[j][qp] * local_x[j];
      var_grad_qp_values[qp] += shape_grad_values[j][qp] * local_x[j];
    }

  //======================================== Compute old time values if needed
  const size_t max_t = AreTimeKernelsActive() ? num_old_blocks_ : 0;
  MatDbl old_var_qp_values(max_t, VecDbl(num_qps, 0.0));
  if (AreTimeKernelsActive())
    for (size_t t = 0; t < max_t; ++t)
    {
      auto& old_var_qp_values_t = old_var_qp_values[t];

      for (size_t j = 0; j < num_nodes; ++j)
        for (uint32_t qp : qp_indices)
          old_var_qp_values_t[qp] += shape_values[j][qp] * old_local_x[t][j];
    }

  //====================================== Set reference data
  auto ref_data_ptr = std::make_shared<FEMKernelRefData>(time_data_,
                                                         cell_qp_data_ptr,
                                                         var_qp_values,
                                                         var_grad_qp_values,
                                                         old_var_qp_values);

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

    auto bndry_condition = GetBoundaryCondition(face.neighbor_id_);

    if (not bndry_condition) continue;

    bndry_conditions.emplace_back(f, bndry_condition);

    std::shared_ptr<FaceQPData> face_qp_data_ptr = nullptr;

    VecDbl var_qp_values;
    VecVec3 var_grad_qp_values;

    //====================================== For non-dirichlet BCs we compute
    //                                       variable values and gradients on
    //                                       each qp
    if (not bndry_condition->IsDirichlet())
    {
      face_qp_data_ptr = cell_mapping.NewFaceQuadraturePointData(f);
      const auto& face_qp_data = *face_qp_data_ptr;
      const auto& qp_indices = face_qp_data.QuadraturePointIndices();
      const size_t num_qps = qp_indices.size();

      const auto& shape_values = face_qp_data.ShapeValues();
      const auto& shape_grad_values = face_qp_data.ShapeGradValues();

      var_qp_values.assign(num_qps, 0.0);
      var_grad_qp_values.assign(num_qps, {});

      for (size_t j = 0; j < num_nodes; ++j)
        for (uint32_t qp : qp_indices)
        {
          var_qp_values[qp] += shape_values[j][qp] * local_x[j];
          var_grad_qp_values[qp] += shape_grad_values[j][qp] * local_x[j];
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
      std::make_shared<FEMBCRefData>(time_data_,
                                     face_qp_data_ptr,
                                     var_qp_values,
                                     var_grad_qp_values,
                                     cur_cell_data.local_x_,
                                     cur_cell_data.node_locations_);

    bndry_condition->SetFaceReferenceData(f, ref_bc_data_ptr);
  } // for face

  return bndry_conditions;
}

/**Returns a set of dirichlet nodes by looking at the BCs applied on
  * faces. Does not get filtered by time status.*/
std::set<uint32_t>
FEMKernelSystem::IdentifyLocalDirichletNodes(const chi_mesh::Cell& cell) const
{
  const auto& cell_mapping = *cur_cell_data.cell_mapping_ptr_;

  std::set<uint32_t> dirichlet_nodes;
  const size_t num_faces = cell.faces_.size();
  for (size_t f = 0; f < num_faces; ++f)
  {
    const auto& face = cell.faces_[f];
    if (face.has_neighbor_) continue;

    auto bc_it = bid_2_BCKernel_map_.find(face.neighbor_id_);
    if (bc_it == bid_2_BCKernel_map_.end()) continue;

    const auto& bndry_condition = *bc_it->second;
    if (bndry_condition.IsDirichlet())
    {
      const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);

      for (size_t fi=0; fi<num_face_nodes; ++fi)
        dirichlet_nodes.insert(cell_mapping.MapFaceNode(f, fi));
    }//if dirichlet
  }// for face

  return dirichlet_nodes;
}

} // namespace chi_math