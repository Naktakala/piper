#include "FEMKernelSystem.h"

#include "math/SpatialDiscretization/spatial_discretization.h"
#include "math/KernelSystem/FEMKernels/FEMKernel.h"
#include "math/KernelSystem/FEMBCs/FEMDirichletBC.h"
#include "math/ParallelMatrix/ParallelMatrix.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"

#define scint64_t(x) static_cast<int64_t>(x)
#define cint64_t const int64_t

namespace chi_math
{

// ##################################################################
/**Initializes cell data prior kernel and BC setup.*/
void FEMKernelSystem::InitCellData(const ParallelVector& x,
                                   const chi_mesh::Cell& cell)
{
  const auto& field_info = field_block_info_.at(current_field_index_);
  const auto& field = field_info.field_;
  const auto& sdm = field->SDM();
  const auto& uk_man = field->UnkManager();

  const auto& cell_mapping = sdm.GetCellMapping(cell);
  const size_t num_nodes = cell_mapping.NumNodes();

  cur_cell_data.cell_mapping_ptr_ = &cell_mapping;
  cur_cell_data.node_locations_ = std::move(cell_mapping.GetNodeLocations());

  auto& dof_map = cur_cell_data.dof_map_;
  dof_map.assign(num_nodes, 0);
  for (size_t i = 0; i < num_nodes; ++i)
  {
    cint64_t block_dof_id = sdm.MapDOF(cell, i, uk_man, 0, 0);
    dof_map[i] = MapBlockGlobalIDToSystem(field_info, block_dof_id);
  }

  const size_t max_t = EquationTermsScope() & EqTermScope::TIME_TERMS
                         ? num_solution_histories_
                         : 0;

  VecDbl local_x(num_nodes, 0.0);
  VecDbl local_x_dot(num_nodes, 0.0);
  const double dt = time_data_.dt_;

  for (size_t i = 0; i < num_nodes; ++i)
  {
    cint64_t block_dof_id = sdm.MapDOFLocal(cell, i, uk_man, 0, 0);
    cint64_t dof_id = MapBlockLocalIDToSystem(field_info, block_dof_id);

    local_x[i] = x[dof_id];

    if (max_t > 0)
      local_x_dot[i] = (local_x[i] - (*old_solution_vectors_[0])[dof_id]) / dt;
  }

  cur_cell_data.local_x_ = std::move(local_x);
  cur_cell_data.local_x_dot_ = std::move(local_x_dot);
}

// ##################################################################
/**Prepares all the necessary data for internal kernels.*/
std::vector<std::shared_ptr<FEMKernel>>
FEMKernelSystem::SetupAndGetCellInternalKernels(const chi_mesh::Cell& cell)
{
  const auto& cell_mapping = *cur_cell_data.cell_mapping_ptr_;

  cur_cell_data.qp_data_ =
    std::move(cell_mapping.MakeVolumeQuadraturePointData());

  const auto& kernels = GetMaterialKernels(cell.material_id_);

  const auto& local_x = cur_cell_data.local_x_;
  const auto& local_x_dot_ = cur_cell_data.local_x_dot_;
  const size_t num_nodes = cell_mapping.NumNodes();

  //====================================== Compute variable values on each qp
  const auto& qp_data = cur_cell_data.qp_data_;
  const auto& qp_indices = qp_data.QuadraturePointIndices();
  const size_t num_qps = qp_indices.size();

  const auto& shape_values = qp_data.ShapeValues();
  const auto& shape_grad_values = qp_data.ShapeGradValues();

  auto& var_qp_values = cur_cell_data.var_qp_values_;
  auto& var_grad_qp_values = cur_cell_data.var_grad_qp_values_;

  var_qp_values.assign(num_qps, 0.0);
  var_grad_qp_values.assign(num_qps, {});

  for (size_t j = 0; j < num_nodes; ++j)
    for (uint32_t qp : qp_indices)
    {
      var_qp_values[qp] += shape_values[j][qp] * local_x[j];
      var_grad_qp_values[qp] += shape_grad_values[j][qp] * local_x[j];
    }

  //======================================== Compute old time values if needed
  const bool time_kernels_active =
    EquationTermsScope() & EqTermScope::TIME_TERMS;

  const size_t max_t = time_kernels_active ? num_solution_histories_ : 0;

  auto& var_dot_qp_values = cur_cell_data.var_dot_qp_values_;
  var_dot_qp_values.assign(num_qps, 0.0);
  if (time_kernels_active and max_t > 0)
    for (size_t j = 0; j < num_nodes; ++j)
      for (uint32_t qp : qp_indices)
        var_dot_qp_values[qp] += shape_values[j][qp] * local_x_dot_[j];

  return kernels;
}

// ##################################################################
/**Prepares all the necessary data for boundary kernels.*/
std::vector<std::pair<size_t, FEMBoundaryConditionPtr>>
FEMKernelSystem::GetCellBCKernels(const chi_mesh::Cell& cell)
{
  std::vector<std::pair<size_t, FEMBoundaryConditionPtr>> bndry_conditions;

  const size_t num_faces = cell.faces_.size();
  for (size_t f = 0; f < num_faces; ++f)
  {
    const auto& face = cell.faces_[f];
    if (face.has_neighbor_) continue;

    auto bndry_condition = GetBoundaryCondition(face.neighbor_id_);

    if (not bndry_condition) continue;

    bndry_conditions.emplace_back(f, bndry_condition);
  } // for face

  return bndry_conditions;
}

void FEMKernelSystem::SetupFaceIntegralBCKernel(size_t face_index)
{
  const auto& cell_mapping = *cur_cell_data.cell_mapping_ptr_;

  const auto& local_x = cur_cell_data.local_x_;

  const size_t num_nodes = cell_mapping.NumNodes();

  cur_face_data.qp_data_ =
    std::move(cell_mapping.MakeFaceQuadraturePointData(face_index));

  const auto& face_qp_data = cur_face_data.qp_data_;
  const auto& qp_indices = face_qp_data.QuadraturePointIndices();
  const size_t num_qps = qp_indices.size();

  const auto& shape_values = face_qp_data.ShapeValues();
  const auto& shape_grad_values = face_qp_data.ShapeGradValues();

  auto& var_qp_values = cur_face_data.var_qp_values_;
  auto& var_grad_qp_values = cur_face_data.var_grad_qp_values_;

  var_qp_values.assign(num_qps, 0.0);
  var_grad_qp_values.assign(num_qps, {});

  for (size_t j = 0; j < num_nodes; ++j)
    for (uint32_t qp : qp_indices)
    {
      var_qp_values[qp] += shape_values[j][qp] * local_x[j];
      var_grad_qp_values[qp] += shape_grad_values[j][qp] * local_x[j];
    }
}

/**Returns a set of dirichlet nodes by looking at the BCs applied on
 * faces. Does not get filtered by time status.*/
std::set<uint32_t>
FEMKernelSystem::IdentifyLocalDirichletNodes(const chi_mesh::Cell& cell) const
{
  const auto& current_field = field_block_info_.at(current_field_index_).field_;
  const auto& field_name = current_field->TextName();

  const auto& cell_mapping = *cur_cell_data.cell_mapping_ptr_;

  const auto& bid2bc_map =
    varname_comp_2_bid2bc_map_.at({field_name, current_field_component_});

  std::set<uint32_t> dirichlet_nodes;
  const size_t num_faces = cell.faces_.size();
  for (size_t f = 0; f < num_faces; ++f)
  {
    const auto& face = cell.faces_[f];
    if (face.has_neighbor_) continue;

    auto bc_it = bid2bc_map.find(face.neighbor_id_);
    if (bc_it == bid2bc_map.end()) continue;

    const auto& bndry_condition = *bc_it->second;
    if (bndry_condition.IsDirichlet())
    {
      const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);

      for (size_t fi = 0; fi < num_face_nodes; ++fi)
        dirichlet_nodes.insert(cell_mapping.MapFaceNode(f, fi));
    } // if dirichlet
  }   // for face

  return dirichlet_nodes;
}

} // namespace chi_math