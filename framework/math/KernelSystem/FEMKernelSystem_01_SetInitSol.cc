#include "FEMKernelSystem.h"

#include "math/SpatialDiscretization/spatial_discretization.h"
#include "math/KernelSystem/FEMBCs/FEMDirichletBC.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#define cint64_t const int64_t

namespace chi_math
{
void FEMKernelSystem::SetInitialSolution()
{
  const auto& grid = sdm_.Grid();

  auto count_vector = solution_vector_;
  solution_vector_.Set(0.0);
  count_vector.Set(0.0);

  for (const auto& cell : grid.local_cells)
  {
    const size_t num_faces = cell.faces_.size();
    const auto& cell_mapping = sdm_.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto node_locations = cell_mapping.GetNodeLocations();

    //====================================== Pack cell-local solution
    std::vector<double> local_x(num_nodes, 0.0);
    for (size_t i = 0; i < num_nodes; ++i)
    {
      cint64_t dof_id = sdm_.MapDOFLocal(cell, i, uk_man_, 0, 0);
      local_x[i] = solution_vector_[dof_id];
    }

    //====================================== Compute cell-local solution
    std::vector<double> local_x_out(num_nodes, 0.0);
    std::vector<double> local_c_out(num_nodes, 0.0);
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces_[f];
      if (face.has_neighbor_) continue;

      auto bc_it = bid_2_boundary_conditions_map_.find(face.neighbor_id_);

      if (bc_it == bid_2_boundary_conditions_map_.end())
        continue; // Default natural BC

      auto bndry_condition = bc_it->second;

      const size_t face_num_nodes = cell_mapping.NumFaceNodes(f);
      if (bndry_condition->IsDirichlet())
      {
        auto dirichlet_condition =
          std::static_pointer_cast<FEMDirichletBC>(
            bndry_condition);

        if (dirichlet_condition->AllowApplyBeforeSolve())
        {
          for (size_t fi = 0; fi < face_num_nodes; ++fi)
          {
            const size_t i = cell_mapping.MapFaceNode(f, fi);

            local_x_out[i] += dirichlet_condition->BCValue();
            local_c_out[i] += 1.0;
          } // for fi
        }   // if preset allowed
      }     // if dirichlet
    }       // for f

    //====================================== Contribute to main residual
    for (size_t i = 0; i < num_nodes; ++i)
    {
      cint64_t dof_id = sdm_.MapDOF(cell, i, uk_man_, 0, 0);

      solution_vector_.SetValue(dof_id, local_x_out[i], VecOpType::ADD_VALUE);
      count_vector.SetValue(dof_id, local_c_out[i], VecOpType::ADD_VALUE);
    }
  } // for cell

  solution_vector_.Assemble();
  count_vector.Assemble();

  auto& sol_raw = solution_vector_.RawValues();
  auto& cnt_raw = count_vector.RawValues();

  for (size_t i = 0; i < solution_vector_.LocalSize(); ++i)
    if (cnt_raw[i] > 0.0) sol_raw[i] /= cnt_raw[i];

  solution_vector_.CommunicateGhostEntries();
}
} // namespace chi_math