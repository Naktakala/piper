#include "FEMKernelSystem.h"

#include "math/SpatialDiscretization/spatial_discretization.h"
#include "math/KernelSystem/FEMBCs/FEMDirichletBC.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"

#define cint64_t const int64_t

namespace chi_math
{
/**Used to create an initial guess. Mostly applies Dirichlet BCs.*/
void FEMKernelSystem::SetInitialSolution()
{
  const auto& grid = sdm_.Grid();

  auto dirichlet_vector_ptr = main_solution_vector_->MakeNewVector();
  auto count_vector_ptr = main_solution_vector_->MakeNewVector();
  auto& dirichlet_vector = *dirichlet_vector_ptr;
  auto& count_vector = *count_vector_ptr;
  dirichlet_vector.Set(0.0);
  count_vector.Set(0.0);

  for (const auto& cell : grid.local_cells)
  {
    const size_t num_faces = cell.faces_.size();
    const auto& cell_mapping = sdm_.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto node_locations = cell_mapping.GetNodeLocations();

    cur_cell_data.cell_mapping_ptr_ = &cell_mapping;

    VecDbl local_x(num_nodes, 0.0);

    for (size_t i = 0; i < num_nodes; ++i)
    {
      cint64_t dof_id = sdm_.MapDOFLocal(cell, i, uk_man_, 0, 0);
      local_x[i] = main_solution_vector_->operator[](dof_id);
    }

    const std::set<uint32_t> dirichlet_nodes =
      IdentifyLocalDirichletNodes(cell);

    //====================================== Compute cell-local solution
    std::vector<double> local_x_out(num_nodes, 0.0);
    std::vector<double> local_c_out(num_nodes, 0.0);
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces_[f];
      if (face.has_neighbor_) continue;

      auto bc_it = bid_2_BCKernel_map_.find(face.neighbor_id_);

      if (bc_it == bid_2_BCKernel_map_.end()) continue; // Default natural BC

      auto bndry_condition = bc_it->second;

      const size_t face_num_nodes = cell_mapping.NumFaceNodes(f);
      if (bndry_condition->IsDirichlet())
      {
        auto dirichlet_condition =
          std::static_pointer_cast<FEMDirichletBC>(bndry_condition);

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

      if (dirichlet_nodes.find(i) != dirichlet_nodes.end())
      {
        dirichlet_vector.SetValue(dof_id, local_x_out[i], VecOpType::ADD_VALUE);
        count_vector.SetValue(dof_id, local_c_out[i], VecOpType::ADD_VALUE);
      }
      else
      {
        dirichlet_vector.SetValue(dof_id, local_x[i], VecOpType::ADD_VALUE);
        count_vector.SetValue(dof_id, 1.0, VecOpType::ADD_VALUE);
      }
    }
  } // for cell

  dirichlet_vector.Assemble();
  count_vector.Assemble();

  auto& sol_raw = dirichlet_vector.RawValues();
  auto& cnt_raw = count_vector.RawValues();

  for (size_t i = 0; i < dirichlet_vector.LocalSize(); ++i)
    if (std::fabs(cnt_raw[i]) > 1.0e-12) sol_raw[i] /= cnt_raw[i];

  dirichlet_vector.CommunicateGhostEntries();

  main_solution_vector_->CopyValues(dirichlet_vector);
}
} // namespace chi_math