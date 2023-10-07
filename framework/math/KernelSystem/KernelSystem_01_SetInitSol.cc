#include "KernelSystem.h"

#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/KernelSystem/FEMBCs/FEMDirichletBC.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"

#define cint64_t const int64_t

namespace chi_math
{
/**Used to create an initial guess. Mostly applies Dirichlet BCs.*/
void KernelSystem::SetInitialSolution()
{
  if (verbosity_ >= 2) Chi::log.LogAll() << "SetInitialSolution " << std::endl;

  auto& x = *main_solution_vector_;

  for (auto& field_info : *primary_fields_container_)
    x.BlockSet(field_info.field_->FieldVector(),
                           field_info.local_offset_,
                           field_info.num_local_dofs_);

  x.CommunicateGhostEntries();

  auto dirichlet_vector_ptr = main_solution_vector_->MakeClone();
  auto count_vector_ptr = main_solution_vector_->MakeClone();
  auto& dirichlet_vector = *dirichlet_vector_ptr;
  auto& count_vector = *count_vector_ptr;

  current_field_index_ = 0;
  for (const auto& field_info : *primary_fields_container_)
  {
    const auto& field = field_info.field_;
    const auto& sdm = field->GetSpatialDiscretization();

    const auto& bid2bc_map = varname_comp_2_bid2bc_map_.at(
      {field->TextName(), current_field_component_});

    const auto& grid = sdm.Grid();
    for (const auto& cell : grid.local_cells)
    {
      const size_t num_faces = cell.faces_.size();
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();
      const auto node_locations = cell_mapping.GetNodeLocations();

      cur_cell_data.cell_mapping_ptr_ = &cell_mapping;

      VecDbl local_x(num_nodes, 0.0);

      for (size_t i = 0; i < num_nodes; ++i)
      {
        cint64_t dof_id =
          primary_fields_container_->MapDOFLocal(cell, i, field_info, 0);

        local_x[i] = x[dof_id];
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

        auto bc_it = bid2bc_map.find(face.neighbor_id_);

        if (bc_it == bid2bc_map.end()) continue; // Default natural BC

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
        cint64_t dof_id =
          primary_fields_container_->MapDOF(cell, i, field_info, 0);

        if (dirichlet_nodes.find(i) != dirichlet_nodes.end())
        {
          dirichlet_vector.SetValue(
            dof_id, local_x_out[i], VecOpType::ADD_VALUE);
          count_vector.SetValue(dof_id, local_c_out[i], VecOpType::ADD_VALUE);
        }
        else
        {
          dirichlet_vector.SetValue(dof_id, local_x[i], VecOpType::ADD_VALUE);
          count_vector.SetValue(dof_id, 1.0, VecOpType::ADD_VALUE);
        }
      }
    } // for cell
    ++current_field_index_;
  } // for field

  dirichlet_vector.Assemble();
  count_vector.Assemble();

  double* sol_raw = dirichlet_vector.Data();
  const double* cnt_raw = count_vector.Data();

  for (size_t i = 0; i < dirichlet_vector.LocalSize(); ++i)
    if (std::fabs(cnt_raw[i]) > 1.0e-12) sol_raw[i] /= cnt_raw[i];

  main_solution_vector_->CopyLocalValues(dirichlet_vector);

  main_solution_vector_->CommunicateGhostEntries();
  if (num_solution_histories_ > 0)
    for (size_t t = 0; t < num_solution_histories_; ++t)
      old_solution_vectors_[t] = main_solution_vector_->MakeCopy();
}
} // namespace chi_math