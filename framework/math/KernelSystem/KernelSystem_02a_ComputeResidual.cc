#include "KernelSystem.h"

#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/KernelSystem/FEMKernels/FEMKernel.h"
#include "math/KernelSystem/FEMBCs/FEMDirichletBC.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"

#define scint64_t(x) static_cast<int64_t>(x)
#define cint64_t const int64_t

namespace chi_math
{

/**Collective method for computing the system residual.*/
void KernelSystem::ComputeResidual(const ParallelVector& x, ParallelVector& r)
{
  if (verbosity_ >= 2) Chi::log.LogAll() << "Compute Residual " << std::endl;

  const auto& grid = primary_fields_container_->GetSystemCommonGrid();

  current_field_index_ = 0;
  for (const auto& field_info : *primary_fields_container_)
  {
    current_sdm_ = &field_info.field_->GetSpatialDiscretization();
    current_spatial_weight_function_ =
      current_sdm_->GetSpatialWeightingFunction();

    for (const auto& cell : grid.local_cells)
    {
      InitCellData(x, cell);

      auto kernels = SetupAndGetCellInternalKernels();
      auto bndry_conditions = GetCellBCKernels();

      const auto& cell_mapping = *cur_cell_data.cell_mapping_ptr_;
      const size_t num_nodes = cell_mapping.NumNodes();
      const auto& dof_map = cur_cell_data.dof_map_;

      const std::set<uint32_t> dirichlet_nodes = IdentifyLocalDirichletNodes();

      std::vector<double> cell_local_r(num_nodes, 0.0);

      // Apply domain kernels
      for (const auto& kernel : kernels)
        for (size_t i = 0; i < num_nodes; ++i)
          if (dirichlet_nodes.find(i) == dirichlet_nodes.end())
            cell_local_r[i] += kernel->ComputeLocalResidual(i);

      // Apply nodal BCs
      for (const auto& [f, bndry_condition] : bndry_conditions)
      {
        if (not bndry_condition->IsDirichlet()) continue;
        bndry_condition->PreComputeValues();
        const size_t face_num_nodes = cell_mapping.NumFaceNodes(f);
        for (size_t fi = 0; fi < face_num_nodes; ++fi)
        {
          const size_t i = cell_mapping.MapFaceNode(f, fi);

          cell_local_r[i] += bndry_condition->ComputeLocalResidual(f, i);
        }
      }

      // Apply integral BCs
      for (const auto& [f, bndry_condition] : bndry_conditions)
      {
        if (bndry_condition->IsDirichlet()) continue;
        SetupFaceIntegralBCKernel(f);

        bndry_condition->PreComputeValues();

        const size_t face_num_nodes = cell_mapping.NumFaceNodes(f);
        for (size_t fi = 0; fi < face_num_nodes; ++fi)
        {
          const size_t i = cell_mapping.MapFaceNode(f, fi);

          if (dirichlet_nodes.find(i) == dirichlet_nodes.end())
            cell_local_r[i] += bndry_condition->ComputeLocalResidual(f, i);
        }
      }

      // Contribute to main residual
      for (size_t i = 0; i < num_nodes; ++i)
        r.SetValue(dof_map[i], cell_local_r[i], VecOpType::ADD_VALUE);

    } // for cell

    ++current_field_index_;
  } // for field

  if (verbosity_ >= 2)
    Chi::log.LogAll() << "Compute Residual Done" << std::endl;

  r.Assemble();
}

} // namespace chi_math