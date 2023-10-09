#include "KernelSystem.h"

#include "math/SpatialDiscretization/SpatialDiscretization.h"
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

/**Collective method for computing the system Jacobian-matrix.*/
void KernelSystem::ComputeJacobian(const ParallelVector& x, ParallelMatrix& J)
{
  if (verbosity_ >= 2)
    Chi::log.LogAll() << "Compute ComputeJacobian " << std::endl;

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

      MatDbl cell_local_J(num_nodes, VecDbl(num_nodes, 0.0));

      // Apply domain kernels
      for (const auto& kernel : kernels)
        for (size_t i = 0; i < num_nodes; ++i)
          if (dirichlet_nodes.find(i) == dirichlet_nodes.end())
            for (size_t j = 0; j < num_nodes; ++j)
              cell_local_J[i][j] += kernel->ComputeLocalJacobian(i, j);

      // Apply nodal BCs
      for (size_t i : dirichlet_nodes)
        cell_local_J[i][i] += 1.0;

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
            for (size_t j = 0; j < num_nodes; ++j)
              cell_local_J[i][j] +=
                bndry_condition->ComputeLocalJacobian(f, i, j);
        }
      }

      // Contribute to main jacobian
      for (size_t i = 0; i < num_nodes; ++i)
        for (size_t j = 0; j < num_nodes; ++j)
          J.AddValue(dof_map[i], dof_map[j], cell_local_J[i][j]);
    } // for cell
    ++current_field_index_;
  } // for field

  if (verbosity_ >= 2)
    Chi::log.LogAll() << "Compute ComputeJacobian Done" << std::endl;

  J.Assemble(/*final=*/true);
}

} // namespace chi_math