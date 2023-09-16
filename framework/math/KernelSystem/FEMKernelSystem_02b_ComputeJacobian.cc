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

/**Collective method for computing the system Jacobian-matrix.*/
void FEMKernelSystem::ComputeJacobian(const GhostedParallelVector& x,
                                      ParallelMatrix& J)
{
  SetActiveKernels(TIME_KERNELS | STD_KERNELS | BNDRY_KERNELS);
  const auto& grid = sdm_.Grid();

  for (const auto& cell : grid.local_cells)
  {
    InitCellData(x, cell);

    auto kernels = SetupCellInternalKernels(cell);
    auto bndry_conditions = SetupCellBCKernels(cell);

    const auto& cell_mapping = *cur_cell_data.cell_mapping_ptr_;
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto& dof_map = cur_cell_data.dof_map_;

    MatDbl cell_local_J(num_nodes, VecDbl(num_nodes, 0.0));

    // Apply nodal BCs
    const std::set<uint32_t> dirichlet_nodes =
      IdentifyLocalDirichletNodes(cell);

    for (size_t i : dirichlet_nodes)
      cell_local_J[i][i] += 1.0;

    // Apply integral BCs
    for (const auto& [f, bndry_condition] : bndry_conditions)
    {
      if (bndry_condition->IsDirichlet()) continue;
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

    for (const auto& kernel : kernels)
      for (size_t i = 0; i < num_nodes; ++i)
        if (dirichlet_nodes.find(i) == dirichlet_nodes.end())
          for (size_t j = 0; j < num_nodes; ++j)
            cell_local_J[i][j] += kernel->ComputeLocalJacobian(i, j);

    // Contribute to main jacobian
    for (size_t i = 0; i < num_nodes; ++i)
      for (size_t j = 0; j < num_nodes; ++j)
        J.AddValue(dof_map[i], dof_map[j], cell_local_J[i][j]);
  } // for cell

  SetActiveKernels(STD_KERNELS | BNDRY_KERNELS);
}

} // namespace chi_math