#include "FEMKernelSystem.h"

#include "math/SpatialDiscretization/spatial_discretization.h"
#include "math/KernelSystem/FEMKernels/FEMKernel.h"
#include "math/KernelSystem/FEMBCs/FEMDirichletBC.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"

#define scint64_t(x) static_cast<int64_t>(x)
#define cint64_t const int64_t

namespace chi_math
{

void FEMKernelSystem::ComputeResidual(ParallelVector& r)
{
  const auto& grid = sdm_.Grid();

  for (const auto& cell : grid.local_cells)
  {
    InitCellKernelData(cell);

    auto kernels = SetupCellInternalKernels(cell);
    auto bndry_conditions = SetupCellBCKernels(cell);

    const auto& internal_nodes = cur_cell_data.node_id_sets_.first;
    const auto& cell_mapping = *cur_cell_data.cell_mapping_ptr_;
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto& dof_map = cur_cell_data.dof_map_;

    std::vector<double> cell_local_r(num_nodes, 0.0);

    for (const auto& kernel : kernels)
      for (size_t i : internal_nodes)
        cell_local_r[i] += kernel->ComputeLocalResidual(i);

    for (const auto& [f, bndry_condition] : bndry_conditions)
    {
      const size_t face_num_nodes = cell_mapping.NumFaceNodes(f);
      for (size_t fi = 0; fi < face_num_nodes; ++fi)
      {
        const size_t i = cell_mapping.MapFaceNode(f, fi);

        cell_local_r[i] += bndry_condition->ComputeLocalResidual(i);
      }
    }

    // Contribute to main residual
    for (size_t i = 0; i < num_nodes; ++i)
      r.SetValue(dof_map[i], cell_local_r[i], VecOpType::ADD_VALUE);
  } // for cell

  r.Assemble();
}

} // namespace chi_math