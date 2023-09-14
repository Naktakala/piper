#include "FEMKernelSystem.h"

#include "math/SpatialDiscretization/spatial_discretization.h"
#include "math/SpatialDiscretization/FiniteElement/finite_element.h"
#include "math/KernelSystem/FEMKernels/FEMKernel.h"
#include "math/KernelSystem/FEMBCs/FEMBoundaryCondition.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"

#define scint64_t(x) static_cast<int64_t>(x)
#define cint64_t const int64_t

namespace chi_math
{
/**\brief Basic constructor for a FEMKernelSystem.
* This constructor also sorts kernels and BCs into convenient maps.*/
FEMKernelSystem::FEMKernelSystem(
  const SpatialDiscretization& sdm,
  const UnknownManager& uk_man,
  std::vector<chi_math::FEMKernelPtr>& volume_kernels,
  std::vector<chi_math::FEMBoundaryConditionPtr>& boundary_conditions)
  : KernelSystem(scint64_t(sdm.GetNumLocalDOFs(uk_man)),
                 scint64_t(sdm.GetNumGlobalDOFs(uk_man)),
                 sdm.GetGhostDOFIndices(uk_man)),
    sdm_(sdm),
    uk_man_(uk_man)
{
  const auto& grid = sdm_.Grid();

  //======================================== Map mat-ids to kernels
  std::set<int> material_ids;
  for (const auto& cell : grid.local_cells)
    material_ids.insert(cell.material_id_);

  for (uint64_t global_id : grid.cells.GetGhostGlobalIDs())
    material_ids.insert(grid.cells[global_id].material_id_);

  for (int mat_id : material_ids)
  {
    auto& material_kernels = matid_2_volume_kernels_map_[mat_id];
    for (auto& kernel : volume_kernels)
    {
      const auto& mat_scope = kernel->GetMaterialIDScope();

      if (mat_scope.empty() or
          std::find(mat_scope.begin(), mat_scope.end(), mat_id) !=
            mat_scope.end())
        material_kernels.push_back(kernel);
    }

    ChiLogicalErrorIf(material_kernels.empty(),
                      "Material " + std::to_string(mat_id) +
                        " does not have a kernel assigned.");
  }

  //======================================== Map boundary ids to BCs
  const auto& bndry_id_map = grid.GetBoundaryIDMap();

  for (const auto& [bid, bname] : bndry_id_map)
  {
    bool assigned = false;
    for (auto& bc : boundary_conditions)
    {
      const auto& bndry_scope = bc->GetBoundaryScope();
      if (std::find(bndry_scope.begin(), bndry_scope.end(), bname) !=
          bndry_scope.end())
      {
        ChiLogicalErrorIf(
          bid_2_BCKernel_map_.count(bid) != 0,
          "More than one boundary condition specified on boundary \"" + bname +
            "\".");

        if (not assigned)
        {
          bid_2_BCKernel_map_[bid] = bc;
          assigned = true;
        }
        else
          ChiLogicalError("Boundary found with multiple BCs assignment");
      }
    } // for bc
  } // for item in bndry map
}

/**Returns the Spatial Discretization Method (SDM).*/
const SpatialDiscretization& FEMKernelSystem::SDM() const { return sdm_; }
/**Returns the nodal unknown structure of each variable set.*/
const UnknownManager& FEMKernelSystem::UnknownStructure() const
{
  return uk_man_;
}

/**Obtains the kernels associated with a material.*/
const std::vector<FEMKernelPtr>& FEMKernelSystem::GetMaterialKernels(int mat_id)
{
  auto iter = matid_2_volume_kernels_map_.find(mat_id);

  ChiLogicalErrorIf(iter == matid_2_volume_kernels_map_.end(),
                    "No kernel for material id " + std::to_string(mat_id));

  const auto& kernels = iter->second;

  ChiLogicalErrorIf(kernels.empty(),
                    "No kernel for material id " + std::to_string(mat_id));

  return kernels;
}

} // namespace chi_math