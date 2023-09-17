#include "FEMKernelSystem.h"

#include "math/SpatialDiscretization/spatial_discretization.h"
#include "math/SpatialDiscretization/FiniteElement/finite_element.h"
#include "math/KernelSystem/FEMKernels/FEMKernel.h"
#include "math/KernelSystem/FEMBCs/FEMBoundaryCondition.h"
#include "math/ParallelMatrix/ParallelMatrix.h"

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
  std::vector<chi_math::FEMBoundaryConditionPtr>& boundary_conditions,
  TimeID oldest_time_id /*=TimeID::T_PLUS_1*/)
  : KernelSystem(scint64_t(sdm.GetNumLocalDOFs(uk_man)),
                 scint64_t(sdm.GetNumGlobalDOFs(uk_man)),
                 sdm.GetGhostDOFIndices(uk_man),
                 oldest_time_id),
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
  }   // for item in bndry map
}

/**Returns the Spatial Discretization Method (SDM).*/
const SpatialDiscretization& FEMKernelSystem::SDM() const { return sdm_; }
/**Returns the nodal unknown structure of each variable set.*/
const UnknownManager& FEMKernelSystem::UnknownStructure() const
{
  return uk_man_;
}

/**Obtains the kernels associated with a material.*/
std::vector<FEMKernelPtr> FEMKernelSystem::GetMaterialKernels(int mat_id)
{
  auto iter = matid_2_volume_kernels_map_.find(mat_id);

  ChiLogicalErrorIf(iter == matid_2_volume_kernels_map_.end(),
                    "No kernel for material id " + std::to_string(mat_id));

  const auto& kernels = iter->second;

  ChiLogicalErrorIf(kernels.empty(),
                    "No kernel for material id " + std::to_string(mat_id));

  std::vector<std::shared_ptr<FEMKernel>> filtered_kernels;
  const bool time_terms_active = QueryTermsActive(EqTermScope::TIME_TERMS);
  const bool domain_terms_active = QueryTermsActive(EqTermScope::DOMAIN_TERMS);

  for (const auto& kernel : kernels)
  {
    if (kernel->IsTimeKernel() and time_terms_active)
      filtered_kernels.push_back(kernel);
    if (not kernel->IsTimeKernel() and domain_terms_active)
      filtered_kernels.push_back(kernel);
  }

  if (time_terms_active and filtered_kernels.empty())
    ChiLogicalError("No time kernels in system");

  return filtered_kernels;
}

/**Obtains the boundary kernel associated with a boundary id.*/
FEMBoundaryConditionPtr
FEMKernelSystem::GetBoundaryCondition(uint64_t boundary_id)
{
  if (not(EquationTermsScope() & EqTermScope::BOUNDARY_TERMS)) return nullptr;

  auto bc_it = bid_2_BCKernel_map_.find(boundary_id);

  if (bc_it == bid_2_BCKernel_map_.end()) return nullptr;

  return bc_it->second;
}

ParallelMatrixSparsityPattern
FEMKernelSystem::BuildMatrixSparsityPattern() const
{
  auto& sdm = SDM();
  auto& uk_man = UnknownStructure();

  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm.BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, uk_man);

  return {nodal_nnz_in_diag, nodal_nnz_off_diag};
}

} // namespace chi_math