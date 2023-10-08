#include "KernelSystem.h"

#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/KernelSystem/FEMKernels/FEMKernel.h"
#include "math/KernelSystem/FEMBCs/FEMBoundaryCondition.h"
#include "math/KernelSystem/Coupling/FEMMaterialProperty.h"
#include "math/ParallelMatrix/ParallelMatrix.h"

#include "materials/MaterialPropertiesData.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

#define scint64_t(x) static_cast<int64_t>(x)
#define cint64_t const int64_t

namespace chi_math
{

RegisterChiObject(chi_math, KernelSystem);

chi::InputParameters KernelSystem::GetInputParameters()
{
  chi::InputParameters params = FieldEquationSystem::GetInputParameters();

  params.SetGeneralDescription("An equation system based on kernels.");
  params.SetDocGroup("doc_KernelSystem");

  params.AddRequiredParameterArray("kernels",
                                   "A list of kernel in parameter form.");
  params.AddRequiredParameterArray("bcs",
                                   "A list of boundary condition handles.");

  return params;
}

/**\brief Constructor for a KernelSystem.
 * This constructor also sorts kernels and BCs into convenient maps.*/
KernelSystem::KernelSystem(const chi::InputParameters& params)
  : FieldEquationSystem(params)
{
  const auto& volume_kernels_inputs = params.GetParam("kernels");
  const auto boundary_condition_inputs = params.GetParam("bcs");

  //======================================== Make reference data
  auto reference_data =
    std::make_shared<FEMKernelSystemData>(material_properties_data_,
                                          time_data_,
                                          *main_solution_vector_,
                                          cur_cell_data,
                                          cur_face_data);

  Chi::object_stack.push_back(reference_data);
  const size_t data_handle = Chi::object_stack.size() - 1;

  //======================================== Make all kernels
  const std::vector<FEMKernelPtr> kernels =
    MakeFEMKernels(volume_kernels_inputs, data_handle);

  //======================================== Make all BC-Kernels
  const std::vector<FEMBoundaryConditionPtr> boundary_conditions =
    MakeBCs(boundary_condition_inputs, data_handle);

  //======================================== Map mat-ids to kernels
  const auto& grid = primary_fields_container_->GetSystemCommonGrid();
  std::set<int> material_ids;
  for (const auto& cell : grid.local_cells)
    material_ids.insert(cell.material_id_);

  for (uint64_t global_id : grid.cells.GetGhostGlobalIDs())
    material_ids.insert(grid.cells[global_id].material_id_);

  for (int mat_id : material_ids)
  {
    auto& material_kernels = matid_2_volume_kernels_map_[mat_id];
    for (const auto& kernel_ptr : kernels) // making a copy
    {
      const auto& mat_scope = kernel_ptr->GetMaterialIDScope();

      if (mat_scope.empty() or
          std::find(mat_scope.begin(), mat_scope.end(), mat_id) !=
            mat_scope.end())
        material_kernels.push_back(kernel_ptr);
    }

    ChiLogicalErrorIf(material_kernels.empty(),
                      "Material " + std::to_string(mat_id) +
                        " does not have a kernel assigned.");
  }

  //======================================== Map boundary ids to BCs
  const auto& bndry_id_map = grid.GetBoundaryIDMap();

  for (const auto& [bid, bname] : bndry_id_map)
  {
    for (auto& bc : boundary_conditions)
    {
      const auto& bc_varname_comp = bc->ActiveVariableAndComponent();
      const auto& bndry_scope = bc->GetBoundaryScope();
      if (std::find(bndry_scope.begin(), bndry_scope.end(), bname) !=
          bndry_scope.end())
      {
        auto& bid2bc_map = varname_comp_2_bid2bc_map_[bc_varname_comp];
        ChiLogicalErrorIf(
          bid2bc_map.count(bid) != 0,
          "More than one boundary condition specified on boundary \"" + bname +
            "\" for variable \"" + bc_varname_comp.first + " component " +
            std::to_string(bc_varname_comp.second));

        bid2bc_map[bid] = bc;
      }
    } // for bc
  }   // for item in bndry map

  //======================================== Precompute QPData
  for (auto& field_info : *primary_fields_container_)
  {
    auto& field = field_info.field_;
    auto& sdm = field->GetSpatialDiscretization();
    if (sdm_stored_cell_qp_data_.count(&sdm) != 0) continue;

    auto& data_for_cells = sdm_stored_cell_qp_data_[&sdm];

    for (const auto& cell : grid.local_cells)
    {
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      CellQPData cell_qp_data =
        cell_mapping.MakeVolumetricQuadraturePointData();

      std::map<size_t, FaceQPData> faces_qp_data;
      size_t f = 0;
      for (const auto& face : cell.faces_)
      {
        if (not face.has_neighbor_)
          faces_qp_data[f] = cell_mapping.MakeSurfaceQuadraturePointData(f);
        ++f;
      }

      data_for_cells.push_back(
        {std::move(cell_qp_data), std::move(faces_qp_data)});
    }
  }
}

} // namespace chi_math