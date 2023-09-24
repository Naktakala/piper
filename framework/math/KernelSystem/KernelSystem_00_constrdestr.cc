#include "KernelSystem.h"

#include "math/SpatialDiscretization/spatial_discretization.h"
#include "math/SpatialDiscretization/FiniteElement/finite_element.h"
#include "math/KernelSystem/FEMKernels/FEMKernel.h"
#include "math/KernelSystem/FEMBCs/FEMBoundaryCondition.h"
#include "math/ParallelMatrix/ParallelMatrix.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

#define scint64_t(x) static_cast<int64_t>(x)
#define cint64_t const int64_t

namespace chi_math
{

FEMKernelSystemData::FEMKernelSystemData(
  const std::vector<std::shared_ptr<FEMMaterialProperty>>&
    fem_material_properties,
  const EquationSystemTimeData& time_data,
  const CellQPData& qp_data,
  const VecDbl& var_qp_values,
  const VecVec3& var_grad_qp_values,
  const VecDbl& var_dot_qp_values,
  const VecDbl& nodal_var_values,
  const VecVec3& node_locations,

  const FaceQPData& face_qp_data,
  const VecDbl& face_var_qp_values,
  const VecVec3& face_var_grad_qp_values)

  : ChiObject(chi::InputParameters{}),
    fem_material_properties_(fem_material_properties),
    time_data_(time_data),
    qp_data_(qp_data),
    var_qp_values_(var_qp_values),
    var_grad_qp_values_(var_grad_qp_values),
    var_dot_qp_values_(var_dot_qp_values),
    nodal_var_values_(nodal_var_values),
    node_locations_(node_locations),

    face_qp_data_(face_qp_data),
    face_var_qp_values_(face_var_qp_values),
    face_var_grad_qp_values_(face_var_grad_qp_values)
{
}

RegisterChiObject(chi_math, KernelSystem);

chi::InputParameters KernelSystem::GetInputParameters()
{
  chi::InputParameters params = EquationSystem::GetInputParameters();

  params.SetGeneralDescription("An equation system based on kernels.");
  params.SetDocGroup("doc_KernelSystem");

  params.AddRequiredParameterArray("kernels",
                                   "A list of kernel in parameter form.");
  params.AddRequiredParameterArray("bcs",
                                   "A list of boundary condition handles.");

  return params;
}

/**\brief Basic constructor for a FEMKernelSystem.
 * This constructor also sorts kernels and BCs into convenient maps.*/
KernelSystem::KernelSystem(const chi::InputParameters& params)
  : EquationSystem(params)
{
  const auto& volume_kernels_inputs = params.GetParam("kernels");
  const auto boundary_condition_inputs = params.GetParam("bcs");

  PopulateMaterialProperties();

  //======================================== Make reference data
  auto reference_data =
    std::make_shared<FEMKernelSystemData>(fem_material_properties_,
                                          time_data_,
                                          cur_cell_data.qp_data_,
                                          cur_cell_data.var_qp_values_,
                                          cur_cell_data.var_grad_qp_values_,
                                          cur_cell_data.var_dot_qp_values_,
                                          cur_cell_data.local_x_,
                                          cur_cell_data.node_locations_,

                                          cur_face_data.qp_data_,
                                          cur_face_data.var_qp_values_,
                                          cur_face_data.var_grad_qp_values_);

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

  //====== EXPERIMENTAL
  for (auto& field_info : *primary_fields_container_)
    field_info.field_->SDM().InitializeQPData(/*internal_faces=*/false,
                                              /*bndry_faces=*/true);
}

} // namespace chi_math