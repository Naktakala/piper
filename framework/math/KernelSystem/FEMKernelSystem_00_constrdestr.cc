#include "FEMKernelSystem.h"

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
  const EquationSystemTimeData& time_data,
  const CellQPData& qp_data,
  const VecDbl& var_qp_values,
  const VecVec3& var_grad_qp_values,
  const MatDbl& old_var_qp_values,
  const VecDbl& nodal_var_values,
  const VecVec3& node_locations,

  const FaceQPData& face_qp_data,
  const VecDbl& face_var_qp_values,
  const VecVec3& face_var_grad_qp_values)

  : ChiObject(chi::InputParameters{}),
    time_data_(time_data),
    qp_data_(qp_data),
    var_qp_values_(var_qp_values),
    var_grad_qp_values_(var_grad_qp_values),
    old_var_qp_values_(old_var_qp_values),
    nodal_var_values_(nodal_var_values),
    node_locations_(node_locations),

    face_qp_data_(face_qp_data),
    face_var_qp_values_(face_var_qp_values),
    face_var_grad_qp_values_(face_var_grad_qp_values)
{
}

RegisterChiObject(chi_math, FEMKernelSystem);

chi::InputParameters FEMKernelSystem::GetInputParameters()
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
FEMKernelSystem::FEMKernelSystem(const chi::InputParameters& params)
  : EquationSystem(params)
{
  const auto& volume_kernels_inputs = params.GetParam("kernels");
  const auto boundary_condition_inputs = params.GetParam("bcs");

  //======================================== Make reference data
  auto reference_data =
    std::make_shared<FEMKernelSystemData>(time_data_,
                                          cur_cell_data.qp_data_,
                                          cur_cell_data.var_qp_values_,
                                          cur_cell_data.var_grad_qp_values_,
                                          cur_cell_data.old_var_qp_values_,
                                          cur_cell_data.local_x_,
                                          cur_cell_data.node_locations_,

                                          cur_face_data.qp_data_,
                                          cur_face_data.var_qp_values_,
                                          cur_face_data.var_grad_qp_values_);

  Chi::object_stack.push_back(reference_data);
  const size_t data_handle = Chi::object_stack.size() - 1;

  //======================================== Make all kernels
  auto& object_factory = ChiObjectFactory::GetInstance();

  std::vector<FEMKernelPtr> kernels;
  for (auto kernel_input : volume_kernels_inputs) // making a copy
  {
    ChiInvalidArgumentIf(not kernel_input.Has("type"),
                         "A kernel input is missing the \"type\" parameter");

    kernel_input.AddParameter("fem_data_handle", data_handle);

    const auto obj_type = kernel_input.GetParamValue<std::string>("type");
    chi::InputParameters in_params =
      object_factory.GetRegisteredObjectParameters(obj_type);
    in_params.AssignParameters(kernel_input);

    const size_t kernel_handle =
      object_factory.MakeRegisteredObjectOfType(obj_type, in_params);

    auto kernel = Chi::GetStackItemPtrAsType<FEMKernel>(
      Chi::object_stack, kernel_handle, __FUNCTION__);

    kernels.push_back(kernel);
  }

  //======================================== Make all BC-Kernels
  std::vector<FEMBoundaryConditionPtr> boundary_conditions;
  for (auto bc_input : boundary_condition_inputs) // making a copy
  {
    ChiInvalidArgumentIf(
      not bc_input.Has("type"),
      "A BoundaryCondition input is missing the \"type\" parameter");

    bc_input.AddParameter("fem_data_handle", data_handle);

    const auto obj_type = bc_input.GetParamValue<std::string>("type");
    chi::InputParameters in_params =
      object_factory.GetRegisteredObjectParameters(obj_type);
    in_params.AssignParameters(bc_input);

    const size_t bc_handle =
      object_factory.MakeRegisteredObjectOfType(obj_type, in_params);

    auto boundary_condition = Chi::GetStackItemPtrAsType<FEMBoundaryCondition>(
      Chi::object_stack, bc_handle, __FUNCTION__);

    boundary_conditions.push_back(boundary_condition);
  }

  //======================================== Map mat-ids to kernels
  const auto& grid = field_block_info_.front().field_->SDM().Grid();
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
  const auto& current_field = field_block_info_.at(current_field_index_).field_;

  for (const auto& kernel : kernels)
  {
    bool kernel_active = false;
    if (kernel->IsTimeKernel() and time_terms_active) kernel_active = true;
    if (not kernel->IsTimeKernel() and domain_terms_active)
      kernel_active = true;

    if (kernel->ActiveVariableAndComponent().first != current_field->TextName())
      kernel_active = false;
    if (kernel->ActiveVariableAndComponent().second != current_field_component_)
      kernel_active = false;

    if (kernel_active) filtered_kernels.push_back(kernel);
  }

  if (time_terms_active and filtered_kernels.empty())
    ChiLogicalError("No time kernels in system");

  if (filtered_kernels.empty())
    ChiLogicalError("No kernels for material id " + std::to_string(mat_id) +
                    " on field \"" + current_field->TextName() + "\"");

  return filtered_kernels;
}

/**Obtains the boundary kernel associated with a boundary id.*/
FEMBoundaryConditionPtr
FEMKernelSystem::GetBoundaryCondition(uint64_t boundary_id)
{
  if (not(EquationTermsScope() & EqTermScope::BOUNDARY_TERMS)) return nullptr;

  const auto& current_field = field_block_info_.at(current_field_index_).field_;
  const auto& field_name = current_field->TextName();

  auto& bid2bc_map =
    varname_comp_2_bid2bc_map_.at({field_name, current_field_component_});

  auto bc_it = bid2bc_map.find(boundary_id);

  if (bc_it == bid2bc_map.end()) return nullptr;

  const auto& bndry_condition = bc_it->second;

  const auto& bc_var_comp = bndry_condition->ActiveVariableAndComponent();

  if (bc_var_comp.first != current_field->TextName() or
      bc_var_comp.second != current_field_component_)
    return nullptr;

  return bc_it->second;
}

ParallelMatrixSparsityPattern
FEMKernelSystem::BuildMatrixSparsityPattern() const
{
  std::vector<int64_t> master_row_nnz_in_diag;
  std::vector<int64_t> mastero_row_nnz_off_diag;

  auto& a1 = master_row_nnz_in_diag;
  auto& a2 = mastero_row_nnz_off_diag;
  for (const auto& field_info : field_block_info_)
  {
    auto& field = field_info.field_;
    auto& sdm = field->SDM();
    auto& uk_man = field->UnkManager();

    std::vector<int64_t> nodal_nnz_in_diag;
    std::vector<int64_t> nodal_nnz_off_diag;
    sdm.BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, uk_man);

    // Now append these to the master list
    auto& b1 = nodal_nnz_in_diag;
    auto& b2 = nodal_nnz_off_diag;

    using std::begin, std::end;
    a1.insert(end(a1), begin(b1), end(b1));
    a2.insert(end(a2), begin(b2), end(b2));
  }

  return {master_row_nnz_in_diag, mastero_row_nnz_off_diag};
}

} // namespace chi_math