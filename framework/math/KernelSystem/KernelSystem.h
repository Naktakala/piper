#ifndef CHITECH_KERNELSYSTEM_H
#define CHITECH_KERNELSYSTEM_H

#include "math/Systems/EquationSystem.h"
#include "math/UnknownManager/unknown_manager.h"
#include "math/SpatialDiscretization/FiniteElement/finite_element.h"
#include "math/SpatialDiscretization/CellMappings/cell_mapping_base.h"

namespace chi
{
class MaterialPropertiesData;
}

namespace chi_mesh
{
class Cell;
}

namespace chi_math
{
class SpatialDiscretization;
class FEMKernel;
class FEMBoundaryCondition;
class ParallelMatrix;
class FEMMaterialProperty;

typedef finite_element::InternalQuadraturePointData CellQPData;
typedef finite_element::FaceQuadraturePointData FaceQPData;
typedef std::shared_ptr<FEMKernel> FEMKernelPtr;
typedef std::shared_ptr<FEMBoundaryCondition> FEMBoundaryConditionPtr;

/**Data pack that can be placed on the object stack for kernels and bc's to
 * use when being constructed.*/
class FEMKernelSystemData : public ChiObject
{
public:
  typedef chi_mesh::Vector3 Vec3;
  typedef std::vector<Vec3> VecVec3;

  FEMKernelSystemData(const std::vector<std::shared_ptr<FEMMaterialProperty>>&
                        fem_material_properties,
                      const EquationSystemTimeData& time_data,
                      const ParallelVector& main_solution_vector,
                      const chi_mesh::Cell*& cell_ptr,
                      const chi_math::CellMapping*& cell_mapping_ptr,
                      const CellQPData& qp_data,
                      const VecDbl& var_qp_values,
                      const VecVec3& var_grad_qp_values,
                      const VecDbl& var_dot_qp_values,
                      const VecDbl& nodal_var_values,
                      const VecVec3& node_locations,

                      const FaceQPData& face_qp_data,
                      const VecDbl& face_var_qp_values,
                      const VecVec3& face_var_grad_qp_values);

  const std::vector<std::shared_ptr<FEMMaterialProperty>>&
    fem_material_properties_;
  const EquationSystemTimeData& time_data_;

  const ParallelVector& main_solution_vector_;

  const chi_mesh::Cell*& cell_ptr_;
  const chi_math::CellMapping*& cell_mapping_ptr_;

  const CellQPData& qp_data_;
  const VecDbl& var_qp_values_;
  const VecVec3& var_grad_qp_values_;
  const VecDbl& var_dot_qp_values_;
  const VecDbl& nodal_var_values_;
  const VecVec3& node_locations_;

  const FaceQPData& face_qp_data_;
  const VecDbl& face_var_qp_values_;
  const VecVec3& face_var_grad_qp_values_;
};

/**General Finite Element Method Kernel system.*/
class KernelSystem : public EquationSystem
{
public:
  // 00_constrdestr
  static chi::InputParameters GetInputParameters();
  /**\brief Basic constructor for a FEMKernelSystem.*/
  explicit KernelSystem(const chi::InputParameters& params);

  /**Obtains the kernels associated with a material.*/
  std::vector<FEMKernelPtr> GetMaterialKernels(int mat_id);

  /**Obtains the boundary kernel associated with a boundary id.*/
  FEMBoundaryConditionPtr GetBoundaryCondition(uint64_t boundary_id);

  // 01_SetInitSol
  /**Used to create an initial guess. Mostly applies Dirichlet BCs.*/
  void SetInitialSolution() override;

  // 02a
  /**Collective method for computing the system residual.*/
  void ComputeResidual(const ParallelVector& x, ParallelVector& r) override;
  // 02b
  /**Collective method for computing the system Jacobian-matrix.*/
  void ComputeJacobian(const ParallelVector& x, ParallelMatrix& J) override;

protected:
  // 00a_PopulateMaterialProperties
  void PopulateMaterialProperties();

  // 00b_MakeKernels
  std::vector<FEMKernelPtr>
  MakeFEMKernels(const chi::ParameterBlock& volume_kernels_inputs,
                 size_t fem_data_handle);

  // 00c_MakeBCs
  std::vector<FEMBoundaryConditionPtr>
  MakeBCs(const chi::ParameterBlock& boundary_condition_inputs,
          size_t fem_data_handle);

  // 03_kernel_setup
  /**Initializes cell data prior kernel and BC setup.*/
  void InitCellData(const ParallelVector& x, const chi_mesh::Cell& cell);

  /**Prepares all the necessary data for internal kernels.*/
  std::vector<std::shared_ptr<FEMKernel>>
  SetupAndGetCellInternalKernels(const chi_mesh::Cell& cell);

  std::vector<std::pair<size_t, FEMBoundaryConditionPtr>>
  GetCellBCKernels(const chi_mesh::Cell& cell);

  void PrecomputeMaterialProperties(const chi_mesh::Cell& cell);

  void SetupFaceIntegralBCKernel(const chi_mesh::Cell& cell, size_t face_index);

  /**Returns a set of dirichlet nodes by looking at the BCs applied on
   * faces. Does not get filtered by time status.*/
  std::set<uint32_t>
  IdentifyLocalDirichletNodes(const chi_mesh::Cell& cell) const;

  std::map<int, std::vector<FEMKernelPtr>> matid_2_volume_kernels_map_;
  typedef std::map<uint64_t, FEMBoundaryConditionPtr> BID2BCMap;
  typedef std::pair<std::string, uint32_t> VarNameComp;
  std::map<VarNameComp, BID2BCMap> varname_comp_2_bid2bc_map_;

  std::vector<std::shared_ptr<FEMMaterialProperty>> fem_material_properties_;

  size_t current_field_index_ = 0;
  uint32_t current_field_component_ = 0;

  /**Utility data structure to store data ahead of executing the kernels*/
  struct CurrentCellData
  {
    chi_mesh::Cell const* cell_ptr_ = nullptr;
    const chi_math::CellMapping* cell_mapping_ptr_ = nullptr;

    std::vector<chi_mesh::Vector3> node_locations_;
    std::vector<int64_t> dof_map_;
    VecDbl local_x_;
    VecDbl local_x_dot_;

    CellQPData qp_data_;
    VecDbl var_qp_values_;
    VecVec3 var_grad_qp_values_;
    VecDbl var_dot_qp_values_;
  } cur_cell_data;

  struct CurrentFaceData
  {
    FaceQPData qp_data_;
    VecDbl var_qp_values_;
    VecVec3 var_grad_qp_values_;
  } cur_face_data;
};

} // namespace chi_math

#endif // CHITECH_KERNELSYSTEM_H
