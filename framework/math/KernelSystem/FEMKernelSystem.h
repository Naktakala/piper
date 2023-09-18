#ifndef CHI_FEMKERNELSYSTEM_H
#define CHI_FEMKERNELSYSTEM_H

#include "KernelSystem.h"
#include "math/UnknownManager/unknown_manager.h"
#include "math/SpatialDiscretization/FiniteElement/finite_element.h"
#include "math/SpatialDiscretization/CellMappings/cell_mapping_base.h"

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

  FEMKernelSystemData(const EquationSystemTimeData& time_data,
                      const CellQPData& qp_data,
                      const VecDbl& var_qp_values,
                      const VecVec3& var_grad_qp_values,
                      const MatDbl& old_var_qp_values,
                      const VecDbl& nodal_var_values,
                      const VecVec3& node_locations,

                      const FaceQPData& face_qp_data,
                      const VecDbl& face_var_qp_values,
                      const VecVec3& face_var_grad_qp_values);
  const EquationSystemTimeData& time_data_;
  const CellQPData& qp_data_;
  const VecDbl& var_qp_values_;
  const VecVec3& var_grad_qp_values_;
  const MatDbl& old_var_qp_values_;
  const VecDbl& nodal_var_values_;
  const VecVec3& node_locations_;

  const FaceQPData& face_qp_data_;
  const VecDbl& face_var_qp_values_;
  const VecVec3& face_var_grad_qp_values_;
};

/**General Finite Element Method Kernel system.*/
class FEMKernelSystem : public KernelSystem
{
public:
  // 00_constrdestr
  /**\brief Basic constructor for a FEMKernelSystem.*/
  FEMKernelSystem(
    const SpatialDiscretization& sdm,
    const UnknownManager& uk_man,
    const std::vector<chi::ParameterBlock>& volume_kernels_inputs,
    const std::vector<chi::ParameterBlock>& boundary_condition_inputs,
    TimeID oldest_time_id);

  /**Returns the Spatial Discretization Method (SDM).*/
  const SpatialDiscretization& SDM() const;

  /**Returns the nodal unknown structure of each variable set.*/
  const UnknownManager& UnknownStructure() const;

  /**Obtains the kernels associated with a material.*/
  std::vector<FEMKernelPtr> GetMaterialKernels(int mat_id);

  /**Obtains the boundary kernel associated with a boundary id.*/
  FEMBoundaryConditionPtr GetBoundaryCondition(uint64_t boundary_id);

  // 01_SetInitSol
  /**Used to create an initial guess. Mostly applies Dirichlet BCs.*/
  void SetInitialSolution() override;

  // 02a
  /**Collective method for computing the system residual.*/
  void ComputeResidual(const GhostedParallelVector& x,
                       ParallelVector& r) override;
  // 02b
  /**Collective method for computing the system Jacobian-matrix.*/
  void ComputeJacobian(const GhostedParallelVector& x,
                       ParallelMatrix& J) override;

  ParallelMatrixSparsityPattern BuildMatrixSparsityPattern() const override;

protected:
  // 03_kernel_setup
  /**Initializes cell data prior kernel and BC setup.*/
  void InitCellData(const GhostedParallelVector& x, const chi_mesh::Cell& cell);

  /**Prepares all the necessary data for internal kernels.*/
  std::vector<std::shared_ptr<FEMKernel>>
  SetupAndGetCellInternalKernels(const chi_mesh::Cell& cell);

  ///**Prepares all the necessary data for boundary kernels.*/
  // std::vector<std::pair<size_t, FEMBoundaryConditionPtr>>
  // SetupCellBCKernels(const chi_mesh::Cell& cell);

  std::vector<std::pair<size_t, FEMBoundaryConditionPtr>>
  GetCellBCKernels(const chi_mesh::Cell& cell);

  void SetupFaceIntegralBCKernel(size_t face_index);

  /**Returns a set of dirichlet nodes by looking at the BCs applied on
   * faces. Does not get filtered by time status.*/
  std::set<uint32_t>
  IdentifyLocalDirichletNodes(const chi_mesh::Cell& cell) const;

  const SpatialDiscretization& sdm_;
  const UnknownManager uk_man_;

  std::map<int, std::vector<FEMKernelPtr>> matid_2_volume_kernels_map_;
  std::map<uint64_t, FEMBoundaryConditionPtr> bid_2_BCKernel_map_;

  /**Utility data structure to store data ahead of executing the kernels*/
  struct CurrentCellData
  {
    const chi_math::CellMapping* cell_mapping_ptr_ = nullptr;

    std::vector<chi_mesh::Vector3> node_locations_;
    std::vector<int64_t> dof_map_;
    VecDbl local_x_;
    MatDbl old_local_x_;

    CellQPData qp_data_;
    VecDbl var_qp_values_;
    VecVec3 var_grad_qp_values_;
    MatDbl old_var_qp_values_;
  } cur_cell_data;

  struct CurrentFaceData
  {
    FaceQPData qp_data_;
    VecDbl var_qp_values_;
    VecVec3 var_grad_qp_values_;
  } cur_face_data;
};

} // namespace chi_math

#endif // CHI_FEMKERNELSYSTEM_H
