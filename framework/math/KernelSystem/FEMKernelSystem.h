#ifndef CHI_FEMKERNELSYSTEM_H
#define CHI_FEMKERNELSYSTEM_H

#include "KernelSystem.h"
#include "math/UnknownManager/unknown_manager.h"
#include "math/SpatialDiscretization/FiniteElement/finite_element.h"
#include "math/SpatialDiscretization/CellMappings/cell_mapping_base.h"
#include "mesh/chi_mesh.h"

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

/**General Finite Element Method Kernel system.*/
class FEMKernelSystem : public KernelSystem
{
public:
  // 00_constrdestr
  /**\brief Basic constructor for a FEMKernelSystem.*/
  FEMKernelSystem(
    const SpatialDiscretization& sdm,
    const UnknownManager& uk_man,
    std::vector<chi_math::FEMKernelPtr>& volume_kernels,
    std::vector<chi_math::FEMBoundaryConditionPtr>& boundary_conditions);

  /**Returns the Spatial Discretization Method (SDM).*/
  const SpatialDiscretization& SDM() const;

  /**Returns the nodal unknown structure of each variable set.*/
  const UnknownManager& UnknownStructure() const;

  /**Obtains the kernels associated with a material.*/
  virtual const std::vector<FEMKernelPtr>& GetMaterialKernels(int mat_id);

  // 01_SetInitSol
  /**Used to create an initial guess. Mostly applies Dirichlet BCs.*/
  void SetInitialSolution() override;

  // 02a
  /**Collective method for computing the system residual.*/
  void ComputeResidual(ParallelVector& r) override;
  // 02b
  /**Collective method for computing the system Jacobian-matrix.*/
  void ComputeJacobian(ParallelMatrix& J) override;

protected:
  // 03_kernel_setup
  /**Initializes cell data prior kernel and BC setup.*/
  void InitCellData(const chi_mesh::Cell& cell);

  /**Prepares all the necessary data for internal kernels.*/
  std::vector<std::shared_ptr<FEMKernel>>
  SetupCellInternalKernels(const chi_mesh::Cell& cell);

  /**Prepares all the necessary data for boundary kernels.*/
  std::vector<std::pair<size_t, FEMBoundaryConditionPtr>>
  SetupCellBCKernels(const chi_mesh::Cell& cell);

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
    std::vector<double> local_x_;
  } cur_cell_data;
};

} // namespace chi_math

#endif // CHI_FEMKERNELSYSTEM_H
