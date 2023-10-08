#ifndef CHITECH_KERNELSYSTEM_H
#define CHITECH_KERNELSYSTEM_H

#include "math/Systems/FieldEquationSystem.h"
#include "FEMKernelSystemData.h"

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

typedef std::shared_ptr<FEMKernel> FEMKernelPtr;
typedef std::shared_ptr<FEMBoundaryCondition> FEMBoundaryConditionPtr;

/**General Finite Element Method Kernel system.*/
class KernelSystem : public FieldEquationSystem
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

  void SetupFaceIntegralBCKernel(const chi_mesh::Cell& cell, size_t face_index);

  /**Returns a set of dirichlet nodes by looking at the BCs applied on
   * faces. Does not get filtered by time status.*/
  std::set<uint32_t>
  IdentifyLocalDirichletNodes(const chi_mesh::Cell& cell) const;

  // Members requiring initialization
  std::map<int, std::vector<FEMKernelPtr>> matid_2_volume_kernels_map_;
  typedef std::map<uint64_t, FEMBoundaryConditionPtr> BID2BCMap;
  typedef std::pair<std::string, uint32_t> VarNameComp;
  std::map<VarNameComp, BID2BCMap> varname_comp_2_bid2bc_map_;

  struct CellQPDataContainer
  {
    CellQPData cell_qp_data_;
    std::map<size_t, FaceQPData> faces_qp_data_;
  };
  std::map<const SpatialDiscretization*,
           std::vector<CellQPDataContainer>> sdm_stored_cell_qp_data_;

  // Runtime data
  size_t current_field_index_ = 0;
  uint32_t current_field_component_ = 0;
  const SpatialDiscretization* current_sdm_ = nullptr;
  const CellQPDataContainer* current_cell_qp_container_ = nullptr;

  typedef std::function<double(const chi_mesh::Vector3&)> SpatialWeightFunction;
  SpatialWeightFunction current_spatial_weight_function_;

  FEMKernelSystemData::CurrentCellData cur_cell_data;

  FEMKernelSystemData::CurrentFaceData cur_face_data;
};

} // namespace chi_math

#endif // CHITECH_KERNELSYSTEM_H
