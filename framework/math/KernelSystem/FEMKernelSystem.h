#ifndef PIPER_FEMKERNELSYSTEM_H
#define PIPER_FEMKERNELSYSTEM_H

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

typedef finite_element::InternalQuadraturePointData CellQPData;
typedef finite_element::FaceQuadraturePointData FaceQPData;

class SpatialDiscretization;
class FEMKernel;
class FEMBoundaryCondition;
class ParallelMatrix;

typedef std::shared_ptr<FEMKernel> FEMKernelPtr;
typedef std::shared_ptr<FEMBoundaryCondition> FEMBoundaryConditionPtr;

class FEMKernelSystem : public KernelSystem
{
public:
  FEMKernelSystem(
    const SpatialDiscretization& sdm,
    const UnknownManager& uk_man,
    std::vector<chi_math::FEMKernelPtr>& volume_kernels,
    std::vector<chi_math::FEMBoundaryConditionPtr>& boundary_conditions);

  void SetInitialSolution() override;

  void ComputeResidual(ParallelVector& r) override;
  void ComputeJacobian(ParallelMatrix& J) override;

  const SpatialDiscretization& SDM() const;
  const UnknownManager& UnknownStructure() const;

  virtual const std::vector<FEMKernelPtr>& GetMaterialKernels(int mat_id);

protected:
  void InitCellKernelData(const chi_mesh::Cell& cell);
  std::vector<std::shared_ptr<FEMKernel>> SetupCellInternalKernels(const chi_mesh::Cell& cell);
  std::vector<std::pair<size_t, FEMBoundaryConditionPtr>> SetupCellBCKernels(const chi_mesh::Cell& cell);
  void SetupInternalKernelsRefData(const std::vector<FEMKernelPtr>& kernels);
  void SetupBoundaryConditionRefData(chi_math::FEMBoundaryCondition& bndry_condition);

  const SpatialDiscretization& sdm_;
  const UnknownManager uk_man_;

  std::map<int, std::vector<FEMKernelPtr>> matid_2_volume_kernels_map_;
  std::map<uint64_t, FEMBoundaryConditionPtr> bid_2_boundary_conditions_map_;

  struct CurrentCellData
  {
    const chi_mesh::Cell* cur_cell_ptr_ = nullptr;
    const chi_math::CellMapping* cell_mapping_ptr_ = nullptr;
    /**Node sets. First = internal nodes, Second = Bndry nodes.*/
    std::pair<std::set<uint32_t>, std::set<uint32_t>> node_id_sets_;

    std::vector<chi_mesh::Vector3> node_locations_;
    std::vector<int64_t> dof_map_;
    std::vector<double> local_x_;

    std::unique_ptr<CellQPData> cell_qp_data_ = nullptr;
    std::vector<double> cell_var_values_;
    std::vector<chi_mesh::Vector3> cell_var_grad_values_;

    std::unique_ptr<FaceQPData> face_qp_data_ = nullptr;
    std::vector<double> face_var_values_;
    std::vector<chi_mesh::Vector3> face_var_grad_values_;

    size_t current_face_index_ = 0;
  } cur_cell_data;
};

} // namespace chi_math

#endif // PIPER_FEMKERNELSYSTEM_H
