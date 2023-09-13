#ifndef PIPER_HEATCONDUCTIONSYSTEM_H
#define PIPER_HEATCONDUCTIONSYSTEM_H

#include "ChiObject.h"

namespace chi_mesh
{
class MeshContinuum;
}

namespace chi_math
{
class SpatialDiscretization;
class FEMKernel;
class FEMBoundaryCondition;

typedef std::shared_ptr<FEMKernel> FEMKernelPtr;
typedef std::shared_ptr<FEMBoundaryCondition> FEMBoundaryConditionPtr;
} // namespace chi_math

namespace hcm
{

class HeatConductionSystem : public ChiObject
{
public:
  static chi::InputParameters GetInputParameters();
  explicit HeatConductionSystem(const chi::InputParameters& params);

  const std::string& GetTemperatureFFName() const;

  const chi_mesh::MeshContinuum& Grid() const;
  const chi_math::SpatialDiscretization& SDM() const;
  chi_math::SpatialDiscretizationPtr& SDMPtr();

  std::vector<chi_math::FEMKernelPtr>& VolumeKernels();
  std::vector<chi_math::FEMBoundaryConditionPtr>& BoundaryConditions();


protected:
  const std::string temperature_ff_name_;

  std::shared_ptr<const chi_mesh::MeshContinuum> grid_ptr_ = nullptr;
  std::shared_ptr<chi_math::SpatialDiscretization> sdm_ptr_ = nullptr;

  std::vector<chi_math::FEMKernelPtr> volume_kernels_;
  std::vector<chi_math::FEMBoundaryConditionPtr> boundary_conditions_;
};

} // namespace hcm

#endif // PIPER_HEATCONDUCTIONSYSTEM_H
