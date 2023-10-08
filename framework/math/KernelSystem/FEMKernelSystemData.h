#ifndef CHITECH_FEMKERNELDATA_H
#define CHITECH_FEMKERNELDATA_H

#include "ChiObject.h"
#include "math/SpatialDiscretization/CellMappings/CellMapping.h"
#include "mesh/chi_mesh.h"
#include "math/chi_math.h"
#include "math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"

namespace chi
{
class MaterialPropertiesData;
}

namespace chi_math
{

class EquationSystemTimeData;
class ParallelVector;

typedef finite_element::VolumetricQuadraturePointData CellQPData;
typedef finite_element::SurfaceQuadraturePointData FaceQPData;

/**Data pack that can be placed on the object stack for kernels and bc's to
 * use when being constructed.*/
class FEMKernelSystemData : public ChiObject
{
public:
  typedef chi_mesh::Vector3 Vec3;
  typedef std::vector<Vec3> VecVec3;

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
    VecDbl coord_qp_values_;
    VecDbl var_dot_qp_values_;
  };

  struct CurrentFaceData
  {
    FaceQPData qp_data_;
    VecDbl var_qp_values_;
    VecVec3 var_grad_qp_values_;
    VecDbl coord_qp_values_;
  };

  FEMKernelSystemData(const chi::MaterialPropertiesData& mat_props_data,
                      const EquationSystemTimeData& time_data,
                      const ParallelVector& main_solution_vector,
                      const CurrentCellData& cell_data,
                      const CurrentFaceData& face_data);

  const chi::MaterialPropertiesData& mat_props_data_;

  const EquationSystemTimeData& time_data_;

  const ParallelVector& main_solution_vector_;

  const CurrentCellData& cell_data_;
  const CurrentFaceData& face_data_;
};

}

#endif // CHITECH_FEMKERNELDATA_H
