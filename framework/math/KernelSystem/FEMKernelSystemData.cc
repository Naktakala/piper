#include "FEMKernelSystemData.h"

namespace chi_math
{

FEMKernelSystemData::FEMKernelSystemData(
  const chi::MaterialPropertiesData& mat_props_data,
  const EquationSystemTimeData& time_data,
  const ParallelVector& main_solution_vector,
  const CurrentCellData& cell_data,
  const CurrentFaceData& face_data)

  : ChiObject(chi::InputParameters{}),
    mat_props_data_(mat_props_data),
    time_data_(time_data),
    main_solution_vector_(main_solution_vector),
    cell_data_(cell_data),
    face_data_(face_data)
{
}

}