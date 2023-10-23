#include "BoundaryAverage.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

namespace chi::derived_object
{

RegisterChiObject(chi::derived_object, BoundaryAverage);

InputParameters BoundaryAverage::GetInputParameters()
{
  InputParameters params = DerivedObject::GetInputParameters();
  params += chi_physics::GridBasedFieldFunctionInterface::GetInputParameters();

  params.AddOptionalParameter(
    "reference_component", 0, "Reference unknown component on parent field.");

  params.AddOptionalParameter(
    "initial_value", 0.0, "Initial value to assign to the layers.");

  params.AddOptionalParameterArray("boundaries",
                                   std::vector<std::string>{},
                                   "A list of boundary names on which this "
                                   "boundary condition mush be supplied.");

  return params;
}

BoundaryAverage::BoundaryAverage(const InputParameters& params)
  : DerivedObject(params),
    chi_physics::GridBasedFieldFunctionInterface(params),
    reference_component_(params.GetParamValue<uint32_t>("reference_component")),
    boundary_scope_(params.GetParamVectorValue<std::string>("boundaries")),
    value_(params.GetParamValue<double>("initial_value"))
{
}

void BoundaryAverage::Update(const Event& event)
{
  auto gridbased_ff_ptr = GetGridBasedFieldFunction();

  ChiLogicalErrorIf(not gridbased_ff_ptr, "Invalid grid based field function.");

  const auto& sdm = gridbased_ff_ptr->GetSpatialDiscretization();
  const auto& grid = sdm.Grid();
  const auto& uk_man = gridbased_ff_ptr->GetUnknownManager();
  const auto& local_data = gridbased_ff_ptr->FieldVectorRead();

  auto coord = sdm.GetSpatialWeightingFunction();

  const auto& bndry_id_map = grid.GetBoundaryIDMap();

  std::map<std::string, uint64_t> inverse_map;
  for (const auto& [bid, bname] : bndry_id_map)
    inverse_map[bname] = bid;

  std::vector<uint64_t> wanted_boundaries;
  for (const auto& bname : boundary_scope_)
  {
    auto it = inverse_map.find(bname);
    ChiLogicalErrorIf(it == inverse_map.end(),
                      "Boundary name \"" + bname + "\" not found.");
    wanted_boundaries.push_back(it->second);
  }

  const uint32_t c = reference_component_;

  double local_integral = 0.0;
  double local_area = 0.0;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    const size_t num_faces = cell.faces_.size();
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces_[f];
      if (face.has_neighbor_) continue;

      const bool do_boundary =
        wanted_boundaries.empty() or
        std::find(wanted_boundaries.begin(),
                  wanted_boundaries.end(),
                  face.neighbor_id_) != wanted_boundaries.end();

      if (not do_boundary) continue;

      const auto qp_data = cell_mapping.MakeSurfaceQuadraturePointData(f);

      const auto& qp_indices = qp_data.QuadraturePointIndices();
      const auto& shape_values = qp_data.ShapeValues();
      const auto& JxW_values = qp_data.JxW_Values();

      double element_intgl = 0.0;
      double element_area = 0.0;
      for (uint32_t qp : qp_indices)
      {
        double phi_qp = 0.0;
        for (size_t j = 0; j < num_nodes; ++j)
        {
          const int64_t dof_map = sdm.MapDOFLocal(cell, j, uk_man, 0, c);
          phi_qp += shape_values[j][qp] * local_data[dof_map];
          element_area +=
            coord(qp_data.QPointXYZ(qp)) * shape_values[j][qp] * JxW_values[qp];
        }

        element_intgl += coord(qp_data.QPointXYZ(qp)) * JxW_values[qp] * phi_qp;
      } // for qp

      local_integral += element_intgl;
      local_area += element_area;
    } // for f

  } // for cell

  double global_integral;
  double global_area;
  MPI_Allreduce(&local_integral,  // sendbuf
                &global_integral, // recvbuf
                1,                // count
                MPI_DOUBLE,       // datatype
                MPI_SUM,          // operation
                Chi::mpi.comm);   // communicator

  MPI_Allreduce(&local_area,    // sendbuf
                &global_area,   // recvbuf
                1,              // count
                MPI_DOUBLE,     // datatype
                MPI_SUM,        // operation
                Chi::mpi.comm); // communicator

  value_ = global_integral / global_area;
}

double BoundaryAverage::SpatialValue(const chi_mesh::Vector3&) const
{
  return value_;
}

} // namespace chi::derived_object