#include "LayeredVolumeAverage.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

namespace chi::derived_object
{

RegisterChiObject(chi::derived_object, LayeredVolumeAverage);

InputParameters LayeredVolumeAverage::GetInputParameters()
{
  InputParameters params = DerivedObject::GetInputParameters();
  params += chi_physics::GridBasedFieldFunctionInterface::GetInputParameters();

  params.AddRequiredParameterArray("nodes",
                                   "An array of 1 dimensional nodes (minimum 2 "
                                   "required) defining the layers.");

  params.AddOptionalParameterArray(
    "root", std::vector<double>{0.0, 0.0, 0.0}, "Root of the layer system.");
  params.AddOptionalParameterArray("direction",
                                   std::vector<double>{0.0, 0.0, 1.0},
                                   "A direction vector indicating the axis of "
                                   "the layer system. Default 0.0,0.0,1.0");

  params.AddOptionalParameter(
    "initial_value", 0.0, "Initial value to assign to the layers.");

  params.AddOptionalParameter(
    "reference_component", 0, "Reference unknown component on parent field.");

  params.AddOptionalParameterArray("parent_up_direction",
                                   std::vector<double>{0.0, 0.0, 1.0},
                                   "Default upward z-direction.");
  params.AddOptionalParameterArray("parent_forward_direction",
                                   std::vector<double>{0.0, 1.0, 0.0},
                                   "Default x-direction.");

  return params;
}

LayeredVolumeAverage::LayeredVolumeAverage(const InputParameters& params)
  : DerivedObject(params),
    chi_physics::GridBasedFieldFunctionInterface(params),
    nodes_(params.GetParamVectorValue<double>("nodes")),
    root_(chi_mesh::Vector3(params.GetParamVectorValue<double>("root"))),
    direction_(
      chi_mesh::Vector3(params.GetParamVectorValue<double>("direction"))),
    reference_component_(params.GetParamValue<uint32_t>("reference_component")),
    parent_rotation_matrix_(MakeRotationMatrix(params))
{
  ChiInvalidArgumentIf(nodes_.size() < 2,
                       "Node set \"nodes\" must have at least 2 entries.");
  {
    double prev_node = nodes_.front();
    bool monotonic = true;
    for (size_t k = 1; k < nodes_.size(); ++k)
    {
      if (nodes_[k] <= prev_node)
      {
        monotonic = false;
        break;
      }
      prev_node = nodes_[k];
    }
    ChiInvalidArgumentIf(
      not monotonic, "Node set \"nodes\" must be monotonically increasing.");
  }

  const double dir_norm = direction_.Norm();
  ChiInvalidArgumentIf(dir_norm < 1.0e-8,
                       "The norm of the direction vector is too small, " +
                         std::to_string(dir_norm) +
                         ". We do allow unnormalized direction vectors but"
                         " a norm this small normally indicates a problem.");
  direction_ /= dir_norm;

  values_.assign(nodes_.size() - 1, 0.0);
}

/**Finds the layer index in which this point lies. If not within the layers,
 * returns the number of layers (end).*/
size_t LayeredVolumeAverage::FindLayer(const chi_mesh::Vector3& position) const
{
  const auto ppo = position - root_;
  const double oneD_pos = ppo.Dot(direction_);

  for (size_t k = 0; k < (nodes_.size() - 1); ++k)
  {
    if (k == 0)
    {
      if (oneD_pos >= nodes_[k] and oneD_pos <= nodes_[k + 1]) return k;
    }
    else
    {
      if (oneD_pos > nodes_[k] and oneD_pos <= nodes_[k + 1]) return k;
    }
  }

  return values_.size();
}

chi_mesh::Matrix3x3
LayeredVolumeAverage::MakeRotationMatrix(const InputParameters& params)
{
  auto normal = chi_mesh::Vector3(
    params.GetParamVectorValue<double>("parent_up_direction"));
  auto binorm = chi_mesh::Vector3(
    params.GetParamVectorValue<double>("parent_forward_direction"));

  const double normal_norm = normal.Norm();
  const double binorm_norm = binorm.Norm();

  ChiInvalidArgumentIf(
    normal_norm < 1.0e-10,
    "Vector supplied in parameter \"parent_up_direction\" has a norm that is "
    "too small. Perfect normalization is not required but at least try.");

  ChiInvalidArgumentIf(
    binorm_norm < 1.0e-10,
    "Vector supplied in parameter \"parent_forward_direction\" has a norm that "
    "is too small. Perfect normalization is not required but at least try.");

  normal /= normal_norm;
  binorm /= binorm_norm;

  auto tangent = binorm.Cross(normal);

  const double tangent_norm = tangent.Norm();

  ChiLogicalErrorIf(std::fabs(tangent_norm) < 1.0e-10,
                    "Tangent norm very small. This can "
                    "happen if the vectors in \"parent_up_direction\" and "
                    "\"parent_forward_direction\" are in the same or "
                    "opposite directions.");
  tangent /= tangent_norm;

  chi_mesh::Matrix3x3 R;
  R.SetColJVec(0, tangent);
  R.SetColJVec(1, binorm);
  R.SetColJVec(2, normal);

  return R;
}

void LayeredVolumeAverage::Update(const Event& event)
{
  auto gridbased_ff_ptr = GetGridBasedFieldFunction();

  ChiLogicalErrorIf(not gridbased_ff_ptr, "Invalid grid based field function.");

  const auto& sdm = gridbased_ff_ptr->GetSpatialDiscretization();
  const auto& grid = sdm.Grid();
  const auto& uk_man = gridbased_ff_ptr->GetUnknownManager();
  const auto& local_data = gridbased_ff_ptr->FieldVectorRead();

  auto coord = sdm.GetSpatialWeightingFunction();

  const uint32_t c = reference_component_;

  const size_t num_layers = values_.size();
  std::vector<double> local_layer_intgls(num_layers, 0.0);
  std::vector<double> local_layer_volumes(num_layers, 0.0);
  for (const auto& cell : grid.local_cells)
  {
    const size_t layer_index =
      FindLayer(parent_rotation_matrix_ * cell.centroid_);
    if (layer_index == num_layers) continue;

    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto qp_data = cell_mapping.MakeVolumetricQuadraturePointData();
    const size_t num_nodes = cell_mapping.NumNodes();

    const auto& qp_indices = qp_data.QuadraturePointIndices();
    const auto& shape_values = qp_data.ShapeValues();
    const auto& JxW_values = qp_data.JxW_Values();

    double element_intgl = 0.0;
    double element_volume = 0.0;
    for (uint32_t qp : qp_indices)
    {
      double phi_qp = 0.0;
      for (size_t j = 0; j < num_nodes; ++j)
      {
        const int64_t dof_map = sdm.MapDOFLocal(cell, j, uk_man, 0, c);
        phi_qp += shape_values[j][qp] * local_data[dof_map];
        element_volume +=
          coord(qp_data.QPointXYZ(qp)) * shape_values[j][qp] * JxW_values[qp];
      }

      element_intgl += coord(qp_data.QPointXYZ(qp)) * JxW_values[qp] * phi_qp;
    } // for qp

    local_layer_intgls[layer_index] += element_intgl;
    local_layer_volumes[layer_index] += element_volume;
  } // for cell

  std::vector<double> global_layer_intgls(num_layers, 0.0);
  std::vector<double> global_layer_volumes(num_layers, 0.0);
  MPI_Allreduce(local_layer_intgls.data(),  // sendbuf
                global_layer_intgls.data(), // recvbuf
                num_layers,                 // count
                MPI_DOUBLE,                 // datatype
                MPI_SUM,                    // operation
                Chi::mpi.comm);             // communicator

  MPI_Allreduce(local_layer_volumes.data(),  // sendbuf
                global_layer_volumes.data(), // recvbuf
                num_layers,                  // count
                MPI_DOUBLE,                  // datatype
                MPI_SUM,                     // operation
                Chi::mpi.comm);              // communicator

  for (size_t k = 0; k < num_layers; ++k)
    if (global_layer_volumes[k] > 1.0e-10)
      values_[k] = global_layer_intgls[k] / global_layer_volumes[k];

  //std::stringstream outstr;
  //for (auto val : values_)
  //  outstr << val << " ";
  //
  //Chi::log.Log() << Name() << " " << outstr.str();
}

double
LayeredVolumeAverage::SpatialValue(const chi_mesh::Vector3& position) const
{
  const size_t layer_index = FindLayer(position);

  if (layer_index >= values_.size()) return 0.0;
  else
    return values_[layer_index];
}

} // namespace chi::derived_object