#include "PiperMeshGenerator.h"

#include "piper/Piper.h"
#include "piper/components/HardwareComponent.h"

#include "chi_log.h"

namespace piper
{

chi::InputParameters PiperMeshGenerator::GetInputParameters()
{
  chi::InputParameters params = chi_mesh::MeshGenerator::GetInputParameters();

  return params;
}

PiperMeshGenerator::PiperMeshGenerator(const chi::InputParameters& params)
  : chi_mesh::MeshGenerator(params)
{
}

void PiperMeshGenerator::SetPipeSystem(const Piper& pipe_system)
{
  pipe_system_ptr_ = &pipe_system;
}

std::unique_ptr<chi_mesh::UnpartitionedMesh>
PiperMeshGenerator::GenerateUnpartitionedMesh(
  std::unique_ptr<chi_mesh::UnpartitionedMesh> input_umesh)
{
  Chi::log.Log0Verbose1()
    << "PiperMeshGenerator::GenerateUnparitionedMesh begin";
  ChiLogicalErrorIf(not pipe_system_ptr_, "Called without pipe_system_set");
  const auto& pipe_system = *pipe_system_ptr_;

  auto umesh = std::make_unique<chi_mesh::UnpartitionedMesh>();

  size_t vertex_counter = 0;
  size_t cell_counter = 0;

  const auto& components = pipe_system.HardwareComponents();
  const auto& volume_comp_IDs = pipe_system.VolumeComponentIDs();
  for (const auto& volume_comp_id : volume_comp_IDs)
  {
    const auto& component = components.at(volume_comp_id);
    const auto& con_points = component->ConnectionPoints();

    const auto& v0 = con_points.front().position_;
    const auto& v1 = con_points.back().position_;

    umesh->GetVertices().push_back(v0);
    ++vertex_counter;
    umesh->GetVertices().push_back(v1);
    ++vertex_counter;

    auto new_cell = new chi_mesh::UnpartitionedMesh::LightWeightCell(
      chi_mesh::CellType::SLAB, chi_mesh::CellType::SLAB);

    new_cell->vertex_ids = {vertex_counter - 1, vertex_counter - 2};

    new_cell->faces.emplace_back(std::vector<uint64_t>{vertex_counter - 1});
    new_cell->faces.emplace_back(std::vector<uint64_t>{vertex_counter - 2});

    umesh->GetRawCells().push_back(new_cell);

    volume_comp_name_2_cell_id_map_[component->Name()] = cell_counter;
    ++cell_counter;
  } // for volume_names

  umesh->ComputeCentroidsAndCheckQuality();
  umesh->BuildMeshConnectivity();

  Chi::log.Log0Verbose1()
    << "PiperMeshGenerator::GenerateUnparitionedMesh done";
  return umesh;
}

const std::map<std::string, uint64_t>&
PiperMeshGenerator::GetVolumeComponent2CellGIDMap() const
{
  return volume_comp_name_2_cell_id_map_;
}

} // namespace piper
