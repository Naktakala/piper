#ifndef PIPER_PIPERMESHGENERATOR_H
#define PIPER_PIPERMESHGENERATOR_H

#include "mesh/MeshGenerator/MeshGenerator.h"

namespace piper
{

class Piper;

class PiperMeshGenerator : public chi_mesh::MeshGenerator
{
public:
  static chi::InputParameters GetInputParameters();
  explicit PiperMeshGenerator(const chi::InputParameters& params);

  void SetPipeSystem(const Piper& pipe_system);

  std::unique_ptr<chi_mesh::UnpartitionedMesh> GenerateUnpartitionedMesh(
    std::unique_ptr<chi_mesh::UnpartitionedMesh> input_umesh) override;

  std::map<std::string, uint64_t> GetVolumeComponent2CellGIDMap();

protected:
  const Piper* pipe_system_ptr_ = nullptr;
  std::map<std::string, uint64_t> volume_comp_name_2_cell_id_map_;
};

} // namespace piper

#endif // PIPER_PIPERMESHGENERATOR_H
