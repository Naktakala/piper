#ifndef PIPER_FLUIDPHYSICS_H
#define PIPER_FLUIDPHYSICS_H

#include "ChiObject.h"

#include "piper/models/ComponentModel.h"
#include "mesh/chi_mesh.h"

namespace chi_mesh
{
class MeshContinuum;
}

namespace piper
{
class Piper;
class PiperMeshGenerator;

class FluidPhysics : public ChiObject
{
public:
  static chi::InputParameters GetInputParameters();
  explicit FluidPhysics(const chi::InputParameters& params);

  void SetPipeSystem(const Piper& pipe_system);
  const Piper& PipeSystem() const;

  const chi_mesh::Vector3& GravityVector() const;

  const chi_mesh::MeshContinuum& Grid() const;

  virtual void MakeMesh();
  virtual void InitializeUnknowns() = 0;
  /**Makes a list of variable names that should be added for a given
  * component category.*/
  virtual std::vector<std::string>
  MakeVariableNamesList(ComponentCategory hw_comp_category) = 0;

  virtual void Step() = 0;

protected:
  const size_t min_cells_per_processor_;
  const chi::ParameterBlock initializer_param_block_;

  const Piper* pipe_system_ptr_;
  std::unique_ptr<PiperMeshGenerator> mesh_generator_;
  std::shared_ptr<const chi_mesh::MeshContinuum> grid_ptr_;

  std::vector<std::unique_ptr<ComponentModel>> component_models_;

  chi_mesh::Vector3 gravity_;
};

} // namespace piper

#endif // PIPER_FLUIDPHYSICS_H
