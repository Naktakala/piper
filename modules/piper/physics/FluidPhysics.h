#ifndef PIPER_FLUIDPHYSICS_H
#define PIPER_FLUIDPHYSICS_H

#include "ChiObject.h"

#include "piper/models/ComponentModel.h"
#include "mesh/chi_mesh.h"

#include <functional>

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

  double DeltaT() const;
  void SetTimeStep(double dt);
  double Time() const;
  void SetTime(double time);
  double EndTime() const;
  void SetEndTime(double end_time);

  std::function<double(const ComponentModel&)> FrictionFactorFuncion() const;

  virtual void Initialize() = 0;
  virtual void MakeMesh();
  virtual void InitializeUnknowns() = 0;
  /**Makes a list of variable names that should be added for a given
   * component category.*/
  virtual std::vector<std::string>
  MakeVariableNamesList(ComponentCategory hw_comp_category) = 0;

  virtual void Step() = 0;
  virtual void Advance() = 0;
  virtual void Execute();

protected:
  const size_t min_cells_per_processor_;
  const chi::ParameterBlock initializer_param_block_;

  const Piper* pipe_system_ptr_;
  std::unique_ptr<PiperMeshGenerator> mesh_generator_;
  std::shared_ptr<const chi_mesh::MeshContinuum> grid_ptr_;

  chi_mesh::Vector3 gravity_;

  std::vector<std::unique_ptr<ComponentModel>> component_models_;

  std::function<double(const ComponentModel&)> friction_factor_function_;

  std::map<std::string, chi::ParameterBlock> compononent_model_parameters_;

private:
  static std::map<std::string, chi::ParameterBlock>
  ConstructComponentModelParametersMap(
    const chi::ParameterBlock& params);
  double dt_ = 0.01;
  double time_ = 0.0;
  double end_time_ = 1.0;
  int max_time_steps_ = -1;
  size_t t_index_ = 0;
};

} // namespace piper

#endif // PIPER_FLUIDPHYSICS_H
