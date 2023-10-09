#ifndef PIPER_FLUIDPHYSICS_H
#define PIPER_FLUIDPHYSICS_H

#include "physics/SolverBase/chi_solver.h"

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

class FluidPhysics : public chi_physics::Solver
{
public:
  static chi::InputParameters GetInputParameters();
  explicit FluidPhysics(const chi::InputParameters& params);

  const Piper& PipeSystem() const;

  const chi_mesh::Vector3& GravityVector() const;

  const chi_mesh::MeshContinuum& Grid() const;

  std::function<double(const ComponentModel&)> FrictionFactorFuncion() const;

  virtual void MakeMesh();
  virtual void InitializeUnknowns() = 0;
  /**Makes a list of variable names that should be added for a given
   * component category.*/
  virtual std::vector<std::string>
  MakeVariableNamesList(ComponentCategory hw_comp_category) = 0;

  void Execute() override;
  chi::ParameterBlock GetInfo(const chi::ParameterBlock& params) const override;

  static void BroadcastStateMap(const std::vector<std::string>& map_keys,
                                std::map<std::string, double>& state_map,
                                uint64_t root);
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

  double step_time_ = 0.0;
  double intgl_step_time_ = 0.0;
  double step_counter_ = 0.0;

private:
  static std::map<std::string, chi::ParameterBlock>
  ConstructComponentModelParametersMap(const chi::ParameterBlock& params);
};

} // namespace piper

#endif // PIPER_FLUIDPHYSICS_H
