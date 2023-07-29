#ifndef PIPER_PIPER_H
#define PIPER_PIPER_H

#include "physics/SolverBase/chi_solver.h"
#include "mesh/chi_mesh.h"

namespace piper
{
class Component;

class Piper : public chi_physics::Solver
{
public:
  static chi::InputParameters GetInputParameters();
  explicit Piper(const chi::InputParameters& params);

  const std::string& SystemName() const;

  virtual void Initialize() override;

protected:
  std::map<std::string, std::shared_ptr<Component>> components_;

  std::vector<std::string> boundary_component_names_;
  std::vector<std::string> volume_component_names_;
  std::vector<std::string> junction_component_names_;

private:
  const std::string system_name_;
  std::string root_component_;
  const chi_mesh::Vector3 datum_;
  const bool print_nodalization_;
};

}

#endif // PIPER_PIPER_H
