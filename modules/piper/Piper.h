#ifndef PIPER_PIPER_H
#define PIPER_PIPER_H

#include "physics/SolverBase/chi_solver.h"
#include "mesh/chi_mesh.h"

namespace piper
{
class HardwareComponent;
class FluidPhysics;

class Piper : public chi_physics::Solver
{
public:
  static chi::InputParameters GetInputParameters();
  explicit Piper(const chi::InputParameters& params);

  const std::string& SystemName() const;
  const std::vector<std::shared_ptr<HardwareComponent>>& HardwareComponents() const;

  const std::vector<size_t>& BoundaryComponentIDs() const;
  const std::vector<size_t>& VolumeComponentIDs() const;
  const std::vector<size_t>& JunctionComponentIDs() const;

  /**Returns the hardware-component id for a component with a given name.*/
  size_t MapHWCompName2ID(const std::string& name) const;

  size_t RootComponentID() const;

  virtual void Initialize() override;
  virtual void Step() override;
private:
  void ConnectComponents();

protected:
  std::map<std::string, size_t> hw_comp_name_2_id_map_;
  std::vector<std::shared_ptr<HardwareComponent>> hardware_components_;

  std::vector<size_t> boundary_component_ids_;
  std::vector<size_t> volume_component_ids_;
  std::vector<size_t> junction_component_ids_;

  std::shared_ptr<FluidPhysics> physics_;

private:
  const std::string system_name_;
  const std::string root_component_name_;
  size_t root_component_id_;
  const chi_mesh::Vector3 datum_;
  const bool print_nodalization_;
};

}

#endif // PIPER_PIPER_H
