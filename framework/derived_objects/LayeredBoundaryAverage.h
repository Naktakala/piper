#ifndef CHITECH_LAYEREDBOUNDARYAVERAGE_H
#define CHITECH_LAYEREDBOUNDARYAVERAGE_H

#include "DerivedObject.h"
#include "physics/FieldFunction/GridBasedFieldFunctionInterface.h"
#include "mesh/chi_mesh.h"

namespace chi::derived_object
{

class LayeredBoundaryAverage
  : public DerivedObject,
    public chi_physics::GridBasedFieldFunctionInterface
{
public:
  static InputParameters GetInputParameters();
  explicit LayeredBoundaryAverage(const InputParameters& params);

  void Update(const Event& event) override;

  double SpatialValue(const chi_mesh::Vector3& position) const override;

protected:
  /**Finds the layer index in which this point lies. If not within the layers,
  * returns the number of layers (end).*/
  size_t FindLayer(const chi_mesh::Vector3& position) const;

  static chi_mesh::Matrix3x3 MakeRotationMatrix(const InputParameters& params);

  const std::vector<double> nodes_;
  const chi_mesh::Vector3 root_;
  chi_mesh::Vector3 direction_;
  uint32_t reference_component_;

  chi_mesh::Matrix3x3 parent_rotation_matrix_;

  const std::vector<std::string> boundary_scope_;

  std::vector<double> values_;
};

} // namespace chi::derived_object

#endif // CHITECH_LAYEREDBOUNDARYAVERAGE_H
