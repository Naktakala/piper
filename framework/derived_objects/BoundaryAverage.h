#ifndef CHITECH_BOUNDARYAVERAGE_H
#define CHITECH_BOUNDARYAVERAGE_H

#include "DerivedObject.h"
#include "physics/FieldFunction/GridBasedFieldFunctionInterface.h"

namespace chi::derived_object
{

class BoundaryAverage : public DerivedObject, public chi_physics::GridBasedFieldFunctionInterface
{
public:
  static InputParameters GetInputParameters();
  explicit BoundaryAverage(const InputParameters& params);

  void Update(const Event& event) override;

  double SpatialValue(const chi_mesh::Vector3& position) const override;

protected:
  uint32_t reference_component_;
  const std::vector<std::string> boundary_scope_;

  double value_;
};

}

#endif // CHITECH_BOUNDARYAVERAGE_H
