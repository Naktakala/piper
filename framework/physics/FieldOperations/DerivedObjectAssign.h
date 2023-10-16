#ifndef CHITECH_DERIVEDOBJECTASSIGN_H
#define CHITECH_DERIVEDOBJECTASSIGN_H

#include "physics/FieldOperations/field_operation.h"
#include "physics/FieldFunction/GridBasedFieldFunctionInterface.h"

namespace chi
{
class DerivedObject;
}

namespace chi_physics::field_operations
{

class DerivedObjectAssign : public FieldOperation,
                            public GridBasedFieldFunctionInterface
{
public:
  static chi::InputParameters GetInputParameters();
  explicit DerivedObjectAssign(const chi::InputParameters& params);

  void Execute() override;

protected:
  const chi::DerivedObject& derived_object_;
};

} // namespace chi_physics::field_operations

#endif // CHITECH_DERIVEDOBJECTASSIGN_H
