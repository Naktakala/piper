#ifndef CHITECH_MATERIALIDSCOPEINTERFACE_H
#define CHITECH_MATERIALIDSCOPEINTERFACE_H

#include "parameters/input_parameters.h"

namespace chi
{

class MaterialIDScopeInterface
{
public:
  /**Returns a list of material ids to which this object is applied.*/
  const std::vector<int>& GetMaterialIDScope() const;

protected:
  static InputParameters GetInputParameters();
  explicit MaterialIDScopeInterface(const InputParameters& params);

private:
  std::vector<int> mat_ids_;
};

}

#endif // CHITECH_MATERIALIDSCOPEINTERFACE_H
