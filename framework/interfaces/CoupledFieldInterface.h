#ifndef CHITECH_COUPLEDFIELDINTERFACE_H
#define CHITECH_COUPLEDFIELDINTERFACE_H

#include "parameters/input_parameters.h"
#include "math/KernelSystem/Coupling/FEMCoupledField.h"

namespace chi_math
{

class MultiFieldContainer;
class FEMCoupledField;

class CoupledFieldInterface
{
public:
  static chi::InputParameters GetInputParameters();
  explicit CoupledFieldInterface(const chi::InputParameters& params);

  void PreComputeInternalCoupledFields();
  void PreComputeFaceCoupledFields();

protected:
  const FEMCoupledField& GetCoupledField(const std::string& field_name);

private:
  const std::vector<size_t> multifieldcontainer_handles_;
  const size_t fem_data_handle_;
  std::vector<std::unique_ptr<FEMCoupledField>> coupled_field_ptrs_;
};

}

#endif // CHITECH_COUPLEDFIELDINTERFACE_H
