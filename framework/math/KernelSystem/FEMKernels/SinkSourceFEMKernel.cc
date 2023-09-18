#include "SinkSourceFEMKernel.h"

#include "ChiObjectFactory.h"

namespace chi_math
{

RegisterChiObject(chi_math, SinkSourceFEMKernel);

chi::InputParameters SinkSourceFEMKernel::GetInputParameters()
{
  chi::InputParameters params = FEMKernel::GetInputParameters();

  params.SetGeneralDescription(
    "A kernel simulating a distributed, cell-constant, sink or source.");
  params.SetDocGroup("doc_Kernels");

  params.AddRequiredParameter<double>(
    "value", "The sink/source value. A negative value would be a sink.");

  return params;
}

SinkSourceFEMKernel::SinkSourceFEMKernel(const chi::InputParameters& params)
  : FEMKernel(params), value_(params.GetParamValue<double>("value"))
{
}

double SinkSourceFEMKernel::ResidualEntryAtQP()
{
  return -value_ * test_values_[i_][qp_];
}

} // namespace chi_math