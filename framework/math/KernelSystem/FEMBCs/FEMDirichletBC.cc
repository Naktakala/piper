#include "FEMDirichletBC.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

namespace chi_math
{

RegisterChiObject(chi_math, FEMDirichletBC);

chi::InputParameters FEMDirichletBC::GetInputParameters()
{
  chi::InputParameters params = FEMBoundaryCondition::GetInputParameters();

  params.SetGeneralDescription(
    "A basic, zero-valued Dirichlet boundary condtion.");
  params.SetDocGroup("doc_KernelSystem");

  params.AddOptionalParameter(
    "bc_value", 0.0, "Value of the boundary condition.");

  params.AddOptionalParameter("apply_before_solve",
                              true,
                              "Flag, when true, inserts this value into the"
                              "solution vector before the solve begins.");

  return params;
}

FEMDirichletBC::FEMDirichletBC(const chi::InputParameters& params)
  : FEMBoundaryCondition(params),
    bc_value_(params.GetParamValue<double>("bc_value")),
    apply_before_solve_(params.GetParamValue<bool>("apply_before_solve"))
{
}

bool FEMDirichletBC::IsDirichlet() const { return true; }


double FEMDirichletBC::ComputeLocalResidual(size_t f, uint32_t i)
{
  auto& ref_data = *face_id_2_ref_data_ptr_map_.at(f);
  i_ = i;
  node_xyz_ = ref_data.node_locations_[i];
  var_node_value_ = ref_data.var_nodal_values_[i];

  return ResidualEntryAtQP();
}

double FEMDirichletBC::ComputeLocalJacobian(size_t f, uint32_t i, uint32_t j)
{
  ChiLogicalError("If this function was called, it must be a mistake!");
}

double FEMDirichletBC::ResidualEntryAtQP()
{
  return var_node_value_ - bc_value_;
}

double FEMDirichletBC::BCValue() const { return bc_value_; }

bool FEMDirichletBC::AllowApplyBeforeSolve() const
{
  return apply_before_solve_;
}

} // namespace chi_math