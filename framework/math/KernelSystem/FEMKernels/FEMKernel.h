#ifndef CHITECH_FEMKERNEL_H
#define CHITECH_FEMKERNEL_H

#include "ChiObject.h"
#include "materials/MaterialIDScopeInterface.h"
#include "interfaces/CoupledFieldInterface.h"
#include "math/KernelSystem/Coupling/FEMCoupledField.h"
#include "mesh/chi_mesh.h"
#include "math/chi_math.h"

namespace chi_math
{

class FEMKernelSystemData;
class FEMMaterialProperty;

/**The abstract base class of a Finite Element Method Kernel.*/
class FEMKernel : public ChiObject,
                  public chi::MaterialIDScopeInterface,
                  public chi_math::CoupledFieldInterface
{
public:
  static chi::InputParameters GetInputParameters();
  explicit FEMKernel(const chi::InputParameters& params);

  /**Computes the current cell's contribution to the residual at cell node i.*/
  virtual double ComputeLocalResidual(uint32_t i);

  /**Computes the current cell's contribution to the jacobian at cell node i,
  * referenced to the node associated with cell node j.*/
  virtual double ComputeLocalJacobian(uint32_t i, uint32_t j);

  /**Returns the residual entry at the quadrature point.*/
  virtual double ResidualEntryAtQP() = 0;

  /**Returns the jacobian entry at the quadrature point.*/
  virtual double JacobianEntryAtQP() = 0;

  /**True if this kernel is derived from a Time kernel*/
  virtual bool IsTimeKernel() const;

  const std::pair<std::string, uint32_t>& ActiveVariableAndComponent() const;

protected:
  const chi_math::FEMMaterialProperty& GetFEMMaterialProperty(const std::string& name);
  const std::pair<std::string, uint32_t> var_name_component_;
  const FEMKernelSystemData& fem_data_;

  typedef chi_mesh::Vector3 Vec3;
  typedef std::vector<Vec3> VecVec3;
  typedef std::vector<VecVec3> MatVec3;

  const double& dt_;
  const double& time_;
  const double& var_dot_dvar_;

  const VecDbl& JxW_values_;
  const MatDbl& test_values_;
  const MatVec3& test_grad_values_;
  const MatDbl& shape_values_;
  const MatVec3& shape_grad_values_;
  const VecDbl& var_value_;
  const VecVec3& var_grad_value_;
  const VecDbl& var_dot_value_;
  const VecVec3& qp_xyz_;

  size_t i_ = 0;
  size_t j_ = 0;
  size_t qp_ = 0;
};

} // namespace chi_math

#endif // CHITECH_FEMKERNEL_H
