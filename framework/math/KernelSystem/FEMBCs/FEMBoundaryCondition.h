#ifndef CHITECH_FEMBOUNDARYCONDITION_H
#define CHITECH_FEMBOUNDARYCONDITION_H

#include "ChiObject.h"
#include "interfaces/CoupledFieldInterface.h"
#include "math/KernelSystem/Coupling/FEMMaterialPropertyInterface.h"
#include "math/KernelSystem/Coupling/FEMCoupledField.h"
#include "mesh/chi_mesh.h"
#include "math/chi_math.h"

namespace chi_math
{

namespace finite_element
{
class FaceQuadraturePointData;
}

class FEMKernelSystemData;

/**The abstract base class of a Finite Element Method Boundary Condition
 * Kernel.*/
class FEMBoundaryCondition : public ChiObject,
                             public chi_math::CoupledFieldInterface

{
public:
  static chi::InputParameters GetInputParameters();
  explicit FEMBoundaryCondition(const chi::InputParameters& params);

  virtual bool IsDirichlet() const;

  void PreComputeValues();

  /**Computes the current face's contribution to the residual at cell node i.*/
  virtual double ComputeLocalResidual(size_t f, uint32_t i);
  /**Computes the current face's contribution to the jacobian at cell node i,
   * referenced to the node associated with cell node j.*/
  virtual double ComputeLocalJacobian(size_t f, uint32_t i, uint32_t j);

  /**Returns the residual entry at the quadrature point.*/
  virtual double ResidualEntryAtQP();
  /**Returns the jacobian entry at the quadrature point.*/
  virtual double JacobianEntryAtQP();

  const std::vector<std::string>& GetBoundaryScope() const;

  const std::pair<std::string, uint32_t>& ActiveVariableAndComponent() const;

protected:
  const std::pair<std::string, uint32_t> var_name_component_;
  const std::vector<std::string> boundary_scope_;

  const FEMKernelSystemData& fem_data_;

  typedef chi_mesh::Vector3 Vec3;
  typedef std::vector<Vec3> VecVec3;
  typedef std::vector<VecVec3> MatVec3;

  const double& dt_;
  const double& time_;

  const VecDbl& JxW_values_;
  const MatDbl& test_values_;
  const MatVec3& test_grad_values_;
  const MatDbl& shape_values_;
  const MatVec3& shape_grad_values_;
  const VecDbl& var_value_;
  const VecVec3& var_grad_value_;
  const VecVec3& qp_xyz_;
  const VecVec3& normal_;
  const VecDbl& nodal_var_values_;
  const VecVec3& node_locations_;
  const VecDbl& coord_;

  size_t i_ = 0;
  size_t j_ = 0;
  size_t qp_ = 0;
};

} // namespace chi_math

#endif // CHITECH_FEMBOUNDARYCONDITION_H
