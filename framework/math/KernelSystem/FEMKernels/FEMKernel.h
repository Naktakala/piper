#ifndef CHI_FEMKERNEL_H
#define CHI_FEMKERNEL_H

#include "ChiObject.h"
#include "mesh/chi_mesh.h"

namespace chi_math
{

namespace finite_element
{
class InternalQuadraturePointData;
}

/**A data structure to hold reference data for FEM Kernels.*/
struct FEMKernelRefData
{
  typedef std::vector<chi_mesh::Vector3> VecVec3;
  typedef std::vector<VecVec3> VecVecVec3;
  typedef std::vector<double> VecDbl;
  typedef std::vector<VecDbl> VecVecDbl;

  FEMKernelRefData(
    const std::shared_ptr<const finite_element::InternalQuadraturePointData>&
      qp_data_ptr_,
    VecDbl var_qp_values,
    VecVec3 var_grad_qp_values);

  const std::shared_ptr<const finite_element::InternalQuadraturePointData>
    qp_data_ptr_;
  const std::vector<unsigned int>& qp_indices_;
  const VecVec3& qpoints_xyz_;
  const VecVecDbl& shape_values_;
  const VecVecVec3& shape_grad_values_;
  const VecDbl& JxW_values_;
  const VecDbl var_qp_values_;
  const VecVec3 var_grad_qp_values_;
};

/**The abstract base class of a Finite Element Method Kernel.*/
class FEMKernel : public ChiObject
{
public:
  static chi::InputParameters GetInputParameters();
  explicit FEMKernel(const chi::InputParameters& params);

  /**Stores reference data as a pointer so that it can be reassigned
  * and automatically garbage collected.*/
  void SetReferenceData(std::shared_ptr<FEMKernelRefData>& ref_data_ptr);

  /**Computes the current cell's contribution to the residual at cell node i.*/
  virtual double ComputeLocalResidual(uint32_t i);

  /**Computes the current cell's contribution to the jacobian at cell node i,
  * referenced to the node associated with cell node j.*/
  virtual double ComputeLocalJacobian(uint32_t i, uint32_t j);

  /**Returns the residual entry at the quadrature point.*/
  virtual double ResidualEntryAtQP(); // default 0.0

  /**Returns the jacobian entry at the quadrature point.*/
  virtual double JacobianEntryAtQP(); // default 0.0

  /**Returns a list of material ids to which this kernel is applied.*/
  const std::vector<int>& GetMaterialIDScope() const;

protected:
  std::shared_ptr<FEMKernelRefData> ref_data_ptr_ = nullptr;

  double test_i_qp_ = 0.0;
  chi_mesh::Vector3 test_grad_i_qp_;
  double var_qp_value_ = 0.0;
  chi_mesh::Vector3 var_grad_qp_value_;
  chi_mesh::Vector3 qp_xyz_;

  double shape_j_qp_ = 0.0;
  chi_mesh::Vector3 shape_grad_j_qp_;

private:
  std::vector<int> mat_ids_;
};

} // namespace chi_math

#endif // CHI_FEMKERNEL_H
