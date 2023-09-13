#ifndef PIPER_FEMKERNEL_H
#define PIPER_FEMKERNEL_H

#include "ChiObject.h"
#include "mesh/chi_mesh.h"

namespace chi_math
{

namespace finite_element
{
class InternalQuadraturePointData;
}

struct FEMKernelRefData
{
  typedef std::vector<chi_mesh::Vector3> VecVec3;
  typedef std::vector<VecVec3> VecVecVec3;
  typedef std::vector<double> VecDbl;
  typedef std::vector<VecDbl> VecVecDbl;

  FEMKernelRefData(const finite_element::InternalQuadraturePointData& qp_data,
                   const VecDbl& var_values,
                   const VecVec3& var_grad_values);

  const std::vector<unsigned int>& qp_indices_;
  const VecVec3& qpoints_xyz_;
  const VecVecDbl& shape_values_;
  const VecVecVec3& shape_grad_values_;
  const VecDbl& JxW_values_;
  const VecDbl& var_values_;
  const VecVec3& var_grad_values_;
};

class FEMKernel : public ChiObject
{
public:
  static chi::InputParameters GetInputParameters();
  explicit FEMKernel(const chi::InputParameters& params);

  void SetReferenceData(std::shared_ptr<FEMKernelRefData>& ref_data_ptr);

  virtual double ComputeLocalResidual(uint32_t i);
  virtual double ComputeLocalJacobian(uint32_t i, uint32_t j);

  virtual double ResidualEntryAtQP();
  virtual double JacobianEntryAtQP();

  const std::vector<int>& GetMaterialIDScope() const;

protected:
  std::shared_ptr<FEMKernelRefData> ref_data_ptr_ = nullptr;

  double test_i_qp_ = 0.0;
  chi_mesh::Vector3 test_grad_i_qp_;
  double var_value_qp_ = 0.0;
  chi_mesh::Vector3 var_grad_value_qp_;
  chi_mesh::Vector3 qp_xyz_;

  double shape_j_qp_ = 0.0;
  chi_mesh::Vector3 shape_grad_j_qp_;

private:
  std::vector<int> mat_ids_;
};

} // namespace chi_math

#endif // PIPER_FEMKERNEL_H
