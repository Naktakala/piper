#ifndef CHI_FEMBOUNDARYCONDITION_H
#define CHI_FEMBOUNDARYCONDITION_H

#include "ChiObject.h"
#include "mesh/chi_mesh.h"

namespace chi_math
{

namespace finite_element
{
class FaceQuadraturePointData;
}

struct EquationSystemTimeData;

/**A data structure to hold reference data for FEM Boundary Condition Kernels.*/
struct FEMBCRefData
{
  typedef std::vector<chi_mesh::Vector3> VecVec3;
  typedef std::vector<VecVec3> VecVecVec3;
  typedef std::vector<double> VecDbl;
  typedef std::vector<VecDbl> VecVecDbl;

  FEMBCRefData(const EquationSystemTimeData& time_data,
               const std::shared_ptr<
                 const finite_element::FaceQuadraturePointData>& qp_data,
               VecDbl var_qp_values,
               VecVec3 var_grad_qp_values,
               const VecDbl& var_nodal_values,
               const VecVec3& node_locations);

  const EquationSystemTimeData& time_data_;
  const std::shared_ptr<const finite_element::FaceQuadraturePointData> qp_data_;
  const std::vector<unsigned int>& qp_indices_;
  const VecVec3& qpoints_xyz_;
  const VecVecDbl& shape_values_;
  const VecVecVec3& shape_grad_values_;
  const VecDbl& JxW_values_;
  const VecVec3& normals_;
  const VecDbl var_qp_values_;
  const VecVec3 var_grad_qp_values_;
  const VecDbl& var_nodal_values_;
  const VecVec3& node_locations_;
};

/**The abstract base class of a Finite Element Method Boundary Condition
 * Kernel.*/
class FEMBoundaryCondition : public ChiObject
{
public:
  typedef size_t FaceID;
  typedef std::shared_ptr<FEMBCRefData> FEMBCRefDataPtr;

  static chi::InputParameters GetInputParameters();
  explicit FEMBoundaryCondition(const chi::InputParameters& params);

  virtual bool IsDirichlet() const;

  /**Stores reference data as a pointer so that it can be reassigned
   * and automatically garbage collected.*/
  void SetFaceReferenceData(FaceID face_index,
                            std::shared_ptr<FEMBCRefData>& ref_data_ptr);

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

protected:
  const std::vector<std::string> boundary_scope_;

  // std::shared_ptr<FEMBCRefData> ref_data_ptr_ = nullptr;

  std::map<size_t, FEMBCRefDataPtr> face_id_2_ref_data_ptr_map_;

  size_t f_ = 0;
  size_t i_ = 0;
  double test_i_qp_ = 0.0;
  chi_mesh::Vector3 test_grad_i_qp_;
  double var_value_qp_ = 0.0;
  chi_mesh::Vector3 var_grad_value_qp_;
  chi_mesh::Vector3 qp_xyz_;
  chi_mesh::Vector3 normal_qp_;

  double dt_ = 1.0;
  double time_ = 0.0;

  double shape_j_qp_ = 0.0;
  chi_mesh::Vector3 shape_grad_j_qp_;
};

} // namespace chi_math

#endif // CHI_FEMBOUNDARYCONDITION_H
