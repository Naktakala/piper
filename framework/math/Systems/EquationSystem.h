#ifndef CHITECH_EQUATIONSYSTEM_H
#define CHITECH_EQUATIONSYSTEM_H

#include "ChiObject.h"

#include "EquationSystemTimeData.h"
#include "math/ParallelVector/ParallelVector.h"
#include "math/Containers/MultiFieldContainer.h"
#include "mesh/chi_mesh.h"

namespace chi
{
class MaterialPropertiesData;
}

namespace chi_physics
{
class FieldFunctionGridBased;
}

namespace chi_math
{

class MultiFieldContainer;

typedef std::vector<std::shared_ptr<chi_physics::FieldFunctionGridBased>>
  FieldList;

class ParallelMatrix;
struct ParallelMatrixSparsityPattern;
class TimeIntegrator;

/**Enum for controlling which terms go into a specific residual/jacobian
 * computation. Mostly set by Executioners.*/
enum class EqTermScope : int
{
  NONE = 0,
  TIME_TERMS = (1 << 0),
  DOMAIN_TERMS = (1 << 2),
  BOUNDARY_TERMS = (1 << 3)
};

/**Defining the bit-wise "or" operator in order to easily add flags.*/
inline EqTermScope operator|(const EqTermScope f1, const EqTermScope f2)
{
  return static_cast<EqTermScope>(static_cast<int>(f1) | static_cast<int>(f2));
}
/**Defining the bit-wise "and" operator to allow us to easily check
 * if a term is active.*/
inline bool operator&(const EqTermScope f1, const EqTermScope f2)
{
  return static_cast<int>(f1) & static_cast<int>(f2);
}

/**Base abstract class for a system of equations.*/
class EquationSystem : public ChiObject
{
public:
  typedef chi_mesh::Vector3 Vec3;
  typedef std::vector<Vec3> VecVec3;

  /**Returns the number of local DOFs across all unknowns.*/
  int64_t NumLocalDOFs() const;
  /**Returns the number of local DOFs across all unknowns.*/
  int64_t NumGlobalDOFs() const;

  /**Returns a reference to the current solution vector.*/
  ParallelVector& SolutionVector(TimeID time_id = TimeID::T_PLUS_1);
  /**Returns a reference to the current residual vector.*/
  ParallelVector& ResidualVector(TimeID time_id = TimeID::T_PLUS_1);

  /**Returns a reference to the system's time integrator.*/
  TimeIntegrator& GetTimeIntegrator();

  /**Sets the current solution vector.*/
  virtual void SetInitialSolution(){};

  /**Computes the residual vector \p r given a solution vector \p x.
   * This method is generally only called by executioners.*/
  virtual void ComputeResidual(const ParallelVector& x, ParallelVector& r) = 0;

  /**Computes the Jacobian matrix \p J given a solution vector \p x.
   * This method is generally only called by executioners.*/
  virtual void ComputeJacobian(const ParallelVector& x, ParallelMatrix& J) = 0;

  /**Sets the current time data.*/
  void SetTimeData(EquationSystemTimeData time_data);

  /**Returns the current equation terms-scope.*/
  EqTermScope EquationTermsScope() const;
  /**Sets the scope of current equations.*/
  void SetEquationTermsScope(EqTermScope eq_term_scope);

  /**Determines if a particular equation term is active.*/
  bool QueryTermsActive(EqTermScope query_scope) const;

  /**Returns a reference to the system current time data.*/
  const EquationSystemTimeData& GetTimeData() const;

  /**Uses the underlying system to build a sparsity pattern.*/
  virtual ParallelMatrixSparsityPattern BuildMatrixSparsityPattern() const;

  /**Updates the fields.*/
  void UpdateFields();

  /**Output fields to VTK. The filename passed via the options will be used
   * plus a time index (if transient).*/
  void OutputFields(int time_index = -1);

  /**Advances the system in time.*/
  void Advance(EquationSystemTimeData time_data,
               std::map<TimeID, const ParallelVector*>& latest_std_residuals);

protected:
  static chi::InputParameters GetInputParameters();
  explicit EquationSystem(const chi::InputParameters& params);

  struct FieldBlockInfo;

  const chi::MaterialPropertiesData& material_properties_data_;
  const int verbosity_;
  const std::string output_file_base_;

  std::shared_ptr<MultiFieldContainer> primary_fields_container_;

  const int64_t num_local_dofs_;
  const int64_t num_globl_dofs_;

  std::shared_ptr<TimeIntegrator> time_integrator_;
  const size_t num_solution_histories_;
  const size_t num_residual_histories_;

  std::unique_ptr<ParallelVector> main_solution_vector_;

  std::vector<std::unique_ptr<ParallelVector>> old_solution_vectors_;
  std::vector<std::unique_ptr<ParallelVector>> old_residual_vectors_;

  EquationSystemTimeData time_data_;

  std::vector<size_t> t_tags_;

private:
  static const chi::MaterialPropertiesData&
  GetOrMakeMaterialPropertiesData(const chi::InputParameters& params);

  static std::shared_ptr<MultiFieldContainer>
  MakeMultifieldContainer(const chi::ParameterBlock& params);

  std::unique_ptr<ParallelVector> MakeSolutionVector();
  static std::shared_ptr<TimeIntegrator>
  InitTimeIntegrator(const chi::InputParameters& params);
  std::vector<std::unique_ptr<ParallelVector>> InitSolutionHistory();
  std::vector<std::unique_ptr<ParallelVector>> InitResidualHistory();

  EqTermScope eq_term_scope_;
};

} // namespace chi_math

#endif // CHITECH_EQUATIONSYSTEM_H
