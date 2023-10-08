#ifndef CHITECH_EQUATIONSYSTEMBASE_H
#define CHITECH_EQUATIONSYSTEMBASE_H

#include "ChiObject.h"

#include "EquationSystemTimeData.h"

#include "math/ParallelVector/ParallelVector.h"

#include <functional>

namespace chi_math
{
class ParallelMatrix;
class TimeIntegrator;
struct ParallelMatrixSparsityPattern;

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

class EquationSystem : public ChiObject
{
public:
  /**Returns the number of local DOFs across all unknowns.*/
  virtual int64_t NumLocalDOFs() const = 0;
  /**Returns the number of local DOFs across all unknowns.*/
  virtual int64_t NumGlobalDOFs() const = 0;

  /**Returns a reference to the system's time integrator.*/
  TimeIntegrator& GetTimeIntegrator();

  /**Returns a reference to the system current time data.*/
  const EquationSystemTimeData& GetTimeData() const;
  /**Sets the current time data.*/
  void SetTimeData(EquationSystemTimeData time_data);

  /**Returns the current equation terms-scope.*/
  EqTermScope EquationTermsScope() const;
  /**Sets the scope of current equations.*/
  void SetEquationTermsScope(EqTermScope eq_term_scope);
  /**Determines if a particular equation term is active.*/
  bool QueryTermsActive(EqTermScope query_scope) const;

  /**Returns a reference to the current solution vector.*/
  ParallelVector& SolutionVector(TimeID time_id = TimeID::T_PLUS_1);
  /**Returns a const reference to the current solution vector.*/
  const ParallelVector& SolutionVector(TimeID time_id = TimeID::T_PLUS_1) const;
  /**Returns a reference to the current residual vector.*/
  ParallelVector& ResidualVector(TimeID time_id = TimeID::T_PLUS_1);

  /**Sets the current solution vector.*/
  virtual void SetInitialSolution(){};

  /**Computes the residual vector \p r given a solution vector \p x.
   * This method is generally only called by executioners.*/
  virtual void ComputeResidual(const ParallelVector& x, ParallelVector& r) = 0;

  /**Computes the Jacobian matrix \p J given a solution vector \p x.
   * This method is generally only called by executioners.*/
  virtual void ComputeJacobian(const ParallelVector& x, ParallelMatrix& J) = 0;

  /**Uses the underlying system to build a sparsity pattern.*/
  virtual ParallelMatrixSparsityPattern BuildMatrixSparsityPattern() const = 0;

  /**Updates the fields.*/
  virtual void UpdateFields() {}

  /**Output fields to VTK. The filename passed via the options will be used
   * plus a time index (if transient).*/
  virtual void OutputFields(int time_index){};

  /**Advances the system in time.*/
  void Advance(EquationSystemTimeData time_data,
               std::map<TimeID, const ParallelVector*>& latest_std_residuals);

  /**General query for things like post-processors.*/
  virtual chi::ParameterBlock GetInfo(const chi::ParameterBlock& params) const;

protected:
  static chi::InputParameters GetInputParameters();
  explicit EquationSystem(const chi::InputParameters& params);

  const int verbosity_;
  std::shared_ptr<TimeIntegrator> time_integrator_;
  const size_t num_solution_histories_;
  const size_t num_residual_histories_;

  EquationSystemTimeData time_data_;

  std::unique_ptr<ParallelVector> main_solution_vector_;

  std::vector<std::unique_ptr<ParallelVector>> old_solution_vectors_;
  std::vector<std::unique_ptr<ParallelVector>> old_residual_vectors_;

  std::vector<std::function<void()>> callbacks_on_pre_advance_;

private:
  static std::shared_ptr<TimeIntegrator>
  InitTimeIntegrator(const chi::InputParameters& params);

  EqTermScope eq_term_scope_;
};

} // namespace chi_math

#endif // CHITECH_EQUATIONSYSTEMBASE_H
