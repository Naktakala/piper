#ifndef CHITECH_FIELDEQUATIONSYSTEM_H
#define CHITECH_FIELDEQUATIONSYSTEM_H

#include "EquationSystem.h"

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

/**Base abstract class for a system of equations.*/
class FieldEquationSystem : public EquationSystem
{
public:
  typedef chi_mesh::Vector3 Vec3;
  typedef std::vector<Vec3> VecVec3;

  /**Returns the number of local DOFs across all unknowns.*/
  int64_t NumLocalDOFs() const override;
  /**Returns the number of local DOFs across all unknowns.*/
  int64_t NumGlobalDOFs() const override;

  /**Uses the underlying system to build a sparsity pattern.*/
  ParallelMatrixSparsityPattern BuildMatrixSparsityPattern() const override;

  /**Updates the fields.*/
  void UpdateFields() override;

  /**Output fields to VTK. The filename passed via the options will be used
   * plus a time index (if transient).*/
  void OutputFields(int time_index) override;

protected:
  static chi::InputParameters GetInputParameters();
  explicit FieldEquationSystem(const chi::InputParameters& params);

  const chi::MaterialPropertiesData& material_properties_data_;
  const std::string output_file_base_;

  std::shared_ptr<MultiFieldContainer> primary_fields_container_;

  const int64_t num_local_dofs_;
  const int64_t num_globl_dofs_;

  std::vector<size_t> t_tags_;

private:
  static const chi::MaterialPropertiesData&
  GetOrMakeMaterialPropertiesData(const chi::InputParameters& params);

  static std::shared_ptr<MultiFieldContainer>
  MakeMultifieldContainer(const chi::ParameterBlock& params);

  std::unique_ptr<ParallelVector> MakeSolutionVector();

  std::vector<std::unique_ptr<ParallelVector>> InitSolutionHistory();
  std::vector<std::unique_ptr<ParallelVector>> InitResidualHistory();


};

} // namespace chi_math

#endif // CHITECH_FIELDEQUATIONSYSTEM_H
