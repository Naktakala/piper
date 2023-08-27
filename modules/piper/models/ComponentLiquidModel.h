#ifndef PIPER_COMPONENTLIQUIDMODEL_H
#define PIPER_COMPONENTLIQUIDMODEL_H

#include "ComponentModel.h"
#include "parameters/input_parameters.h"

namespace piper
{

class LiquidPhysics;

class ComponentLiquidModel : public ComponentModel
{
public:
  /**This structure is used to assemble a local equation.*/
  struct EqCoeffs
  {
    double rhs_;
    std::vector<std::vector<double>> coeff_sets_;
    std::vector<std::vector<size_t>> id_maps_;
  };

  /**Assembles an momentum equation if defined.*/
  virtual void AssembleEquations();

  const EqCoeffs& GetEquationCoefficients(size_t id) const;
  EqCoeffs& GetEquationCoefficients(size_t id);
protected:
  static chi::InputParameters GetInputParameters();
  ComponentLiquidModel(const chi::InputParameters& params,
                       std::vector<std::unique_ptr<ComponentModel>>& family,
                       LiquidPhysics& liquid_physics,
                       const HardwareComponent& hardware_component,
                       const chi_mesh::Cell* cell,
                       const std::vector<std::string>& variable_names);

  LiquidPhysics& physics_;
  std::vector<EqCoeffs> equation_coefficients_;
};

} // namespace piper

#endif // PIPER_COMPONENTLIQUIDMODEL_H
