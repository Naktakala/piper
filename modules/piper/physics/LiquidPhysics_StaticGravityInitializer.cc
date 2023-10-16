#include "LiquidPhysics.h"

#include "piper/Piper.h"
#include "piper/components/HardwareComponent.h"

#include "chi_log.h"

namespace piper
{

// ###################################################################
void LiquidPhysics::StaticGravityInitializer(const chi::ParameterBlock& params)
{
  typedef chi_mesh::Vector3 Vec3;

  //======================================== Ensure state has been set
  ChiInvalidArgumentIf(not params.Has("state"),
                       "Parameters supplied must have sub-parameter"
                       " \"state\".");
  const chi::ParameterBlock& state_params = params.GetParam("state");

  //======================================== Get pressure specification
  if (not state_params.Has("p"))
    ChiInvalidArgument("This initializer requires pressure, \"p\", to be "
                       "specified.");

  const double root_comp_p = state_params.GetParamValue<double>("p");

  Chi::log.Log() << "StaticGravityInitializer root pressure = " << root_comp_p;

  //======================================== Collect independent state
  // specifications
  std::vector<StateVal> auxiliary_state_specs;
  for (const auto& param : state_params)
  {
    if (param.Name() != "p")
    {
      auxiliary_state_specs.emplace_back(param.Name(),
                                         param.GetValue<double>());
      Chi::log.Log() << "Auxiliary state parameter " << param.Name() << "="
                     << param.GetValue<double>();
    }
  }

  //======================================== Construct root comp state spec
  auto root_comp_state_spec = auxiliary_state_specs;
  root_comp_state_spec.emplace_back("p", root_comp_p);

  //======================================== Get references
  const size_t root_comp_id = pipe_system_ptr_->RootComponentID();
  const ComponentModel& root_comp_model = *component_models_.at(root_comp_id);
  const Vec3& root_node_position = root_comp_model.GetRootNodePosition();

  //======================================== Compute root density
  // The fundamental assumption here is that
  // the density will not change
  const auto root_comp_state =
    EvaluateState({"rho", "e", "T", "p", "k", "mu", "beta"}, root_comp_state_spec);
  const double root_comp_rho = root_comp_state.at("rho");

  //======================================== Print the root state
  {
    std::stringstream outstr;
    for (const auto& val : root_comp_state)
      outstr << "  " << val.first << " " << val.second << "\n";
    Chi::log.Log0Verbose1()
      << "Root state:\n"
      << "  position " << root_node_position.PrintStr() << "\n"
      << outstr.str();
  }

  for (const auto& model_ptr : component_models_)
  {
    auto& model = *model_ptr;
    std::stringstream outstr;

    if (model.Category() == ComponentCategory::BoundaryLike or
        model.Category() == ComponentCategory::Volumetric)
    {
      const Vec3 hw_comp_centroid = model.MakeCentroid();
      const Vec3 delta_h_vec = hw_comp_centroid - root_node_position;

      const double delta_z = delta_h_vec.Dot(gravity_.Normalized());
      const double g = gravity_.Norm();

      double p = 0.0;
      double avg_rho = root_comp_rho;
      for (size_t k = 0; k < 100; ++k)
      {
        p = root_comp_p + avg_rho * delta_z * g;

        auto perturbed_state = auxiliary_state_specs;
        perturbed_state.emplace_back("p", p);

        const double end_rho =
          EvaluateState({"rho"}, perturbed_state).at("rho");

        const double new_avg_rho = 0.5 * (end_rho + root_comp_rho);
        const double drho_rel_rho = std::fabs(new_avg_rho - avg_rho);

        avg_rho = new_avg_rho;

        if (drho_rel_rho < 1.0e-12) break;
      }

      const double adjusted_p = p;

      auto aux_specs = auxiliary_state_specs;
      aux_specs.emplace_back("p", adjusted_p);

      outstr << "Initialized pressure in \""
             << model.GetHardwareComponent().Name() << "\" with ";
      for (const auto& [spec_name, val] : aux_specs)
        outstr << "\"" << spec_name << "\"=" << val << " ";

      const std::vector<std::string> var_names = {
        "rho", "e", "T", "p", "k", "mu", "beta"};
      const auto state = EvaluateState(var_names, aux_specs);

      for (const auto& var_name : var_names)
      {
        const double state_var_val = state.at(var_name);
        model.VarOld(var_name) = state_var_val;
        model.VarNew(var_name) = state_var_val;
      }
    }
    else if (model.Category() == ComponentCategory::JunctionLike)
    {
      model.VarOld("u") = 0.0;
    }
    else
      ChiLogicalError("Unhandled component category");

    Chi::log.Log0Verbose1() << outstr.str();
  }

  Chi::mpi.Barrier();
}

} // namespace piper
