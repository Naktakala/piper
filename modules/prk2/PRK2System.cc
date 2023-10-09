#include "PRK2System.h"

#include "math/ParallelVector/ParallelSTLVector.h"
#include "math/ParallelMatrix/ParallelMatrix.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

#include <numeric>

#define scint64_t static_cast<int64_t>

namespace prk2
{

RegisterChiObject(prk2, PRK2System);

chi::InputParameters PRK2System::GetInputParameters()
{
  chi::InputParameters params = chi_math::EquationSystem::GetInputParameters();

  std::vector<double> default_lambdas = {
    0.0124, 0.0304, 0.111, 0.301, 1.14, 3.01};
  std::vector<double> default_betas = {
    0.00021, 0.00142, 0.00127, 0.00257, 0.00075, 0.00027};

  params.AddOptionalParameterArray(
    "precursor_lambdas", default_lambdas, "An array of decay constants");
  params.AddOptionalParameterArray(
    "precursor_betas",
    default_betas,
    "An array of fractional delayed neutron fractions");

  params.AddOptionalParameter(
    "gen_time", 1.0e-5, "Neutron generation time [s]");
  params.AddOptionalParameter("initial_rho", 0.0, "Initial reactivity [$]");
  params.AddOptionalParameter(
    "initial_source", 1.0, "Initial source strength [/s]");

  params.AddOptionalParameter(
    "initial_population", 1.0, "Initial neutron population");

  using namespace chi_data_types;
  params.ConstrainParameterRange("gen_time",
                                 AllowableRangeLowLimit::New(1.0e-12));
  params.ConstrainParameterRange("initial_source",
                                 AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("initial_population",
                                 AllowableRangeLowLimit::New(0.0));

  return params;
}

PRK2System::PRK2System(const chi::InputParameters& params)
  : chi_math::EquationSystem(params),
    lambdas_(params.GetParamVectorValue<double>("precursor_lambdas")),
    betas_(params.GetParamVectorValue<double>("precursor_betas")),
    beta_(std::accumulate(betas_.begin(), betas_.end(), /*init_val=*/0.0)),
    gen_time_(params.GetParamValue<double>("gen_time")),
    rho_(params.GetParamValue<double>("initial_rho")),
    source_strength_(params.GetParamValue<double>("initial_source")),
    num_precursors_(lambdas_.size()),
    num_local_dofs_(Chi::mpi.location_id == 0 ? scint64_t(num_precursors_ + 1)
                                              : 0),
    num_global_dofs_(scint64_t(num_precursors_ + 1)),
    A_(num_local_dofs_, num_local_dofs_, 0.0),
    x_initial_(num_local_dofs_, 0.0),
    q_(num_local_dofs_, 0.0)
{
  main_solution_vector_ = std::make_unique<chi_math::ParallelSTLVector>(
    num_local_dofs_, num_global_dofs_, Chi::mpi.comm);

  for (size_t t = 0; t < num_solution_histories_; ++t)
  {
    auto vec = main_solution_vector_->MakeClone();
    old_solution_vectors_.push_back(std::move(vec));
  }

  for (size_t t = 0; t < num_residual_histories_; ++t)
  {
    auto vec = std::make_unique<chi_math::ParallelSTLVector>(
      num_local_dofs_, num_global_dofs_, Chi::mpi.comm);
    old_residual_vectors_.push_back(std::move(vec));
  }

  callbacks_on_pre_advance_.emplace_back(
    std::bind(&PRK2System::PreAdvanceCallback, this));

  // Precompute most of the system elements
  if (Chi::mpi.location_id == 0) InitializeSystem();

  MPI_Bcast(&population_,   // send/rcv buffer
            1,              // count
            MPI_DOUBLE,     // datatype
            0,              // root
            Chi::mpi.comm); // communicator
}

void PRK2System::InitializeSystem()
{
  const size_t J = num_precursors_;

  A_[0][0] = beta_ * (rho_ - 1.0) / gen_time_;

  for (size_t j = 1; j <= J; ++j)
  {
    A_[0][j] = lambdas_[j - 1];
    A_[j][j] = -lambdas_[j - 1];
    A_[j][0] = betas_[j - 1] / gen_time_;
  }

  q_.resize(J + 1, 0.0);
  q_[0] = source_strength_;

  // ============================= Initializing x
  // If there is a source and the reactivity is < 0 then
  // there exists a unique solution.
  if (source_strength_ > 0.0 and rho_ < 0.0)
  {
    const auto b_theta = -1.0 * q_;

    x_initial_ = A_.Inverse() * b_theta;
  }
  // Otherwise we initialize the system as a critical system with
  // no source.
  else
  {
    auto A_temp = A_;
    auto b_temp = x_initial_;
    b_temp.Set(0.0);
    b_temp[0] = 1.0;
    for (auto& val : A_temp[0])
      val = 0.0;
    A_temp[0][0] = 1.0;

    x_initial_ = A_temp.Inverse() * b_temp;
  }

  Chi::log.Log() << "PRK Initial solution: " << x_initial_.PrintStr();
}

chi_math::ParallelMatrixSparsityPattern
PRK2System::BuildMatrixSparsityPattern() const
{
  chi_math::ParallelMatrixSparsityPattern pattern;
  if (Chi::mpi.location_id == 0)
  {
    pattern.nodal_nnz_in_diag_ =
      std::vector<int64_t>(num_precursors_ + 1, scint64_t(num_precursors_ + 1));
    pattern.nodal_nnz_off_diag_ = std::vector<int64_t>(num_precursors_ + 1);
  }

  return pattern;
}

void PRK2System::SetInitialSolution()
{
  main_solution_vector_->Set(x_initial_.elements_);

  if (num_solution_histories_ > 0)
    for (size_t t = 0; t < num_solution_histories_; ++t)
      old_solution_vectors_[t] = main_solution_vector_->MakeCopy();
}

void PRK2System::ComputeResidual(const chi_math::ParallelVector& x,
                                 chi_math::ParallelVector& r)
{

  if (Chi::mpi.location_id == 0)
  {
    if (QueryTermsActive(chi_math::EqTermScope::DOMAIN_TERMS))
    {
      A_[0][0] = beta_ * (rho_ - 1.0) / gen_time_;
      for (int64_t i = 0; i < num_local_dofs_; ++i)
      {
        double r_i = -q_[i];
        for (int64_t j = 0; j < num_local_dofs_; ++j)
          r_i -= A_[i][j] * x[j];

        r[i] = r_i;
      }
    }

    if (QueryTermsActive(chi_math::EqTermScope::TIME_TERMS))
    {
      const auto& x_old = SolutionVector(chi_math::TimeID::T);
      for (int64_t i = 0; i < num_local_dofs_; ++i)
        r[i] += (x[i] - x_old[i]) * time_data_.var_dot_dvar_;
    }
  }

  r.Assemble();
}

void PRK2System::ComputeJacobian(const chi_math::ParallelVector& x,
                                 chi_math::ParallelMatrix& J)
{
  if (Chi::mpi.location_id == 0)
  {
    MatDbl cell_A(num_local_dofs_, VecDbl(num_local_dofs_, 0.0));

    if (QueryTermsActive(chi_math::EqTermScope::DOMAIN_TERMS))
    {
      A_[0][0] = beta_ * (rho_ - 1.0) / gen_time_;
      for (int64_t i = 0; i < num_local_dofs_; ++i)
        for (int64_t j = 0; j < num_local_dofs_; ++j)
          cell_A[i][j] -= A_[i][j];
    }

    if (QueryTermsActive(chi_math::EqTermScope::TIME_TERMS))
      for (int64_t i = 0; i < num_local_dofs_; ++i)
        cell_A[i][i] += time_data_.var_dot_dvar_;

    for (int64_t i = 0; i < num_local_dofs_; ++i)
      for (int64_t j = 0; j < num_local_dofs_; ++j)
        J.AddValue(i, j, cell_A[i][j]);
  }

  J.Assemble(/*final=*/true);
}

void PRK2System::PreAdvanceCallback()
{
  std::array<double, 2> population_period = {0.0, 0.0};
  if (Chi::mpi.location_id == 0)
  {
    const auto& x_tp1 = SolutionVector(chi_math::TimeID::T_PLUS_1);
    const auto& x_t = SolutionVector(chi_math::TimeID::T);
    const double dt = time_data_.dt_;

    period_tph_ = dt / log(x_tp1[0] / x_t[0]);

    if (period_tph_ > 0.0 and period_tph_ > 1.0e6) period_tph_ = 1.0e6;
    if (period_tph_ < 0.0 and period_tph_ < -1.0e6) period_tph_ = -1.0e6;

    population_period[0] = x_tp1[0];
    population_period[1] = period_tph_;
  }

  MPI_Bcast(population_period.data(), // send/rcv buffer
            2,                        // count
            MPI_DOUBLE,               // datatype
            0,                        // root
            Chi::mpi.comm);           // communicator

  population_ = population_period[0];
  period_tph_ = population_period[1];
}

void PRK2System::SetProperties(const chi::ParameterBlock& params)
{
  for (const auto& param : params)
  {
    if (param.Name() == "rho") rho_ = param.GetValue<double>();
  }
}

chi::ParameterBlock PRK2System::GetInfo(const chi::ParameterBlock& params) const
{
  const auto param_name = params.GetParamValue<std::string>("name");

  if (param_name == "period") return chi::ParameterBlock("", period_tph_);
  else if (param_name == "population")
    return chi::ParameterBlock("", population_);
  else
    ChiInvalidArgument("Unsupported info name \"" + param_name + "\".");
}

} // namespace prk2