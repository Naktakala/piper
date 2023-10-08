-- Cylindrical coordinate diffusion with convection boundary condition
-- Analytical solution Tmax = 0.25 * q/k * R^2 + Tb + 0.5 * R^2 * q/h
-- With k=2, q=1e4, Tb=300, R=1.35, h=100.0
-- Tmax = 2669.25

rnodes = {}
znodes = {}

Nz = 1
zmin = 0.0
L = 100

Nr = 400
rmin = 0.0
R = 1.35

dz = L / Nz
for i = 0, Nz do
  znodes[i + 1] = zmin + dz * i
end

dr = R / Nr
for i = 0, Nr do
  rnodes[i + 1] = rmin + dr * i
end

meshgen1 = chi_mesh.OrthogonalMeshGenerator.Create({ node_sets = { rnodes, znodes } })
chi_mesh.MeshGenerator.Execute(meshgen1)

system1 = chi_math.KernelSystem.Create
({
  fields = {
    chi_physics.FieldFunctionGridBased.Create({
      name = "T",
      sdm_type = "LagrangeC",
      initial_value = 300.0,
      coordinate_system = "rz"
    })
  },
  kernels = {
    {type = hcm.ThermalConductionKernel.type, var="T", k = 2.0 },
    {type = chi_math.SinkSourceFEMKernel.type, var="T", value = 100.0e2 },
    {type = chi_math.FEMTimeDerivativeKernel.type, var="T",}
  },
  bcs = {
    {
      type = hcm.ConvectiveHeatFluxBC.type,
      boundaries = { "XMAX" },
      T_bulk = 300.0,
      convection_coefficient = 100.0,
      var="T"
    }
    --{
    --  type = chi_math.FEMDirichletBC.type,
    --  boundaries = {"XMAX"},
    --  bc_value = 300.0,
    --  var = "T"
    --}
  },
  time_integrator = chi_math.CrankNicolsonTimeIntegrator.Create({}),
  --verbosity = 2
  --output_filename_base = "transient_cylTest1"
})

phys1 = chi_math.TransientNonLinearExecutioner.Create
--phys1 = chi_math.SteadyNonLinearExecutioner.Create
({
  name = "phys1",
  system = system1,
  solver_params =
  {
    nl_method = "NEWTON",
    l_rel_tol = 1.0e-5,
  },

  dt = 0.1,
  end_time = 10.0,

  print_header = false,
  print_footer = false,
  print_nl_residual = false,
  print_l_residual = false,
})

chi.AggregateNodalValuePostProcessor.Create
({
  name = "maxval",
  field_function = "T",
  operation = "max",
  print_on = {"ProgramExecuted"}
})

chi.PostProcessorPrinterSetOptions
({
  print_scalar_time_history = false,
  --csv_filename = "transient_cylTest1.csv"
})

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

