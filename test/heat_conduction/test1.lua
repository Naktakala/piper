znodes = {}
N = 10
zmin = 0.0
L = 1.0
znodes = {}

dx = L / N
for i = 0, N do
  znodes[i + 1] = zmin + dx * i
end
meshgen1 = chi_mesh.OrthogonalMeshGenerator.Create({ node_sets = { znodes } })
chi_mesh.MeshGenerator.Execute(meshgen1)

system1 = chi_math.KernelSystem.Create
({
  fields = {
    chi_physics.FieldFunctionGridBased.Create({
      name = "T",
      sdm_type = "PWLC"
    })
  },
  kernels = {
    { type = hcm.ThermalConductionKernel.type, var = "T", k = 16.0 }
  },
  bcs = {
    { type = chi_math.FEMDirichletBC.type, var = "T", boundaries = { "ZMAX" }, bc_value = -1.0 },
    { type = chi_math.FEMDirichletBC.type, var = "T", boundaries = { "ZMIN" }, bc_value = 10.0 }
  }
  --verbosity = 2
})

phys1 = chi_math.SteadyNonLinearExecutioner.Create
({
  name = "phys1",
  system = system1,
  solver_params =
  {
    nl_method = "PJFNK",
    l_rel_tol = 1.0e-5,
  }
})

chi.AggregateNodalValuePostProcessor.Create
({
  name = "maxval",
  field_function = "T",
  operation = "max",
  print_on = {"ProgramExecuted"}
})

chi.PostProcessorPrinterSetOptions
(
  {
    print_scalar_time_history = false
  }
)

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

if (master_export == nil) then
  chiExportMultiFieldFunctionToVTK({ "T" }, "test1")
end

chiLogPrintTimingGraph()