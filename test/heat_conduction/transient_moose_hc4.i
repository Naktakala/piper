[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 100
  []
  parallel_type = DISTRIBUTED
[]

[Variables]
  [T]
    order = FIRST
    family = LAGRANGE
    initial_condition = 300.0
  []
[]

[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = T
  []
  [heat_generation]
    type = BodyForce
    variable = T
    value = 100e4
  []
  [time_derivative]
    type = SpecificHeatConductionTimeDerivative
    variable = T
  []
[]

[BCs]
  [bc1]
    type = DirichletBC
    variable = T
    value = 293.15
    boundary = 'left bottom top right'
  []
[]

[Functions]
  [tc_zrh]
    type = PiecewiseLinear
    x = '0.0 200	250	300	350	400	450	500	550	600	650	700'
    y = '0.01 4.63	10.10	12.68	13.96	14.60	14.89	14.98	14.95	14.84	14.69	14.52'
    # x = '0.0 1000.0'
    # y = '12.0 12.0'
  []
  [cp_zrh]
    type = PiecewiseLinear
    x = '0.0  100.0  300.0 800.0'
    y = '1.0 0.18e3 0.32e3 0.65e3'
    # x = '0.0 1000.0'
    # y = '500.0 500.0'
  []
[]

[Materials]
  [k]
    type = HeatConductionMaterial
    temp = T
    thermal_conductivity_temperature_function = tc_zrh
    specific_heat_temperature_function = cp_zrh
  []
  [density]
    type = Density
    density = 5586.0
  []
[]

[Postprocessors]
  [maxval]
    type = NodalExtremeValue
    variable = T
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  # type = Steady

  # solve_type = PJFNK
  solve_type = NEWTON

  num_steps = 20
  dt = 20.0
  end_time = 20000.0

  [TimeIntegrator]
    type = ImplicitEuler
    # type = CrankNicolson
  []

  # l_tol = 1e-03
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
  # line_search = l2
  automatic_scaling = false
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 50 0.7'
  # petsc_options_iname = '-pc_type'
  # petsc_options_value = 'none'
  auto_preconditioning=false
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  perf_graph = true
[]
