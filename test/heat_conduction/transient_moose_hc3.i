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
    initial_condition = 0.0
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
    value = 100e2
  []
  [time_derivative]
    type = TimeDerivative
    variable = T
  []
[]

[BCs]
  [bc1]
    type = DirichletBC
    variable = T
    value = 0
    boundary = 'left bottom top right'
  []
  # [bc2]
  #   type = ConvectiveHeatFluxBC
  #   variable = T
  #   T_infinity = 100
  #   heat_transfer_coefficient = 10000
  #   boundary = 'right'
  # []
[]

[Materials]
  [k]
    type = HeatConductionMaterial
    temp = T
    thermal_conductivity = 16.0
  []
[]

[Executioner]
  type = Transient

  solve_type = PJFNK
  # solve_type = NEWTON

  num_steps = 2
  dt = 0.001

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
[]
