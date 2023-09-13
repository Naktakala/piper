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
[]

[BCs]
  [bc1]
    type = DirichletBC
    variable = T
    value = -1
    boundary = 'left'
    # preset = False
  []
  [bc2]
    type = DirichletBC
    variable = T
    value = 10
    boundary = 'right'
    # preset = False
  []
  [bc3]
    type = DirichletBC
    variable = T
    value = 20
    boundary = 'bottom'
    # preset = False
  []
  [bc4]
    type = DirichletBC
    variable = T
    value = 30
    boundary = 'top'
    # preset = False
  []
[]

[Materials]
  [k]
    type = HeatConductionMaterial
    temp = T
    thermal_conductivity = 16.0
  []
[]

[Executioner]
  type = Steady

  solve_type = PJFNK
  # solve_type = NEWTON

  # l_tol = 1e-03
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
  # line_search = l2
  automatic_scaling = false
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 50 0.7'
  # petsc_options_iname = '-pc_type'
  # petsc_options_value = 'none'
  # auto_preconditioning=false
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]
