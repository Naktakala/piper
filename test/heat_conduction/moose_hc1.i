[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 10
  []
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

  # l_tol = 1e-03
  # nl_abs_tol = 1e-4
  # nl_rel_tol = 1e-4
  # line_search = l2
  automatic_scaling = false
  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  # petsc_options_value = 'hypre boomeramg 50 0.7'
  auto_preconditioning=false
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]
