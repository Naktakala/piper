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
  [heat_generation]
    type = BodyForce
    variable = T
    value = 100e2
  []
[]

[BCs]
  [bc1]
    type = DirichletBC
    variable = T
    value = 1
    boundary = 'left'
    # preset = False
  []
  [bc2]
    type = DirichletBC
    variable = T
    value = 1
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

[Postprocessors]
  [max_val]
    type = NodalExtremeValue
    variable = T
  []
[]

[Executioner]
  type = Steady

  solve_type = JFNK

  l_tol = 1e-05
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
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
