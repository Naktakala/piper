[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 30
    ny = 30
    nz = 30
    xmax = 2.0
    ymax = 2.0
    zmax = 2.0
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
    type = Diffusion
    variable = T
  []
  [heat_generation]
    type = BodyForce
    variable = T
    value = 1.0
  []
[]

[BCs]
  [bc1]
    type = DirichletBC
    variable = T
    value = 0
    boundary = 'left right bottom top front back'
    # preset = False
  []
[]

[Postprocessors]
  [maxval]
    type = NodalExtremeValue
    variable = T
  []
[]

[Executioner]
  type = Steady

  solve_type = LINEAR
  # solve_type = NEWTON

  # l_tol = 1e-03
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
  # line_search = l2
  automatic_scaling = false
  petsc_options_iname = '-pc_type -pc_hypre_type
  -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold
  -pc_mg_galerkin_mat_product_algorithm'
  petsc_options_value = 'hypre boomeramg 50 0.7 hypre'
  # petsc_options_iname = '-pc_type'
  # petsc_options_value = 'none'
  # auto_preconditioning=false
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]
