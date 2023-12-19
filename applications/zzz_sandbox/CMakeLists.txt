CMAKE_MINIMUM_REQUIRED(VERSION 3.5)
  

  ADD_SUBDIRECTORY(ex_petscDebug_test_solvers_for_laplacian/)
  ADD_SUBDIRECTORY(ex40_other/) #erdi
# ADD_SUBDIRECTORY(ex06c_navier_stokes_poiseuille_finite_volume_petsc/)  #copied from petsc tutorial, should remove it

  ADD_SUBDIRECTORY(./class_2021_2022/)
  
  ADD_SUBDIRECTORY(./Karthik/)

  ADD_SUBDIRECTORY(./shankar/)
  ADD_SUBDIRECTORY(./shankar_biharmonic_coupled_systems/)
  ADD_SUBDIRECTORY(./shankar_biharmonic_uncoupled_systems/)
