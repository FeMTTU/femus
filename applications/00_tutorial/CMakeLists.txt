CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

get_filename_component(APP_FOLDER_NAME_PARENT ${CMAKE_CURRENT_LIST_DIR} NAME)

ADD_SUBDIRECTORY(00_laplacian_all_dims_dir_neu/)
ADD_SUBDIRECTORY(00_laplacian_all_dims_dir_neu_class_2021_2022/)


ADD_SUBDIRECTORY(ex01_function_approximation_no_eqn/)
ADD_SUBDIRECTORY(ex02_poisson_a_manufactured_sol/)
ADD_SUBDIRECTORY(ex02_poisson_b_manufactured_sol_incremental_conv_rate/)
ADD_SUBDIRECTORY(ex02_poisson_convergence/)

ADD_SUBDIRECTORY(ex03_poisson_plus_nonlin_advection/)
ADD_SUBDIRECTORY(ex03b_div_grad_nonlinear_diffusivity/)

ADD_SUBDIRECTORY(ex04_biharmonic_coupled_systems/)
ADD_SUBDIRECTORY(ex05_biharmonic_uncoupled_systems/)

ADD_SUBDIRECTORY(ex06_navier_stokes/)
ADD_SUBDIRECTORY(ex06b_navier_stokes/)
# ADD_SUBDIRECTORY(ex06c_navier_stokes_poiseuille/)  #copied from petsc tutorial, should remove it

ADD_SUBDIRECTORY(ex07_boussinesq/)
ADD_SUBDIRECTORY(ex08_boussinesq/)

ADD_SUBDIRECTORY(ex20_time_dependent/)

ADD_SUBDIRECTORY(ex40_other/)



