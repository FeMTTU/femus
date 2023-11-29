CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

get_filename_component(APP_FOLDER_NAME_PARENT ${CMAKE_CURRENT_LIST_DIR} NAME)



# Steady-state - BEGIN

# Steady-state Poisson, No Field Split - BEGIN
# all these are very very sim
ADD_SUBDIRECTORY(00_poisson_ex40_1/)  #@todo, some runtime error
ADD_SUBDIRECTORY(00_poisson_ex40_2/)  #@todo, some runtime error
ADD_SUBDIRECTORY(00_poisson_ex41/)
ADD_SUBDIRECTORY(00_poisson_ex42/)
ADD_SUBDIRECTORY(00_poisson_ex43/)
# Steady-state Poisson, No Field Split - END


# Steady-state NS - BEGIN
# all these are very very sim

# No Field Split - BEGIN
ADD_SUBDIRECTORY(ex12/)  # square  #@todo, some runtime error
ADD_SUBDIRECTORY(ex44/)  # with pressure stabilization   #@todo, some runtime error
ADD_SUBDIRECTORY(ex44_ns_pcd/)  # only initialized
# No Field Split - END

# Field Split - BEGIN
ADD_SUBDIRECTORY(ex13/)  # square
ADD_SUBDIRECTORY(ex14/)  # square
# Field Split - END

# Steady-state NS - END




# Steady-state Boussinesq  - BEGIN

# Steady-state Boussinesq, Field Split  - BEGIN
ADD_SUBDIRECTORY(ex1/)     # square
ADD_SUBDIRECTORY(ex10/)    # square
ADD_SUBDIRECTORY(ex11/)    # square
ADD_SUBDIRECTORY(ex11_1/)  # square
ADD_SUBDIRECTORY(ex11_2/)  # square16
# Steady-state Boussinesq, Field Split  - END

# Steady-state Boussinesq, No Field Split - BEGIN
ADD_SUBDIRECTORY(ex46/)    # cube   # extension of ex44
# Steady-state Boussinesq, No Field Split - END

# Steady-state Boussinesq  - END



# Steady-state MHD, Field Split  - BEGIN
ADD_SUBDIRECTORY(ex50_mhd/)   # square
ADD_SUBDIRECTORY(ex50_1_mhd/)  #very very sim, with manufactured sol
# Steady-state MHD, Field Split  - END



# Steady-state Optimal control, Field Split  - BEGIN
ADD_SUBDIRECTORY(ex60_opt_control/)
# Steady-state Optimal control, Field Split  - END

# Steady-state - END




# Time dependent - BEGIN

# Time dependent, Boussinesq, Field Split - BEGIN
ADD_SUBDIRECTORY(ex20_1/)  # square
ADD_SUBDIRECTORY(ex20_2/)  # rectangle
ADD_SUBDIRECTORY(ex21/)    # rectangle8x64
ADD_SUBDIRECTORY(ex21_1/)  # rectangle8x64
ADD_SUBDIRECTORY(ex21_2/)  # rectangle8x64
ADD_SUBDIRECTORY(ex22_1/)  # rectangle_w4_h1
ADD_SUBDIRECTORY(ex23/)    # rectangle_w4_h1
ADD_SUBDIRECTORY(ex24/)    # rectangle_w1_h8
ADD_SUBDIRECTORY(ex25/)    # box_4_1_1

# Time dependent, Boussinesq, Field Split - END


# Time dependent - END
