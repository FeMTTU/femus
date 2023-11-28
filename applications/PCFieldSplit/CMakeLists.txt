CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

get_filename_component(APP_FOLDER_NAME_PARENT ${CMAKE_CURRENT_LIST_DIR} NAME)



# Steady-state - BEGIN

# Steady-state Poisson, No Field Split - BEGIN
ADD_SUBDIRECTORY(ex40_1/) 
ADD_SUBDIRECTORY(ex40_2/) 
ADD_SUBDIRECTORY(ex41/)   
ADD_SUBDIRECTORY(ex42/)   
ADD_SUBDIRECTORY(ex43/)   
# Steady-state Poisson, No Field Split - END


# Steady-state NS - BEGIN

# Steady-state NS, Field Split - BEGIN
ADD_SUBDIRECTORY(ex12/)
ADD_SUBDIRECTORY(ex13/)
ADD_SUBDIRECTORY(ex14/)
ADD_SUBDIRECTORY(ex_ns_pcd/)  #wip
# Steady-state NS, Field Split - END

# Steady-state NS, No Field Split - BEGIN
ADD_SUBDIRECTORY(ex44/)
# Steady-state NS, No Field Split - END

# Steady-state NS - END


# Steady-state Boussinesq  - BEGIN

# Steady-state Boussinesq, Field Split  - BEGIN
ADD_SUBDIRECTORY(ex1/)

ADD_SUBDIRECTORY(ex10/)
ADD_SUBDIRECTORY(ex11/)
ADD_SUBDIRECTORY(ex11_1/)
ADD_SUBDIRECTORY(ex11_2/) ## @todo see why it is all different, must be EOF or similar 
# Steady-state Boussinesq, Field Split  - END

# Steady-state Boussinesq, No Field Split - BEGIN
ADD_SUBDIRECTORY(ex46/)  # extension of ex44
# Steady-state Boussinesq, No Field Split - END

# Steady-state Boussinesq  - END



# Steady-state MHD, Field Split  - BEGIN
ADD_SUBDIRECTORY(ex50_mhd/)
ADD_SUBDIRECTORY(ex50_1_mhd/)  #very very sim, with manufactured sol
# Steady-state MHD, Field Split  - END



# Steady-state Optimal control, Field Split  - BEGIN
ADD_SUBDIRECTORY(ex60_opt_control/)
# Steady-state Optimal control, Field Split  - END

# Steady-state - END




# Time dependent - BEGIN

# Time dependent, Boussinesq, Field Split - BEGIN
ADD_SUBDIRECTORY(ex20_1/)
ADD_SUBDIRECTORY(ex20_2/)
ADD_SUBDIRECTORY(ex21/)
ADD_SUBDIRECTORY(ex21_1/)
ADD_SUBDIRECTORY(ex21_2/)
ADD_SUBDIRECTORY(ex22_1/)
ADD_SUBDIRECTORY(ex23/)
ADD_SUBDIRECTORY(ex24/)
ADD_SUBDIRECTORY(ex25/)
# Time dependent, Boussinesq, Field Split - END


# Time dependent - END

