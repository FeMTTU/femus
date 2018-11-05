
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#ns_poiseuille - Poiseuille flow in a 2d box
ADD_SUBDIRECTORY(./ns_poiseuille) 

#ns_fulldir - Full Dirichlet NS problem
ADD_SUBDIRECTORY(./ns_fulldir) 

# NS opt with lifting approach
ADD_SUBDIRECTORY(./opt_lift)

# NS opt with standard approach with Lag mult for compatibility constraint
ADD_SUBDIRECTORY(./opt_bdry_ctrl)
