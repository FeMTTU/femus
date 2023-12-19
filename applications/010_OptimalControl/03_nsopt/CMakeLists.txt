
CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

ADD_SUBDIRECTORY(./00_state) 

# NS opt with standard approach with Lag mult for compatibility constraint
ADD_SUBDIRECTORY(./01_opt_bdry_ctrl)

# NS opt with lifting approach
ADD_SUBDIRECTORY(./02_opt_lift)
