
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

ADD_SUBDIRECTORY(./ns_state) 

# NS opt with lifting approach
ADD_SUBDIRECTORY(./opt_lift)

# NS opt with standard approach with Lag mult for compatibility constraint
ADD_SUBDIRECTORY(./opt_bdry_ctrl)
