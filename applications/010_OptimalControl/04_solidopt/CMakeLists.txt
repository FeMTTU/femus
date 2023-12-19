CMAKE_MINIMUM_REQUIRED(VERSION 3.5)


# Full Dirichlet solid mechanics problem
ADD_SUBDIRECTORY(./00_solid_state) 

# boundary control solid mechanics problem
ADD_SUBDIRECTORY(./01_solid_control_boundary) 

# distributed control solid mechanics problem
ADD_SUBDIRECTORY(./02_solid_control_lifting) 
