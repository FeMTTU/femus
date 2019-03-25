CMAKE_MINIMUM_REQUIRED(VERSION 2.8)


# Full Dirichlet solid mechanics problem
ADD_SUBDIRECTORY(./state_fulldir) 

# distributed control solid mechanics problem
ADD_SUBDIRECTORY(./solid_dist_ctrl) 
