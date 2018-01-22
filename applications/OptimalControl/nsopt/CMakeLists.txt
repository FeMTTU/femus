
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#ns_poiseuille - Poiseuille flow in a 2d box
ADD_SUBDIRECTORY(./ns_poiseuille) 

#ns_fulldir - Full Dirichlet NS problem
ADD_SUBDIRECTORY(./ns_fulldir) 

#Stokes_dist_adj - Stokes equation just with 2x2 block
ADD_SUBDIRECTORY(./Stokes_dist_adj)

#testOPT - Stokes equation with distributed control - lifting - 3x3 block
ADD_SUBDIRECTORY(./opt_lift)

#control is with p
ADD_SUBDIRECTORY(./control_p)

#compatibility condition on boundary _boundary integral
ADD_SUBDIRECTORY(./opt_bdry_ctrl)