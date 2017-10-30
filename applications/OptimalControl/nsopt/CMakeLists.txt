
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#testNS - copy of Ex6 with 2D Box b.c.
ADD_SUBDIRECTORY(./testNS) 

#Stokes_dist_adj - Stokes equation just with 2x2 block
ADD_SUBDIRECTORY(./Stokes_dist_adj)

#testOPT - Stokes equation with distributed control - lifting - 3x3 block
ADD_SUBDIRECTORY(./opt_lift)

#control is with p
ADD_SUBDIRECTORY(./control_p)
