CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(fsi_steady_cylinder_beam_3d_bdc)

SET(bdc_src cylinder_beam_3d_benchmark.cpp)

ADD_LIBRARY(${PROJECT_NAME} SHARED ${bdc_src})
