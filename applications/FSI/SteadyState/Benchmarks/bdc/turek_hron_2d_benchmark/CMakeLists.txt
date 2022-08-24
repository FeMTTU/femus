CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(fsi_steady_2d_turek_hron_benchmark_bdc)

SET(bdc_src 
fsi_2d_turek_hron_benchmark.cpp)

ADD_LIBRARY(${PROJECT_NAME} SHARED ${bdc_src})