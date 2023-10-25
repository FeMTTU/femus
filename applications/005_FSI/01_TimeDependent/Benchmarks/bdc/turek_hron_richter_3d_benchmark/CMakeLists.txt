CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(fsi_td_3d_turek_hron_richter_benchmark_bdc)

SET(bdc_src 
fsi_td_3d_turek_hron_richter_benchmark.cpp)

ADD_LIBRARY(${PROJECT_NAME} SHARED ${bdc_src})