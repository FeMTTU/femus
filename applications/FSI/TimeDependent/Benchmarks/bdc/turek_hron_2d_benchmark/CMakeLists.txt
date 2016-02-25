CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(fsi_td_2d_turek_hron_benchmark_bdc)

SET(bdc_src 
fsi_td_2d_turek_hron_benchmark.cpp)

ADD_LIBRARY(${PROJECT_NAME} SHARED ${bdc_src})