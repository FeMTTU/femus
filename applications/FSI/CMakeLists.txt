
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
  
# FSI Time Dependent Benchmarks
  ADD_SUBDIRECTORY(./TimeDependent/)
# FSI SteadyState Benchmarks
  ADD_SUBDIRECTORY(./SteadyState/)
# FSI SteadyState Benchmark Turek3
  ADD_SUBDIRECTORY(./FSI_TUREK_3/)
# FSI SteadyState in porous media
 ADD_SUBDIRECTORY(./fsi_porous/)
