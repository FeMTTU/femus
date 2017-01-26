CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
  
# FSI Time Dependent Benchmarks
  ADD_SUBDIRECTORY(./TimeDependent/)
# FSI SteadyState Benchmarks
  ADD_SUBDIRECTORY(./SteadyState/)
 # FSI aneurysm
 ADD_SUBDIRECTORY(./aneurysm/)

 # FSI fsiOPT
  ADD_SUBDIRECTORY(./fsiOPT/)