CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(TimedepNavierStokes)

# Build the executable
ADD_EXECUTABLE(backeuler ${CMAKE_PROJECT_DIR}main.cpp)

# Link the executable to the petsc anf femttu libs
TARGET_LINK_LIBRARIES(backeuler femlib)
TARGET_LINK_LIBRARIES(backeuler ${PETSC_LIBRARIES})
TARGET_LINK_LIBRARIES(backeuler ${CMAKE_SOURCE_DIR}/ExternalLibraries/b64-1.4.2/lib/libb64.gcc47.a)

SET(CMAKE_BUILD_TYPE Debug)

IF(MPI_FOUND)
  #TARGET_LINK_LIBRARIES(steadyns ${MPI_LIBRARY})
  #TARGET_LINK_LIBRARIES(steadyns ${MPI_CXX_LINK_FLAGS})
  TARGET_LINK_LIBRARIES(backeuler ${MPI_EXTRA_LIBRARY})
ENDIF(MPI_FOUND)