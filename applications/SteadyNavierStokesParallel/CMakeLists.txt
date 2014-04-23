CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(SteadyNavierStokes)

# Build the executable
ADD_EXECUTABLE(steadyns ${CMAKE_PROJECT_DIR}main.cpp)

# Link the executable to the petsc anf femttu libs
TARGET_LINK_LIBRARIES(steadyns femlib)
TARGET_LINK_LIBRARIES(steadyns ${PETSC_LIBRARIES})
TARGET_LINK_LIBRARIES(steadyns b64)

IF(MPI_FOUND)
  TARGET_LINK_LIBRARIES(steadyns ${MPI_EXTRA_LIBRARY})
ENDIF(MPI_FOUND)