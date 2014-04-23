CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(TimedepNavierStokes)

# Build the executable
ADD_EXECUTABLE(backeuler ${CMAKE_PROJECT_DIR}main.cpp)

# Link the executable to the petsc anf femttu libs
TARGET_LINK_LIBRARIES(backeuler femus)
TARGET_LINK_LIBRARIES(backeuler ${PETSC_LIBRARIES})
TARGET_LINK_LIBRARIES(backeuler b64)

IF(MPI_FOUND)
  TARGET_LINK_LIBRARIES(backeuler ${MPI_EXTRA_LIBRARY})
ENDIF(MPI_FOUND)