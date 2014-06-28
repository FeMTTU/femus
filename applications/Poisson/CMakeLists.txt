CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(Poisson)

# Build the executable
ADD_EXECUTABLE(Poisson ${CMAKE_PROJECT_DIR}main.cpp)

# Link the executable to the petsc anf femttu libs
TARGET_LINK_LIBRARIES(Poisson femus)
TARGET_LINK_LIBRARIES(Poisson ${PETSC_LIBRARIES})
TARGET_LINK_LIBRARIES(Poisson b64)
TARGET_LINK_LIBRARIES(Poisson /home/simone/Software/jsoncpp-src-0.5.0/libs/linux-gcc-4.8/libjson_linux-gcc-4.8_libmt.so)

IF(MPI_FOUND)
  TARGET_LINK_LIBRARIES(Poisson ${MPI_EXTRA_LIBRARY})
ENDIF(MPI_FOUND)