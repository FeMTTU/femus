CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

set(THIS_APPLICATION gencase)

PROJECT(${THIS_APPLICATION})

set(CMAKE_VERBOSE_MAKEFILE on)

ADD_EXECUTABLE(${THIS_APPLICATION}
${CMAKE_PROJECT_DIR}main.cpp 
${CMAKE_PROJECT_DIR}GenCase.cpp
)

# have to reproduce this in CMake syntax
# libmesh_INCLUDE  := $(shell ${FM_LM_DIR}/bin/libmesh-config --include)
# libmesh_CPPFLAGS := $(shell ${FM_LM_DIR}/bin/libmesh-config --cppflags)
# libmesh_LIBS     := $(shell ${FM_LM_DIR}/bin/libmesh-config --libs)

execute_process(COMMAND $ENV{FM_BASEPATH_TO_LM}/$ENV{FM_LM_FOLDER}/bin/libmesh-config --cppflags
                COMMAND tr "\n" " " OUTPUT_VARIABLE LIBMESH_CPP)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBMESH_CPP}")
execute_process(COMMAND $ENV{FM_BASEPATH_TO_LM}/$ENV{FM_LM_FOLDER}/bin/libmesh-config --include
                COMMAND tr "\n" " " OUTPUT_VARIABLE LIBMESH_INC)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBMESH_INC}")
execute_process(COMMAND $ENV{FM_BASEPATH_TO_LM}/$ENV{FM_LM_FOLDER}/bin/libmesh-config --libs
                COMMAND tr "\n" " " OUTPUT_VARIABLE LIBMESH_LIBS)
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LIBMESH_LIBS}")


TARGET_LINK_LIBRARIES(${THIS_APPLICATION} femlib)
TARGET_LINK_LIBRARIES(${THIS_APPLICATION} ${PETSC_LIBRARIES})
TARGET_LINK_LIBRARIES(${THIS_APPLICATION} ${CMAKE_SOURCE_DIR}/ExternalLibraries/b64-1.4.2/lib/libb64.gcc47.a)

IF(MPI_FOUND)
  TARGET_LINK_LIBRARIES(${THIS_APPLICATION} ${MPI_EXTRA_LIBRARY})
ENDIF(MPI_FOUND)