CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(femttu)

# The FEMTTU library version number.
SET(FEMTTU_VERSION_MAJOR 1)
SET(FEMTTU_VERSION_MINOR 0)
SET(FEMTTU_BUILD_VERSION 0)
SET(FEMTTU_VERSION
    "${FEMTTU_MAJOR_VERSION}.${FEMTTU_MINOR_VERSION}.${FEMTTU_BUILD_VERSION}")

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")

# If the user specifies -DCMAKE_BUILD_TYPE on the command line, take their definition
# and dump it in the cache along with proper documentation, otherwise set CMAKE_BUILD_TYPE
# to Debug prior to calling PROJECT()
IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ENDIF()

# where to place output binaries
SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib64/ CACHE INTERNAL
     "Single output directory for building all libraries.")
# SET (EXECUTABLE_OUTPUT_PATH ${MAF_BINARY_DIR}/bin/ CACHE INTERNAL 
#      "Single output directory for building all executables.")

MARK_AS_ADVANCED (
  LIBRARY_OUTPUT_PATH 
  #EXECUTABLE_OUTPUT_PATH
)

# Find Petsc Library 
FIND_PACKAGE(PETSc REQUIRED)

IF(PETSC_FOUND)
  SET(HAVE_PETSC 1)
ENDIF(PETSC_FOUND) 

# Find MPI (optional) 
FIND_PACKAGE(MPI)

IF(MPI_FOUND)
  SET(HAVE_MPI 1)
ENDIF(MPI_FOUND) 

# configure a header file to pass some of the CMake settings
# to the source code
CONFIGURE_FILE(
  "${PROJECT_SOURCE_DIR}/config/FEMTTUConfig.h.in"
  "${PROJECT_BINARY_DIR}/config/FEMTTUConfig.h" @ONLY IMMEDIATE
  )

#-- Add an Option to toggle the generation of the API documentation
OPTION(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
if(BUILD_DOCUMENTATION)
  FIND_PACKAGE(Doxygen)
  if (NOT DOXYGEN_FOUND)
    message(FATAL_ERROR 
      "Doxygen is needed to build the documentation. Please install it correctly")
  endif()
  #-- Configure the Template Doxyfile for our specific project
  CONFIGURE_FILE(
    "${PROJECT_SOURCE_DIR}/Doxygen/DoxyfileFemTTUlib.in"
    "${PROJECT_BINARY_DIR}/DoxygenDocumentation/DoxyfileFemTTUlib" @ONLY IMMEDIATE
  ) 

  #-- Add a custom target to run Doxygen when ever the project is built
  add_custom_target (Docs ALL 
    COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}//DoxygenDocumentation/DoxyfileFemTTUlib
    SOURCES ${PROJECT_BINARY_DIR}//DoxygenDocumentation/DoxyfileFemTTUlib)
  # IF you do NOT want the documentation to be generated EVERY time you build the project
  # then leave out the 'ALL' keyword from the above command.
endif()
 

# Include petsc files
INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
ADD_DEFINITIONS(${PETSC_DEFINITIONS})

#include femttu files
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/config)

#include b64 files
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/ExternalLibraries/b64-1.4.2/include/)

# Build the FEMTTUlib library
SUBDIRS(src)

# Ask the user which executable have to be built
OPTION(BUILD_STEADYNAVIERSTOKES "Build the Steady NavierStokes application" ON)

# Build Applications
IF(BUILD_STEADYNAVIERSTOKES)
    SUBDIRS(Applications/SteadyNavierStokesParallel)
ENDIF(BUILD_STEADYNAVIERSTOKES)  

OPTION(BUILD_TIMEDEPNAVIERSTOKES "Build the Time dependent NavierStokes application" ON)

IF(BUILD_TIMEDEPNAVIERSTOKES)
    SUBDIRS(Applications/BackwardEulerNavierStokesParallel/)
ENDIF(BUILD_TIMEDEPNAVIERSTOKES)   

OPTION(BUILD_TAYLORVORTEX "Build the Time dependent Taylor vortex application" ON)

IF(BUILD_TAYLORVORTEX)
    SUBDIRS(Applications/TaylorGreenVortexFlowGeneralisedalpha/)
ENDIF(BUILD_TAYLORVORTEX)   

OPTION(BUILD_FSI_CMS_NWM "Build the FSI COMSOL Benchmark" ON)

IF(BUILD_FSI_CMS_NWM)
    SUBDIRS(Applications/FSI_CMS_NWM/)
ENDIF(BUILD_FSI_CMS_NWM)

OPTION(BUILD_FSI_TUREK_3 "Build the FSI Turek 3 Benchmark" ON)

IF(BUILD_FSI_TUREK_3)
    SUBDIRS(Applications/FSI_TUREK_3/)
ENDIF(BUILD_FSI_TUREK_3)

OPTION(BUILD_MHDOPT "Build the FSI Turek 3 Benchmark" ON)

IF(BUILD_MHDOPT)
    SUBDIRS(Applications/mhdopt/)
    LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib/)
ENDIF(BUILD_MHDOPT)