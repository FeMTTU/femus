CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(femttu)

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")

#############################################################################################
### Femus library version
#############################################################################################

SET(FEMTTU_VERSION_MAJOR 1)
SET(FEMTTU_VERSION_MINOR 0)
SET(FEMTTU_BUILD_VERSION 0)
SET(FEMTTU_VERSION "${FEMTTU_MAJOR_VERSION}.${FEMTTU_MINOR_VERSION}.${FEMTTU_BUILD_VERSION}")



#############################################################################################
### Compilers
#############################################################################################

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



#############################################################################################
### Output libraries folder
#############################################################################################

# where to place output binaries
SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib64/ CACHE INTERNAL
     "Single output directory for building all libraries.")
# SET (EXECUTABLE_OUTPUT_PATH ${MAF_BINARY_DIR}/bin/ CACHE INTERNAL 
#      "Single output directory for building all executables.")

MARK_AS_ADVANCED (
  LIBRARY_OUTPUT_PATH 
  #EXECUTABLE_OUTPUT_PATH
)



#############################################################################################
### Find packages
#############################################################################################

# Find Petsc Library 
FIND_PACKAGE(PETSc REQUIRED)

FIND_PACKAGE(FParser REQUIRED)
MESSAGE(STATUS "FPARSER_FOUND = ${FPARSER_FOUND}")
SET (HAVE_FPARSER 0)
IF(FPARSER_FOUND)
  SET(HAVE_FPARSER 1)
ENDIF(FPARSER_FOUND)

SET(HAVE_PETSC 0)
IF(PETSC_FOUND)
  SET(HAVE_PETSC 1)
ENDIF(PETSC_FOUND) 

# Find MPI (optional) 
FIND_PACKAGE(MPI)

SET(HAVE_MPI 0)
IF(MPI_FOUND)
  SET(HAVE_MPI 1)
ENDIF(MPI_FOUND) 

# Find HDF5 (optional) 
#FIND_PACKAGE (HDF5 NAMES hdf5)

#SET(HAVE_HDF5 0)
#IF(HDF5_FOUND)
  SET(HAVE_HDF5 1)
#ENDIF(HDF5_FOUND)

# Find Libmesh (optional)
FIND_PACKAGE(LIBMESH)
MESSAGE(STATUS "LIBMESH_FOUND = ${LIBMESH_FOUND}")
SET (HAVE_LIBMESH 0)
IF(LIBMESH_FOUND)
  SET(HAVE_LIBMESH 1)
ENDIF(LIBMESH_FOUND)



#############################################################################################
### Configuration project file
#############################################################################################

# configure a header file to pass some of the CMake settings
# to the source code
CONFIGURE_FILE(
  "${PROJECT_SOURCE_DIR}/config/FEMTTUConfig.h.in"
  "${PROJECT_BINARY_DIR}/config/FEMTTUConfig.h" @ONLY IMMEDIATE
  )

  
  
#############################################################################################
### Build options
#############################################################################################

#-- build Adept package
OPTION(BUILD_ADEPT "Build Adept for automatic differentiation library" ON)

#-- Add an Option to toggle the generation of the API documentation
OPTION(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)

# Ask the user which executable have to be built
OPTION(BUILD_NAVIERSTOKESAPPS "Build the NavierStokes-line applications" ON)

# Ask the user which executable have to be built
OPTION(BUILD_FSIAPPS "Build the FSI applications" ON)

# --
OPTION(BUILD_GENCASE "Build the case generator" ON)

# --
OPTION(BUILD_MHDOPT "Build the MHD optimal control problem" ON)

# --
OPTION(BUILD_TEMPER "Build the temperature optimal control problem" ON)

# --
OPTION(BUILD_POISSON "Build the Poisson problem" ON)

# --
OPTION(BUILD_TEMPERATURE "Build the Temperature problem" ON)

# --
OPTION(BUILD_POISSON_AMR "Build the Poisson_AMR problem" ON)



#############################################################################################
### Build documentation (Doxygen)
#############################################################################################

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



#############################################################################################
### Include files
#############################################################################################

# Include petsc files
INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
ADD_DEFINITIONS(${PETSC_DEFINITIONS})

#include femus include files
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/algebra)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/mesh)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/utils)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/quadrature)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/parallel)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/equations)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/solution)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/enums)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/fe)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/physics)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/config)

#include b64 include files
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/ExternalLibraries/b64-1.4.2/include/)

#include jsoncpp include files
INCLUDE_DIRECTORIES(ExternalLibraries/jsoncpp-src-0.5.0/include/)

# Include Fparser files
IF(FPARSER_FOUND)
  INCLUDE_DIRECTORIES(${FPARSER_INCLUDE_DIR})
ENDIF(FPARSER_FOUND)



#############################################################################################
### Build femus and external libraries 
#############################################################################################

# Build the b64 library
ADD_SUBDIRECTORY(ExternalLibraries/b64-1.4.2)

# Build the jsoncpp library
ADD_SUBDIRECTORY(ExternalLibraries/jsoncpp-src-0.5.0)

IF(BUILD_ADEPT)
  INCLUDE(${PROJECT_SOURCE_DIR}/ExternalLibraries/Adept/External_Adept.cmake)
ENDIF(BUILD_ADEPT)

# Build the femus library
ADD_SUBDIRECTORY(src)

#############################################################################################
### Build Applications 
#############################################################################################

IF(BUILD_NAVIERSTOKESAPPS)
  
  # Steady NavierStokes
  ADD_SUBDIRECTORY(applications/NavierStokesApps/SteadyNavierStokesParallel)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/SteadyNavierStokesParallel/output/)
  FILE(MAKE_DIRECTORY applications/NavierStokesApps/SteadyNavierStokesParallel/input/)
  FILE(COPY           applications/NavierStokesApps/SteadyNavierStokesParallel/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/SteadyNavierStokesParallel/input/)

  # Time-dependent BE NavierStokes
  ADD_SUBDIRECTORY(applications/NavierStokesApps/BackwardEulerNavierStokesParallel/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/BackwardEulerNavierStokesParallel/output/)
  FILE(MAKE_DIRECTORY applications/NavierStokesApps/BackwardEulerNavierStokesParallel/input/)
  FILE(COPY           applications/NavierStokesApps/BackwardEulerNavierStokesParallel/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/BackwardEulerNavierStokesParallel/input/)
 
  # Time-dependent Generalized-alpha NavierStokes
  ADD_SUBDIRECTORY(applications/NavierStokesApps/TaylorGreenVortexFlowGeneralisedalpha/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/TaylorGreenVortexFlowGeneralisedalpha/output/)
  FILE(MAKE_DIRECTORY applications/NavierStokesApps/TaylorGreenVortexFlowGeneralisedalpha/input/)
  FILE(COPY           applications/NavierStokesApps/TaylorGreenVortexFlowGeneralisedalpha/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/TaylorGreenVortexFlowGeneralisedalpha/input/)
 
  # Steady Stokes
  ADD_SUBDIRECTORY(applications/NavierStokesApps/SteadyStokes/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/SteadyStokes/output/)
  FILE(MAKE_DIRECTORY applications/NavierStokesApps/SteadyStokes/input/)
  FILE(COPY           applications/NavierStokesApps/SteadyStokes/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/SteadyStokes/input/)

  # Steady Navier Stokes steady
  ADD_SUBDIRECTORY(applications/NavierStokesApps/NavierStokesTemp_AMR/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/NavierStokesTemp_AMR/output/)
  FILE(MAKE_DIRECTORY applications/NavierStokesApps/NavierStokesTemp_AMR/input/)
  FILE(COPY           applications/NavierStokesApps/NavierStokesTemp_AMR/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/NavierStokesTemp_AMR/input/)
 
ENDIF(BUILD_NAVIERSTOKESAPPS)


IF(BUILD_FSIAPPS)

# FSI SteadyState Benchmark Comsol
    ADD_SUBDIRECTORY(applications/FSIApps/FSI_CMS_NWM/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/FSIApps/FSI_CMS_NWM/output/)
    FILE(MAKE_DIRECTORY applications/FSIApps/FSI_CMS_NWM/input/)
    FILE(COPY           applications/FSIApps/FSI_CMS_NWM/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/FSIApps/FSI_CMS_NWM/input/)

# FSI SteadyState Benchmark Turek1
    ADD_SUBDIRECTORY(applications/FSIApps/FSI_TUREK_1/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/FSIApps/FSI_TUREK_1/output/)
    FILE(MAKE_DIRECTORY applications/FSIApps/FSI_TUREK_1/input/)
    FILE(COPY           applications/FSIApps/FSI_TUREK_1/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/FSIApps/FSI_TUREK_1/input/)

# FSI SteadyState Benchmark Turek3
    ADD_SUBDIRECTORY(applications/FSIApps/FSI_TUREK_3/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/FSIApps/FSI_TUREK_3/output/)
    FILE(MAKE_DIRECTORY applications/FSIApps/FSI_TUREK_3/input/)
    FILE(COPY           applications/FSIApps/FSI_TUREK_3/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/FSIApps/FSI_TUREK_3/input/)

# FSI SteadyState in porous media
    ADD_SUBDIRECTORY(applications/FSIApps/fsi_porous/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/FSIApps/fsi_porous/output/)
    FILE(MAKE_DIRECTORY applications/FSIApps/fsi_porous/input/)
    FILE(COPY           applications/FSIApps/fsi_porous/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/FSIApps/fsi_porous/input/)

ENDIF(BUILD_FSIAPPS)


IF(BUILD_GENCASE)
    ADD_SUBDIRECTORY(applications/gencase/)
ENDIF(BUILD_GENCASE)


IF(BUILD_MHDOPT)
       ADD_SUBDIRECTORY(applications/mhdopt/)
    FILE(MAKE_DIRECTORY applications/mhdopt/input/)
    FILE(COPY applications/mhdopt/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/mhdopt/input/)
ENDIF(BUILD_MHDOPT)


IF(BUILD_TEMPER)
       ADD_SUBDIRECTORY(applications/temper/)
    FILE(MAKE_DIRECTORY applications/temper/input/)
    FILE(COPY applications/temper/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/temper/input/)
ENDIF(BUILD_TEMPER)


IF(BUILD_POISSON)
       ADD_SUBDIRECTORY(applications/Poisson/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/Poisson/output/)
    FILE(MAKE_DIRECTORY applications/Poisson/input/)
    FILE(COPY           applications/Poisson/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/Poisson/input/)
ENDIF(BUILD_POISSON)


IF(BUILD_TEMPERATURE)
       ADD_SUBDIRECTORY(applications/temperature/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/temperature/output/)
    FILE(MAKE_DIRECTORY applications/temperature/input/)
    FILE(COPY           applications/temperature/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/temperature/input/)
ENDIF(BUILD_TEMPERATURE)


IF(BUILD_POISSON_AMR)
       ADD_SUBDIRECTORY(applications/Poisson_AMR/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/Poisson_AMR/output/)
    FILE(MAKE_DIRECTORY applications/Poisson_AMR/input/)
    FILE(COPY           applications/Poisson_AMR/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/Poisson_AMR/input/)
ENDIF(BUILD_POISSON_AMR)



#############################################################################################
### Unit tests
#############################################################################################

INCLUDE(CTest)
ENABLE_TESTING()

ADD_SUBDIRECTORY(unittests/testNSSteadyDD/)
    FILE(MAKE_DIRECTORY unittests/testNSSteadyDD/input/)
    FILE(COPY           unittests/testNSSteadyDD/input/ DESTINATION  ${PROJECT_BINARY_DIR}/unittests/testNSSteadyDD/input/)
    
ADD_SUBDIRECTORY(unittests/testFSISteady/)
    FILE(MAKE_DIRECTORY unittests/testFSISteady/input/)
    FILE(COPY           unittests/testFSISteady/input/ DESTINATION  ${PROJECT_BINARY_DIR}/unittests/testFSISteady/input/)

    
    