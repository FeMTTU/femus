CMAKE_MINIMUM_REQUIRED(VERSION 3.5)


PROJECT(femus)


#############################################################################################
### Femus library version - BEGIN
#############################################################################################

SET(FEMUS_VERSION_MAJOR 1)
SET(FEMUS_VERSION_MINOR 0)
SET(FEMUS_BUILD_VERSION 0)
SET(FEMUS_VERSION "${FEMTTU_MAJOR_VERSION}.${FEMTTU_MINOR_VERSION}.${FEMTTU_BUILD_VERSION}")

#############################################################################################
### Femus library version - END
#############################################################################################



#############################################################################################
### CMake options - BEGIN
#############################################################################################
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}     -O3 -march=native -mtune=native")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -mtune=native")
set(CMAKE_MACOSX_RPATH 1)
#############################################################################################
### CMake options - END
#############################################################################################




   
#############################################################################################
### Compiling - BEGIN
#############################################################################################

# If the user specifies -DCMAKE_BUILD_TYPE on the command line, take their definition
# and dump it in the cache along with proper documentation, otherwise set CMAKE_BUILD_TYPE
# to Debug prior to calling PROJECT()
#SET(CMAKE_BUILD_TYPE "release")

IF((CMAKE_BUILD_TYPE MATCHES "debug") OR (CMAKE_BUILD_TYPE MATCHES "Debug"))
  SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE})
ELSE()
  SET(CMAKE_BUILD_TYPE "release")
ENDIF()

#message("CMAKE_BUILD_TYPE = " ${CMAKE_BUILD_TYPE})

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-parameter")



set(CMAKE_CXX_STANDARD 11)


#############################################################################################
### Compiling - END
#############################################################################################



#############################################################################################
### Compiling, Output libraries folder - BEGIN
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
### Compiling, Output libraries folder - END
#############################################################################################



### Library extension - BEGIN (@todo where is it used?)
IF(UNIX)
SET (DYLIB so)
ENDIF(UNIX)

IF(APPLE)
SET (DYLIB dylib)
ENDIF(APPLE)
### Library extension - END



#############################################################################################
### Packages, Find - BEGIN
#############################################################################################

# Find Petsc Library - BEGIN
FIND_PACKAGE(PETSc REQUIRED)
MESSAGE(STATUS "PETSC_FOUND = ${PETSC_FOUND}")

SET(HAVE_PETSC 0)
IF(PETSC_FOUND)
  SET(HAVE_PETSC 1)
ENDIF(PETSC_FOUND)
# Find Petsc Library - END


# Find Slepc Library (optional) - BEGIN
FIND_PACKAGE(SLEPc)
MESSAGE(STATUS "SLEPC_FOUND = ${SLEPC_FOUND}")
SET(HAVE_SLEPC 0)
IF(SLEPC_FOUND)
  SET(HAVE_SLEPC 1)
ENDIF(SLEPC_FOUND)
# Find Slepc Library (optional) - END


# Find MPI (optional) - BEGIN
if(NOT DEFINED CMAKE_SUPPRESS_DEVELOPER_WARNINGS)
     set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "No dev warnings")
endif()

FIND_PACKAGE(MPI)
MESSAGE(STATUS "MPI_FOUND = ${MPI_FOUND}")
SET(HAVE_MPI 0)
IF(MPI_FOUND)
SET(HAVE_MPI 1)
ENDIF(MPI_FOUND)
# Find MPI (optional) - END


# Find HDF5 (optional) - BEGIN

# FIND_PACKAGE (HDF5)
# MESSAGE(STATUS "HDF5_FOUND = ${HDF5_FOUND}")
#
# SET(HAVE_HDF5 0)
# IF(HDF5_FOUND)
SET(HAVE_HDF5 1)
# SET(HDF5_FOUND 1)
# SET (HDF5_FOUND 1)
# ENDIF(HDF5_FOUND)

# Find HDF5 (optional) - END



# Find Metis (optional) - BEGIN

# FIND_PACKAGE(METIS)
# MESSAGE(STATUS "METIS_FOUND = ${METIS_FOUND}")
#
# SET (HAVE_METIS 0)
# IF(METIS_FOUND)
SET(HAVE_METIS 1)
# SET(METIS_FOUND 1)
# ENDIF(METIS_FOUND)

# Find Metis (optional) - END



# Find FParser - BEGIN
FIND_PACKAGE(FParser)
MESSAGE(STATUS "FPARSER_FOUND = ${FPARSER_FOUND}")
SET (HAVE_FPARSER 0)
IF(FPARSER_FOUND)
  SET(HAVE_FPARSER 1)
ENDIF(FPARSER_FOUND)
# Find FParser - END






# Find Libmesh (optional) - BEGIN
FIND_PACKAGE(LIBMESH)
MESSAGE(STATUS "LIBMESH_FOUND = ${LIBMESH_FOUND}")
SET (HAVE_LIBMESH 0)
IF(LIBMESH_FOUND)
  SET(HAVE_LIBMESH 1)
ENDIF(LIBMESH_FOUND)
# Find Libmesh (optional) - END



#############################################################################################
### Packages, Find - END
#############################################################################################




#############################################################################################
### External libraries, Add  - BEGIN
#############################################################################################


### External Libraries, Build options - BEGIN

#-- build Adept library
OPTION(BUILD_B64 "Build Adept library for automatic differentiation" ON)

#-- build JsonCPP library
OPTION(BUILD_JSONCPP "Build JsonCPP library for json input parsing" ON)

#-- build b64 library
OPTION(BUILD_ADEPT "Build b64 library for 64-encoding" ON)

### External Libraries, Build options - END 



# Build the b64 library - BEGIN
SET(HAVE_B64 0)
IF(BUILD_B64)
  SET(HAVE_B64 1)
  INCLUDE(external/External_b64.cmake)
ENDIF(BUILD_B64)
# Build the b64 library - END 


# Build the jsoncpp library - BEGIN
SET(HAVE_JSONCPP 0)
IF(BUILD_JSONCPP)
  SET(HAVE_JSONCPP 1)
  INCLUDE(external/External_jsoncpp.cmake)
ENDIF(BUILD_JSONCPP)
# Build the jsoncpp library - END 


# Build the Adept external library - BEGIN
SET(HAVE_ADEPT 0)
IF(BUILD_ADEPT)
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DADEPT_RECORDING_PAUSABLE")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DADEPT_RECORDING_PAUSABLE")
  SET(HAVE_ADEPT 1)
  INCLUDE(external/External_adept.cmake)
ENDIF(BUILD_ADEPT)
# Build the Adept external library - END 



#message(${CMAKE_CXX_FLAGS_RELEASE})
#message(${CMAKE_CXX_FLAGS_DEBUG})

#############################################################################################
### External libraries, Add - END
#############################################################################################


#############################################################################################
### Femus, Include files, Configuration project file - BEGIN
#############################################################################################

# configure a header file to pass some of the CMake settings to the source code
CONFIGURE_FILE(
  "${PROJECT_SOURCE_DIR}/src/00_utils/FemusConfig.hpp.in"
  "${PROJECT_BINARY_DIR}/include/FemusConfig.hpp" @ONLY IMMEDIATE
  )

CONFIGURE_FILE(
  "${PROJECT_SOURCE_DIR}/UseFemus.CMake.in"
  "${PROJECT_BINARY_DIR}/UseFemus.CMake" @ONLY IMMEDIATE
  )

#############################################################################################
### Femus, Include files, Configuration project file - END
#############################################################################################


#############################################################################################
### Femus, Include files - BEGIN
#############################################################################################

#include femus include files - BEGIN
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/02_reference_geom_elements/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/02_reference_geom_elements/00_definition/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/02_reference_geom_elements/01_fe/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/06_solution/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/06_solution/00_single_level/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/07_mesh_or_solution/01_multiple_levels/01_output/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/08_equations/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/08_equations/01_time_dependent/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/08_equations/boundary_conditions/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/08_fsi/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/algebra/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/algebra/01_matrices/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/algebra/02_preconditioners/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/algebra/03_solvers_with_preconditioner/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/ism/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_enums/uq/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_utils/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_utils/00_application_initialization/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_utils/00_parallel/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_utils/data_types/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_file_handling/application_file_handling/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_file_handling/runtime_input_parsing/file/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/00_file_handling/runtime_input_parsing/command_line/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/01_materials/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/01_parameters/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_calculus/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_calculus/function_parser/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_reference_geom_elements/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_reference_geom_elements/00_definition/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_reference_geom_elements/00_definition/1d/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_reference_geom_elements/00_definition/2d/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_reference_geom_elements/00_definition/3d/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_reference_geom_elements/01_fe/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_reference_geom_elements/01_fe/1d/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_reference_geom_elements/01_fe/2d/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_reference_geom_elements/01_fe/3d/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_reference_geom_elements/02_quadrature/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_reference_geom_elements/02_quadrature/0d/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_reference_geom_elements/02_quadrature/1d/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_reference_geom_elements/02_quadrature/2d/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_reference_geom_elements/02_quadrature/3d/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/02_reference_geom_elements/03_fe_evaluations_at_quadrature/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/03_algebra_dense/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/03_algebra_dense/00_vectors/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/03_algebra_dense/01_matrices/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/03_algebra_only_mpi/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/03_algebra_only_mpi/00_vectors/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/03_algebra_only_mpi/01_matrices/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/03_algebra/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/03_algebra/00_vectors/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/03_algebra/01_matrices/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/03_algebra/02_preconditioners/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/05_problem/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_mesh/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_mesh/00_single_level/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_mesh/00_single_level/00_definition/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_mesh/00_single_level/01_input/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_mesh/00_single_level/01_input/01_from_external_file/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_mesh/00_single_level/01_input/02_from_implemented_code/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_mesh/00_single_level/02_partitioning/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_mesh/00_single_level/03_refinement/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_mesh/01_multiple_levels/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_mesh/01_multiple_levels/00_definition/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_mesh/gencase_deprecated/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_solution/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_solution/00_single_level/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_solution/00_single_level/00_definition/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_solution/01_multiple_levels/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_solution/01_multiple_levels/00_definition/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/06_solution/01_multiple_levels/02_FE_convergence/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/07_mesh_or_solution/00_single_level/01_output/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/07_mesh_or_solution/01_multiple_levels/01_output/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/08_algebra_dependent_on_Mesh_and_Solution_but_independent_of_Systems/03_solvers_with_preconditioner/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/08_algebra_dependent_on_Mesh_and_Solution_but_independent_of_Systems/03_solvers_with_preconditioner/petsc_asm/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/08_algebra_dependent_on_Mesh_and_Solution_but_independent_of_Systems/03_solvers_with_preconditioner/petsc_field_split/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/08_equations/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/08_equations/00_stationary/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/08_equations/01_time_dependent/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/08_equations/assemble/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/08_fsi/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/08_fsi/00_stationary/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/08_fsi/01_time_dependent/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/09_optimal_control/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/ism/)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/uq/)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include/)
#include femus include files - END


#############################################################################################
### Femus, Include files - END
#############################################################################################



#############################################################################################
### Packages, Include files - BEGIN
#############################################################################################



# Include petsc files - BEGIN
IF(PETSC_FOUND)
INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
ADD_DEFINITIONS(${PETSC_DEFINITIONS})
ENDIF(PETSC_FOUND)
# Include petsc files - END

# Include slepc files - BEGIN
IF(SLEPC_FOUND)
 INCLUDE_DIRECTORIES(${SLEPC_INCLUDE_DIRS})
ENDIF(SLEPC_FOUND)
# Include slepc files - END


# Include Fparser files - BEGIN
IF(FPARSER_FOUND)
  INCLUDE_DIRECTORIES(${FPARSER_INCLUDE_DIR})
ENDIF(FPARSER_FOUND)
# Include Fparser files - END


# Libmesh includes and flags - BEGIN
IF(LIBMESH_FOUND)
execute_process(COMMAND $ENV{FM_BASEPATH_TO_LM}/$ENV{FM_LM_FOLDER}/bin/libmesh-config --cppflags
                COMMAND tr "\n" " " OUTPUT_VARIABLE LIBMESH_CPP)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBMESH_CPP}")
execute_process(COMMAND $ENV{FM_BASEPATH_TO_LM}/$ENV{FM_LM_FOLDER}/bin/libmesh-config --include
                COMMAND tr "\n" " " OUTPUT_VARIABLE LIBMESH_INC)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBMESH_INC}")
execute_process(COMMAND $ENV{FM_BASEPATH_TO_LM}/$ENV{FM_LM_FOLDER}/bin/libmesh-config --libs
                COMMAND tr "\n" " " OUTPUT_VARIABLE LIBMESH_LIBS)
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LIBMESH_LIBS}")
ENDIF(LIBMESH_FOUND)
# old Make syntax:
# libmesh_INCLUDE  := $(shell ${FM_LM_DIR}/bin/libmesh-config --include)
# libmesh_CPPFLAGS := $(shell ${FM_LM_DIR}/bin/libmesh-config --cppflags)
# libmesh_LIBS     := $(shell ${FM_LM_DIR}/bin/libmesh-config --libs)
# Libmesh includes and flags - END


# HPCC TTU
INCLUDE_DIRECTORIES(/opt/ohpc/pub/mpi/openmpi-gnu/1.10.6/include/)
INCLUDE_DIRECTORIES(/opt/ohpc/pub/libs/gnu/openmpi/boost/1.63.0/include/)


#############################################################################################
### Packages, Include files - END
#############################################################################################





#############################################################################################
### Femus, Add sources - BEGIN
#############################################################################################

# Add the femus library sources
ADD_SUBDIRECTORY(src)

#############################################################################################
### Femus, Add sources - END
#############################################################################################



#############################################################################################
### Femus, Applications, Add  - BEGIN
#############################################################################################
# add femus macro
INCLUDE(${CMAKE_MODULE_PATH}/femusMacroBuildApplication.cmake)


# Add applications
ADD_SUBDIRECTORY(applications)

#############################################################################################
### Femus, Applications, Add - END
#############################################################################################



#############################################################################################
### Femus, Add Unit tests - BEGIN
#############################################################################################

INCLUDE(CTest)
# Add unit tests
IF(BUILD_TESTING)
  ENABLE_TESTING()
  ADD_SUBDIRECTORY(unittests)
ENDIF(BUILD_TESTING)

#############################################################################################
### Femus, Add Unit tests - END
#############################################################################################




#############################################################################################
### Femus, Mesh files - BEGIN
#############################################################################################

FILE(COPY           ${PROJECT_SOURCE_DIR}/src/06_mesh/00_single_level/01_input/00_mesh_files/ DESTINATION ${PROJECT_BINARY_DIR}/src/06_mesh/00_single_level/01_input/00_mesh_files/)

#############################################################################################
### Femus, Mesh files - END
#############################################################################################



#############################################################################################
### Femus, Build documentation (Doxygen) - BEGIN
#############################################################################################

### Build options, Documentation - BEGIN

#-- Add an Option to toggle the generation of the API documentation
OPTION(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)

### Build options, Documentation - END


IF(BUILD_DOCUMENTATION)
  FIND_PACKAGE(Doxygen)
  IF(NOT DOXYGEN_FOUND)
    MESSAGE(FATAL_ERROR
      "Doxygen is needed to build the documentation. Please install it correctly")
  ENDIF(NOT DOXYGEN_FOUND)
  #-- Configure the Template Doxyfile for our specific project
  CONFIGURE_FILE(
    "${PROJECT_SOURCE_DIR}/doc/doxygen/DoxyfileFemuslib.in"
    "${PROJECT_BINARY_DIR}/doc/doxygen/DoxyfileFemuslib" @ONLY IMMEDIATE
  )

  #-- Add a custom target to run Doxygen when ever the project is built
  ADD_CUSTOM_TARGET (Docs ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doc/doxygen/DoxyfileFemuslib
    SOURCES ${PROJECT_BINARY_DIR}/doc/doxygen/DoxyfileFemuslib)
  # IF you do NOT want the documentation to be generated EVERY time you build the project
  # then leave out the 'ALL' keyword from the above command.
ENDIF(BUILD_DOCUMENTATION)

#############################################################################################
### Femus, Build documentation (Doxygen) - END
#############################################################################################



# message(petsc library = ${PETSC_LIBRARIES})
# message(b64 library = ${B64_LIBRARIES})
# message(json library = ${JSONCPP_LIBRARIES})
# message(adept library = ${ADEPT_LIBRARIES})
# message(fparser library = ${FPARSER_LIBRARY})
# message(mpi extra library = ${MPI_EXTRA_LIBRARY})
# message(hdf5 extra library = ${HDF5_LIBRARIES})
# message(slepc library = ${SLEPC_LIBARIES})
# message(slepc library = SLEPC::slepc)
