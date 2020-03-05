#############################################################
# Try to find SLEPc                                         #
#                                                           #
# Once done this will define:                               #
#  SLEPC_FOUND     - system has SLEPc                       #
#  SLEPC_DIR       - SLEPc directory                        #
#  SLEPC_INC       - SLEPc include directory                #
#  SLEPC_LIB       - SLEPc library (static or dynamic)      #
#                                                           #
# Usage:                                                    #
#  find_package(SLEPc)                                      #
#                                                           #
# Setting these changes the behavior of the search          #
#  SLEPC_DIR       - SLEPc directory                        #
#############################################################

######## MADE BY RAMAH SHARAF FOR GMSH #############

## Try to set SLEPC_DIR ##
##########################
if(NOT DEFINED SLEPC_DIR)
  set(SLEPC_DIR $ENV{SLEPC_DIR})
endif()

## Includes ##
##############
if(EXISTS "${SLEPC_DIR}/include" AND
   EXISTS "${SLEPC_DIR}/${PETSC_ARCH}/include")
 set(SLEPC_INC "${SLEPC_DIR}/include" "${SLEPC_DIR}/${PETSC_ARCH}/include")
else()
  message(SEND_ERROR "SLEPc includes not found")
endif()

## Library ##
#############
if(EXISTS "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.so")
  set(SLEPC_LIB "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.so")
elseif(EXISTS "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.a")
  set(SLEPC_LIB "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.a")
elseif(EXISTS "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.dylib")
  set(SLEPC_LIB "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.dylib")
else()
  message(SEND_ERROR "SLEPc library not found")
endif()

####### PUT BACK THE OLD FILE ###########

# Configure SLEPc IMPORT (this involves creating an 'imported' target
# and attaching 'properties')
if (SLEPC_FOUND AND NOT TARGET SLEPC::slepc)
  add_library(SLEPC::slepc INTERFACE IMPORTED)

  # Add include paths
  set_property(TARGET SLEPC::slepc PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${SLEPC_INCLUDE_DIRS})

  # Add libraries
  unset(_libs)
  foreach (lib ${SLEPC_LIBRARIES})
    find_library(LIB_${lib} NAMES ${lib} PATHS ${SLEPC_LIBRARY_DIRS} NO_DEFAULT_PATH)
    list(APPEND _libs ${LIB_${lib}})
  endforeach()
  set_property(TARGET SLEPC::slepc PROPERTY INTERFACE_LINK_LIBRARIES "${_libs}")
endif()

if (SLEPC_FOUND AND NOT TARGET SLEPC::slepc_static)
  add_library(SLEPC::slepc_static INTERFACE IMPORTED)

  # Add libraries (static)
  unset(_libs)
  foreach (lib ${SLEPC_STATIC_LIBRARIES})
    find_library(LIB_${lib} ${lib} HINTS ${SLEPC_STATIC_LIBRARY_DIRS})
    list(APPEND _libs ${LIB_${lib}})
  endforeach()
  set_property(TARGET SLEPC::slepc_static PROPERTY INTERFACE_LINK_LIBRARIES "${_libs}")

endif()

# Compile and run test
if (DOLFIN_SKIP_BUILD_TESTS)

  # Assume SLEPc works
  set(SLEPC_TEST_RUNS TRUE)

elseif (SLEPC_FOUND)

  # Create SLEPc test program
  set(SLEPC_TEST_LIB_CPP
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/slepc_test_lib.cpp")
  file(WRITE ${SLEPC_TEST_LIB_CPP} "
#include \"petsc.h\"
#include \"slepceps.h\"
int main()
{
  PetscErrorCode ierr;
  int argc = 0;
  char** argv = NULL;
  ierr = SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
  EPS eps;
  ierr = EPSCreate(PETSC_COMM_SELF, &eps); CHKERRQ(ierr);
  //ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);
  ierr = EPSDestroy(&eps); CHKERRQ(ierr);
  ierr = SlepcFinalize(); CHKERRQ(ierr);
  return 0;
}
")

  # Add MPI variables if MPI has been found
  if (MPI_C_FOUND)
    set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${MPI_C_INCLUDE_PATH})
    set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${MPI_C_LIBRARIES})
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${MPI_C_COMPILE_FLAGS}")
  endif()

  # Try to run test program (shared linking)
  try_run(
    SLEPC_TEST_LIB_EXITCODE
    SLEPC_TEST_LIB_COMPILED
    ${CMAKE_CURRENT_BINARY_DIR}
    ${SLEPC_TEST_LIB_CPP}
    CMAKE_FLAGS
    "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
    "-DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES}"
    #LINK_LIBRARIES PETSC::petsc SLEPC::slepc
    LINK_LIBRARIES PETSC_LIBRARIES
    LINK_LIBRARIES SLEPC::slepc
    COMPILE_OUTPUT_VARIABLE SLEPC_TEST_LIB_COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE SLEPC_TEST_LIB_OUTPUT
    )

  if (SLEPC_TEST_LIB_COMPILED AND SLEPC_TEST_LIB_EXITCODE EQUAL 0)

    message(STATUS "Test SLEPC_TEST_RUNS with shared library linking - Success")
    set(SLEPC_TEST_RUNS TRUE)

    # Static libraries not required, so unset
    set_property(TARGET SLEPC::slepc_static PROPERTY INTERFACE_LINK_LIBRARIES)

  else()

    message(STATUS "Test SLEPC_TEST_RUNS with shared library linking - Failed")

    # Add MPI variables if MPI has been found
    if (MPI_C_FOUND)
      set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${MPI_C_INCLUDE_PATH})
      set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${MPI_C_LIBRARIES})
      set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${MPI_C_COMPILE_FLAGS}")
    endif()

    # Try to run test program (static linking)
    try_run(
      SLEPC_TEST_STATIC_LIBS_EXITCODE
      SLEPC_TEST_STATIC_LIBS_COMPILED
      ${CMAKE_CURRENT_BINARY_DIR}
      ${SLEPC_TEST_LIB_CPP}
      CMAKE_FLAGS
      "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
      "-DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES}"
      LINK_LIBRARIES SLEPC::slepc SLEPC::slepc_static
      COMPILE_OUTPUT_VARIABLE SLEPC_TEST_STATIC_LIBS_COMPILE_OUTPUT
      RUN_OUTPUT_VARIABLE SLEPC_TEST_STATIC_LIBS_OUTPUT
      )

    if (SLEPC_TEST_STATIC_LIBS_COMPILED AND SLEPC_TEST_STATIC_LIBS_EXITCODE EQUAL 0)
      message(STATUS "Test SLEPC_TEST__RUNS with static linking - Success")
      set(SLEPC_TEST_RUNS TRUE)
    else()
      message(STATUS "Test SLEPC_TETS_RUNS with static linking - Failed")
      set(SLEPC_TEST_RUNS FALSE)
    endif()
  endif()
endif()

# Standard package handling
include(FindPackageHandleStandardArgs)
if (SLEPC_FOUND)
  find_package_handle_standard_args(SLEPc
    REQUIRED_VARS SLEPC_FOUND SLEPC_TEST_RUNS
    VERSION_VAR SLEPC_VERSION
    FAIL_MESSAGE "SLEPc could not be configured.")
else()
  find_package_handle_standard_args(SLEPc
    REQUIRED_VARS SLEPC_FOUND
    FAIL_MESSAGE "SLEPc could not be found. Be sure to set SLEPC_DIR.")
endif()

######## END OLD FILE #########

## CMake check and done ##
##########################
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SLEPc
  "SLEPc could not be found: be sure to set SLEPC_DIR in your environment variables"
  SLEPC_LIB SLEPC_INC SLEPC_DIR)
