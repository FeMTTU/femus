#
#  femusMacroBuildApplication.cmake
#
#  Created by Simone Bnà 
#
#


MACRO(femusMacroBuildApplication mainname appname)

# Build the executable
ADD_EXECUTABLE(${appname} ${PROJECT_SOURCE_DIR}/${mainname}.cpp)
set_target_properties(${appname} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")


# Link the executable to the petsc anf femttu libs
TARGET_LINK_LIBRARIES(${appname} femus)
TARGET_LINK_LIBRARIES(${appname} ${PETSC_LIBRARIES})
TARGET_LINK_LIBRARIES(${appname} ${B64_LIBRARIES})
TARGET_LINK_LIBRARIES(${appname} ${JSONCPP_LIBRARIES})
TARGET_LINK_LIBRARIES(${appname} ${ADEPT_LIBRARIES})

IF(SLEPC_FOUND)
  #TARGET_LINK_LIBRARIES(${appname} ${SLEPC_LIBARIES})
  TARGET_LINK_LIBRARIES(${appname} SLEPC::slepc SLEPC::slepc_static)
ENDIF(SLEPC_FOUND)

IF(FPARSER_FOUND)
  TARGET_LINK_LIBRARIES(${appname} ${FPARSER_LIBRARY})
ENDIF(FPARSER_FOUND)

IF(MPI_FOUND) 
  TARGET_LINK_LIBRARIES(${appname} ${MPI_CXX_LIBRARIES})
  #TARGET_LINK_LIBRARIES(${appname} ${MPI_EXTRA_LIBRARY})
ENDIF(MPI_FOUND)

IF(HDF5_FOUND)
  TARGET_LINK_LIBRARIES(${appname} ${HDF5_LIBRARIES})
ENDIF(HDF5_FOUND)

FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output/)
FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/input/)
FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/save/)
FILE(COPY           ${PROJECT_SOURCE_DIR}/input/ DESTINATION ${PROJECT_BINARY_DIR}/input/)
# TODO this file copy does not generate a dependency rule in the makefiles, maybe we should think of how to obtain that,
# to avoid re-running cmake when new input files are added in the applications

ENDMACRO(femusMacroBuildApplication)
