#-----------------------------------------------------------------------------
# Adpet
#-----------------------------------------------------------------------------

SET(proj Adept)

INCLUDE(ExternalProject)
INCLUDE(${PROJECT_SOURCE_DIR}/ExternalLibraries/modules/PackagesMacro.cmake)

SET(ADEPT_PACKAGE_DIR ${PROJECT_SOURCE_DIR}/ExternalLibraries/Adept/Package)
SET(ADEPT_VERSION 1.0)

#look if there is a package for VTK
IF (EXISTS "${ADEPT_PACKAGE_DIR}")   
  MESSAGE(STATUS "ADEPT: Found ADEPT packages")
  # here we should unpack into "${VTK_BINARY_DIR}/Sources"
  SET (ADEPT_PACKAGE_PATH "${ADEPT_PACKAGE_DIR}")
  SET (ADEPT_UNPACK_PATH "${PROJECT_BINARY_DIR}/ExternalLibraries/Adept/Sources")

  FIND_AND_UNPACK_PACKAGE (adept-${ADEPT_VERSION} ${ADEPT_PACKAGE_PATH} ${ADEPT_UNPACK_PATH})
ENDIF (EXISTS "${ADEPT_PACKAGE_DIR}")

ExternalProject_Add(${proj}
  SOURCE_DIR         ${ADEPT_UNPACK_PATH}/adept-${ADEPT_VERSION}
  CONFIGURE_COMMAND  ""
  BUILD_COMMAND make libadept.a
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND ""
) 

SET (ADEPT_INCLUDE_DIRS  ${ADEPT_UNPACK_PATH}/adept-${ADEPT_VERSION}/include)
SET (ADEPT_LIBRARIES     ${ADEPT_UNPACK_PATH}/adept-${ADEPT_VERSION}/lib/libadept.a)