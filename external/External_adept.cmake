#-----------------------------------------------------------------------------
# Adept
#-----------------------------------------------------------------------------

SET(ADEPT_VERSION 1.1)
SET(ADEPT_SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/adept/adept-${ADEPT_VERSION})

#look if there is a package for adept
IF (EXISTS "${ADEPT_SOURCE_DIR}")   
  MESSAGE(STATUS "ADEPT_FOUND = TRUE")
ENDIF (EXISTS "${ADEPT_SOURCE_DIR}")

SET (ADEPT_INCLUDE_DIRS  ${ADEPT_SOURCE_DIR}/include)
SET (ADEPT_LIBRARIES     ${PROJECT_BINARY_DIR}/lib64/libadept.a)

INCLUDE_DIRECTORIES(${ADEPT_INCLUDE_DIRS})

#ExternalProject_Add(${proj}
#  SOURCE_DIR         ${ADEPT_UNPACK_PATH}/adept-${ADEPT_VERSION}
#  CONFIGURE_COMMAND  ""
#  BUILD_COMMAND make libadept.a
#  BUILD_IN_SOURCE 1
#  INSTALL_COMMAND ""
#) 
ADD_SUBDIRECTORY(${ADEPT_SOURCE_DIR})