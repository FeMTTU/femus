#-----------------------------------------------------------------------------
# b64
#-----------------------------------------------------------------------------

SET(B64_VERSION 1.4.2)
SET(B64_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ExternalLibraries/b64/b64-${B64_VERSION})

#look if there is a package for b64
IF (EXISTS "${B64_SOURCE_DIR}")   
  MESSAGE (STATUS "B64_FOUND = TRUE")
ENDIF (EXISTS "${B64_SOURCE_DIR}")

SET (B64_INCLUDE_DIRS  ${B64_SOURCE_DIR}/include/)
SET (B64_LIBRARIES     ${PROJECT_BINARY_DIR}/lib64/libb64.a)

INCLUDE_DIRECTORIES(${B64_INCLUDE_DIRS})

#ExternalProject_Add(${proj}
#  SOURCE_DIR         ${B64_UNPACK_PATH}/b64-${B64_VERSION}
#  CONFIGURE_COMMAND  ""
#  BUILD_COMMAND make libb64.a
#  BUILD_IN_SOURCE 1
#  INSTALL_COMMAND ""
#) 
ADD_SUBDIRECTORY(${B64_SOURCE_DIR})