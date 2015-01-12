#-----------------------------------------------------------------------------
# jsoncpp
#-----------------------------------------------------------------------------

SET(JSONCPP_VERSION 0.5.0)
SET(JSONCPP_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ExternalLibraries/jsoncpp/jsoncpp-src-${JSONCPP_VERSION})

#look if there is a package for jsoncpp
IF (EXISTS "${JSONCPP_SOURCE_DIR}")   
  MESSAGE (STATUS "JSONCPP_FOUND = TRUE")
ENDIF (EXISTS "${JSONCPP_SOURCE_DIR}")

SET (JSONCPP_INCLUDE_DIRS  ${JSONCPP_SOURCE_DIR}/include)
SET (JSONCPP_LIBRARIES     ${PROJECT_BINARY_DIR}/lib64/libjsoncpp.so)

INCLUDE_DIRECTORIES(${JSONCPP_INCLUDE_DIRS})

#ExternalProject_Add(${proj}
#  SOURCE_DIR         ${JSONCPP_UNPACK_PATH}/jsoncpp-${JSONCPP_VERSION}
#  CONFIGURE_COMMAND  ""
#  BUILD_COMMAND make libjsoncpp.so
#  BUILD_IN_SOURCE 1
#  INSTALL_COMMAND ""
#) 
ADD_SUBDIRECTORY(${JSONCPP_SOURCE_DIR})