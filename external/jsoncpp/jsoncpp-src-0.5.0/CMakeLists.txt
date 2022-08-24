CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(jsoncpp)
 
SET(jsoncpplib_src
src/lib_json/json_reader.cpp 
src/lib_json/json_value.cpp 
src/lib_json/json_writer.cpp
) 

IF (NOT LIBRARY_OUTPUT_PATH)
  SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib64 CACHE PATH
       "Single output directory for building all libraries.")
  MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
ENDIF (NOT LIBRARY_OUTPUT_PATH)
 
ADD_LIBRARY(${PROJECT_NAME} SHARED ${jsoncpplib_src})