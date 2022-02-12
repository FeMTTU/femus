PROJECT(adept)

SET(adeptlib_src
adept/adept.cpp adept/adept_openmp.cpp 
) 

# SET(adeptlib_src
# adept/adept_openmp.cpp 
# ) 

IF (NOT LIBRARY_OUTPUT_PATH)
  SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib64 CACHE PATH
       "Single output directory for building all libraries.")
  MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
ENDIF (NOT LIBRARY_OUTPUT_PATH)
 
ADD_LIBRARY(${PROJECT_NAME} STATIC ${adeptlib_src})