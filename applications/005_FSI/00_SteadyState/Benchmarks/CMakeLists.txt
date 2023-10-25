CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(FSI_STEADY)

# add bdc directory
ADD_SUBDIRECTORY("./bdc")

SET(MAIN_FILE "fsiSteady")
SET(EXEC_FILE "fsisteady")

femusMacroBuildApplication(${MAIN_FILE} ${EXEC_FILE})

FILE(COPY ./run_fsi_2d.sh DESTINATION ${PROJECT_BINARY_DIR})
FILE(COPY ./run_fsi_3d.sh DESTINATION ${PROJECT_BINARY_DIR})
FILE(COPY ./run_from_kdevelop_with_redirect_or_tee.sh DESTINATION ${PROJECT_BINARY_DIR})


