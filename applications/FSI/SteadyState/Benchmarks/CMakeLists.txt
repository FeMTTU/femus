CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(FSI_STEADY)

# add bdc directory
ADD_SUBDIRECTORY("./bdc")

SET(MAIN_FILE "fsiSteady")
SET(EXEC_FILE "fsisteady")

femusMacroBuildApplication(${MAIN_FILE} ${EXEC_FILE})

FILE(COPY ./run_fsi_turek.sh DESTINATION ${PROJECT_BINARY_DIR})
FILE(COPY ./run_fsi_richter.sh DESTINATION ${PROJECT_BINARY_DIR})


