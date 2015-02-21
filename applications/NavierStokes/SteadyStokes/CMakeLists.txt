CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(stokes)


SET(DIR_NAME "NavierStokes/SteadyStokes")
SET(MAIN_FILE "main")
SET(EXEC_FILE "steadyst")

femusMacroBuildApplication( ${DIR_NAME} ${MAIN_FILE} ${EXEC_FILE})

