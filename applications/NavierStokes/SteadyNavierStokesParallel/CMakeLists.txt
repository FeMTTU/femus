CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(SteadyNavierStokes)

SET(DIR_NAME "NavierStokes/SteadyNavierStokesParallel")
SET(MAIN_FILE "main")
SET(EXEC_FILE "steadyns")

femusMacroBuildApplication( ${DIR_NAME} ${MAIN_FILE} ${EXEC_FILE})

