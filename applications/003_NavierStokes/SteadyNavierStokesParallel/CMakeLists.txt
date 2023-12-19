CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

PROJECT(SteadyNavierStokes)

SET(MAIN_FILE "main")
SET(EXEC_FILE "steadyns")

femusMacroBuildApplication(${MAIN_FILE} ${EXEC_FILE})
