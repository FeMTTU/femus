CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(TimedepNavierStokes)

SET(DIR_NAME "NavierStokes/BackwardEulerNavierStokesParallel")
SET(MAIN_FILE "main")
SET(EXEC_FILE "backeuler")

femusMacroBuildApplication( ${DIR_NAME} ${MAIN_FILE} ${EXEC_FILE})

