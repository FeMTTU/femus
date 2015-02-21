CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(NavierStokesTemp_AMR)


SET(DIR_NAME "NavierStokes/NavierStokesTemp_AMR")
SET(MAIN_FILE "main")
SET(EXEC_FILE "nst_amr")

femusMacroBuildApplication( ${DIR_NAME} ${MAIN_FILE} ${EXEC_FILE})


