CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(NavierStokesTemp_AMR)


SET(MAIN_FILE "main")
SET(EXEC_FILE "nst_amr")

femusMacroBuildApplication(${MAIN_FILE} ${EXEC_FILE})


