CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(stokes)


SET(MAIN_FILE "main")
SET(EXEC_FILE "steadyst")

femusMacroBuildApplication(${MAIN_FILE} ${EXEC_FILE})

