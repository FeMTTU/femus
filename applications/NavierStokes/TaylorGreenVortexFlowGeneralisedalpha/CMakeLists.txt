CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(TaylorVortex)

SET(MAIN_FILE "main")
SET(EXEC_FILE "taylorvortex")

femusMacroBuildApplication(${MAIN_FILE} ${EXEC_FILE})

