CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(TaylorVortex)

SET(DIR_NAME "NavierStokes/TaylorGreenVortexFlowGeneralisedalpha")
SET(MAIN_FILE "main")
SET(EXEC_FILE "taylorvortex")

femusMacroBuildApplication( ${DIR_NAME} ${MAIN_FILE} ${EXEC_FILE})

