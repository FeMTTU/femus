CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(Poisson)

SET(DIR_NAME "Poisson")
SET(MAIN_FILE "main")
SET(EXEC_FILE "Poisson")

femusMacroBuildApplication(${DIR_NAME} ${MAIN_FILE} ${EXEC_FILE})

