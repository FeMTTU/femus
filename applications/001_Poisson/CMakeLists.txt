CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(Poisson)

SET(MAIN_FILE "main")
SET(EXEC_FILE "Poisson")

femusMacroBuildApplication(${MAIN_FILE} ${EXEC_FILE})

