CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(Poisson)

SET(MAIN_FILE "main")
SET(EXEC_FILE "Poisson")

femusMacroBuildApplication(Poisson ${MAIN_FILE} ${EXEC_FILE})

