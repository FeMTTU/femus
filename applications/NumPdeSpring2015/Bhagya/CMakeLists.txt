CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(NPdeClassBhagya)

SET(DIR_NAME "NumPdeSpring2015/Bhagya")
SET(MAIN_FILE "main")
SET(EXEC_FILE "npde_bhagya")

femusMacroBuildApplication(${DIR_NAME} ${MAIN_FILE} ${EXEC_FILE})
