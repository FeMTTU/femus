CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(NPdeClassGiacomo)

SET(DIR_NAME "NumPdeSpring2015/Giacomo")
SET(MAIN_FILE "main")
SET(EXEC_FILE "npde_giacomo")

femusMacroBuildApplication(${DIR_NAME} ${MAIN_FILE} ${EXEC_FILE})


