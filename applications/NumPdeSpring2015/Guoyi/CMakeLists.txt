CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(NPdeClassGuoyi)

SET(DIR_NAME "NumPdeSpring2015/Guoyi")
SET(MAIN_FILE "main")
SET(EXEC_FILE "npde_guoyi")

femusMacroBuildApplication(${DIR_NAME} ${MAIN_FILE} ${EXEC_FILE})
