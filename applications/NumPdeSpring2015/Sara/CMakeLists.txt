CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(NPdeClassSara)

SET(DIR_NAME "NumPdeSpring2015/Sara")
SET(MAIN_FILE "main")
SET(EXEC_FILE "npde_sara")

femusMacroBuildApplication(${DIR_NAME} ${MAIN_FILE} ${EXEC_FILE})