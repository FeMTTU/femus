CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(NPdeClassThanuja)

SET(DIR_NAME "NumPdeSpring2015/Thanuja")
SET(MAIN_FILE "main")
SET(EXEC_FILE "npde_thanuja")

femusMacroBuildApplication(NumPdeSpring2015/Thanuja ${MAIN_FILE} ${EXEC_FILE})