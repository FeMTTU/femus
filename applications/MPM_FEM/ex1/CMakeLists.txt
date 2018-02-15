CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

get_filename_component(APP_FOLDER_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

set(THIS_APPLICATION ${APP_FOLDER_NAME})

PROJECT("${APP_FOLDER_NAME_PARENT}_${THIS_APPLICATION}")


SET(MAIN_FILE "${THIS_APPLICATION}") # the name of the main file with no extension
SET(EXEC_FILE "${APP_FOLDER_NAME_PARENT}_${MAIN_FILE}") # the name of the executable file

femusMacroBuildApplication(${MAIN_FILE} ${EXEC_FILE})
