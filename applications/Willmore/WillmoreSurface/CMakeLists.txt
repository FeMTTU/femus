CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

get_filename_component(APP_FOLDER_NAME_PARENT ${CMAKE_CURRENT_LIST_DIR} NAME)

ADD_SUBDIRECTORY(ex1/)
ADD_SUBDIRECTORY(ex2/)
ADD_SUBDIRECTORY(ex3/)

FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/Willmore/WillmoreSurface/input/)
FILE(COPY           ${PROJECT_SOURCE_DIR}/applications/Willmore/WillmoreSurface/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/Willmore/WillmoreSurface/input/)


