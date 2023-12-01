CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

get_filename_component(APP_FOLDER_NAME_PARENT ${CMAKE_CURRENT_LIST_DIR} NAME)

ADD_SUBDIRECTORY(ex1/) #no eqn  (move to tutorial)
ADD_SUBDIRECTORY(ex2/) #Boussinesq  (assemble from library now)
ADD_SUBDIRECTORY(ex3/) #NS (assemble from library now)
ADD_SUBDIRECTORY(ex4/) #Poisson
ADD_SUBDIRECTORY(ex5/) #Poisson
ADD_SUBDIRECTORY(ex6/) #Poisson
ADD_SUBDIRECTORY(ex7/) #Poisson
ADD_SUBDIRECTORY(ex8/) #Poisson

ADD_SUBDIRECTORY(ex9/) #Poisson
