#!/bin/bash

 
 APP_NAME=$1
 
 
 cp -r  create_application_template/  ../../applications/
 
 cd ../../applications/
 
 mv  create_application_template  $APP_NAME
 
 cd $APP_NAME
 
 mv empty.cpp $APP_NAME.cpp
 
 
 echo "You have to manually add the application to the Cmakelists file with ADD_SUBDIRECTORY($APP_NAME)"
 
  echo "Then run Cmake to generate stuff in the build"
 

