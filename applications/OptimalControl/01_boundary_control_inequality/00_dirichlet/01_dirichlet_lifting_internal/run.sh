#!/bin/bash


EXECUTABLE_NAME="./01_dirichlet_lifting_internal"
CONF_FILE="../../../../../../femus/applications/OptimalControl/param.hpp"



for myval in 0.1 0.01 0.001 0.0001 0.00001 0.000001  0.0000001 0.00000001
do
   echo $myval
     sed '/#define ALPHA_CTRL_VOL/ c #define ALPHA_CTRL_VOL '${myval}' ' -i  ${CONF_FILE}
     make
     $EXECUTABLE_NAME
     sleep 1
done
