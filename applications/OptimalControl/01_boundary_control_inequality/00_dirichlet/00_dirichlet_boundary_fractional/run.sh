#!/bin/bash


EXECUTABLE_NAME="./00_dirichlet_boundary_fractional"
CONF_FILE="../../../../../../femus/applications/OptimalControl/param.hpp"


for level in 2 3 4 5
do
     sed '/#define N_UNIFORM_LEVELS/ c #define N_UNIFORM_LEVELS '${level}' ' -i  ${CONF_FILE}

for myval in 0.1 0.01 0.001 0.0001 0.00001 0.000001  0.0000001 0.00000001
do
   echo $myval
     sed '/#define ALPHA_CTRL_BDRY/ c #define ALPHA_CTRL_BDRY '${myval}' ' -i  ${CONF_FILE}
     make
     $EXECUTABLE_NAME
     sleep 1
done

done
