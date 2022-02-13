#!/bin/bash


EXECUTABLE_NAME="./00_dirichlet_boundary_fractional"
CONF_FILE="../../../../../../femus/applications/OptimalControl/param.hpp"
MAIN_FILE="../../../../../../femus/applications/OptimalControl/01_boundary_control_inequality/00_dirichlet/00_dirichlet_boundary_fractional/00_dirichlet_boundary_fractional.cpp"



for bdry_norm_flag in 0  1
do

     sed '/#define IS_CTRL_FRACTIONAL_SOBOLEV/ c #define IS_CTRL_FRACTIONAL_SOBOLEV '${bdry_norm_flag}' ' -i  ${MAIN_FILE}

for level in 1 2 3 4 5
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


done
