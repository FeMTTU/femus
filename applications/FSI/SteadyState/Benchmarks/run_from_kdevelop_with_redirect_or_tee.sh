#!/bin/bash

# This script can be loaded in Kdevelop in Run -> Configure Launches
# It allows to use redirect ">" or "| tee" when doing Execute from Kdevelop

EXECUTABLE_NAME="./fsisteady"

 ${EXECUTABLE_NAME} -input "turek_FSI1.med" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nlevel 1 -nrefinement 1 -std_output run_output.txt -output_time output/kdev | tee output/kdev/run_output.txt
############################################
