#!/bin/bash

# cmake -L # to view all cache variables

MAX_REFS=5 \
EXECUTABLE_NAME="./fsisteady"

POISSON_RATIO=0.5
############################################

# Direct solver for time info
for ref in $(seq 1 1 ${MAX_REFS});
do
 output_time=2d_direct_time_`date +%Y`-`date +%m`-`date +%d`_`date +%H`-`date +%M`-`date +%S`;  echo $output_time;  mkdir  output/$output_time; \
 ${EXECUTABLE_NAME} -input "turek_FSI1.med" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni ${POISSON_RATIO} -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nlevel 1 -nrefinement $ref -std_output run_output.txt -output_time output/$output_time |& tee output/$output_time/run_output.txt; \
 sleep 2;  \
done
#to be 100% sure two folders with the same name are not generated

# Direct solver for convergence and MUMPS info
for ref in $(seq 1 1 ${MAX_REFS});
do
 output_time=2d_direct_conv_`date +%Y`-`date +%m`-`date +%d`_`date +%H`-`date +%M`-`date +%S`;  echo $output_time;  mkdir  output/$output_time; 
 ${EXECUTABLE_NAME} -input "turek_FSI1.med" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni ${POISSON_RATIO} -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -ksp_view -mat_mumps_icntl_11 1 -ksp_monitor_true_residual -nlevel 1 -nrefinement $ref -std_output run_output.txt -output_time output/$output_time |& tee output/$output_time/run_output.txt;
 sleep 2; #to be 100% sure two folders with the same name are not generated
done


############################################

# Multigrid solver for time info
${EXECUTABLE_NAME} -input "turek_FSI1.med" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni ${POISSON_RATIO} -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -nrefinement 5  -nlevel 5  -std_output FSI1_mg.txt  > FSI1_mg.txt

# Multigrid solver for convergence and memory info
for ref in $(seq 1 1 ${MAX_REFS});
do
 output_time=2d_mg_`date +%Y`-`date +%m`-`date +%d`_`date +%H`-`date +%M`-`date +%S`;  echo $output_time;  mkdir  output/$output_time;  
 ${EXECUTABLE_NAME} -input "turek_FSI1.med" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni ${POISSON_RATIO} -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -ksp_monitor_true_residual -nrefinement $ref  -nlevel $ref  -mem_infos 1  -std_output run_output.txt -output_time output/$output_time |& tee output/$output_time/run_output.txt;
 sleep 2; #to be 100% sure two folders with the same name are not generated
done
