#!/bin/bash

############################################

MAX_REFS=3
EXECUTABLE_NAME="./fsisteady"

POISSON_RATIO=0.5
############################################

# Direct solver for time info
for ref in $(seq 1 1 ${MAX_REFS});
do
 output_time=3d_direct_time_`date +%Y`-`date +%m`-`date +%d`_`date +%H`-`date +%M`-`date +%S`;  echo $output_time;  mkdir  output/$output_time; 
 ${EXECUTABLE_NAME} -input "turek_FSI1_3d.med" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni ${POISSON_RATIO}  -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_cylinder_beam_3d_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nlevel 1 -nrefinement $ref -std_output run_output.txt -output_time output/$output_time | tee output/$output_time/run_output.txt;
 sleep 2; #to be 100% sure two folders with the same name are not generated
done





# Direct solver for time info
./fsisteady -input "./input/richter3d.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni ${POISSON_RATIO} -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nlevel 1 -nrefinement 1 -std_output richter3d_1_preonly.txt > richter3d_1_preonly.txt
./fsisteady -input "./input/richter3d.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni ${POISSON_RATIO} -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nlevel 1 -nrefinement 2 -std_output richter3d_2_preonly.txt > richter3d_2_preonly.txt
./fsisteady -input "./input/richter3d.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni ${POISSON_RATIO} -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nlevel 1 -nrefinement 3 -std_output richter3d_3_preonly.txt > richter3d_3_preonly.txt

# Direct solver for convergence and MUMPS info
./fsisteady -input "./input/richter3d.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni ${POISSON_RATIO} -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -ksp_view -mat_mumps_icntl_11 1 -ksp_monitor_true_residual -nlevel 1 -nrefinement 1 -std_output richter3d_1_gmres.txt  > richter3d_1_gmres.txt
./fsisteady -input "./input/richter3d.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni ${POISSON_RATIO} -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -ksp_view -mat_mumps_icntl_11 1 -ksp_monitor_true_residual -nlevel 1 -nrefinement 2 -std_output richter3d_2_gmres.txt  > richter3d_2_gmres.txt
./fsisteady -input "./input/richter3d.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni ${POISSON_RATIO} -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -ksp_view -mat_mumps_icntl_11 1 -ksp_monitor_true_residual -nlevel 1 -nrefinement 3 -std_output richter3d_3_gmres.txt  > richter3d_3_gmres.txt

############################################

# Multigrid solver for time info
./fsisteady -input "./input/richter3d.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni ${POISSON_RATIO} -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -nrefinement 3  -nlevel 3  -std_output richter3d_mg.txt  > richter3d_mg.txt

# Multigrid solver for convergence and memory info
./fsisteady -input "./input/richter3d.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni ${POISSON_RATIO} -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -ksp_monitor_true_residual -mem_infos 1 -nrefinement 3  -nlevel 3  -std_output richter3d_mg_meminfo_3.txt  > richter3d_mg_meminfo_3.txt
./fsisteady -input "./input/richter3d.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni ${POISSON_RATIO} -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -ksp_monitor_true_residual -mem_infos 1 -nrefinement 2  -nlevel 2  -std_output richter3d_mg_meminfo_2.txt  > richter3d_mg_meminfo_2.txt
./fsisteady -input "./input/richter3d.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni ${POISSON_RATIO} -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -ksp_monitor_true_residual -mem_infos 1 -nrefinement 1  -nlevel 1  -std_output richter3d_mg_meminfo_1.txt  > richter3d_mg_meminfo_1.txt
