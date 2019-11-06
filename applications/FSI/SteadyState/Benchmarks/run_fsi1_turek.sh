#!/bin/bash

# output_time=`date +%Y`-`date +%m`-`date +%d`_`date +%H`-`date +%M`-`date +%S`
# echo $output_time
# mkdir  output/$output_time #gives error if the directory exists

# output_time=`date +%Y`-`date +%m`-`date +%d`_`date +%H`-`date +%M`-`date +%S`; echo $output_time; mkdir  output/$output_time; ./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -ksp_view -mat_mumps_icntl_11 1 -ksp_monitor_true_residual -nlevel 1 -nrefinement 1 -std_output run_output.txt -output_time output/$output_time | tee output/$output_time/run_output.txt

# output_time=`date +%Y`-`date +%m`-`date +%d`_`date +%H`-`date +%M`-`date +%S`; echo $output_time; mkdir  output/$output_time; ./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -ksp_view -mat_mumps_icntl_11 1 -ksp_monitor_true_residual -nlevel 1 -nrefinement 1 -std_output FSI1_1_gmres.txt -output_time output/$output_time > output/$output_time/FSI1_1_gmres.txt
############################################

# Direct solver for time info
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nlevel 1 -nrefinement 1 -std_output FSI1_1_preonly.txt > FSI1_1_preonly.txt
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nlevel 1 -nrefinement 2 -std_output FSI1_2_preonly.txt > FSI1_2_preonly.txt
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nlevel 1 -nrefinement 3 -std_output FSI1_3_preonly.txt > FSI1_3_preonly.txt
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nlevel 1 -nrefinement 4 -std_output FSI1_4_preonly.txt > FSI1_4_preonly.txt
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nlevel 1 -nrefinement 5 -std_output FSI1_5_preonly.txt > FSI1_5_preonly.txt

# Direct solver for convergence and MUMPS info
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -ksp_view -mat_mumps_icntl_11 1 -ksp_monitor_true_residual -nlevel 1 -nrefinement 1 -std_output FSI1_1_gmres.txt  > FSI1_1_gmres.txt
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -ksp_view -mat_mumps_icntl_11 1 -ksp_monitor_true_residual -nlevel 1 -nrefinement 2 -std_output FSI1_2_gmres.txt  > FSI1_2_gmres.txt
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -ksp_view -mat_mumps_icntl_11 1 -ksp_monitor_true_residual -nlevel 1 -nrefinement 3 -std_output FSI1_3_gmres.txt  > FSI1_3_gmres.txt
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -ksp_view -mat_mumps_icntl_11 1 -ksp_monitor_true_residual -nlevel 1 -nrefinement 4 -std_output FSI1_4_gmres.txt  > FSI1_4_gmres.txt
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -ksp_view -mat_mumps_icntl_11 1 -ksp_monitor_true_residual -nlevel 1 -nrefinement 5 -std_output FSI1_5_gmres.txt  > FSI1_5_gmres.txt

############################################

# Multigrid solver for time info
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -nrefinement 5  -nlevel 5  -std_output FSI1_mg.txt  > FSI1_mg.txt

# Multigrid solver for convergence and memory info
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -ksp_monitor_true_residual -nrefinement 5  -nlevel 5  -mem_infos 1 -std_output FSI1_mg_mem_5.txt  > FSI1_mg_mem_5.txt
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -ksp_monitor_true_residual -nrefinement 4  -nlevel 4  -mem_infos 1 -std_output FSI1_mg_mem_4.txt  > FSI1_mg_mem_4.txt
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -ksp_monitor_true_residual -nrefinement 3  -nlevel 3  -mem_infos 1 -std_output FSI1_mg_mem_3.txt  > FSI1_mg_mem_3.txt
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -ksp_monitor_true_residual -nrefinement 2  -nlevel 2  -mem_infos 1 -std_output FSI1_mg_mem_2.txt  > FSI1_mg_mem_2.txt
./fsisteady -input "./input/turek_FSI1.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -ksp_monitor_true_residual -nrefinement 1  -nlevel 1  -mem_infos 1 -std_output FSI1_mg_mem_1.txt  > FSI1_mg_mem_1.txt
