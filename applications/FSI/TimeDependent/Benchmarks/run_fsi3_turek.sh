#!/bin/bash

# add -ksp_view to view the configurations options of the solver and to view the memory consumption of the direct solver
# add -ksp_monitor_true_residual to plot the residual of the outer gmres
# add -log_summary to view a summary of the run



# Direct solver
# modify -nrefinement to increase or decrease the number of refinements

# for FSI3: store the solution after 120 iterations
./fsitimedependent -input "./input/turek_FSI3.neu" -n_timesteps 120 -autosave_time_interval 120 -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 6000000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi3_td_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nrefinement 4


# Load solution files:
# Modify the file name turek_FSI3_level_time6.420000 according to the number of -nrefinements: level = nrefinements.

# Load solution files directly from the server http://www.math.ttu.edu/~eaulisa/Benchmarks/FSI3/save/
# For time info
#./fsitimedependent -input "./input/turek_FSI3.neu" -restart_file_name http://www.math.ttu.edu/~eaulisa/Benchmarks/runs/save/FSI3/turek_FSI3_4_time6.420000 -n_timesteps 20 -autosave_time_interval 40 -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 6000000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi3_td_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nrefinement 4 -std_output stdOutput.txt > stdOutput.txt
# For MUMPS and convergence info
#./fsitimedependent -input "./input/turek_FSI3.neu" -restart_file_name http://www.math.ttu.edu/~eaulisa/Benchmarks/runs/save/FSI3/turek_FSI3_4_time6.420000 -n_timesteps 20 -autosave_time_interval 40 -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 6000000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi3_td_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -ksp_monitor_true_residual -ksp_view -mat_mumps_icntl_11 1 -nrefinement 4 -std_output stdOutput.txt > stdOutput.txt


# Load solution files from the local folder ./save/
# Be sure you already copied into the local folder ./save/ the solution files from http://www.math.ttu.edu/~eaulisa/Benchmarks/FSI3/save/
# For time info
#./fsitimedependent -input "./input/turek_FSI3.neu" -restart_file_name ./save/turek_FSI3_4_time6.420000 -n_timesteps 20 -autosave_time_interval 40 -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 6000000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi3_td_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nrefinement 4 -std_output stdOutput.txt > stdOutput.txt
# For MUMPS and convergence info
#./fsitimedependent -input "./input/turek_FSI3.neu" -restart_file_name ./save/turek_FSI3_4_time6.420000 -n_timesteps 20 -autosave_time_interval 40 -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 6000000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi3_td_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -ksp_monitor_true_residual -ksp_view -mat_mumps_icntl_11 1 -nrefinement 4 -std_output stdOutput.txt > stdOutput.txt


# GMRES-MG-Richardson-ASM-ILU-MUMPS solver
#./fsitimedependent -input "./input/turek.neu" -n_timesteps 500 -autosave_time_interval 50 -rhof 1000 -muf 1 -rhos 1000 -E 6000000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi3_td_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -ksp_monitor_true_residual -nrefinement 2  -nlevel 2  -std_output turek_mg.txt > turek_mg.txt
