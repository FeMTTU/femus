#!/bin/bash

# add -ksp_view to view the configurations options of the solver and to view the memory consumption of the direct solver
# add -ksp_monitor_true_residual to plot the residual of the outer gmres
# add -log_summary to view a summary of the run

# Direct solver
# modify -nrefinement to increase or decrease the number of refinement
./fsisteady -input "./input/turek.neu"  -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nrefinement 1 -std_output stdOutput.txt > stdOutput.txt

# Direct solver with ksp_view and condition number evaluation
# modify -nrefinement to increase or decrease the number of refinement
# notice that the std output should be redirected in the stdOutput.txt file
# ./fsisteady -input "./input/turek.neu"  -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -ksp_view -mat_mumps_icntl_11 1 -nrefinement 1 -std_output stdOutput.txt  > stdOutput.txt

# Direct solver with ksp_view and residual evaluation
# modify -nrefinement to increase or decrease the number of refinement
# ./fsisteady -ksp_view -mat_mumps_icntl_11 1 -input "./input/turek.neu" -nrefinement 6 -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nrefinement 1 -std_output stdOutput.txt > stdOutput.txt


# Direct solver with with ksp_monitor_true_residual
# modify -nrefinement to increase or decrease the number of refinement
# ./fsisteady -ksp_monitor_true_residual -input "./input/turek.neu" -nrefinement 6 -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1

# To automatically print the convergence info use "gmres" and specify the std::output redirection
# notice that the std output should be redirected in the stdOutput.txt file
# ./fsisteady -ksp_monitor_true_residual -input "./input/turek.neu" -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -nrefinement 1 -std_output stdOutput.txt > stdOutput.txt

# GMRES-MG-Richardson-ASM-ILU-MUMPS solver
./fsisteady -input "./input/turek.neu" -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -ksp_monitor_true_residual -nrefinement 2  -nlevel 2  -std_output turek_mg.txt > turek_mg.txt
