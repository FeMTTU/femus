#!/bin/bash

# add -ksp_view to view the configurations options of the solver and to view the memory consumption of the direct solver
# add -ksp_monitor_true_residual to plot the residual of the outer gmres
# add -log_summary to view a summary of the run

# Direct solver
# modify -nrefinement to increase or decrease the number of refinement
./fsisteady -input "./input/richter3d.neu"  -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nrefinement 1

# Direct solver with ksp_view and condition number evaluation
# modify -nrefinement to increase or decrease the number of refinement
# notice that the std output should be redirected in the stdOutput.txt file
#./fsisteady -input "./input/richter3d.neu"  -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -ksp_view -mat_mumps_icntl_11 1 -nrefinement 1 -std_output stdOutput.txt > stdOutput.txt



# Direct solver with ksp_view and residual evaluation
# modify -nrefinement to increase or decrease the number of refinement
# ./fsisteady -ksp_view -mat_mumps_icntl_11 1 -input "./input/richter3d.neu" -nrefinement 1 -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1

# Direct solver with ksp_monitor_true_residual
# modify -nrefinement to increase or decrease the number of refinement
# ./fsisteady -ksp_monitor_true_residual -input "./input/richter3d.neu" -nrefinement 1 -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1

# To automatically print the convergence info use "gmres" and specify the std::output redirection
# Notice that the std output should be redirected in the stdOutput.txt file
# ./fsisteady -ksp_monitor_true_residual -input "./input/richter3d.neu" -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_steady_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -nrefinement 1 -std_output stdOutput.txt > stdOutput.txt

# GMRES-MG-Richardson-ASM-ILU-MUMPS solver
#./fsisteady -input "./input/turek.neu" -nrefinement 3 -nlevel 3 -rhof 1000 -muf 1 -rhos 1000 -E 1400000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 40
