#!/bin/bash

# add -ksp_view to view the configurations options of the solver and to view the memory consumption of the direct solver
# add -ksp_monitor_true_residual to plot the residual of the outer gmres
# add -log_summary to view a summary of the run

# Direct solver
# modify -nrefinement to increase or decrease the number of refinement
./fsisteady -input "./input/turek.neu" -dim 2 -nrefinement 5 -nlevel 1 -sim 1 -rhof 1000 -muf 1 -rhos 1000 -E 1400000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1


# GMRES-MG-Richardson-ASM-ILU-MUMPS solver
#./fsisteady -input "./input/turek.neu" -dim 2 -nrefinement 3 -nlevel 3 -sim 1 -rhof 1000 -muf 1 -rhos 1000 -E 1400000 -ni 0.5 -ic_bdc "../../../../lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 40
