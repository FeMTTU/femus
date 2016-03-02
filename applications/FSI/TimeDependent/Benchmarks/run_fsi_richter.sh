#!/bin/bash

# add -ksp_view to view the configurations options of the solver and to view the memory consumption of the direct solver
# add -ksp_monitor_true_residual to plot the residual of the outer gmres
# add -log_summary to view a summary of the run

# Direct solver
# modify -nrefinement to increase or decrease the number of refinement
./fsitimedependent -input "./input/richter3d.neu" -time_step 0.1 -n_timesteps 100 -autosave_time_interval 20 -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_td_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nrefinement 2 -std_output richter3d.txt > richter3d.txt


# GMRES-MG-Richardson-ASM-ILU-MUMPS solver
#./fsitimedependent -input "./input/richter3d.neu" -time_step 0.1 -n_timesteps 100 -autosave_time_interval 20 -nlevel 2 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_td_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 80 -nrefinement 2 -std_output richter3dmg.txt > richter3dmg.txt

