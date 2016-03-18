#!/bin/bash

# Direct solver to generate initial data
./fsitimedependent -input "./input/richter3d.neu" -time_step 0.1 -n_timesteps 90 -autosave_time_interval 90 -nlevel 1 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_td_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nrefinement 2 -std_output richter3d.txt > richter3d.txt

###########################################

# Direct solver for time info
./fsitimedependent -input "./input/richter3d.neu" -restart_file_name http://www.math.ttu.edu/~eaulisa/Benchmarks/runs/save/richter3D/richter3D_2_time3.600000 -time_step 0.1 -n_timesteps 10 -autosave_time_interval 20 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_td_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nlevel 1 -nrefinement 2 -std_output richter3D_2_preonly.txt > richter3D_2_preonly.txt
./fsitimedependent -input "./input/richter3d.neu" -restart_file_name http://www.math.ttu.edu/~eaulisa/Benchmarks/runs/save/richter3D/richter3D_3_time3.600000 -time_step 0.1 -n_timesteps 10 -autosave_time_interval 20 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_td_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1 -nlevel 1 -nrefinement 3 -std_output richter3D_3_preonly.txt > richter3D_3_preonly.txt

# Direct solver for convergence and MUMPS info
./fsitimedependent -input "./input/richter3d.neu" -restart_file_name http://www.math.ttu.edu/~eaulisa/Benchmarks/runs/save/richter3D/richter3D_2_time3.600000 -time_step 0.1 -n_timesteps 10 -autosave_time_interval 20 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_td_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -ksp_monitor_true_residual -ksp_view -mat_mumps_icntl_11 1 -nlevel 1 -nrefinement 2 -std_output richter3D_2_gmres.txt > richter3D_2_gmres.txt
./fsitimedependent -input "./input/richter3d.neu" -restart_file_name http://www.math.ttu.edu/~eaulisa/Benchmarks/runs/save/richter3D/richter3D_3_time3.600000 -time_step 0.1 -n_timesteps 10 -autosave_time_interval 20 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_td_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 1 -ksp_monitor_true_residual -ksp_view -mat_mumps_icntl_11 1 -nlevel 1 -nrefinement 3 -std_output richter3D_3_gmres.txt > richter3D_3_gmres.txt

############################################

# Multigrid solver for time info
./fsitimedependent -input "./input/richter3d.neu" -restart_file_name http://www.math.ttu.edu/~eaulisa/Benchmarks/runs/save/richter3D/richter3D_2_time3.600000 -time_step 0.1 -n_timesteps 10 -autosave_time_interval 20 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_td_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -nlevel 2 -nrefinement 2 -std_output richter3D_2_gmres.txt > richter3D_2_gmres.txt
./fsitimedependent -input "./input/richter3d.neu" -restart_file_name http://www.math.ttu.edu/~eaulisa/Benchmarks/runs/save/richter3D/richter3D_3_time3.600000 -time_step 0.1 -n_timesteps 10 -autosave_time_interval 20 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_td_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -nlevel 3 -nrefinement 3 -std_output richter3D_3_gmres.txt > richter3D_3_gmres.txt

# Multigrid solver for convergence and memory info
./fsitimedependent -input "./input/richter3d.neu" -restart_file_name http://www.math.ttu.edu/~eaulisa/Benchmarks/runs/save/richter3D/richter3D_2_time3.600000 -time_step 0.1 -n_timesteps 10 -autosave_time_interval 20 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_td_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -ksp_monitor_true_residual -mem_infos 1 -nlevel 2 -nrefinement 2 -std_output richter3D_2_gmres_meminfo.txt > richter3D_2_gmres_meminfo.txt
./fsitimedependent -input "./input/richter3d.neu" -restart_file_name http://www.math.ttu.edu/~eaulisa/Benchmarks/runs/save/richter3D/richter3D_3_time3.600000 -time_step 0.1 -n_timesteps 10 -autosave_time_interval 20 -rhof 1000 -muf 1 -rhos 1000 -E 1500000 -ni 0.5 -nnonlin_iter 15 -ic_bdc "../../../../lib64/libfsi_td_3d_turek_hron_richter_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 60 -ksp_monitor_true_residual -mem_infos 1 -nlevel 3 -nrefinement 3 -std_output richter3D_3_gmres_meminfo.txt > richter3D_3_gmres_meminfo.txt

