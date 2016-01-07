#!/bin/bash
# Direct solver
#./fsisteady -ksp_view -input "./input/turek.neu" -dim 2 -nrefinement 4 -nlevel 1 -sim 1 -rhof 1000 -muf 1 -rhos 1000 -E 1400000 -ni 0.5 -ic_bdc "/home/simone/Software/Devel/femusopt/lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "preonly" -max_outer_solver_iter 1

# MG solver
./fsisteady -input "./input/turek.neu" -dim 2 -nrefinement 3 -nlevel 3 -sim 1 -rhof 1000 -muf 1 -rhos 1000 -E 1400000 -ni 0.5 -ic_bdc "/home/simone/Software/Devel/femusopt/lib64/libfsi_steady_2d_turek_hron_benchmark_bdc.so" -outer_ksp_solver "gmres" -max_outer_solver_iter 40
