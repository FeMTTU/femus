#!/bin/bash

nprocs=6
dt=(0.001 0.005 0.01)
alphaVel=(600 700 750 800 900)
udes=1.4
Bref=5.e-3

# for t in $(seq 0 1 2); do
  for a in $(seq 0 1 4); do
      mpiexec -n ${nprocs} -mca btl ^sm $PWD/applications/main/main-${FEMUS_METHOD} --dt 1. --alphaVel ${alphaVel[a]}  --udes ${udes} --Bref ${Bref}
  done
# done