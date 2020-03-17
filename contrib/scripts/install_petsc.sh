#!/bin/bash

# Here the philosophy is that mpi and hdf5 are downloaded with petsc,
# so we are not going to use the SYSTEM MPI, or the SYSTEM HDF5

# we might install the code in other environments where petsc and mpi are already there,
# so in that case the PETSC_DIR and PETSC_ARCH will be set by the system administrator


echo Install petsc

if test "$1" != "--prefix-external"; then
echo "The first argument must be --prefix-external"; exit;
fi

if test "$2" = ""; then
echo "The second argument must be the directory (either relative or absolute) where you want to install petsc"; exit;
fi


SOFTWARE_DIR=`readlink -f $2`
echo "=========" $SOFTWARE_DIR

mkdir -p $SOFTWARE_DIR
cd $SOFTWARE_DIR



#######################################################################
echo Download, extract, compile PETSC
FM_PETSC_DIR_REL=petsc
FM_PETSC_DIR_ABS=$SOFTWARE_DIR/$FM_PETSC_DIR_REL
export PETSC_DIR=$FM_PETSC_DIR_ABS

myarchs=(arch-linux2-cxx-opt arch-linux2-cxx-dbg)
debugflag=(0 1)


echo =========== Remove previous installations
rm -rf $FM_PETSC_DIR_REL/
echo =========== Clone

git clone https://gitlab.com/petsc/petsc $FM_PETSC_DIR_REL

git checkout v3.12.2 -b petsc_current


cd $FM_PETSC_DIR_REL

for i in 0 #1 let us only install the optimized version, to speed up the installation
do
export PETSC_ARCH=${myarchs[i]}
echo =========== Configure
./configure  --with-debugging=${debugflag[i]} --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --with-clanguage=cxx  COPTFLAGS='-O3 -march=native -mtune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native' --with-shared-libraries=1 --download-openmpi=1 --download-fblaslapack=1 --download-blacs=1 --download-scalapack=1 --download-metis=1 --download-parmetis=1 --download-mumps=1 --download-hdf5=1
#   --with-x=1
#   --with-mpi-dir=$FM_MPI_DIR_ABS 
#   --with-hdf5-dir=$FM_HDF5_DIR_ABS
echo =========== Make
make all #don't put -j here because it gives error! it's automatically -j!
make test
done
