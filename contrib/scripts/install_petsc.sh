#!/bin/bash

# Here the philosophy is that mpi and hdf5 are downloaded with petsc,
# so we are not going to use the SYSTEM MPI, or the SYSTEM HDF5

# we might install the code in other environments where petsc and mpi are already there,
# so in that case the PETSC_DIR and PETSC_ARCH will be set by the system administrator



# these variables must be in sync with the configure_femus script
PETSC_VERSION_GIT_TAG=v3.20.2
#v3.16.4
# this is the tag in the petsc repo. You can directly find the one you want here: https://gitlab.com/petsc/petsc/-/tags
COMMON_PETSC_DIRNAME=petsc$PETSC_VERSION_GIT_TAG



if test "$1" = "--help" || test "$1" = "-h"; then
echo "The first argument must be --prefix-external";
echo "The second argument must be the directory (either relative or absolute) that will contain the Petsc folder";
exit;
fi

if test "$1" != "--prefix-external"; then
echo "The first argument must be --prefix-external";
exit;
fi

if test "$2" = ""; then
echo "The second argument must be the directory (either relative or absolute) that will contain the Petsc folder";
exit;
fi


echo Install petsc

SOFTWARE_DIR=`readlink -f $2`
echo "=========" $SOFTWARE_DIR

mkdir -p $SOFTWARE_DIR
cd $SOFTWARE_DIR



#######################################################################
echo Download, extract, compile PETSC
FM_PETSC_DIR_REL=$COMMON_PETSC_DIRNAME
FM_PETSC_DIR_ABS=$SOFTWARE_DIR/$FM_PETSC_DIR_REL
export PETSC_DIR=$FM_PETSC_DIR_ABS

myarchs=(arch-linux2-cxx-opt arch-linux2-cxx-dbg)
debugflag=(0 1)


echo =========== Remove previous installations
rm -rf $FM_PETSC_DIR_REL/


echo =========== Clone

git clone -b release https://gitlab.com/petsc/petsc.git $FM_PETSC_DIR_REL


cd $FM_PETSC_DIR_REL

git checkout $PETSC_VERSION_GIT_TAG -b petsc_currently_adopted



for i in 0 #1 let us only install the optimized version, to speed up the installation
do
export PETSC_ARCH=${myarchs[i]}
echo =========== Configure
./configure  --with-debugging=${debugflag[i]} --with-cc=gcc --with-cxx=g++ --with-fc=gfortran COPTFLAGS='-O3 -march=native -mtune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native' --with-shared-libraries=1 --download-openmpi=1 --download-fblaslapack=1 --download-blacs=1 --download-scalapack=1 --download-metis=1 --download-parmetis=1 --download-mumps=1 --download-hdf5=1
#   --with-x=1
#   --with-mpi-dir=$FM_MPI_DIR_ABS 
#   --with-hdf5-dir=$FM_HDF5_DIR_ABS
echo =========== Make
make all #don't put -j here because it gives error! it's automatically -j!
make test
done
