#!/bin/bash

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





###################################################################
echo Setup PETSC
FM_PETSC_DIR_REL=petsc
FM_PETSC_DIR_ABS=$SOFTWARE_DIR/$FM_PETSC_DIR_REL
export PETSC_DIR=$FM_PETSC_DIR_ABS

export PETSC_ARCH=arch-linux2-cxx-opt  #let us only install libmesh with optimized petsc


#######################################################################
echo Download, extract, compile LIBMESH


# What to do to prepare a libmesh installation
# It's assumed that mpi and petsc (and what petsc needs... hdf5?!) are properly installed,
# So we only need to setup the environment for them

# Now, the question is: basically  libmesh needs for me mpi and petsc.

# Now, the alternatives when you install libmesh in the configure script are:
# Tell libmesh directly with --with-thelibraryiwant-dir=,
# or modify PATH and   LD_LIBRARY_PATH...
# The best thing is always to use the options of the configure script...
FM_LIBMESH_DIR_REL=libmesh
FM_LIBMESH_BUILD=build
FM_LIBMESH_INSTALL=install

echo =========== Remove previous installations
rm   -rf $FM_LIBMESH_DIR_REL/
mkdir -p $FM_LIBMESH_DIR_REL/${FM_LIBMESH_BUILD}
echo =========== Download from git repository, my fork in github
git clone https://github.com/giorgiobornia/libmesh.git  ${FM_LIBMESH_DIR_REL}/${FM_LIBMESH_BUILD}
# git clone git://github.com/libMesh/libmesh.git ${FM_LIBMESH_DIR_REL}/${FM_LIBMESH_BUILD}




mkdir -p ${FM_LIBMESH_DIR_REL}/${FM_LIBMESH_INSTALL}-petsc-${PETSC_ARCH}
cd  ${FM_LIBMESH_DIR_REL}/${FM_LIBMESH_BUILD}
# git checkout gambit_format
echo =========== Compile
./configure --with-methods="opt dbg" --prefix=$SOFTWARE_DIR/${FM_LIBMESH_DIR_REL}/${FM_LIBMESH_INSTALL}-petsc-${PETSC_ARCH}
# this script should be clever enough to find the external packages in petsc
make -j $(( $(nproc)-1 ))
make install

