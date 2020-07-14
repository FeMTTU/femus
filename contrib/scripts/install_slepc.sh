#!/bin/bash

 
#a release version of slepc cannot be installed with a development version of petsc

SLEPC_VERSION_GIT_TAG=v3.13
 

echo Install slepc

if test "$1" != "--prefix-external"; then
echo "The first argument must be --prefix-external"; exit;
fi 
 
 
echo "The PETSC directory is " $PETSC_DIR
echo "The PETSC architecture is " $PETSC_ARCH

SLEPC_NAME=slepc

#  SOFTWARE_DIR=$HOME/software/
SOFTWARE_DIR=`readlink -f $2`
echo "=========" $SOFTWARE_DIR


cd $SOFTWARE_DIR


echo =========== Remove previous installations
rm -rf $SLEPC_NAME
echo =========== Clone


git clone https://gitlab.com/slepc/slepc


cd $SLEPC_NAME

export SLEPC_DIR=$SOFTWARE_DIR/$SLEPC_NAME

git checkout $SLEPC_VERSION_GIT_TAG -b slepc_current

./configure
make
make check
