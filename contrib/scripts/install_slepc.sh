#!/bin/bash

 
 #a development version of slepc cannot be installed with a release version of slepc
 

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

git clone https://gitlab.com/slepc/slepc

cd $SLEPC_NAME

export SLEPC_DIR=$SOFTWARE_DIR/$SLEPC_NAME

# git checkout v3.12.2 -b slepc_current

./configure
make
make check
