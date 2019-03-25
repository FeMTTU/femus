#!/bin/bash

 
 #a development version of slepc cannot be installed with a release version of slepc
 
echo "The PETSC directory is " $PETSC_DIR
echo "The PETSC architecture is " $PETSC_ARCH

SLEPC_NAME=slepc

REF_DIR=$HOME/software/

cd $REF_DIR

git clone https://bitbucket.org/slepc/slepc

cd $SLEPC_NAME

export SLEPC_DIR=$REF_DIR/$SLEPC_NAME

git checkout v3.10.2 -b slepc_current

./configure
make
make check
