MyFEMuS
======

Welcome to the MyFEMuS project! MyFEMuS is a fork of the FEMuS project administered mainly by Eugenio Aulisa.
The manual intallation, will also work with the FEMuS project.

For the FEMuS project automatic installation see below.


<!-- ![alt tag](https://github.com/FeMTTU/femus/blob/master/doc/images/logo.jpg?raw=true) -->
<!-- ![alt tag](https://github.com/FeMTTU/femus/blob/master/doc/images/FSI.jpg?raw=true) -->

Step by step MyFEMuS manual setup  (largely tested on OpenSuse)
======

Install PETSC

From the directory $INSTALLATION_DIR clone petsc

    git clone -b maint https://bitbucket.org/petsc/petsc petsc
    
    cd petsc
    
Configure, compile, test PETSC with the following options
    
    ./configure --with-debugging=0 --with-x=1 COPTFLAGS="-O3 -march=native -mtune=native" CXXOPTFLAGS="-O3 -march=native -mtune=native" FOPTFLAGS="-O3 -march=native -mtune=native" --download-openmpi=1 --download-fblaslapack=1 --download-hdf5=1 --download-metis=1 --download-parmetis=1 --with-shared-libraries=1 --download-blacs=1 --download-scalapack=1 --download-mumps=1 --download-suitesparse

    make PETSC_DIR=$INSTALLATION_DIR/petsc PETSC_ARCH=arch-linux2-c-opt all

    make PETSC_DIR=$INSTALLATION_DIR/petsc PETSC_ARCH=arch-linux2-c-opt test
 
======

Install SLEPC

From the directory $INSTALLATION_DIR clone slepc

    git clone -b maint https://bitbucket.org/slepc/slepc slepc

    cd slepc
    
Configure, compile, test SLEPC with the following options
    
    export PETSC_DIR=$INSTALLATION_DIR/petsc 
    
    export PETSC_ARCH=arch-linux2-c-opt
    
    ./configure
    
    make SLEPC_DIR=$PWD all
    
    make SLEPC_DIR=$PWD test

======

Install MyFEMuS 

Be sure you have installed al least gcc 7, cmake, cmake-gui. Fparser may be handy for some applications but it is not required. 

Clone the MyFEMuS source code from the github repository

From the directory $INSTALLATION_DIR clone MyFEMuS

    https://github.com/eaulisa/MyFEMuS.git
    
    cd MyFEMuS
    
I generally export the following variables in the ./bashrc file in my user home, so that are available everywhere, otherwise you will need to export them all the times. 
    
    export PETSC_DIR=$INSTALLATION_DIR/petsc 
    
    export PETSC_ARCH=arch-linux2-c-opt
    
    export SLEPC_DIR=$INSTALLATION_DIR/slepc 
    
Configure MyFEMuS using cmake-gui. 

    cmake-gui 

    Where is the source code: $INSTALLATION_DIR/MyFEMuS
    
    Where to build the binaries: $INSTALLATION_DIR/feumsbin
    
    CMAKE_BUILD_TYPE choose between release (default) or debug
    
    Press Configure button
    
    Press Generate button

Compile
    
    cd $INSTALLATION_DIR/femusbin
    
    make
    
Run. All applications are built in the folder $INSTALLATION_DIR/femusbin/applications/..
    
======
    

FEMuS automatic configuration, contact Giorgio Bornia for support.
======

Welcome to the FEMuS project! FEMuS is an open-source Finite Element C++ library 
built on top of PETSc, which allows scientists to build and solve multiphysics 
problems with multigrid and domain decomposition techniques.


<!-- ![alt tag](https://github.com/FeMTTU/femus/blob/master/doc/images/logo.jpg?raw=true) -->
<!-- ![alt tag](https://github.com/FeMTTU/femus/blob/master/doc/images/FSI.jpg?raw=true) -->

Setup
=====


Clone the FEMuS source code from the github repository:


    git clone https://github.com/FeMTTU/femus.git

   
You need PETSc for FEMuS to work.
If PETSc is not already installed in your machine, the script "install_petsc.sh" in contrib/scripts/ will install it automatically,
with the following syntax:
  
    ./femus/contrib/scripts/install_petsc.sh --prefix-external my_dir 
  

where "my_dir" is the directory, either absolute or relative, in which you want PETSc to be installed 
(please put it outside of the femus repo directory, to prevent from potential git tracking).

 Source the "configure_femus.sh" script and execute the function "fm_set_femus" in order to set some environment variables:

    source femus/contrib/scripts/configure_femus.sh

    fm_set_femus  --prefix-external my_dir --method-petsc opt
   

  Create the build directory, cd to it and run cmake:
   
    mkdir femus.build

    cd femus.build

    cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE="[Debug Release RelWithDebInfo MinSizeRel None]"  ../femus



Authors
========

Eugenio Aulisa

Simone Bnà

Giorgio Bornia



License
========

FEMuS is an open-source software distributed under the LGPL license, version 2.1

