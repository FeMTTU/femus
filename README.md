FEMuS
======

Welcome to the FEMuS project! FEMuS is an opensource Finite Element C++ library 
built on top of Petsc which allows scientists to build and solve multiphysics 
problems with multigrid and domain decomposition techniques.


<!-- ![alt tag](https://github.com/FeMTTU/femus/blob/master/doc/images/logo.jpg?raw=true) -->
<!-- ![alt tag](https://github.com/FeMTTU/femus/blob/master/doc/images/FSI.jpg?raw=true) -->

Setup
=====


Clone the FEMuS source code from the github repository:


    git clone https://github.com/FeMTTU/femus.git

   
Install petsc and libmesh in some common directory "my_external_directory" (please put it outside of the femus repo directory). 
If they are not installed already, the script "install_external.sh" in contrib/scripts/ will do it automatically, with the following syntax:

  
    ./femus/contrib/scripts/install_petsc_and_libmesh.sh --prefix-external my_external_directory 

  
Source the "configure_femus.sh" script and execute the function "fm_set_femus" in order to set some environment variables:


    source femus/contrib/scripts/configure_femus.sh

    fm_set_femus  --prefix-external my_external_directory --method-petsc [opt,dbg] --method-libmesh [opt,dbg]

   
Create the build directory, cd to it and run cmake:
   
    mkdir femus.build

    cd femus.build

    cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE="[Debug Release RelWithDebInfo MinSizeRel None]" -DLIBMESH_DIR=my_libmesh_directory   ../femus



Authors
========

Eugenio Aulisa

Simone Bnà      

Giorgio Bornia



License
========

Femus is an open-source software distributed under the LGPL license, version 2.1

