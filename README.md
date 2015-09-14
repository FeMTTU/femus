FEMuS
======

Welcome to the FEMuS project! FEMuS is an open-source Finite Element C++ library 
built on top of Petsc which allows scientists to build and solve multiphysics 
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
  

where "my_dir" is the directory in which PETSc will be installed (please put it outside of the femus repo directory).

 Source the "configure_femus.sh" script and execute the function "fm_set_femus" in order to set some environment variables:

    source femus/contrib/scripts/configure_femus.sh

    fm_set_femus  --prefix-external my_dir --method-petsc [opt,dbg]
   

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

Femus is an open-source software distributed under the LGPL license, version 2.1

