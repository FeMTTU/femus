#ifndef _femusextlibconf_h_
#define _femusextlibconf_h_


#include "FEMTTUConfig.h"

//ok, the next step is to SUBSTITUTE THIS FILE WITH the FEMTTU one and REMOVE the -D macros 
// in the COMPILING COMMAND... this can still happen in the OLD MAKEFILE system
// L'idea potrebbe essere quella di INCLUDERE qui il ttuconfig intanto


//this fils concerns the use of external libraries 
//and it should just be automatically generated 
//by the configure, based on the libraries we choose to activate




//One thing is: What libraries are AVAILABLE, and we link to ALL of them
//Another thing is: WHAT LIBRARIES e currently use,
//and we choose this at RUNTIME.
//We might have TWO I/O libraries, but now we want to use ONLY ONE
//In other cases the libraries are EXCLUSIVE, we can use ONLY ONE of them.
//like the solver packages.


//*******************************************************
//*********** PETSC****************************
//*******************************************************



//******PETSC VERSION ***

//what happens here is done with the configure_femus.sh script!

//was: #include "libmesh_config.h"
//Instead of libmesh_config we put only the variables 
// assuming that the rest of the code picks the LIBMESH variables
//DIRECTLY from somewhere else


#endif