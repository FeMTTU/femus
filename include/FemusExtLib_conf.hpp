#ifndef _femusextlibconf_h_
#define _femusextlibconf_h_


//this fils concerns the use of external libraries 
//and it should just be automatically generated 
//by the configure, based on the libraries we choose to activate




//One thing is: What libraries are AVAILABLE, and we link to ALL of them
//Another thing is: WHAT LIBRARIES e currently use,
//and we choose this at RUNTIME.
//We might have TWO I/O libraries, but now we want to use ONLY ONE
//In other cases the libraries are EXCLUSIVE, we can use ONLY ONE of them.
//like the solver packages.



//**********************************************************
//************ LIBMESH *************************************
//**********************************************************
//This include file rules the Libmesh dependencies
// We must tune these variables with $(fm-enable-libmesh)
// if we define these variables but $(fm-enable-libmesh) is NO, you get errors


// #define LM_INIT   //TODO useless now  *********** LIB-TYPE dependency

// #define LM_REFCOUNT  //uncoupled         *********** LIB-TYPE dependency

//#define LM_REAL      //uncoupled          *********** INCLUDE-TYPE dependency


//*****************
#ifdef LM_REFCOUNT

//refcount needs init (and also debug mode, by the way)
    #ifndef LM_INIT   
      #define LM_INIT
    #endif

//also,if refcount starts then LM_REAL must also start, otherwise
                      //you get ambiguous references to Real
    #ifndef LM_REAL 
    #define LM_REAL
    #endif
    
  #include "libmesh/reference_counted_object.h"

#endif
 //******************


#ifdef LM_INIT
  #include "libmesh/libmesh.h"
#endif
 
//you may want to use LM_INIT either directly or because you use LM_REFCOUNT




//*******************************************************
//*********** PETSC****************************
//*******************************************************

// ======== SOLVER LIBRARY ======
//   #define FEMUS_HAVE_PETSC 
//   #define HAVE_MPI   //what if I want to use Petsc without MPI?
 
#ifdef FEMUS_HAVE_PETSC
 #define LSOLVER  PETSC_SOLVERS
#endif



//******PETSC VERSION ***
// #ifdef FEMUS_HAVE_PETSC

//what happens here is done with the configure_femus.sh script!

//was: #include "libmesh_config.h"
//Instead of libmesh_config we put only the variables 
// assuming that the rest of the code picks the LIBMESH variables
//DIRECTLY from somewhere else


/* PETSc's major version number, as detected by LibMesh */
// #ifndef FEMUS_DETECTED_PETSC_VERSION_MAJOR 
// #define FEMUS_DETECTED_PETSC_VERSION_MAJOR  3 
// #endif

/* PETSc's minor version number, as detected by LibMesh */
// #ifndef FEMUS_DETECTED_PETSC_VERSION_MINOR 
// #define FEMUS_DETECTED_PETSC_VERSION_MINOR  1 
// #endif

/* PETSc's subminor version number, as detected by LibMesh */
// #ifndef FEMUS_DETECTED_PETSC_VERSION_SUBMINOR 
// #define FEMUS_DETECTED_PETSC_VERSION_SUBMINOR  0 
// #endif

// #endif


//*******************************************************
//************ MPI *************************************
//*******************************************************
// this should stay here, not in the separate files!
#ifdef HAVE_MPI// Mpi
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "mpi.h"   //For MPI_COMM_WORLD
#pragma GCC diagnostic warning "-Wunused-parameter"
#endif






//*******************************************************
//************ HDF5 *************************************
//*******************************************************
/* Here you define the FORMAT 
 * for FILE INPUT/OUTPUT
 * But, this is used by no one!
 * This is another sort of LIBRARY configuration:
 * we choose the library for the I/O format.
 * Well, we could also choose more than one i/o format.
 
 */
// ======= OUTPUT FORMAT ======
// #define Vtk 1
#define FEMUS_PRINT_XDMF //not used yet
#define FEMUS_HAVE_HDF5  //not used yet

#endif