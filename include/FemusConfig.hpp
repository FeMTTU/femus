#ifndef __femus_FemusConfig_hpp__
#define __femus_FemusConfig_hpp__


// the configured options and settings for Tutorial
#define FEMTTU_VERSION_MAJOR 
#define FEMTTU_VERSION_MINOR 

//MPI library

#define HAVE_MPI

//PETSc solver library

#define HAVE_PETSC

//MPI library

#define HAVE_METIS

//HDF5 library

#define HAVE_HDF5

//b64 library

#define HAVE_B64

//jsoncpp library

#define HAVE_JSONCPP

//adept library

#define HAVE_ADEPT

//FParser library

/* #undef HAVE_FPARSER */

//libMesh library

/* #undef HAVE_LIBMESH */

//SLEPc library

/* #undef HAVE_SLEPC */

#ifdef HAVE_PETSC
  #undef  LSOLVER
  #define LSOLVER  PETSC_SOLVERS
  
  //******PETSC VERSION ***88
  /* PETSc's major version number, as detected by findPETSC.cmake */
  #ifndef FEMTTU_DETECTED_PETSC_VERSION_MAJOR
    #define FEMTTU_DETECTED_PETSC_VERSION_MAJOR  
  #endif

  /* PETSc's minor version number, as detected by findPETSC.cmake */
  #ifndef FEMTTU_DETECTED_PETSC_VERSION_MINOR
    #define FEMTTU_DETECTED_PETSC_VERSION_MINOR  
  #endif

  /* PETSc's subminor version number, as detected by findPETSC.cmake */
  #ifndef FEMTTU_DETECTED_PETSC_VERSION_SUBMINOR
    #define FEMTTU_DETECTED_PETSC_VERSION_SUBMINOR  
  #endif
  
  

#endif



#endif
