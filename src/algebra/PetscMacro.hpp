#ifndef __femus_algebra_PetscMacro_hpp__
#define __femus_algebra_PetscMacro_hpp__

// #include "SolverlibConf.hpp"
#include "FemusConfig.hpp"

#ifdef HAVE_PETSC

// A convenient macro for comparing PETSc versions.  Returns 1 if the
// current PETSc version is < major.minor.subminor and zero otherwise.
//
// This macro does not require petscversion.h to be included for it to work correctly.
// It instead relies on the PETSc version numbers detected during configure.

// #define PETSC_VERSION_LESS_THAN(major,minor,subminor)			                                            \
//   ((FEMTTU_DETECTED_PETSC_VERSION_MAJOR < (major) ||						                    \
//     (FEMTTU_DETECTED_PETSC_VERSION_MAJOR == (major) && (FEMTTU_DETECTED_PETSC_VERSION_MINOR < (minor) ||	    \
// 				  (FEMTTU_DETECTED_PETSC_VERSION_MINOR == (minor) &&		                    \
// 				   FEMTTU_DETECTED_PETSC_VERSION_SUBMINOR < (subminor))))) ? 1 : 0)


#define PETSC_VERSION_LESS_THAN(major,minor,subminor) 0

// Make up for missing extern "C" in old PETSc versions
#if PETSC_VERSION_LESS_THAN(2,3,0)
  #  define EXTERN_C_FOR_PETSC_BEGIN extern "C" {
  #  define EXTERN_C_FOR_PETSC_END }
 #else
   #  define EXTERN_C_FOR_PETSC_BEGIN
   #  define EXTERN_C_FOR_PETSC_END
#endif


// Petsc include files
EXTERN_C_FOR_PETSC_BEGIN
#include <petsc.h>
EXTERN_C_FOR_PETSC_END

#endif // HAVE_PETSC


#endif // __petsc_macro_h__
