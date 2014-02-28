#ifndef __petsc_macro_h__
#define __petsc_macro_h__

#include "FemusExtLib_conf.hpp"


#ifdef FEMUS_HAVE_PETSC 

// A convenient macro for comparing PETSc versions.  Returns 1 if the
// current PETSc version is < major.minor.subminor and zero otherwise.
//
// This macro does not require petscversion.h to be included for it to work correctly.
// It instead relies on the PETSc version numbers detected during configure.  Note that if
// LIBMESH_HAVE_PETSC is not defined, none of the LIBMESH_DETECTED_PETSC_VERSION_* variables will
// be defined either.
#define PETSC_VERSION_LESS_THAN(major,minor,subminor)			                                            \
  ((FEMUS_DETECTED_PETSC_VERSION_MAJOR < (major) ||						                    \
    (FEMUS_DETECTED_PETSC_VERSION_MAJOR == (major) && (FEMUS_DETECTED_PETSC_VERSION_MINOR < (minor) ||	    \
				  (FEMUS_DETECTED_PETSC_VERSION_MINOR == (minor) &&		                    \
				   FEMUS_DETECTED_PETSC_VERSION_SUBMINOR < (subminor))))) ? 1 : 0)



// In case the configure test some day fails, we can fall back on including petscversion.h.
// In order to support PETSc 2.3.1, however, we need to use a few hacks that allow us to
// include petscversion.h without including petsc.h first.  These are explained below...
//
// // We have to jump through some hoops here: in Petsc 2.3.1 you cannot
// // include petscversion.h without including petsc.h beforehand.  This
// // was because petscversion.h relied on the existence of
// // PETSC_EXTERN_CXX_BEGIN/END in 2.3.1.  The problem is, we don't know
// // how to include petsc.h (wrapped or not wrapped in extern "C") until
// // we know the version!
// 
// // First figure out if we need to define PETSC_EXTERN_CXX_BEGIN/END to
// // make petscversion.h happy
// #ifndef PETSC_EXTERN_CXX_BEGIN
// #  define PETSC_EXTERN_CXX_BEGIN
// #  define PETSC_EXTERN_CXX_END
// #  define LIBMESH_PETSC_EXTERN_C_WORKAROUND
// #endif
// 
// // Now actually include it
// #include <petscversion.h>
// 
// // And finally, get rid of our bogus-definitions.  "petsc.h" will set these itself.
// #ifdef LIBMESH_PETSC_EXTERN_C_WORKAROUND
// #  undef PETSC_EXTERN_CXX_BEGIN
// #  undef PETSC_EXTERN_CXX_END
// #endif

// Make up for missing extern "C" in old PETSc versions
#if !defined(LIBMESH_USE_COMPLEX_NUMBERS) && PETSC_VERSION_LESS_THAN(2,3,0)
#  define EXTERN_C_FOR_PETSC_BEGIN extern "C" {
#  define EXTERN_C_FOR_PETSC_END }
#else
#  define EXTERN_C_FOR_PETSC_BEGIN 
#  define EXTERN_C_FOR_PETSC_END
#endif


// Petsc include files
EXTERN_C_FOR_PETSC_BEGIN
#include "petsc.h"
EXTERN_C_FOR_PETSC_END


#endif // LIBMESH_HAVE_PETSC

// #endif

#endif // __petsc_macro_h__
