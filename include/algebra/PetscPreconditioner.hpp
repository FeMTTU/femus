/*=========================================================================

 Program: FEMUS
 Module: PetscPreconditioner
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __petsc_preconditioner_h__
#define __petsc_preconditioner_h__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "FEMTTUConfig.h"

#ifdef HAVE_PETSC

// Local includes
#include "Preconditioner.hpp"
#include "PrecondtypeEnum.hpp"

// Petsc includes
#include "petscpc.h"

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class SparseMatrix;
class NumericVector;

/// This class provides an interface to  preconditioners  from Petsc.
class PetscPreconditioner : public Preconditioner {

protected:
  PC _pc; ///<Preconditioner context
  Mat _mat;///< Petsc Matrix

public:
  // Constructor
  ///  Constructor. Initializes PetscPreconditioner data structures
  PetscPreconditioner ():  Preconditioner() {}
  /// Destructor.
  virtual ~PetscPreconditioner () {
    this->clear ();
  }
  /// Release all memory and clear data structures.
  virtual void clear () {}
  /// Initialize data structures if not done so already.
  virtual void init ();
  // Return
  /// Returns the actual Petsc PC struct.
  PC pc() {
    return _pc;
  }

  // Compute
  /// Computes the preconditioned vector "y" based on input "x".
  /// Usually by solving Py=x to get the action of P^-1 x.
  virtual void apply(const NumericVector & x, NumericVector & y);
  /// Tells PETSC to use the user-specified preconditioner
  static void set_petsc_preconditioner_type
  (const PreconditionerType & preconditioner_type, PC & pc);
private:
  /**
   * Some PETSc preconditioners (ILU, LU) don't work in parallel.  This function
   * is called from set_petsc_preconditioner_type() to set additional options
   * for those so-called sub-preconditioners.  This method ends up being static
   * so that it can be called from set_petsc_preconditioner_type().  Not sure
   * why set_petsc_preconditioner_type() needs to be static though...
   */
// #if PETSC_VERSION_LESS_THAN(3,0,0)
// //  In Petsc 2.3.3, PCType was #define'd as const char*
//   static void set_petsc_subpreconditioner_type(PCType type, PC& pc);
// #else
//  In later versions, PCType is #define'd as char*, so we need the const
  static void set_petsc_subpreconditioner_type(const PCType type, PC& pc);
// #endif

};




/*----------------------- inline functions ----------------------------------*/
/*
inline PetscPreconditioner::PetscPreconditioner () :  Preconditioner(){}
inline PetscPreconditioner::~PetscPreconditioner ()*/

#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // #ifdef __petsc_linear_solver_h__
