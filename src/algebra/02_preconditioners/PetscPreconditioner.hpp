/*=========================================================================

 Program: FEMuS
 Module: PetscPreconditioner
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_algebra_PetscPreconditioner_hpp__
#define __femus_algebra_PetscPreconditioner_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "FemusConfig.hpp"

#ifdef HAVE_PETSC

// Local includes
#include "Preconditioner.hpp"
#include "PrecondtypeEnum.hpp"

// Petsc includes
#include "petscpc.h"


namespace femus {


//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class SparseMatrix;
class NumericVector;

/**
 * This class provides an interface to  preconditioners  from Petsc.
 *
 */

class PetscPreconditioner : public Preconditioner {

public:

    /**  Constructor. Initializes PetscPreconditioner data structures */
    PetscPreconditioner ():  Preconditioner() {}

    /** Destructor. */
    virtual ~PetscPreconditioner () {
        this->clear ();
    }

    /** Release all memory and clear data structures. */
    virtual void clear () {}

    /** Initialize data structures if not done so already. */
    virtual void init ();

    // Return
    /** Returns the actual Petsc PC struct. */
    PC pc() {
        return _pc;
    }

    // Compute
    /** Computes the preconditioned vector "y" based on input "x". */
    /** Usually by solving Py=x to get the action of P^-1 x. */
    virtual void apply(const NumericVector & x, NumericVector & y);

    /** Tells PETSC to use the user-specified preconditioner */
    static void set_petsc_preconditioner_type  (const PreconditionerType & preconditioner_type, PC & pc, const int &parallelOverlapping = 0);

private:

    /**
     * Some PETSc preconditioners (ILU, LU) don't work in parallel.  This function
     * is called from set_petsc_preconditioner_type() to set additional options
     * for those so-called sub-preconditioners.  This method ends up being static
     * so that it can be called from set_petsc_preconditioner_type().  Not sure
     * why set_petsc_preconditioner_type() needs to be static though...
     */
    static void set_petsc_subpreconditioner_type(const PCType type, PC& pc);

protected:

    PC _pc;   ///<Preconditioner context
    Mat _mat; ///< Petsc Matrix

};

} //end namespace femus


#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // #ifdef __petsc_linear_solver_h__
