/*=========================================================================

 Program: FEMUS
 Module: Preconditioner
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_algebra_Preconditioner_hpp__
#define __femus_algebra_Preconditioner_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "FemusConfig.hpp"
#include "SolverPackageEnum.hpp"
#include "PrecondtypeEnum.hpp"

// C++ includes
#include <memory>
#include <iostream>


namespace femus {



//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class SparseMatrix;
class NumericVector;

/**
 * This class provides a uniform interface for preconditioners.  This base
 * class is overloaded to provide linear solvers from different packages
 * like PETSC or Trilinos.
 *
 * In the below comments P is the matrix to be preconditioned with Apply()
 * performing the equivalent of the matrix vector product P^-1 x.  This
 * can also be thought of as (usually approximately) solving for Py=x.
 */

class Preconditioner {

public:

    /** Constructor. Initializes Preconditioner data structures */
    Preconditioner ();

    /** Destructor. */
    virtual ~Preconditioner ();

    /** Builds a \p Preconditioner using the linear solver package specified by \p solver_package */
    static Preconditioner * build(const SolverPackage solver_package = PETSC_SOLVERS);

    /** @returns true if the data structures */
    bool initialized () const {
        return _is_initialized;
    }

    /** Computes the preconditioned vector "y" based on input "x". Usually by solving Py=x to get the action of P^-1 x. */
    virtual void apply(const NumericVector & x, NumericVector & y) = 0;

    /** Release all memory and clear data structures. */
    virtual void clear () {}

    /** Initialize data structures if not done so already. */
    virtual void init () {};

    /** Sets the matrix P to be preconditioned. */
    void set_matrix(SparseMatrix & mat);

    /** Returns the type of preconditioner to use. */
    PreconditionerType type () const {
        return _preconditioner_type;
    }

    /** Sets the type of preconditioner to use. */
    void set_type (const PreconditionerType pct);

protected:

    /** The matrix P... ie the matrix to be preconditioned. This is often the actual syst matrix of a linear sytem.*/
    SparseMatrix * _matrix;

    /** Enum statitng with type of preconditioner to use. */
    PreconditionerType _preconditioner_type;

    /** Flag indicating if the data structures have been initialized. */
    bool _is_initialized;

};

/**
 *----------------------- inline functions ----------------------------------
 */

inline Preconditioner::Preconditioner () :
    _matrix(NULL),
    _preconditioner_type (ILU_PRECOND),
    _is_initialized      (false) {
}

inline Preconditioner::~Preconditioner () {
    this->clear ();
}

inline void Preconditioner::set_matrix(SparseMatrix & mat) {
    //If the matrix is changing then we (probably) need to reinitialize.
    _is_initialized = false;
    _matrix = &mat;
}

inline void Preconditioner::set_type (const PreconditionerType pct) {
    //If the preconditioner type changes we (probably) need to reinitialize.
    _is_initialized = false;
    _preconditioner_type = pct;
}


} //end namespace femus



#endif
