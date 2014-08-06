/*=========================================================================

  Program: FEMUS
  Module: PetscLinearEquationSolver
  Authors: Eugenio Aulisa, Simone Bnà

  Copyright (c) FEMTTU
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/

#ifndef __GmresPetscLinearEquationSolver_hpp__
#define __GmresPetscLinearEquationSolver_hpp__

#include "FEMTTUConfig.h"

#ifdef HAVE_PETSC

#ifdef HAVE_MPI
#include <mpi.h>
#endif

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "LinearEquationSolver.hpp"

namespace femus {

/**
 * This class inherits the abstract class LinearEquationSolver. In this class the solver is implemented using the PETSc package
 */

class GmresPetscLinearEquationSolver : public LinearEquationSolver {

public:

    /**  Constructor. Initializes Petsc data structures */
    GmresPetscLinearEquationSolver (const unsigned &igrid, mesh *other_mesh);

    /// Destructor.
    ~GmresPetscLinearEquationSolver ();

    /// Release all memory and clear data structures.
    void clear ();

    /// Initialize data structures if not done so already plus much more
    void init (Mat& Amat, Mat &Pmat);


    void set_tolerances(const double &rtol, const double &atol,
                        const double &divtol, const unsigned &maxits);


    void SetDirichletBCsHandling(const unsigned int &DirichletBCsHandlingMode) {
        _DirichletBCsHandlingMode = DirichletBCsHandlingMode;
    }

    // Solvers ------------------------------------------------------
    // ========================================================
    /// Call the GMRES smoother-solver using the PetscLibrary.
    std::pair< int, double> solve(const vector <unsigned> &VankaIndex, const bool &ksp_clean);

    // Setting --------------------------------------------
    ///  Set the user-specified solver stored in \p _solver_type
    void set_petsc_solver_type ();

    clock_t BuildIndex(const vector <unsigned> &variable_to_be_solved);

private:

    // member data
    PC _pc;      ///< Preconditioner context
    KSP _ksp;    ///< Krylov subspace context
    PetscReal _rtol;
    PetscReal _abstol;
    PetscReal _dtol;
    PetscInt  _maxits;
    vector< vector <PetscInt> > _indexai;
    bool _indexai_init;
    vector <IS> _isA;
    Vec _Pw;
    bool _Pw_is_initialized;
    VecScatter _scat;
    bool _scat_is_initialized;
    Mat _Pmat;
    bool _Pmat_is_initialized;
    unsigned int _DirichletBCsHandlingMode; //* 0 Penalty method,  1 Elimination method */

};

inline GmresPetscLinearEquationSolver::GmresPetscLinearEquationSolver (const unsigned &igrid, mesh* other_msh)
    : LinearEquationSolver(igrid, other_msh) {

    if(igrid==0) {
        this->_preconditioner_type = MLU_PRECOND;
        this->_solver_type         = PREONLY;
    }
    else {
        if(_msh->n_processors()==1) {
            this->_preconditioner_type = ILU_PRECOND;
        }
        else {
            this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
        }
        unsigned dim = _msh->GetDimension();
        unsigned base = pow(2,dim);
        unsigned exponent = 5 - dim;
    }

    _rtol   = 1.e-8;
    _abstol = 1.e-40;
    _dtol   = 1.e+50;
    _maxits = 4;

    _indexai_init=0;

    _Pmat_is_initialized = false;
    _Pw_is_initialized = false;
    _scat_is_initialized = false;

    _DirichletBCsHandlingMode = 0;

}

// =============================================
inline GmresPetscLinearEquationSolver::~GmresPetscLinearEquationSolver () {
    this->clear ();

    for(unsigned i=0; i<_isA.size(); i++) {
        ISDestroy(&_isA[i]);
    }

}

} //end namespace femus


#endif
#endif
