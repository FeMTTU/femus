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

#ifndef __femus_algebra_VankaPetscLinearEquationSolver_hpp__
#define __femus_algebra_VankaPetscLinearEquationSolver_hpp__

#include "FemusConfig.hpp"

#ifdef HAVE_PETSC

#ifdef HAVE_MPI
#include <mpi.h>
#endif

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "LinearEquationSolver.hpp"
#include "PetscVector.hpp"
#include "PetscMatrix.hpp"
#include "PetscMacro.hpp"


namespace femus {

/**
 * This class inherits the abstract class LinearEquationSolver. In this class the solver is implemented using the PETSc package
 **/

class VankaPetscLinearEquationSolver : public LinearEquationSolver {


public:

    /**  Constructor. Initializes Petsc data structures */
    VankaPetscLinearEquationSolver (const unsigned &igrid, Mesh *other_mesh);

    /** Destructor. */
    ~VankaPetscLinearEquationSolver ();
private:

    /** Release all memory and clear data structures. */
    void clear ();

    /** Initialize the ksp and pc object */
    void init (Mat& Amat, Mat &Pmat, KSP &ksp, PC &pc);

    /** To be Added */
    void set_tolerances(const double &rtol, const double &atol,
                        const double &divtol, const unsigned &maxits);

    /** To be Added */
    void SetElementBlockNumber(const unsigned & block_elemet_number);

    /** Set the number of Schur variables to be considered */
    void SetNumberOfSchurVariables(const unsigned short & NSchurVar) {
        _NSchurVar=NSchurVar;
    }

    /** Call the Vanka smoother-solver using the PetscLibrary */
    void solve(const vector <unsigned> &VankaIndex, const bool &ksp_clean);

    /**  Set the user-specified solver stored in _solver_type */
    void set_petsc_solver_type (KSP &skp);

    /** To be Added */
    clock_t BuildIndex(const vector <unsigned> &VankaIndex);

    // member data

    vector <PC>  _pc;     ///< Preconditioner context
    vector <KSP> _ksp;    ///< Krylov subspace context
    vector <Mat> _A;
    vector <Mat> _B;
    vector <Vec> _r;
    vector <Vec> _w;
    vector <VecScatter> _scatA;
    vector <Vec> _s;
    vector <VecScatter> _scatB;
    PetscReal _rtol;
    PetscReal _abstol;
    PetscReal _dtol;
    PetscInt  _maxits;
    unsigned _block_element_number;
    short unsigned _NSchurVar;
    vector< vector <PetscInt> > _indexai;
    vector< vector <unsigned> > _Psize;
    bool _indexai_init;
    vector <IS> _isA;
    vector <IS> _isB;

};

inline VankaPetscLinearEquationSolver::VankaPetscLinearEquationSolver (const unsigned &igrid, Mesh* other_msh)
    : LinearEquationSolver(igrid, other_msh) {

    if(igrid==0) {
        this->_preconditioner_type = MLU_PRECOND;
        this->_solver_type         = PREONLY;
        _block_element_number = _msh->el->GetElementNumber();
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
        _block_element_number = pow(base,exponent);
    }

    _rtol = 1.e-8;
    _abstol = 1.e-40;
    _dtol = 1.e+50;
    _maxits = 10;
    _indexai_init=0;
    _NSchurVar=1;

}

// =============================================
inline VankaPetscLinearEquationSolver::~VankaPetscLinearEquationSolver () {
    this->clear ();

    for(unsigned i=0; i<_isA.size(); i++) {
        ISDestroy(&_isA[i]);
    }

    for(unsigned i=0; i<_isB.size(); i++) {
        ISDestroy(&_isB[i]);
    }

}

} //end namespace femus


#endif
#endif
