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

#ifndef __femus_algebra_GmresPetscLinearEquationSolver_hpp__
#define __femus_algebra_GmresPetscLinearEquationSolver_hpp__

#include "FemusConfig.hpp"

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
  GmresPetscLinearEquationSolver ( const unsigned &igrid, Mesh *other_mesh );

  /// Destructor.
  ~GmresPetscLinearEquationSolver ();

  /// Release all memory and clear data structures.
  void clear ();

  /// Initialize data structures if not done so already plus much more
  void init ( Mat& Amat, Mat &Pmat );


  void set_tolerances ( const double &rtol, const double &atol,
                        const double &divtol, const unsigned &maxits );


  void SetDirichletBCsHandling ( const unsigned int &DirichletBCsHandlingMode ) {
    _DirichletBCsHandlingMode = DirichletBCsHandlingMode;
  }

  // Solvers ------------------------------------------------------
  // ========================================================
  /// Call the GMRES smoother-solver using the PetscLibrary.
  void solve ( const vector <unsigned> &VankaIndex, const bool &ksp_clean );

  // Setting --------------------------------------------
  ///  Set the user-specified solver stored in \p _solver_type
  void set_petsc_solver_type ( KSP &ksp);

  clock_t BuildIndex ( const vector <unsigned> &variable_to_be_solved );

  /** @deprecated, remove soon */
  std::pair<unsigned int, double> solve ( SparseMatrix&  matrix_in,
                                          SparseMatrix&  precond_in,  NumericVector& solution_in,  NumericVector& rhs_in,
                                          const double tol,   const unsigned int m_its );

  /** @deprecated, remove soon */
  void init ( SparseMatrix* matrix );

  KSP* GetKSP() {
    return &_ksp;
  };

  void MGsetLevels ( LinearEquationSolver *LinSolver, const unsigned &level, const unsigned &maxlevel,
                      const vector <unsigned> &variable_to_be_solved,
                      SparseMatrix* PP, SparseMatrix* RR,
                      const unsigned &npre, const unsigned &npost );

  void MGsolve ( const bool ksp_clean );

  void MGinit( const MgSmootherType & mg_smoother_type, const unsigned &levelMax ){

    KSPCreate(PETSC_COMM_WORLD,&_ksp);

    KSPGetPC(_ksp,&_pc);
    PCSetType(_pc,PCMG);
    PCMGSetLevels(_pc,levelMax,NULL);

    if( mg_smoother_type == FULL ){
      PCMGSetType(_pc, PC_MG_FULL);
    }
    else if( mg_smoother_type == MULTIPLICATIVE ){
      PCMGSetType(_pc, PC_MG_MULTIPLICATIVE);
    }
    else if( mg_smoother_type == ADDITIVE ){
      PCMGSetType(_pc, PC_MG_ADDITIVE);
    }
    else if( mg_smoother_type == KASKADE ){
      PCMGSetType(_pc, PC_MG_KASKADE);
    }
    else{
      std::cout <<"Wrong mg_type for PETSCsolve()"<<std::endl;
      abort();
    }

  };

  void MGclear(){
    KSPDestroy(&_ksp);
  }

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

inline GmresPetscLinearEquationSolver::GmresPetscLinearEquationSolver ( const unsigned &igrid, Mesh* other_msh )
  : LinearEquationSolver ( igrid, other_msh ) {

  if ( igrid == 0 ) {
    this->_preconditioner_type = MLU_PRECOND;
    this->_solver_type         = PREONLY;
  } else {
    if ( _msh->n_processors() == 1 ) {
      this->_preconditioner_type = ILU_PRECOND;
    } else {
      //this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
      this->_preconditioner_type = ASM_PRECOND;
    }
    unsigned dim = _msh->GetDimension();
    unsigned base = pow ( 2, dim );
    unsigned exponent = 5 - dim;
  }

  _rtol   = 1.e-8;
  _abstol = 1.e-40;
  _dtol   = 1.e+50;
  _maxits = 4;

  _indexai_init = 0;

  _Pmat_is_initialized = false;
  _Pw_is_initialized = false;
  _scat_is_initialized = false;

  _DirichletBCsHandlingMode = 0;

}

// =============================================
inline GmresPetscLinearEquationSolver::~GmresPetscLinearEquationSolver () {
  this->clear ();

  for ( unsigned i = 0; i < _isA.size(); i++ ) {
    ISDestroy ( &_isA[i] );
  }

}

} //end namespace femus


#endif
#endif
