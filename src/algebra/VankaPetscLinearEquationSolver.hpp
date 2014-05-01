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

#ifndef __VankaPetscLinearEquationSolver_hpp__
#define __VankaPetscLinearEquationSolver_hpp__

#include "FEMTTUConfig.h"

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
   */

  class VankaPetscLinearEquationSolver : public LinearEquationSolver {

  public:
  
    // Constructor --------------------------------------
    /**  Constructor. Initializes Petsc data structures */
    VankaPetscLinearEquationSolver (const unsigned &igrid, mesh *other_mesh);
  
    /// Destructor.
    ~VankaPetscLinearEquationSolver ();

    /// Release all memory and clear data structures.
    void clear ();

    /// Initialize data structures if not done so already plus much more
    void init (Mat& Amat, Mat &Pmat);
    // void init(Mat& matrix, const bool pc_flag, const bool Schur);
    //  void init_schur(Mat& matrix);

    void set_tolerances(const double &rtol, const double &atol,
			const double &divtol, const unsigned &maxits,const unsigned &index);

    void SetElementBlockNumber(const unsigned & block_elemet_number);
    void SetNumberOfSchurVariables(const unsigned short & NSchurVar){
      _NSchurVar=NSchurVar;
    }
  
    // Solvers ------------------------------------------------------
    // ========================================================
    /// Call the Vanka(Schur) smoother-solver using the PetscLibrary.
    std::pair< int, double> solve(const vector <unsigned> &VankaIndex, const bool &ksp_clean);
  private:
    // data ---------------------------------
    PC _pc;      ///< Preconditioner context
    vector <KSP> _ksp;    ///< Krylov subspace context
    vector< PetscReal > _rtol;
    vector< PetscReal > _abstol;
    vector< PetscReal > _dtol;
    vector< PetscInt >  _maxits;
    unsigned _block_element_number;
    
    short unsigned _NSchurVar;
  
    vector< vector <PetscInt> > _indexai;
  
    vector< vector <unsigned> > _Psize;
    bool _indexai_init;
  
    vector <IS> _isA;
    vector <IS> _isB;
  
    // Setting --------------------------------------------
    ///  Set the user-specified solver stored in \p _solver_type
    void set_petsc_solver_type ();
    void set_petsc_solver_type2 ();

    clock_t BuildIndex(const vector <unsigned> &VankaIndex);
  
  };

  inline VankaPetscLinearEquationSolver::VankaPetscLinearEquationSolver (const unsigned &igrid, mesh* other_msh)
    : LinearEquationSolver(igrid, other_msh) {
        
    if(igrid==0){
      this->_preconditioner_type = MLU_PRECOND;
      this->_solver_type         = PREONLY;
      _block_element_number = _msh->el->GetElementNumber();  
    }
    else{
      if(_msh->_nprocs==1) {  
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

    _ksp.resize(2);  
  
    _rtol.resize(2,1.e-8);
    _abstol.resize(2, 1.e-40);
    _dtol.resize(2, 1.e+50);
    _maxits.resize(2,10);
  
    _indexai_init=0;
    _NSchurVar=1;
   
  }

  // =============================================
  inline VankaPetscLinearEquationSolver::~VankaPetscLinearEquationSolver () {
    this->clear ();
  
    for(unsigned i=0;i<_isA.size();i++){
      ISDestroy(&_isA[i]); 	
    }
  
    for(unsigned i=0;i<_isB.size();i++){
      ISDestroy(&_isB[i]); 	
    }
    
  }

} //end namespace femus


#endif 
#endif 
