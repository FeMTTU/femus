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

#ifndef __AsmPetscLinearEquationSolver_hpp__
#define __AsmPetscLinearEquationSolver_hpp__

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

namespace femus {
  
  /**
   * This class inherits the abstract class LinearEquationSolver. In this class the solver is implemented using the PETSc package
   **/

  class AsmPetscLinearEquationSolver : public LinearEquationSolver {

  private:
    // data ---------------------------------
    PC _pc;      ///< Preconditioner context
    KSP _ksp;    ///< Krylov subspace context
    PetscReal  _rtol;
    PetscReal  _abstol;
    PetscReal  _dtol;
    PetscInt   _maxits;
  
    unsigned _element_block_number;
  
    vector< vector <PetscInt> > _indexai;
    bool _indexai_init;
   
    unsigned short _NSchurVar;
    
    vector< vector <PetscInt> > _is_ovl_idx;
    vector< vector <PetscInt> > _is_loc_idx;
  
    vector <IS> _is_ovl;
    vector <IS> _is_loc;
   
    KSP       *_subksp;  							
    PetscInt  _nlocal,_first;   						
    PC        _subpc;   		
    
    bool _standard_ASM;
    unsigned _overlap;
    
    Mat _Pmat;
    bool _Pmat_is_initialized;
      
  public:
    // Constructor --------------------------------------
    /**  Constructor. Initializes Petsc data structures */
    AsmPetscLinearEquationSolver (const unsigned &igrid, mesh *other_mesh);
  
    /// Destructor.
    ~AsmPetscLinearEquationSolver ();

  private:
    /// Release all memory and clear data structures.
    void clear ();

    /// Initialize the ksp objects plus much more
    void init (Mat& Amat, Mat &Pmat);

    void set_tolerances(const double &rtol, const double &atol,
			const double &divtol, const unsigned &maxits,const unsigned &index);
    
    void SetElementBlockNumber(const unsigned & block_elemet_number);
    void SetElementBlockNumber(const char all[], const unsigned & overlap=1){
      _standard_ASM=1;
      _overlap=overlap;
    }
      
    void SetSchurVariableNumber(const unsigned short & NSchurVar){ _NSchurVar=NSchurVar; };
    // Solvers ------------------------------------------------------
    // ========================================================
    /// Call the ASM smoother-solver using the PetscLibrary.
    std::pair< int, double> solve(const vector <unsigned> &VankaIndex, const bool &ksp_clean);
    
    
    // Setting --------------------------------------------
    ///  Set the user-specified solver stored in \p _solver_type
    void set_petsc_solver_type ();
  
  
    clock_t BuildIndex(const vector <unsigned> &VankaIndex);
    clock_t BuildAMSIndex(const vector <unsigned> &VankaIndex);
  
  };

  // =================================================

  inline AsmPetscLinearEquationSolver::AsmPetscLinearEquationSolver (const unsigned &igrid, mesh* other_msh)
  : LinearEquationSolver(igrid, other_msh) {
        
    if(igrid==0){
      this->_preconditioner_type = MLU_PRECOND;
      this->_solver_type         = PREONLY;
      _element_block_number = _msh->el->GetElementNumber();  
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
      _element_block_number = pow(base,exponent);
    }

    _ksp;  
  
    _rtol = 1.e-8;
    _abstol = 1.e-40;
    _dtol = 1.e+50;
    _maxits = 4;
  
    _indexai_init=0;
  
    _Pmat_is_initialized = false;
    _NSchurVar=1;
    
    _standard_ASM=1;
    _overlap=1;

  }

  // =============================================
  
  inline AsmPetscLinearEquationSolver::~AsmPetscLinearEquationSolver () {
    this->clear ();
    
    for(unsigned i=0;i<_is_loc.size();i++){
      ISDestroy(&_is_loc[i]); 	
    }
    
    for(unsigned i=0;i<_is_ovl.size();i++){
      ISDestroy(&_is_ovl[i]); 	
    }
  
  }

} //end namespace femus


#endif 
#endif 
