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

#ifndef __femus_algebra_AsmPetscLinearEquationSolver_hpp__
#define __femus_algebra_AsmPetscLinearEquationSolver_hpp__

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

namespace femus {

/**
 * This class inherits the abstract class LinearEquationSolver. In this class the solver is implemented using the PETSc package
 **/

class AsmPetscLinearEquationSolver : public LinearEquationSolver {

public:

    /**  Constructor. Initializes Petsc data structures */
    AsmPetscLinearEquationSolver (const unsigned &igrid, Mesh *other_mesh);

    /** Destructor */
    ~AsmPetscLinearEquationSolver ();
    KSP* GetKSP(){ return &_ksp; };

    void MGsetLevels ( LinearEquationSolver *LinSolver, const unsigned &level, const unsigned &maxlevel,
                      const vector <unsigned> &variable_to_be_solved,
                      SparseMatrix* PP, SparseMatrix* RR ,
                      const unsigned &npre, const unsigned &npost);

    void MGsolve ( const bool ksp_clean );

    void MGinit( const MgSmootherType &mg_smoother_type, const unsigned &levelMax ){

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

    /** Release all memory and clear data structures. */
    void clear ();

    /** Initialize the ksp objects plus much more */
    void init (Mat& Amat, Mat &Pmat);

    /** To be Added */
    void set_tolerances(const double &rtol, const double &atol,
                        const double &divtol, const unsigned &maxits);

    /** To be Added */
    void SetElementBlockNumber(const unsigned & block_elemet_number);
    void SetElementBlockNumberSolid(const unsigned & block_elemet_number, const unsigned & overlap);
    void SetElementBlockNumberFluid(const unsigned & block_elemet_number, const unsigned & overlap);

    /** To be Added */
    void SetElementBlockNumber(const char all[], const unsigned & overlap=1);

    /** To be Added */
    void SetNumberOfSchurVariables(const unsigned short & NSchurVar) {
        _NSchurVar=NSchurVar;
    };

    /** Call the ASM smoother-solver using the PetscLibrary */
    void solve(const vector <unsigned> &variable_to_be_solved, const bool &ksp_clean);

    /**  Set the user-specified solver stored in \p _solver_type */
    void set_petsc_solver_type ( KSP &ksp );

    /** To be Added */
    clock_t BuildBDCIndex(const vector <unsigned> &variable_to_be_solved);

    /** To be Added */
    clock_t BuildAMSIndex(const vector <unsigned> &variable_to_be_solved);

    // member data

    PC _pc;      ///< Preconditioner context
    KSP _ksp;    ///< Krylov subspace context
    KSP       *_ksp_asm;
    vector < PC >  _pc_asm;

    PetscReal  _rtol;
    PetscReal  _abstol;
    PetscReal  _dtol;
    PetscInt   _maxits;
    unsigned _element_block_number[2];
    vector< vector <PetscInt> > _indexai;
    bool _indexai_init;
    unsigned short _NSchurVar;
    vector< vector <PetscInt> > _is_ovl_idx;

    vector< vector <PetscInt> > _is_loc_idx;
    vector <IS> _is_ovl;
    vector <IS> _is_loc;
    PetscInt  _nlocal,_first;
    bool _standard_ASM;
    unsigned _overlap;
    Mat _Pmat;
    bool _Pmat_is_initialized;
    vector <unsigned> _block_type_range;


    //vector < KSP*> _ksp_split;
    //vector< vector < PC > >  _pc_split;
    //vector< vector <PetscInt> > _is_ovl_u_idx;
    //vector< vector <PetscInt> > _is_ovl_p_idx;
    //vector <IS> _is_ovl_u;
    //vector <IS> _is_ovl_p;

};

// =================================================

inline AsmPetscLinearEquationSolver::AsmPetscLinearEquationSolver (const unsigned &igrid, Mesh* other_msh)
    : LinearEquationSolver(igrid, other_msh) {

    if(igrid==0) {
        this->_preconditioner_type = MLU_PRECOND;
        this->_solver_type         = PREONLY;
        _element_block_number[0] = _msh->el->GetElementNumber();
	_element_block_number[1] = _msh->el->GetElementNumber();
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
        _element_block_number[0] = pow(base,exponent);
	_element_block_number[1] = pow(base,exponent);
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
    _overlap=0;

}

// =============================================

inline AsmPetscLinearEquationSolver::~AsmPetscLinearEquationSolver () {
    this->clear ();

    for(unsigned i=0; i<_is_loc.size(); i++) {
        ISDestroy(&_is_loc[i]);
    }

    for(unsigned i=0; i<_is_ovl.size(); i++) {
        ISDestroy(&_is_ovl[i]);
    }

}

} //end namespace femus


#endif
#endif
