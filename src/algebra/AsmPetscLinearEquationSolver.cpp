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

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "FEMTTUConfig.h"

#ifdef HAVE_PETSC

// Local Includes
#include "AsmPetscLinearEquationSolver.hpp"
#include "PetscPreconditioner.hpp"
#include "PetscMatrix.hpp"
#include <iomanip>

namespace femus {
  
  

  using namespace std;
  
  // ====================================================
  // ------------------- Class functions ------------
  // ====================================================

  void AsmPetscLinearEquationSolver::set_tolerances(const double &rtol, const double &atol,
						    const double &divtol, const unsigned &maxits) {

    _rtol   = static_cast<PetscReal>(rtol);
    _abstol = static_cast<PetscReal>(atol);
    _dtol   = static_cast<PetscReal>(divtol);
    _maxits = static_cast<PetscInt>(maxits);
						      
  }

  // ==============================================

  void AsmPetscLinearEquationSolver::SetElementBlockNumber(const unsigned& block_elemet_number) {
    _element_block_number = block_elemet_number;
    _indexai_init=0;
    _standard_ASM=0;
  }

  // ==============================================
  clock_t AsmPetscLinearEquationSolver::BuildBDCIndex(const vector <unsigned> &variable_to_be_solved){
    
    clock_t SearchTime = 0;
    clock_t start_time = clock();
    
    unsigned IndexaSize=KKoffset[KKIndex.size()-1][_msh->_iproc] - KKoffset[0][_msh->_iproc];
    _indexai.resize(2);
    
    _indexai[0].clear();
    _indexai[1].clear();
        
    _indexai[0].resize(IndexaSize);
    _indexai[1].resize(IndexaSize);
    
    unsigned count0=0;
    unsigned count1=0;
    
    
    vector <bool> ThisSolutionIsIncluded(_SolPdeIndex.size(),false);
    for (unsigned iind=0; iind<variable_to_be_solved.size(); iind++) {
      unsigned PdeIndexSol=variable_to_be_solved[iind];
      ThisSolutionIsIncluded[PdeIndexSol]=true;
    }
          
    for(int k=0; k < _SolPdeIndex.size(); k++) {
      unsigned indexSol = _SolPdeIndex[k];
      unsigned soltype = _SolType[indexSol];
      for(unsigned inode_mts = _msh->MetisOffset[soltype][_msh->_iproc]; inode_mts < _msh->MetisOffset[soltype][_msh->_iproc+1]; inode_mts++) {
	int local_mts = inode_mts-_msh->MetisOffset[soltype][_msh->_iproc];
	int idof_kk = KKoffset[k][_msh->_iproc] +local_mts; 
	if( !ThisSolutionIsIncluded[k] || (*(*_Bdc)[indexSol])(inode_mts) < 1.9) {
	  _indexai[0][count0] = idof_kk;
	  count0++;
	}
	else{
	  _indexai[1][count1] = idof_kk;
	  count1++;
	}
      }
    } 
    _indexai[0].resize(count0);
    _indexai[1].resize(count1);
    
    std::sort(_indexai[0].begin(), _indexai[0].end());
    std::sort(_indexai[1].begin(), _indexai[1].end());
    
    clock_t end_time=clock();
    SearchTime = (end_time-start_time);
    
    return SearchTime;
  }


  // ==============================================
  
  clock_t AsmPetscLinearEquationSolver::BuildAMSIndex(const vector <unsigned> &variable_to_be_solved){
    clock_t SearchTime=0;
    clock_t start_time=clock();
    
    unsigned nel=_msh->GetElementNumber();
    
    bool FastVankaBlock=true;
    if(_NSchurVar==!0){
      FastVankaBlock=(_SolType[_SolPdeIndex[variable_to_be_solved[variable_to_be_solved.size()-_NSchurVar]]]<3)?false:true;
    }
    
    unsigned iproc=_msh->_iproc;
    
    unsigned DofOffset = KKoffset[0][iproc];
    unsigned DofOffsetSize=KKoffset[KKIndex.size()-1][iproc] - KKoffset[0][iproc];
    vector < unsigned > indexa(DofOffsetSize,DofOffsetSize);
    vector < unsigned > indexb(DofOffsetSize,DofOffsetSize);
    vector <bool> owned(DofOffsetSize,false);
    
    map<int,bool> mymap;
    
    unsigned ElemOffset   = _msh->MetisOffset[3][iproc];
    unsigned ElemOffsetp1 = _msh->MetisOffset[3][iproc+1];
    unsigned ElemOffsetSize= ElemOffsetp1 - ElemOffset;
    vector <PetscInt> indexci(ElemOffsetSize);
    vector < unsigned > indexc(ElemOffsetSize,ElemOffsetSize);
    
    
    vector <bool> ThisVaribaleIsNonSchur(_SolPdeIndex.size(),true);
    for (unsigned iind=variable_to_be_solved.size()-_NSchurVar; iind<variable_to_be_solved.size(); iind++) {
      unsigned PdeIndexSol=variable_to_be_solved[iind];
      ThisVaribaleIsNonSchur[PdeIndexSol]=false;
    }
    
    // *** Start Vanka Block ***
    bool test_end=0;
    int vb_index=0;
    while (test_end==0){
      
      _is_loc_idx.resize(vb_index+1);
      _is_loc_idx[vb_index].resize(DofOffsetSize);
      
      _is_ovl_idx.resize(vb_index+1);	
      _is_ovl_idx[vb_index].resize(DofOffsetSize);
      
      PetscInt PAsize=0;
      PetscInt PBsize=0;
      
      PetscInt Csize=0;
      test_end=1;
      
      int gel=_msh->IS_Mts2Gmt_elem_offset[iproc] + vb_index*_element_block_number;
      // ***************** NODE/ELEMENT SERCH *******************
      
      for (int iel_mts=gel; iel_mts<gel+_element_block_number && iel_mts< _msh->IS_Mts2Gmt_elem_offset[iproc+1]; iel_mts++) {
	unsigned iel = _msh->IS_Mts2Gmt_elem[iel_mts];     
	
	for (unsigned i=0; i<_msh->el->GetElementDofNumber(iel,0); i++) {
	  unsigned inode=_msh->el->GetElementVertexIndex(iel,i)-1u;
	  unsigned nvei=_msh->el->GetVertexElementNumber(inode);
	  const unsigned *pt_jel=_msh->el->GetVertexElementAddress(inode,0);
	  for (unsigned j=0; j<nvei*(!FastVankaBlock)+FastVankaBlock; j++) {
	    unsigned jel=(!FastVankaBlock)?*(pt_jel++)-1u:iel;
	    //add elements for velocity to be solved
	    
	    unsigned jel_Metis = _msh->GetMetisDof(jel,3);
	    
	    if(jel_Metis >= ElemOffset && jel_Metis < ElemOffsetp1){
	      if(indexc[jel_Metis- ElemOffset]==ElemOffsetSize){
		indexci[Csize]=jel_Metis-ElemOffset;
		indexc[jel_Metis-ElemOffset]=Csize++;
		//----------------------------------------------------------------------------------
		//add non-schur node to be solved
		
		//for (unsigned iind=0; iind<variable_to_be_solved.size()-_NSchurVar; iind++) {
		// unsigned indexSol=variable_to_be_solved[iind];
		for(int indexSol=0; indexSol < _SolPdeIndex.size(); indexSol++){ 
		  if(ThisVaribaleIsNonSchur[indexSol]){
		    unsigned SolPdeIndex = _SolPdeIndex[indexSol];
		    unsigned SolType = _SolType[SolPdeIndex];
		    const unsigned *pt_un=_msh->el->GetElementVertexAddress(jel,0);
		    unsigned nvej=_msh->el->GetElementDofNumber(jel,_msh->_END_IND[SolType]);
		    for (unsigned jj=0; jj<nvej; jj++) {
		      unsigned jnode=(SolType<3)?(*(pt_un++)-1u):(jel+jj*nel);
		      unsigned jnode_Metis = _msh->GetMetisDof(jnode,SolType);
		      unsigned kkdof=GetKKDof(SolPdeIndex, indexSol, jnode);
		      if(jnode_Metis >= _msh->MetisOffset[SolType][iproc] &&
			jnode_Metis <  _msh->MetisOffset[SolType][iproc+1]){
			//unsigned kkdof=GetKKDof(SolPdeIndex, indexSol, jnode);
			if(indexa[kkdof- DofOffset]==DofOffsetSize && owned[kkdof- DofOffset]==false) {
			  owned[kkdof- DofOffset]=true;
			  _is_loc_idx[vb_index][PAsize]=kkdof;
			  indexa[kkdof-DofOffset]=PAsize++;
			}
			if(indexb[kkdof- DofOffset]==DofOffsetSize) {
			  _is_ovl_idx[vb_index][PBsize]=kkdof;
			  indexb[kkdof-DofOffset]=PBsize++;
			}
		      }
		      else mymap[kkdof]=true;
		    }
		  }
		}
	      }
	    }
	  }
	}
	//-----------------------------------------------------------------------------------------
	//Add Schur nodes (generally pressure) to be solved
	{
	  //for (unsigned iind=variable_to_be_solved.size()-_NSchurVar; iind<variable_to_be_solved.size(); iind++) {
	  //  unsigned indexSol=variable_to_be_solved[iind];
	  for(int indexSol=0; indexSol < _SolPdeIndex.size(); indexSol++){ 
	    if(!ThisVaribaleIsNonSchur[indexSol]){
	      unsigned SolPdeIndex = _SolPdeIndex[indexSol];
	      unsigned SolType = _SolType[SolPdeIndex];
	      const unsigned *pt_un=_msh->el->GetElementVertexAddress(iel,0);
	      unsigned nvei=_msh->el->GetElementDofNumber(iel,_msh->_END_IND[SolType]);
	      for (unsigned ii=0; ii<nvei; ii++) {
		unsigned inode=(SolType<3)?(*(pt_un++)-1u):(iel+ii*nel);
		unsigned inode_Metis = _msh->GetMetisDof(inode,SolType);
		unsigned kkdof=GetKKDof(SolPdeIndex, indexSol, inode);
		if(inode_Metis >= _msh->MetisOffset[SolType][iproc] &&
		  inode_Metis <  _msh->MetisOffset[SolType][iproc+1]){
		  if(indexa[kkdof- DofOffset]==DofOffsetSize && owned[kkdof- DofOffset]==false) {
		    owned[kkdof- DofOffset]=true;
		    _is_loc_idx[vb_index][PAsize]=kkdof;
		    indexa[kkdof-DofOffset]=PAsize++;
		  }
		  if(indexb[kkdof- DofOffset]==DofOffsetSize ) {
		    _is_ovl_idx[vb_index][PBsize]=kkdof;
		    indexb[kkdof-DofOffset]=PBsize++;
		  }		
		}
		else mymap[kkdof]=true;
	      }
	    }
	  }
	}
	//-----------------------------------------------------------------------------------------
      }
      if(gel+_element_block_number <_msh->IS_Mts2Gmt_elem_offset[iproc+1] ) test_end=0;     
      
      // *** re-initialize indeces(a,c,d)
      for (PetscInt i=0; i<PAsize; i++) {
	indexa[_is_loc_idx[vb_index][i]-DofOffset]=DofOffsetSize;
      }
      for (PetscInt i=0; i<PBsize; i++) {
	indexb[_is_ovl_idx[vb_index][i]-DofOffset]=DofOffsetSize;
      }
      for (PetscInt i=0; i<Csize; i++) {
	indexc[indexci[i]]=ElemOffsetSize;
      }
      
      _is_loc_idx[vb_index].resize(PAsize);
      
      _is_ovl_idx[vb_index].resize(PBsize+mymap.size());
      
      int i=0;
      for (std::map<int,bool>::iterator it=mymap.begin(); it!=mymap.end(); ++it,++i){
	_is_ovl_idx[vb_index][PBsize+i]= it->first;
      }       
      
      std::sort(_is_loc_idx[vb_index].begin(), _is_loc_idx[vb_index].end());
      std::sort(_is_ovl_idx[vb_index].begin(), _is_ovl_idx[vb_index].end());
      
      vb_index++;
    }
    
    //BEGIN Generate std::vector<IS> for vanka solve ***********
    _is_loc.resize(_is_loc_idx.size());
    _is_ovl.resize(_is_ovl_idx.size());
    for(unsigned vb_index=0;vb_index<_is_loc_idx.size();vb_index++){  
      PetscErrorCode ierr;  
      ierr=ISCreateGeneral(MPI_COMM_SELF,_is_loc_idx[vb_index].size(),&_is_loc_idx[vb_index][0],PETSC_USE_POINTER,&_is_loc[vb_index]);
      CHKERRABORT(MPI_COMM_SELF,ierr);
      ierr=ISCreateGeneral(MPI_COMM_SELF,_is_ovl_idx[vb_index].size(),&_is_ovl_idx[vb_index][0],PETSC_USE_POINTER,&_is_ovl[vb_index]);
      CHKERRABORT(MPI_COMM_SELF,ierr);
    }
    //END Generate std::vector<IS> for vanka solve ***********  
    
    
    clock_t end_time=clock();
    SearchTime+=(end_time-start_time);
    return SearchTime;
  }
  
  // =================================================
  
  std::pair< int, double> AsmPetscLinearEquationSolver::solve(const vector <unsigned> &variable_to_be_solved, const bool &ksp_clean) {
    PetscVector* EPSCp=static_cast<PetscVector*> (_EPSC);  
    Vec EPSC=EPSCp->vec(); 
    PetscVector* RESp=static_cast<PetscVector*> (_RES);  
    Vec RES=RESp->vec();
    PetscMatrix* KKp=static_cast<PetscMatrix*>(_KK); 
    Mat KK=KKp->mat(); 
    
    PetscErrorCode ierr;
    clock_t SearchTime, AssemblyTime, SolveTime, UpdateTime;
    int its;
    double final_resid;
    
    // ***************** NODE/ELEMENT SEARCH *******************
    clock_t start_time=clock();
    if(_indexai_init==0) {
      _indexai_init = 1;
      if(!_standard_ASM)
      	BuildAMSIndex(variable_to_be_solved); 
      BuildBDCIndex(variable_to_be_solved);
    }
    SearchTime = start_time - clock();
    // ***************** END NODE/ELEMENT SEARCH *******************  
    
    // ***************** ASSEMBLE matrix to set Dirichlet BCs by penalty *******************
    start_time=clock();
    if(ksp_clean){
      this->clear();
      // initialize Pmat wiwth penaly diagonal on the Dirichlet Nodes
      MatDuplicate(KK,MAT_COPY_VALUES,&_Pmat);
      MatSetOption(_Pmat,MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE);
      MatZeroRows(_Pmat,_indexai[0].size(),&_indexai[0][0],1.e40,0,0);
      _Pmat_is_initialized = true;
      init(KK,_Pmat);
    }
    
    AssemblyTime = clock() - start_time;
    
    // ***************** END ASSEMBLE ***********
    
    // ***************** SOLVE ******************
    start_time = clock();
    
    ierr = KSPSolve(_ksp, RES, EPSC);			CHKERRABORT(MPI_COMM_WORLD,ierr);
    
    // Get the number of iterations required for convergence
    ierr = KSPGetIterationNumber(_ksp, &its);		CHKERRABORT(MPI_COMM_WORLD,ierr);
    
    // Get the norm of the final residual to return to the user.
    ierr = KSPGetResidualNorm(_ksp, &final_resid); 	CHKERRABORT(MPI_COMM_WORLD,ierr);
    
    SolveTime = clock()-start_time;
    // ***************** END SOLVE ******************
    
    // ***************** RES/EPS UPDATE RES ******************
    start_time=clock();
    
    *_EPS += *_EPSC;
    
    _RESC->matrix_mult(*_EPSC,*_KK);
    *_RES -= *_RESC;
    
    UpdateTime = clock() - start_time;
    
    // **************** END RES/EPS UPDATE RES ***************
    
    // *** Computational info ***
#ifndef NDEBUG   
    cout << "ASM Grid: " << _msh->GetGridNumber()<< "        SOLVER TIME:        "  << std::setw(11) << std::setprecision(6) << std::fixed <<
      static_cast<double>( SearchTime + AssemblyTime + SolveTime + UpdateTime)/ CLOCKS_PER_SEC<<
      "  ITS: " << its  << "\t ksp_clean = "<< ksp_clean<<endl;
#endif
    
    return std::make_pair(its,final_resid);
  }
  
  
  void AsmPetscLinearEquationSolver::clear() {
    int ierr=0;
    if(_Pmat_is_initialized){
      _Pmat_is_initialized = false;
      ierr = MatDestroy(&_Pmat);          CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
    
    if (this->initialized()) {
      this->_is_initialized = false;
      ierr = KSPDestroy(&_ksp);		CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
  }


  // ========================================================

  void AsmPetscLinearEquationSolver::init(Mat& Amat, Mat& Pmat) {
  
    // Initialize the data structures if not done so already.
    if (!this->initialized())    {
      this->_is_initialized = true;
      int ierr=0;
      // Create the linear solver context
      ierr = KSPCreate(MPI_COMM_WORLD, &_ksp);					CHKERRABORT(MPI_COMM_WORLD,ierr);
      // Create the preconditioner context
      ierr = KSPGetPC(_ksp, &_pc);						CHKERRABORT(MPI_COMM_WORLD,ierr);
      // Set operators. The input matrix works as the preconditioning matrix
      this->set_petsc_solver_type();
      
      ierr = KSPSetOperators(_ksp, Amat, Pmat, SAME_PRECONDITIONER);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      
      // Set the tolerances for the iterative solver.  Use the user-supplied
      // tolerance for the relative residual & leave the others at default values.
      ierr = KSPSetTolerances(_ksp,_rtol,_abstol,_dtol,_maxits);	CHKERRABORT(MPI_COMM_WORLD,ierr);
      
      if(_msh->GetGridNumber()!=0)
	KSPSetInitialGuessKnoll(_ksp, PETSC_TRUE);
      
      if(_msh->GetGridNumber()!=0)
	KSPSetNormType(_ksp,KSP_NORM_NONE);
      
      // Set the options from user-input
      // Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  KSPSetFromOptions() is called _after_ any other customization  routines.
      ierr = KSPSetFromOptions(_ksp);						CHKERRABORT(MPI_COMM_WORLD,ierr);
      // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
      //ierr = PCSetFromOptions (_pc);	CHKERRABORT(MPI_COMM_WORLD,ierr);
      
      // Notify PETSc of location to store residual history.
      // This needs to be called before any solves, since
      // it sets the residual history length to zero.  The default
      // behavior is for PETSc to allocate (internally) an array
      // of size 1000 to hold the residual norm history.
      ierr = KSPSetResidualHistory(_ksp,
				   PETSC_NULL,   // pointer to the array which holds the history
				   PETSC_DECIDE, // size of the array holding the history
				   PETSC_TRUE);  // Whether or not to reset the history for each solve.
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      
      PetscPreconditioner::set_petsc_preconditioner_type(ASM_PRECOND,_pc);
      PetscReal zero = 1.e-16;
      PCFactorSetZeroPivot(_pc,zero);
      PCFactorSetShiftType(_pc,MAT_SHIFT_NONZERO);
      if(!_standard_ASM){
      	ierr = PCASMSetLocalSubdomains(_pc,_is_loc_idx.size(),&_is_ovl[0],&_is_loc[0]); CHKERRABORT(MPI_COMM_WORLD,ierr);
      }
      else{
	ierr = PCASMSetOverlap(_pc,_overlap); CHKERRABORT(MPI_COMM_WORLD,ierr);
      }
      
      ierr = KSPSetUp(_ksp);							    CHKERRABORT(MPI_COMM_WORLD,ierr);
      
      ierr = PCASMGetSubKSP(_pc,&_nlocal,&_first,&_subksp);			    CHKERRABORT(MPI_COMM_WORLD,ierr);
      
      for (int i=0; i<_nlocal; i++) {
	ierr = KSPGetPC(_subksp[i],&_subpc);					    CHKERRABORT(MPI_COMM_WORLD,ierr);
	
	PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type,_subpc); 
	ierr = KSPSetTolerances(_subksp[i],_rtol,_abstol,_dtol,1); 	    CHKERRABORT(MPI_COMM_WORLD,ierr);          
	
	ierr = KSPSetFromOptions(_subksp[i]);
	PetscReal zero = 1.e-16;
	PCFactorSetZeroPivot(_subpc,zero);
	PCFactorSetShiftType(_subpc,MAT_SHIFT_NONZERO);
      }
    }
  }

  // =================================================
  void AsmPetscLinearEquationSolver::set_petsc_solver_type() {
    int ierr = 0;
    switch (this->_solver_type) {
    case CG:
      ierr = KSPSetType(_ksp, (char*) KSPCG);						CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case CR:
      ierr = KSPSetType(_ksp, (char*) KSPCR);						CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case CGS:
      ierr = KSPSetType(_ksp, (char*) KSPCGS);						CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case BICG:
      ierr = KSPSetType(_ksp, (char*) KSPBICG);					CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case TCQMR:
      ierr = KSPSetType(_ksp, (char*) KSPTCQMR);					CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case TFQMR:
      ierr = KSPSetType(_ksp, (char*) KSPTFQMR);					CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case LSQR:
      ierr = KSPSetType(_ksp, (char*) KSPLSQR);					CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case BICGSTAB:
      ierr = KSPSetType(_ksp, (char*) KSPBCGS);					CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case MINRES:
      ierr = KSPSetType(_ksp, (char*) KSPMINRES);					CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case GMRES:
      ierr = KSPSetType(_ksp, (char*) KSPGMRES);					CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case RICHARDSON:
      ierr = KSPSetType(_ksp, (char*) KSPRICHARDSON);					CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case CHEBYSHEV:
      ierr = KSPSetType(_ksp, (char*) KSPCHEBYSHEV);					CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case PREONLY:
      ierr = KSPSetType(_ksp, (char*) KSPPREONLY);					CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    default:
      std::cerr << "ERROR:  Unsupported PETSC Solver: "
		<< this->_solver_type               << std::endl
		<< "Continuing with PETSC defaults" << std::endl;
    }
  }
} //end namespace femus

#endif 

