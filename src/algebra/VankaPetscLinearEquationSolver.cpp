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
#include "PetscMacro.hpp"
#include "VankaPetscLinearEquationSolver.hpp"
#include "PetscPreconditioner.hpp"
#include "PetscVector.hpp"
#include "PetscMatrix.hpp"
#include <iomanip>


namespace femus {

  using namespace std;

  // ==============================================
  // ----------------------- functions ------
  // ==============================================

  void VankaPetscLinearEquationSolver::set_tolerances(const double &rtol, const double &atol,
						      const double &divtol, const unsigned &maxits) {
    
    _rtol = static_cast<PetscReal>(rtol);
    _abstol = static_cast<PetscReal>(atol);
    _dtol = static_cast<PetscReal>(divtol);
    _maxits = static_cast<PetscInt>(maxits);

  };

  // ========================================================

  void VankaPetscLinearEquationSolver::SetElementBlockNumber(const unsigned & block_elemet_number){
    _block_element_number = block_elemet_number;
    _indexai_init=0;
  };


  // ========================================================

  clock_t VankaPetscLinearEquationSolver::BuildIndex(const vector <unsigned> &VankaIndex){
    clock_t SearchTime=0;
    clock_t start_time=clock();
    _indexai_init=1;
    unsigned nel=_msh->GetNumberOfElements();
    bool FastVankaBlock=true;
    if(_NSchurVar==!0){
      FastVankaBlock=(_SolType[_SolPdeIndex[VankaIndex[VankaIndex.size()-_NSchurVar]]]<3)?false:true;
    }
 
    unsigned IndexaOffset = KKoffset[0][processor_id()];
    unsigned IndexaSize=KKoffset[KKIndex.size()-1][processor_id()] - KKoffset[0][processor_id()];
    vector < unsigned > indexa(IndexaSize,IndexaSize);
    _Psize.resize(3);
  
    unsigned IndexbOffset = IndexaOffset;
    unsigned IndexbSize=IndexaSize;
    vector <PetscInt> indexbi(IndexbSize);
    vector < unsigned > indexb(IndexbSize,IndexbSize);
  
    unsigned IndexdOffset   =_msh->MetisOffset[3][processor_id()];
    unsigned IndexdOffsetp1 =_msh->MetisOffset[3][processor_id()+1];
    unsigned IndexdSize= IndexdOffsetp1 - IndexdOffset;
    vector <PetscInt> indexdi(IndexdSize);
    vector < unsigned > indexd(IndexdSize,IndexdSize);
    for(unsigned i=0;i<indexd.size();i++) indexd[i] = IndexdSize;
  
    unsigned IndexcOffset   = _msh->MetisOffset[3][processor_id()];
    unsigned IndexcOffsetp1 = _msh->MetisOffset[3][processor_id()+1];
    unsigned IndexcSize= IndexcOffsetp1 - IndexcOffset;
    vector <PetscInt> indexci(IndexcSize);
    vector < unsigned > indexc(IndexcSize,IndexcSize);
  
   
    // *** Start Vanka Block ***
    bool test_end=0;
    int vanka_block_index=0;
    while (test_end==0){
      _indexai.resize(vanka_block_index+1);
      _indexai[vanka_block_index].resize(IndexaSize);
      _Psize[0].resize(vanka_block_index+1);
      _Psize[1].resize(vanka_block_index+1);
      _Psize[2].resize(vanka_block_index+1);
            
      PetscInt Asize=0;
      PetscInt counterb=0;
      PetscInt Csize=0;
      PetscInt Dsize=0;
      PetscInt PDsize=0;
    
      test_end=1;
      for(int isdom=0; isdom<_msh->nsubdom; isdom++) {  
	int gel=_msh->IS_Mts2Gmt_elem_offset[isdom] + vanka_block_index*_block_element_number;
	// ***************** NODE/ELEMENT SERCH *******************
      
	for (int iel_mts=gel; iel_mts<gel+_block_element_number && iel_mts< _msh->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
	  unsigned iel = _msh->IS_Mts2Gmt_elem[iel_mts];     
    
	  for (unsigned i=0; i<_msh->el->GetElementDofNumber(iel,0); i++) {
	    unsigned inode=_msh->el->GetElementVertexIndex(iel,i)-1u;
	    unsigned nvei=_msh->el->GetVertexElementNumber(inode);
	    const unsigned *pt_jel=_msh->el->GetVertexElementAddress(inode,0);
	    for (unsigned j=0; j<nvei*(!FastVankaBlock)+FastVankaBlock; j++) {
	      unsigned jel=(!FastVankaBlock)?*(pt_jel++)-1u:iel;
	      //add elements for velocity to be solved

	      unsigned jel_Metis = _msh->GetMetisDof(jel,3);
	     
	      if(jel_Metis >= IndexcOffsetp1 || jel_Metis < IndexcOffset ||
		 indexc[jel_Metis-IndexcOffset] == IndexcSize){
		if(jel_Metis < IndexcOffsetp1 && jel_Metis >= IndexcOffset){
		  indexci[Csize]=jel_Metis-IndexcOffset;
		  indexc[jel_Metis-IndexcOffset]=Csize++;
		}
		//add non-schur node to be solved
		for (unsigned iind=0; iind<VankaIndex.size()-_NSchurVar; iind++) {
		  unsigned indexSol=VankaIndex[iind];
		  unsigned SolPdeIndex = _SolPdeIndex[indexSol];
		  unsigned SolType = _SolType[SolPdeIndex];
		  const unsigned *pt_un=_msh->el->GetElementVertexAddress(jel,0);
		  unsigned nvej=_msh->el->GetElementDofNumber(jel,SolType);
		  for (unsigned jj=0; jj<nvej; jj++) {
		    unsigned jnode=(SolType<3)?(*(pt_un++)-1u):(jel+jj*nel);

		    unsigned jnode_Metis = _msh->GetMetisDof(jnode,SolType);
		    if(jnode_Metis >= _msh->MetisOffset[SolType][processor_id()] &&
		       jnode_Metis <  _msh->MetisOffset[SolType][processor_id()+1]){
		      unsigned kkdof=GetKKDof(SolPdeIndex, indexSol, jnode);
		      if (indexa[kkdof- IndexaOffset]==IndexaSize && 1.1 <(*(*_Bdc)[SolPdeIndex])(jnode_Metis) ) {
			_indexai[vanka_block_index][Asize]=kkdof;
			indexa[kkdof-IndexaOffset]=Asize++;
		      }
		    }
		  }
		}
		for (unsigned jj=0; jj<_msh->el->GetElementDofNumber(jel,0); jj++) {
		  unsigned jnode=_msh->el->GetElementVertexIndex(jel,jj)-1u;
		  unsigned nvej=_msh->el->GetVertexElementNumber(jnode);
		  const unsigned *pt_kel=_msh->el->GetVertexElementAddress(jnode,0);
		  for (unsigned k=0; k<nvej; k++) {
		    unsigned kel=*(pt_kel++)-1u;
		    //add all variables to be updated
		    unsigned kel_Metis = _msh->GetMetisDof(kel,3);
		    if(kel_Metis >= IndexdOffsetp1 || 
		       (kel_Metis >= IndexdOffset && indexd[kel_Metis-IndexdOffset] == IndexdSize)){
		    
		      if(kel_Metis < IndexdOffsetp1){
			indexdi[Dsize]=kel_Metis-IndexdOffset;
			indexd[kel_Metis-IndexdOffset]=Dsize++;
		      }
		  		    
		      for (unsigned int indexSol=0; indexSol<KKIndex.size()-1u; indexSol++) {
			const unsigned *pt_un=_msh->el->GetElementVertexAddress(kel,0);
			unsigned SolPdeIndex = _SolPdeIndex[indexSol];
			unsigned SolType = _SolType[SolPdeIndex];
			unsigned nvek=_msh->el->GetElementDofNumber(kel,SolType);
			for (unsigned kk=0; kk<nvek; kk++) {
			  unsigned knode=(SolType<3)?(*(pt_un++)-1u):(kel+kk*nel);
		
			  unsigned knode_Metis = _msh->GetMetisDof(knode,SolType);
			  if(knode_Metis >= _msh->MetisOffset[SolType][processor_id()] &&
			     knode_Metis <  _msh->MetisOffset[SolType][processor_id()+1]){
			    unsigned kkdof=GetKKDof(SolPdeIndex, indexSol, knode);
			    if(indexb[kkdof- IndexbOffset]==IndexbSize && 0.1<(*(*_Bdc)[SolPdeIndex])(knode_Metis)) {
			      indexbi[counterb]=kkdof;
			      indexb[kkdof-IndexbOffset]=counterb++;
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  //Add Schur nodes (generally pressure) to be solved
	  //if(iel_mts >= _msh->IS_Mts2Gmt_elem_offset[processor_id()] && iel_mts < _msh->IS_Mts2Gmt_elem_offset[processor_id()+1])
	  {
	    for (unsigned iind=VankaIndex.size()-_NSchurVar; iind<VankaIndex.size(); iind++) {
	      unsigned indexSol=VankaIndex[iind];
	      unsigned SolPdeIndex = _SolPdeIndex[indexSol];
	      unsigned SolType = _SolType[SolPdeIndex];
	      const unsigned *pt_un=_msh->el->GetElementVertexAddress(iel,0);
	      unsigned nvei=_msh->el->GetElementDofNumber(iel,SolType);
	      for (unsigned ii=0; ii<nvei; ii++) {
		unsigned inode=(SolType<3)?(*(pt_un++)-1u):(iel+ii*nel);
		unsigned inode_Metis = _msh->GetMetisDof(inode,SolType);
		if(inode_Metis >= _msh->MetisOffset[SolType][processor_id()] &&
		   inode_Metis <  _msh->MetisOffset[SolType][processor_id()+1]){
		  unsigned kkdof=GetKKDof(SolPdeIndex, indexSol, inode);
		  if (indexa[kkdof- IndexaOffset]==IndexaSize && 1.1<(*(*_Bdc)[SolPdeIndex])(inode_Metis) ) {
		    _indexai[vanka_block_index][Asize]=kkdof;
		    indexa[kkdof - IndexaOffset]=Asize++;
		    PDsize++;
		  }
		}
	      }
	    }
	  }
	}
	if(gel+_block_element_number <_msh->IS_Mts2Gmt_elem_offset[isdom+1] ) test_end=0;     
      }
   
    
      PetscInt PBsize=Asize;
      PetscInt PCsize=PBsize-PDsize;
      for (PetscInt i=0; i<counterb; i++) {
	unsigned jnode=indexbi[i];
	if(indexa[jnode- IndexaOffset]==IndexaSize) {
	  _indexai[vanka_block_index][Asize]=jnode;
	  indexa[jnode-IndexaOffset]=Asize++;
	}
	// *** reinitialize indexb
	indexb[jnode- IndexbOffset]=IndexbSize;
      }
      PetscInt PAmBsize=Asize-PBsize;
  
      // *** re-initialize indeces(a,c,d)
      for (PetscInt i=0; i<Asize; i++) {
	indexa[_indexai[vanka_block_index][i]-IndexaOffset]=IndexaSize;
      }
      for (PetscInt i=0; i<Csize; i++) {
	indexc[indexci[i]]=IndexcSize;
      }
      for (PetscInt i=0; i<Dsize; i++) {
	indexd[indexdi[i]]=IndexdSize;
      }
    
      _indexai[vanka_block_index].resize(Asize);
    
    
      std::sort(_indexai[vanka_block_index].begin(), _indexai[vanka_block_index].begin()+PBsize);
      std::sort(_indexai[vanka_block_index].begin() + PBsize, _indexai[vanka_block_index].end());
   
       
      _Psize[0][vanka_block_index]=PBsize;
      _Psize[1][vanka_block_index]=PCsize;
      _Psize[2][vanka_block_index]=PDsize;
      
      vanka_block_index++;
    }
    clock_t end_time=clock();
    SearchTime+=(end_time-start_time);
  
  
    //BEGIN Generate std::vector<IS> for vanka solve ***********
    _isA.resize(_indexai.size());
    _isB.resize(_indexai.size());
    for(unsigned vanka_block_index=0;vanka_block_index<_indexai.size();vanka_block_index++){  
      unsigned PBsize = _Psize[0][vanka_block_index];
      unsigned PAmBsize = _indexai[vanka_block_index].size()-PBsize;
      PetscErrorCode ierr;     
      ierr=ISCreateGeneral(MPI_COMM_WORLD,PBsize,&_indexai[vanka_block_index][0],PETSC_USE_POINTER,&_isA[vanka_block_index]);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr=ISCreateGeneral(MPI_COMM_WORLD,PAmBsize,&_indexai[vanka_block_index][PBsize],PETSC_USE_POINTER,&_isB[vanka_block_index]);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
    //END Generate std::vector<IS> for vanka solve ***********  
    return SearchTime;
  }


  // ========================================================

  void VankaPetscLinearEquationSolver::solve(const vector <unsigned> &VankaIndex,
								const bool &ksp_clean) {
  
    PetscVector* EPSp=static_cast<PetscVector*> (_EPS);  //TODO
    Vec EPS=EPSp->vec(); //TODO
    PetscVector* RESp=static_cast<PetscVector*> (_RES);  //TODO
    Vec RES=RESp->vec(); //TODO
    PetscMatrix* KKp=static_cast<PetscMatrix*>(_KK); //TODO
    Mat KK=KKp->mat(); //TODO
  
    PetscErrorCode ierr;
    clock_t SearchTime=0, AssemblyTime=0, SolveTime=0, UpdateTime=0;
  
 
    // ***************** NODE/ELEMENT SEARCH *******************
    if(_indexai_init==0) {
      SearchTime += BuildIndex(VankaIndex); 
    }
    // ***************** END NODE/ELEMENT SEARCH *******************  
    
    // ***************** INIT *****************
    if(ksp_clean){	
      this->clear();
      _ksp.resize(_indexai.size());
      _pc.resize(_indexai.size());
      _A.resize(_indexai.size());
      _B.resize(_indexai.size());
      
      _w.resize(_indexai.size());
      _r.resize(_indexai.size());
      _scatA.resize(_indexai.size());
            
      _s.resize(_indexai.size());
      _scatB.resize(_indexai.size());
    
      for(unsigned vb_i=0;vb_i<_indexai.size();vb_i++){  
        ierr = MatGetSubMatrix(KK,_isA[vb_i],_isA[vb_i],MAT_INITIAL_MATRIX,&_A[vb_i]);	CHKERRABORT(MPI_COMM_WORLD,ierr);
        ierr = MatGetSubMatrix(KK,_isB[vb_i],_isA[vb_i],MAT_INITIAL_MATRIX,&_B[vb_i]);	CHKERRABORT(MPI_COMM_WORLD,ierr);
        //init ksp object
	
        this->init(_A[vb_i],_A[vb_i],_ksp[vb_i],_pc[vb_i]);
	
	ierr = MatGetVecs(_A[vb_i],&_w[vb_i],&_r[vb_i]);			CHKERRABORT(MPI_COMM_WORLD,ierr);
	ierr = VecScatterCreate(_r[vb_i],NULL,RES,_isA[vb_i],&_scatA[vb_i]);	CHKERRABORT(MPI_COMM_WORLD,ierr);
	
	ierr = MatGetVecs(_B[vb_i],NULL,&_s[vb_i]);				CHKERRABORT(MPI_COMM_WORLD,ierr);
	ierr = VecScatterCreate(_s[vb_i],NULL,RES,_isB[vb_i],&_scatB[vb_i]);	CHKERRABORT(MPI_COMM_WORLD,ierr);
		
      }
      this->_is_initialized = true;
    }
    // ***************** END INIT *****************
    
    for(unsigned vb_i=0;vb_i<_indexai.size();vb_i++){  
      // ***************** ASSEMBLY ******************
      clock_t start_time = clock();
      
      // Get block residual vector and copy it in r[vb_i]
      Vec res;
      ierr = VecGetSubVector(RES,_isA[vb_i],&res); 	CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecCopy(res,_r[vb_i]); 			CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecDestroy(&res); 				CHKERRABORT(MPI_COMM_WORLD,ierr);
 
      AssemblyTime+=( clock() - start_time);
      // ***************** END ASSEMBLY ******************
      // ***************** SOLVE ******************
      start_time=clock(); 
      
      // Solve
      ierr = KSPSolve(_ksp[vb_i],_r[vb_i],_w[vb_i]);  	CHKERRABORT(MPI_COMM_WORLD,ierr);
      
      SolveTime += (clock() - start_time);
      // ***************** END SOLVE ******************
      
      // ***************** UPDATING ******************
      start_time=clock();
      
      // update solution
      ierr = VecScatterBegin(_scatA[vb_i],_w[vb_i],EPS,ADD_VALUES, SCATTER_FORWARD);	CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecScatterEnd(_scatA[vb_i],_w[vb_i],EPS,ADD_VALUES, SCATTER_FORWARD);	CHKERRABORT(MPI_COMM_WORLD,ierr);
            
      // update Residual for the A bock
      ierr = MatMult(_A[vb_i],_w[vb_i],_r[vb_i]);					CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecScale(_r[vb_i], -1.); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecScatterBegin(_scatA[vb_i],_r[vb_i],RES,ADD_VALUES, SCATTER_FORWARD);	CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecScatterEnd(_scatA[vb_i],_r[vb_i],RES,ADD_VALUES, SCATTER_FORWARD);	CHKERRABORT(MPI_COMM_WORLD,ierr);
      
      // update Residual for the B bock
      ierr = MatMult(_B[vb_i],_w[vb_i],_s[vb_i]);					CHKERRABORT(MPI_COMM_WORLD,ierr);							      
      ierr = VecScale (_s[vb_i], -1.); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecScatterBegin(_scatB[vb_i],_s[vb_i],RES,ADD_VALUES, SCATTER_FORWARD);	CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecScatterEnd(_scatB[vb_i],_s[vb_i],RES,ADD_VALUES, SCATTER_FORWARD);	CHKERRABORT(MPI_COMM_WORLD,ierr);
         
      UpdateTime+=(clock()-start_time);
      // ***************** END UPDATING *******************      
    } //end loop over subdomain
    
#ifndef NDEBUG  
    // *** Computational info ***
    cout << "VANKA Grid: "<<_msh->GetGridNumber()<< "      SOLVER TIME:        " << std::setw(11) << std::setprecision(6) << std::fixed <<
    static_cast<double>(SearchTime + AssemblyTime + SolveTime + UpdateTime)/ CLOCKS_PER_SEC<<
      "  ITS: " << _maxits << endl;
#endif

  }

  // ========================================================

  void VankaPetscLinearEquationSolver::clear() {
  
    if (this->initialized()) {
      this->_is_initialized = false;
      int ierr=0;
      for(unsigned i=0;i<_ksp.size();i++){
	ierr = MatDestroy(&_A[i]); 		CHKERRABORT(MPI_COMM_WORLD,ierr);
	ierr = MatDestroy(&_B[i]); 		CHKERRABORT(MPI_COMM_WORLD,ierr);
	ierr = KSPDestroy(&_ksp[i]);		CHKERRABORT(MPI_COMM_WORLD,ierr);
	
	ierr = VecDestroy(&_w[i]);		CHKERRABORT(MPI_COMM_WORLD,ierr);
	ierr = VecDestroy(&_r[i]);		CHKERRABORT(MPI_COMM_WORLD,ierr);
	ierr = VecScatterDestroy(&_scatA[i]);	CHKERRABORT(MPI_COMM_WORLD,ierr);
		
	ierr = VecDestroy(&_s[i]);		CHKERRABORT(MPI_COMM_WORLD,ierr);
	ierr = VecScatterDestroy(&_scatB[i]);	CHKERRABORT(MPI_COMM_WORLD,ierr);
      }
    }
  }

  // ========================================================
  
  void VankaPetscLinearEquationSolver::init(Mat& Amat, Mat& Pmat,  KSP &ksp, PC &pc){
      
      int ierr;
      // Create the linear solver context
      ierr = KSPCreate(MPI_COMM_WORLD, &ksp);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      // Create the preconditioner context
      ierr = KSPGetPC(ksp, &pc);			CHKERRABORT(MPI_COMM_WORLD,ierr);
  
      // Set user-specified  solver and preconditioner types
      this->set_petsc_solver_type(ksp);
        
      //ierr = KSPSetOperators(ksp, Amat, Pmat, SAME_PRECONDITIONER);	
      ierr = KSPSetOperators(ksp, Amat, Pmat); //PETSC3p5	
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      // Set the tolerances for the iterative solver.  Use the user-supplied
      // tolerance for the relative residual & leave the others at default values.
      ierr = KSPSetTolerances(ksp,_rtol,_abstol,_dtol,_maxits);	CHKERRABORT(MPI_COMM_WORLD,ierr);
    
      if(_msh->GetGridNumber()!=0)
	KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);

      if(_msh->GetGridNumber()!=0)
	KSPSetNormType(ksp,KSP_NORM_NONE);
     
      // Set the options from user-input
      // Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  KSPSetFromOptions() is called _after_ any other customization  routines.
      ierr = KSPSetFromOptions(ksp);						CHKERRABORT(MPI_COMM_WORLD,ierr);
   
      // Notify PETSc of location to store residual history.
      // This needs to be called before any solves, since
      // it sets the residual history length to zero.  The default
      // behavior is for PETSc to allocate (internally) an array
      // of size 1000 to hold the residual norm history.
//       ierr = KSPSetResidualHistory(ksp,
// 				   PETSC_NULL,   // pointer to the array which holds the history
// 				   PETSC_DECIDE, // size of the array holding the history
// 				   PETSC_TRUE);  // Whether or not to reset the history for each solve.
//       CHKERRABORT(MPI_COMM_WORLD,ierr);

      PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type,pc);
      PetscReal zero = 1.e-16;
      PCFactorSetZeroPivot(pc,zero);
      PCFactorSetShiftType(pc,MAT_SHIFT_NONZERO);
  
  }

  // =================================================

  void VankaPetscLinearEquationSolver::set_petsc_solver_type(KSP &ksp) {
    int ierr = 0;
    switch (this->_solver_type) {
    case CG:
      ierr = KSPSetType(ksp, (char*) KSPCG);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case CR:
      ierr = KSPSetType(ksp, (char*) KSPCR);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case CGS:
      ierr = KSPSetType(ksp, (char*) KSPCGS);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case BICG:
      ierr = KSPSetType(ksp, (char*) KSPBICG);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case TCQMR:
      ierr = KSPSetType(ksp, (char*) KSPTCQMR);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case TFQMR:
      ierr = KSPSetType(ksp, (char*) KSPTFQMR);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case LSQR:
      ierr = KSPSetType(ksp, (char*) KSPLSQR);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case BICGSTAB:
      ierr = KSPSetType(ksp, (char*) KSPBCGS);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case MINRES:
      ierr = KSPSetType(ksp, (char*) KSPMINRES);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case GMRES:
      ierr = KSPSetType(ksp, (char*) KSPGMRES);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case RICHARDSON:
      ierr = KSPSetType(ksp, (char*) KSPRICHARDSON);	CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case CHEBYSHEV:
      ierr = KSPSetType(ksp, (char*) KSPCHEBYSHEV);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    case PREONLY:
      ierr = KSPSetType(ksp, (char*) KSPPREONLY);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      return;
    default:
      std::cerr << "ERROR:  Unsupported PETSC Solver: "
		<< this->_solver_type               << std::endl
		<< "Continuing with PETSC defaults" << std::endl;
    }
  }


  // =================================================

 


} //end namespace femus

#endif 

// // ========================================================
// void VankaPetscLinearEquationSolver::init(Mat& matrix, const bool pc_flag, const bool Schur) {
//   // Initialize the data structures if not done so already.
//   //I have to initialize many times the data structure therefore the next line must be commented
//   //   if (!this->initialized()) {
//   this->_is_initialized = true;
//   int ierr=0;
//   int grid=_msh->GetGridNumber();
//   
//   // Since the matrix is sequential for us, we use as communicator MPI_COMM_WORLD (serial) 
//   // instead of PETSC_COMM_WORLD (parallel)
//   ierr = KSPCreate(MPI_COMM_WORLD,&_ksp[0]);					CHKERRABORT(MPI_COMM_WORLD,ierr);
//   ierr = KSPSetOperators(_ksp[0],matrix,matrix,SAME_NONZERO_PATTERN);	CHKERRABORT(MPI_COMM_WORLD,ierr);
//   ierr = KSPSetInitialGuessNonzero(_ksp[0], PETSC_TRUE);			CHKERRABORT(MPI_COMM_WORLD,ierr);
//   
//   if (!Schur || grid==0) {
//     ierr = KSPSetTolerances(_ksp[0],_rtol[0],_abstol[0],_dtol[0],_maxits[0]);	CHKERRABORT(MPI_COMM_WORLD,ierr);
//   } else {
//     ierr = KSPSetTolerances(_ksp[0],_rtol[1],_abstol[1],_dtol[1],_maxits[1]);	CHKERRABORT(MPI_COMM_WORLD,ierr);
//   }
//   ierr = KSPSetFromOptions(_ksp[0]);						CHKERRABORT(MPI_COMM_WORLD,ierr);
//   
//   // Set user-specified  solver and preconditioner types
//   this->set_petsc_solver_type();
//   
//   if (!pc_flag) {
//     ierr = KSPGetPC(_ksp[0],&_pc);						CHKERRABORT(MPI_COMM_WORLD,ierr);  
//     PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);
//     PetscReal zero = 1.e-16;
//     ierr = PCFactorSetZeroPivot(_pc,zero); 					CHKERRABORT(MPI_COMM_WORLD,ierr); 
//     ierr = PCFactorSetShiftType(_pc,MAT_SHIFT_NONZERO); 			CHKERRABORT(MPI_COMM_WORLD,ierr); 
//   }
// }
// 
// // ========================================================
// void VankaPetscLinearEquationSolver::init_schur(Mat& matrix) {
//   // Initialize the data structures if not done so already.
//   //I have to initialize many times the data structure therefore the next line must be commented
//   //   if (!this->initialized()) {
//   this->_is_initialized = true;
//   int ierr=0;
//   
//   ierr = MatSchurComplementGetKSP(matrix,&_ksp[1]);				CHKERRABORT(MPI_COMM_WORLD,ierr);
//   ierr = KSPSetInitialGuessNonzero(_ksp[1], PETSC_TRUE);			CHKERRABORT(MPI_COMM_WORLD,ierr);
//   
//   
//   ierr = KSPSetTolerances(_ksp[1],_rtol[0],_abstol[0],_dtol[0],_maxits[0]);	CHKERRABORT(MPI_COMM_WORLD,ierr);
//   ierr = KSPSetFromOptions(_ksp[1]);						CHKERRABORT(MPI_COMM_WORLD,ierr);
//   
//   // Set user-specified  solver and preconditioner types
//   this->set_petsc_solver_type2();
//   
//   ierr = KSPGetPC(_ksp[1], &_pc);						CHKERRABORT(MPI_COMM_WORLD,ierr);
//   PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);
//   PetscReal zero = 1.e-16;
//   PCFactorSetZeroPivot(_pc,zero);
//   PCFactorSetShiftType(_pc,MAT_SHIFT_NONZERO);
//   
//   //   }
// }


// void VankaPetscLinearEquationSolver::set_petsc_solver_type2() {
//   int ierr = 0;
//   switch (this->_solver_type) {
//     case CG:
//       ierr = KSPSetType(_ksp[1], (char*) KSPCG);		CHKERRABORT(MPI_COMM_WORLD,ierr);
//       return;
//     case CR:
//       ierr = KSPSetType(_ksp[1], (char*) KSPCR);		CHKERRABORT(MPI_COMM_WORLD,ierr);
//       return;
//     case CGS:
//       ierr = KSPSetType(_ksp[1], (char*) KSPCGS);		CHKERRABORT(MPI_COMM_WORLD,ierr);
//       return;
//     case BICG:
//       ierr = KSPSetType(_ksp[1], (char*) KSPBICG);		CHKERRABORT(MPI_COMM_WORLD,ierr);
//       return;
//     case TCQMR:
//       ierr = KSPSetType(_ksp[1], (char*) KSPTCQMR);		CHKERRABORT(MPI_COMM_WORLD,ierr);
//       return;
//     case TFQMR:
//       ierr = KSPSetType(_ksp[1], (char*) KSPTFQMR);		CHKERRABORT(MPI_COMM_WORLD,ierr);
//       return;
//     case LSQR:
//       ierr = KSPSetType(_ksp[1], (char*) KSPLSQR);		CHKERRABORT(MPI_COMM_WORLD,ierr);
//       return;
//     case BICGSTAB:
//       ierr = KSPSetType(_ksp[1], (char*) KSPBCGS);		CHKERRABORT(MPI_COMM_WORLD,ierr);
//       return;
//     case MINRES:
//       ierr = KSPSetType(_ksp[1], (char*) KSPMINRES);		CHKERRABORT(MPI_COMM_WORLD,ierr);
//       return;
//     case GMRES:
//       ierr = KSPSetType(_ksp[1], (char*) KSPGMRES);		CHKERRABORT(MPI_COMM_WORLD,ierr);
//       return;
//     case RICHARDSON:
//       ierr = KSPSetType(_ksp[1], (char*) KSPRICHARDSON);	CHKERRABORT(MPI_COMM_WORLD,ierr);
//       return;
//     case CHEBYSHEV:
//       ierr = KSPSetType(_ksp[1], (char*) KSPCHEBYSHEV);		CHKERRABORT(MPI_COMM_WORLD,ierr);
//       return;
//     case PREONLY:
//       ierr = KSPSetType(_ksp[1], (char*) KSPPREONLY);		CHKERRABORT(MPI_COMM_WORLD,ierr);
//       return;
//     default:
//       std::cerr << "ERROR:  Unsupported PETSC Solver: "
//       << this->_solver_type               << std::endl
//       << "Continuing with PETSC defaults" << std::endl;
//   }
// }






//     else { // *********** if SCHUR COMPLEMENT ****************
//       
//             unsigned PCsize = _Psize[1][vanka_block_index];
//             unsigned PDsize = _Psize[2][vanka_block_index];
//             
//             ierr = ISCreateStride(MPI_COMM_WORLD,PBsize,0,1,&isPB); 			CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = ISGetIndices(isPB,&ind2); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
//             
//             IS isCl,isDl;
//             ierr=ISCreateGeneral(MPI_COMM_WORLD,PCsize,ind2,PETSC_USE_POINTER,&isCl); 		CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr=ISCreateGeneral(MPI_COMM_WORLD,PDsize,&ind2[PCsize],PETSC_USE_POINTER,&isDl); 	CHKERRABORT(MPI_COMM_WORLD,ierr);
//       
//             Vec f;
//             ierr = VecGetSubVector(res,isCl,&f); 					CHKERRABORT(MPI_COMM_WORLD,ierr);
//             Vec w;
//             ierr = VecDuplicate(f,&w); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = VecSet(w,0); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//            
//             Vec g;
//             ierr = VecGetSubVector(res,isDl,&g); 					CHKERRABORT(MPI_COMM_WORLD,ierr);
//             Vec z;
//             ierr = VecDuplicate(g,&z); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = VecSet(z,0); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//       	
//             /* linear system matrix Schur Complement*/
//             Mat A;
//             ierr = MatGetSubMatrix(PA,isCl,isCl,MAT_INITIAL_MATRIX,&A);		CHKERRABORT(MPI_COMM_WORLD,ierr);
//             Mat B;
//             ierr = MatGetSubMatrix(PA,isDl,isCl,MAT_INITIAL_MATRIX,&B);		CHKERRABORT(MPI_COMM_WORLD,ierr);
//             Mat Bt;
//             if (GetMatrixProperties()) {
//       	ierr = MatCreateTranspose(B,&Bt); 					CHKERRABORT(MPI_COMM_WORLD,ierr);
//             } 
//             else {
//       	// In the Monolithic Fluid-Structure formulation Bt is not the transpose of B
//       	ierr = MatGetSubMatrix(PA,isCl,isDl,MAT_INITIAL_MATRIX,&Bt);		CHKERRABORT(MPI_COMM_WORLD,ierr);
//             }
//              
//             ierr = KSPMatRegisterAll();						CHKERRABORT(MPI_COMM_WORLD,ierr);
//             Mat C,D;
//             //C=D-B A^{-1} B^t
//             if (!GetStabilization()) {
//       	ierr = MatCreateSchurComplement(A,A,Bt,B,PETSC_NULL,&C); 		CHKERRABORT(MPI_COMM_WORLD,ierr);
//             } 
//             else {
//       	ierr = MatGetSubMatrix(PA,isDl,isDl,MAT_INITIAL_MATRIX,&D); 		CHKERRABORT(MPI_COMM_WORLD,ierr);
//       	ierr = MatCreateSchurComplement(A,A,Bt,B,D,&C); 			CHKERRABORT(MPI_COMM_WORLD,ierr);
//             }
//       
//             ierr = ISDestroy(&isCl); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = ISDestroy(&isDl); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//       	
//             clock_t end_time=clock();
//             AssemblyTime+=(end_time-start_time);
//             // ***************** END ASSEMBLY ******************
//             // ***************** START SCHUR COMPLEMENT SOLVER ******************
//       	
//             start_time=clock(); // START SOLVE 0 TIME
//             //init Schur Complement
//             this->init_schur(C);
//             //init _Ksp (GMRES - no precond)
//             this->init(C,true,Schur);
//             // Prediction step: solve A^{-1} f -> w (temporary)
//             ierr = KSPSolve(_ksp[1],f,w);  						CHKERRABORT(MPI_COMM_WORLD,ierr);   
//             ierr = KSPGetIterationNumber(_ksp[1],&its); 				CHKERRABORT(MPI_COMM_WORLD,ierr);
//             its_A += its;
//             end_time=clock();  // END SOLVE 0 TIME
//             SolveTime0+=(end_time-start_time);
//       
//             start_time=clock(); // START SOLVE 1 TIME
//             // Projection Step
//             // Compute: g - B A^{-1} f -> g (temporary)
//             ierr = VecScale (w, -1.); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = MatMultAdd(B,w,g,g); 						CHKERRABORT(MPI_COMM_WORLD,ierr); 
//             
//             // Solve z = C^{-1} (g - B A^{-1} f)
//             ierr = KSPSolve(_ksp[0],g,z); 						CHKERRABORT(MPI_COMM_WORLD,ierr); 
//             ierr = KSPGetIterationNumber(_ksp[0],&its); 				CHKERRABORT(MPI_COMM_WORLD,ierr);
//             its_C += its;
//             end_time=clock(); // END SOLVE 1 TIME
//             SolveTime1+=(end_time-start_time);
//       
//             start_time=clock(); // START SOLVE 2 TIME
//       	
//             // Correction Step
//             // Compute: f - Bt z -> f
//             ierr =VecCopy(z,g); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = VecScale(g, -1.); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = MatMultAdd(Bt,g,f,f); 						CHKERRABORT(MPI_COMM_WORLD,ierr);  
//       
//             // solve w=A^{-1}(f-Bt z)
//             ierr = KSPSolve(_ksp[1],f,w); 						CHKERRABORT(MPI_COMM_WORLD,ierr); 
//       
//             PetscScalar *Z[1];
//             ierr = VecGetArray(w,Z); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = VecSetValues(Pw,PCsize,ind2,Z[0],INSERT_VALUES); 			CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = VecRestoreArray(w,Z); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
//       
//             ierr = VecGetArray(z,Z); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = VecSetValues(Pw,PDsize,&ind2[PCsize],Z[0],INSERT_VALUES); 		CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = VecRestoreArray(z,Z); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
//       
//             ierr = VecAssemblyBegin(Pw); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = VecAssemblyEnd(Pw); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
//             
//             ierr=ISRestoreIndices(isPB,&ind2); 					CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = ISDestroy(&isPB); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//       
//           
//             
//             ierr = MatDestroy(&A); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = MatDestroy(&B); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = MatDestroy(&Bt); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = MatDestroy(&C); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//             if (GetStabilization()) {
//       	ierr = MatDestroy(&D); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//             }
//             ierr = VecDestroy(&f); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = VecDestroy(&w); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = VecDestroy(&g); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
//             ierr = VecDestroy(&z); 							CHKERRABORT(MPI_COMM_WORLD,ierr);   
//             
//             end_time=clock(); // END SOLVE 2 TIME
//             SolveTime2+=(end_time-start_time);  
// 							      }
// ***************** END VANKA SCHUR COMPLEMENT SOLVER ******************
