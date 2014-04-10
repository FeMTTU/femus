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
#include "PetscLinearEquationSolver.hpp"
#include "PetscPreconditioner.hpp"
#include "PetscVector.hpp"
#include "PetscMatrix.hpp"

using namespace std;

// ==========================================================
extern "C" {
#if PETSC_VERSION_LESS_THAN(2,2,1)
  typedef int PetscErrorCode;
  typedef int PetscInt;
#endif

#if PETSC_VERSION_LESS_THAN(3,0,1) && PETSC_VERSION_RELEASE
  // ------------------------------------------------------------------
  PetscErrorCode __libmesh_petsc_preconditioner_setup(void * ctx) {
    Preconditioner * preconditioner = static_cast<Preconditioner*>(ctx);
    preconditioner->init();
    return 0;
  }
  // ------------------------------------------------------------------------------
  PetscErrorCode __libmesh_petsc_preconditioner_apply(void *ctx, Vec x, Vec y)  {
    Preconditioner * preconditioner = static_cast<Preconditioner*>(ctx);
    PetscVector x_vec(x);
    PetscVector y_vec(y);
    preconditioner->apply(x_vec,y_vec);
    return 0;
  }
#else
  // ----------------------------------------------------------------
  PetscErrorCode __libmesh_petsc_preconditioner_setup(PC pc) {
    void *ctx;
    PetscErrorCode ierr = PCShellGetContext(pc,&ctx);
    CHKERRQ(ierr);
    Preconditioner * preconditioner = static_cast<Preconditioner*>(ctx);
    preconditioner->init();
    return 0;
  }
  // --------------------------------------------------------------------------
  PetscErrorCode __libmesh_petsc_preconditioner_apply(PC pc, Vec x, Vec y) {
    void *ctx;
    PetscErrorCode ierr = PCShellGetContext(pc,&ctx);
    CHKERRQ(ierr);
    Preconditioner * preconditioner = static_cast<Preconditioner*>(ctx);
    PetscVector x_vec(x);
    PetscVector y_vec(y);
    preconditioner->apply(x_vec,y_vec);
    return 0;
  }
#endif
} // end extern "C
// ================================================


// ==============================================
// ----------------------- functions ------
// ==============================================

void PetscLinearEquationSolver::set_tolerances(const double &rtol, const double &atol,
                                       const double &divtol, const unsigned &maxits, const unsigned &index) {

  if(index>=_rtol.size()){
    cout<<"Error the index in the PetscLinearEquationSolver::set_tolerances is out of range"<<endl;
    exit(0);
  }
  
  _rtol[index]   = static_cast<PetscReal>(rtol);
  _abstol[index] = static_cast<PetscReal>(atol);
  _dtol[index]   = static_cast<PetscReal>(divtol);
  _maxits[index] = static_cast<PetscInt>(maxits);

};

// ******************************************************************

void PetscLinearEquationSolver::set_num_elem_vanka_block(const unsigned num_elem_vanka_block) {
  _num_elem_vanka_block = num_elem_vanka_block;
  _indexai_init=0;
};

// *********************** for GMRES *********************************
clock_t PetscLinearEquationSolver::BuildIndex(){
  
  clock_t SearchTime = 0;
  clock_t start_time = clock();
  _indexai_init = 1;
   
  unsigned IndexaSize=KKoffset[KKIndex.size()-1][_msh->_iproc] - KKoffset[0][_msh->_iproc];
  _indexai.resize(2);
  _indexai[0].resize(IndexaSize);
  _indexai[1].resize(IndexaSize);
    
  unsigned count0=0;
  unsigned count1=0;
  for(int k=0; k < _SolPdeIndex.size(); k++) {
    unsigned indexSol = _SolPdeIndex[k];
    unsigned soltype = _SolType[indexSol];
    for(unsigned inode_mts = _msh->MetisOffset[soltype][_msh->_iproc]; 
	inode_mts < _msh->MetisOffset[soltype][_msh->_iproc+1]; inode_mts++) {
      int local_mts = inode_mts-_msh->MetisOffset[soltype][_msh->_iproc];
      int idof_kk = KKoffset[k][_msh->_iproc] +local_mts; 
      if((*(*_Bdc)[indexSol])(inode_mts) < 1.9) {
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
  
  
  
  //BEGIN Generate std::vector<IS> for GMRES solve by elimination ***********
  _isA.resize(1);
  PetscInt Asize=_indexai[1].size();
  int ierr=ISCreateGeneral(MPI_COMM_WORLD,Asize,&_indexai[1][0], PETSC_USE_POINTER ,&_isA[0]);	   
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  //END Generate std::vector<IS> for GMRES solve by elimination ***********  
    
  
  clock_t end_time=clock();
  SearchTime = (end_time-start_time);
  
  return SearchTime;
}

// *********************** for VANKA *********************************

clock_t PetscLinearEquationSolver::BuildIndex(const vector <unsigned> &VankaIndex,
				      const short unsigned &NSchurVar){
  clock_t SearchTime=0;
  clock_t start_time=clock();
  _indexai_init=1;
  unsigned nel=_msh->GetElementNumber();
  bool FastVankaBlock=true;
  if(NSchurVar==!0){
    FastVankaBlock=(_SolType[_SolPdeIndex[VankaIndex[VankaIndex.size()-NSchurVar]]]<3)?false:true;
  }
 
  unsigned IndexaOffset = KKoffset[0][_msh->_iproc];
  unsigned IndexaSize=KKoffset[KKIndex.size()-1][_msh->_iproc] - KKoffset[0][_msh->_iproc];
  vector < unsigned > indexa(IndexaSize,IndexaSize);
  _Psize.resize(3);
  
  unsigned IndexbOffset = IndexaOffset;
  unsigned IndexbSize=IndexaSize;
  vector <PetscInt> indexbi(IndexbSize);
  vector < unsigned > indexb(IndexbSize,IndexbSize);
  
  unsigned IndexdOffset   =_msh->MetisOffset[3][_msh->_iproc];
  unsigned IndexdOffsetp1 =_msh->MetisOffset[3][_msh->_iproc+1];
  unsigned IndexdSize= IndexdOffsetp1 - IndexdOffset;
  vector <PetscInt> indexdi(IndexdSize);
  vector < unsigned > indexd(IndexdSize,IndexdSize);
  for(unsigned i=0;i<indexd.size();i++) indexd[i] = IndexdSize;
  
  unsigned IndexcOffset   = _msh->MetisOffset[3][_msh->_iproc];
  unsigned IndexcOffsetp1 = _msh->MetisOffset[3][_msh->_iproc+1];
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
      int gel=_msh->IS_Mts2Gmt_elem_offset[isdom] + vanka_block_index*_num_elem_vanka_block;
      // ***************** NODE/ELEMENT SERCH *******************
      
      for (int iel_mts=gel; iel_mts<gel+_num_elem_vanka_block && iel_mts< _msh->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
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
	      //----------------------------------------------------------------------------------
	      //add non-schur node to be solved
	      for (unsigned iind=0; iind<VankaIndex.size()-NSchurVar; iind++) {
		unsigned indexSol=VankaIndex[iind];
		unsigned SolPdeIndex = _SolPdeIndex[indexSol];
		unsigned SolType = _SolType[SolPdeIndex];
		const unsigned *pt_un=_msh->el->GetElementVertexAddress(jel,0);
		unsigned nvej=_msh->el->GetElementDofNumber(jel,_msh->_END_IND[SolType]);
		for (unsigned jj=0; jj<nvej; jj++) {
		  unsigned jnode=(SolType<3)?(*(pt_un++)-1u):(jel+jj*nel);

		  unsigned jnode_Metis = _msh->GetMetisDof(jnode,SolType);
		  if(jnode_Metis >= _msh->MetisOffset[SolType][_msh->_iproc] &&
		     jnode_Metis <  _msh->MetisOffset[SolType][_msh->_iproc+1]){
		    unsigned kkdof=GetKKDof(SolPdeIndex, indexSol, jnode);
		    if (indexa[kkdof- IndexaOffset]==IndexaSize && 1.1 <(*(*_Bdc)[SolPdeIndex])(jnode_Metis) ) {
		      _indexai[vanka_block_index][Asize]=kkdof;
		      indexa[kkdof-IndexaOffset]=Asize++;
		    }
		  }
		}
	      }
	      //-----------------------------------------------------------------------------------
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
		      unsigned nvek=_msh->el->GetElementDofNumber(kel,_msh->_END_IND[SolType]);
		      for (unsigned kk=0; kk<nvek; kk++) {
			unsigned knode=(SolType<3)?(*(pt_un++)-1u):(kel+kk*nel);
		
			unsigned knode_Metis = _msh->GetMetisDof(knode,SolType);
			if(knode_Metis >= _msh->MetisOffset[SolType][_msh->_iproc] &&
			   knode_Metis <  _msh->MetisOffset[SolType][_msh->_iproc+1]){
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
	      //------------------------------------------------------------------------
	    }
	  }
	}
	//-----------------------------------------------------------------------------------------
	//Add Schur nodes (generally pressure) to be solved
	//if(iel_mts >= _msh->IS_Mts2Gmt_elem_offset[_msh->_iproc] && iel_mts < _msh->IS_Mts2Gmt_elem_offset[_msh->_iproc+1])
	{
	  for (unsigned iind=VankaIndex.size()-NSchurVar; iind<VankaIndex.size(); iind++) {
	    unsigned indexSol=VankaIndex[iind];
	    unsigned SolPdeIndex = _SolPdeIndex[indexSol];
	    unsigned SolType = _SolType[SolPdeIndex];
	    const unsigned *pt_un=_msh->el->GetElementVertexAddress(iel,0);
	    unsigned nvei=_msh->el->GetElementDofNumber(iel,_msh->_END_IND[SolType]);
	    for (unsigned ii=0; ii<nvei; ii++) {
	      unsigned inode=(SolType<3)?(*(pt_un++)-1u):(iel+ii*nel);
	      unsigned inode_Metis = _msh->GetMetisDof(inode,SolType);
	      if(inode_Metis >= _msh->MetisOffset[SolType][_msh->_iproc] &&
		 inode_Metis <  _msh->MetisOffset[SolType][_msh->_iproc+1]){
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
	//-----------------------------------------------------------------------------------------
      }
      if(gel+_num_elem_vanka_block <_msh->IS_Mts2Gmt_elem_offset[isdom+1] ) test_end=0;     
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
 
// ***********************************************************************

std::pair< int, double> PetscLinearEquationSolver::solve(const vector <unsigned> &VankaIndex,
						 const short unsigned &NSchurVar,const bool &Schur) {
  PetscVector* EPSp=static_cast<PetscVector*> (_EPS);  //TODO
  Vec EPS=EPSp->vec(); //TODO
  PetscVector* RESp=static_cast<PetscVector*> (_RES);  //TODO
  Vec RES=RESp->vec(); //TODO
  PetscMatrix* KKp=static_cast<PetscMatrix*>(_KK); //TODO
  Mat KK=KKp->mat(); //TODO
  
  PetscErrorCode ierr;
  clock_t SearchTime=0, AssemblyTime=0, SolveTime0=0, SolveTime1=0, SolveTime2=0, UpdateTime=0;

  int its_A=0, its_C=0, its=0;
  
  // ***************** NODE/ELEMENT SEARCH *******************
  if(_indexai_init==0) SearchTime += BuildIndex(VankaIndex,NSchurVar); 
  // ***************** END NODE/ELEMENT SEARCH *******************  
  
  for(unsigned vanka_block_index=0;vanka_block_index<_indexai.size();vanka_block_index++){  
    // ***************** ASSEMBLY ******************
    clock_t start_time = clock(); 
    unsigned PBsize = _Psize[0][vanka_block_index];
    unsigned PAmBsize = _indexai[vanka_block_index].size()-PBsize;
    // generate IS
    IS &isPA=_isA[vanka_block_index];
        
    Vec res;
    ierr = VecGetSubVector(RES,isPA,&res); 		CHKERRABORT(MPI_COMM_WORLD,ierr);
  
    Vec Pr;
    ierr = VecDuplicate(res,&Pr); 			CHKERRABORT(MPI_COMM_WORLD,ierr);
      
    // Solution Vector Pw
    Vec Pw;
    ierr = VecDuplicate(Pr,&Pw); 			CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecSet(Pw,0); 				CHKERRABORT(MPI_COMM_WORLD,ierr);
        
    //Matrix PA (PA is the Total matrix over a patch)
    Mat PA;
    ierr = MatGetSubMatrix(KK,isPA,isPA,MAT_INITIAL_MATRIX,&PA); CHKERRABORT(MPI_COMM_WORLD,ierr);
    IS isPB;  
    const PetscInt *ind,*ind2;
    if (!Schur || _msh->GetGridNumber()==0) { // ******** IF NON-SCHUR COMPLEMENT ***************
      ierr = VecCopy(res,Pr); 				CHKERRABORT(MPI_COMM_WORLD,ierr);
      clock_t end_time=clock();
      AssemblyTime+=(end_time-start_time);
      // ***************** END NON-SCHUR ASSEMBLY ******************
      // ***************** START NON-SCHUR COMPLEMENT SOLVER ******************
      start_time=clock(); 
	
      //init ksp object
      this->init(PA,false,false);
      
      // Solve
      ierr = KSPSolve(_ksp[0],Pr,Pw);  			CHKERRABORT(MPI_COMM_WORLD,ierr);
	
      ierr = KSPGetIterationNumber(_ksp[0],&its); 	CHKERRABORT(MPI_COMM_WORLD,ierr);
      its_A += its;

      end_time=clock();
      SolveTime0+=(end_time-start_time);
      // ***************** END NON-SCHUR COMPLEMENT SOLVER ******************
    }
    else { // *********** if SCHUR COMPLEMENT ****************
	
      unsigned PCsize = _Psize[1][vanka_block_index];
      unsigned PDsize = _Psize[2][vanka_block_index];
      
      ierr = ISCreateStride(MPI_COMM_WORLD,PBsize,0,1,&isPB); 			CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = ISGetIndices(isPB,&ind2); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
      
      IS isCl,isDl;
      ierr=ISCreateGeneral(MPI_COMM_WORLD,PCsize,ind2,PETSC_USE_POINTER,&isCl); 		CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr=ISCreateGeneral(MPI_COMM_WORLD,PDsize,&ind2[PCsize],PETSC_USE_POINTER,&isDl); 	CHKERRABORT(MPI_COMM_WORLD,ierr);

      Vec f;
      ierr = VecGetSubVector(res,isCl,&f); 					CHKERRABORT(MPI_COMM_WORLD,ierr);
      Vec w;
      ierr = VecDuplicate(f,&w); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSet(w,0); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
     
      Vec g;
      ierr = VecGetSubVector(res,isDl,&g); 					CHKERRABORT(MPI_COMM_WORLD,ierr);
      Vec z;
      ierr = VecDuplicate(g,&z); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSet(z,0); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
	
      /* linear system matrix Schur Complement*/
      Mat A;
      ierr = MatGetSubMatrix(PA,isCl,isCl,MAT_INITIAL_MATRIX,&A);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      Mat B;
      ierr = MatGetSubMatrix(PA,isDl,isCl,MAT_INITIAL_MATRIX,&B);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      Mat Bt;
      if (GetMatrixProperties()) {
	ierr = MatCreateTranspose(B,&Bt); 					CHKERRABORT(MPI_COMM_WORLD,ierr);
      } 
      else {
	// In the Monolithic Fluid-Structure formulation Bt is not the transpose of B
	ierr = MatGetSubMatrix(PA,isCl,isDl,MAT_INITIAL_MATRIX,&Bt);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      }
       
      ierr = KSPMatRegisterAll();						CHKERRABORT(MPI_COMM_WORLD,ierr);
      Mat C,D;
      //C=D-B A^{-1} B^t
      if (!GetStabilization()) {
	ierr = MatCreateSchurComplement(A,A,Bt,B,PETSC_NULL,&C); 		CHKERRABORT(MPI_COMM_WORLD,ierr);
      } 
      else {
	ierr = MatGetSubMatrix(PA,isDl,isDl,MAT_INITIAL_MATRIX,&D); 		CHKERRABORT(MPI_COMM_WORLD,ierr);
	ierr = MatCreateSchurComplement(A,A,Bt,B,D,&C); 			CHKERRABORT(MPI_COMM_WORLD,ierr);
      }

      ierr = ISDestroy(&isCl); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = ISDestroy(&isDl); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
	
      clock_t end_time=clock();
      AssemblyTime+=(end_time-start_time);
      // ***************** END ASSEMBLY ******************
      // ***************** START SCHUR COMPLEMENT SOLVER ******************
	
      start_time=clock(); // START SOLVE 0 TIME
      //init Schur Complement
      this->init_schur(C);
      //init _Ksp (GMRES - no precond)
      this->init(C,true,Schur);
      // Prediction step: solve A^{-1} f -> w (temporary)
      ierr = KSPSolve(_ksp[1],f,w);  						CHKERRABORT(MPI_COMM_WORLD,ierr);   
      ierr = KSPGetIterationNumber(_ksp[1],&its); 				CHKERRABORT(MPI_COMM_WORLD,ierr);
      its_A += its;
      end_time=clock();  // END SOLVE 0 TIME
      SolveTime0+=(end_time-start_time);

      start_time=clock(); // START SOLVE 1 TIME
      // Projection Step
      // Compute: g - B A^{-1} f -> g (temporary)
      ierr = VecScale (w, -1.); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = MatMultAdd(B,w,g,g); 						CHKERRABORT(MPI_COMM_WORLD,ierr); 
      
      // Solve z = C^{-1} (g - B A^{-1} f)
      ierr = KSPSolve(_ksp[0],g,z); 						CHKERRABORT(MPI_COMM_WORLD,ierr); 
      ierr = KSPGetIterationNumber(_ksp[0],&its); 				CHKERRABORT(MPI_COMM_WORLD,ierr);
      its_C += its;
      end_time=clock(); // END SOLVE 1 TIME
      SolveTime1+=(end_time-start_time);

      start_time=clock(); // START SOLVE 2 TIME
	
      // Correction Step
      // Compute: f - Bt z -> f
      ierr =VecCopy(z,g); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecScale(g, -1.); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = MatMultAdd(Bt,g,f,f); 						CHKERRABORT(MPI_COMM_WORLD,ierr);  

      // solve w=A^{-1}(f-Bt z)
      ierr = KSPSolve(_ksp[1],f,w); 						CHKERRABORT(MPI_COMM_WORLD,ierr); 

      PetscScalar *Z[1];
      ierr = VecGetArray(w,Z); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSetValues(Pw,PCsize,ind2,Z[0],INSERT_VALUES); 			CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecRestoreArray(w,Z); 						CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = VecGetArray(z,Z); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSetValues(Pw,PDsize,&ind2[PCsize],Z[0],INSERT_VALUES); 		CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecRestoreArray(z,Z); 						CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = VecAssemblyBegin(Pw); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecAssemblyEnd(Pw); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
      
      ierr=ISRestoreIndices(isPB,&ind2); 					CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = ISDestroy(&isPB); 							CHKERRABORT(MPI_COMM_WORLD,ierr);

    
      
      ierr = MatDestroy(&A); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = MatDestroy(&B); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = MatDestroy(&Bt); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = MatDestroy(&C); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
      if (GetStabilization()) {
	ierr = MatDestroy(&D); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
      }
      ierr = VecDestroy(&f); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecDestroy(&w); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecDestroy(&g); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecDestroy(&z); 							CHKERRABORT(MPI_COMM_WORLD,ierr);   
      
      end_time=clock(); // END SOLVE 2 TIME
      SolveTime2+=(end_time-start_time);  
    }
    // ***************** END VANKA SCHUR COMPLEMENT SOLVER ******************
    // ***************** START UPDATING AND CLEANING ******************
    start_time=clock();
    //update Residual for PA
    ierr = MatMult(PA,Pw,Pr); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecScale(Pr, -1.); 						CHKERRABORT(MPI_COMM_WORLD,ierr);

    PetscScalar *R[1];
    ierr = VecGetArray(Pr,R); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
    PetscScalar *W[1];
    ierr = VecGetArray(Pw,W); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = ISGetIndices(isPA,&ind); 					CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr=VecSetValues(RES,PBsize,ind,R[0],ADD_VALUES); 			CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr=VecSetValues(EPS,PBsize,ind,W[0],ADD_VALUES); 			CHKERRABORT(MPI_COMM_WORLD,ierr);
    _EPS->close();

    ierr = ISRestoreIndices(isPA,&ind); 				CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecRestoreArray(Pr,R); 					CHKERRABORT(MPI_COMM_WORLD,ierr);

    IS &isB=_isB[vanka_block_index];
    Vec Ps;
    ierr = VecCreateMPI(MPI_COMM_WORLD, PAmBsize, PETSC_DETERMINE, &Ps); 	CHKERRABORT(MPI_COMM_WORLD,ierr);
	
    ierr = VecSetFromOptions(Ps); 						CHKERRABORT(MPI_COMM_WORLD,ierr);

    Mat PB;
    ierr = MatGetSubMatrix(KK,isB,isPA,MAT_INITIAL_MATRIX,&PB); 		CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = MatMult(PB,Pw,Ps); 							CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = VecScale (Ps, -1.); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
    PetscScalar *S[1];
    ierr = VecGetArray(Ps,S); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = ISGetIndices(isB,&ind); 						CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr=VecSetValues(RES,PAmBsize,ind,S[0],ADD_VALUES); 			CHKERRABORT(MPI_COMM_WORLD,ierr);
	
    ierr = ISRestoreIndices(isB,&ind); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecRestoreArray(Ps,S); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecDestroy(&Ps); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = MatDestroy(&PB); 							CHKERRABORT(MPI_COMM_WORLD,ierr);

    _RES->close();
     
    ierr = VecRestoreArray(Pw,W); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecDestroy(&Pw); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecDestroy(&Pr); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecDestroy(&res); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = MatDestroy(&PA); 							CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = KSPDestroy(&_ksp[0]); 						CHKERRABORT(MPI_COMM_WORLD,ierr);
    clock_t end_time=clock();
    UpdateTime+=(end_time-start_time);
    // ***************** END UPDATING AND CLEANING *******************
    //vanka_block_index++;  
    //  } //end loop over vanka block
  } //loop over subdomain

  // *** Computational info ***
  cout << "Grid: "<<_msh->GetGridNumber()<< "      SOLVER TIME:              " <<
    static_cast<double>(SearchTime + AssemblyTime + SolveTime0 + SolveTime1 + SolveTime2 + UpdateTime)/ CLOCKS_PER_SEC<<
    "  ITS: " << its_A + its_C << endl;


  //    cout <<  " SearchTime: " << static_cast<double>(SearchTime)/ CLOCKS_PER_SEC << endl;
  //    cout << " AssemblyTime: " << static_cast<double>(AssemblyTime)/ CLOCKS_PER_SEC << endl;
  //    cout << " SolverTime: " <<  static_cast<double>(SolveTime0 + SolveTime1 + SolveTime2)/ CLOCKS_PER_SEC << endl;
  //    cout << " UpdateTime: " <<  static_cast<double>(UpdateTime)/ CLOCKS_PER_SEC << endl;

  return std::make_pair(its_A+its_C, 
			static_cast<double>(SearchTime + AssemblyTime 
					    + SolveTime0 + SolveTime1 
					    + SolveTime2 + UpdateTime)/ CLOCKS_PER_SEC);

}

// ********************************************************************************

std::pair< int, double> PetscLinearEquationSolver::solve() {
    
  clock_t SearchTime, AssemblyTime, SolveTime, UpdateTime;
  int its;
  double final_resid;
  PetscErrorCode ierr; 
    
  // ***************** NODE/ELEMENT SEARCH *******************
  clock_t start_time=clock();
  if(_indexai_init==0) BuildIndex();
  SearchTime = clock() - start_time;
  // ***************** END NODE/ELEMENT SEARCH *******************  
    
  if(_DirichletBCsHandlingMode==0) // By penalty
    {
      PetscVector* EPSCp=static_cast<PetscVector*> (_EPSC);  //TODO
      Vec EPSC=EPSCp->vec(); //TODO
      PetscVector* RESp=static_cast<PetscVector*> (_RES);  //TODO
      Vec RES=RESp->vec(); //TODO
      PetscMatrix* KKp=static_cast<PetscMatrix*>(_KK); //TODO
      Mat KK=KKp->mat(); //TODO
            
      
      // ***************** ASSEMBLE matrix to set Dirichlet BCs by penalty *******************
      start_time=clock();
     
      Mat Pmat;
      MatDuplicate(KK,MAT_COPY_VALUES,&Pmat);
      MatSetOption(Pmat,MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE);
      MatZeroRows(Pmat,_indexai[0].size(),&_indexai[0][0],1.e40,0,0);
          
      AssemblyTime = clock()-start_time;
            
      //vector <double>  value(_indexai[0].size());
      //_KK->matrix_get_diagonal_values(_indexai[0],value);

      //double penalty=1.0e20;
      //_KK->matrix_set_diagonal_values(_indexai[0],penalty);
      //_KK->close();
               
      //for(int i=0;i<_indexai[0].size();i++){ 
      //  int ierr = MatSetValuesBlocked(Pmat,1,&_indexai[0][i],1,&_indexai[0][i],&penalty,INSERT_VALUES);  CHKERRABORT(MPI_COMM_WORLD,ierr);
      //}
            
      // ***************** END ASSEMBLE ******************

      // ***************** SOLVE ******************
      start_time=clock();

      this->init(KK,Pmat);  //Pmat has penaly diagonal on the Dirichlet Nodes
      _EPSC->zero();

      // Solve the linear system
      ierr = KSPSolve(_ksp[0], RES, EPSC);			CHKERRABORT(MPI_COMM_WORLD,ierr);

      SolveTime = clock()-start_time;
      // ***************** END SOLVE ******************

      // ***************** UPDATE and Restore Matrix ******************
      start_time=clock();

      *_EPS += *_EPSC;

      //_KK->matrix_set_diagonal_values(_indexai[0],value);
      //_KK->close();

      _RESC->matrix_mult(*_EPSC,*_KK);
      *_RES -= *_RESC;
                  
      // Get the number of iterations required for convergence
      ierr = KSPGetIterationNumber(_ksp[0], &its);		CHKERRABORT(MPI_COMM_WORLD,ierr);

      // Get the norm of the final residual to return to the user.
      ierr = KSPGetResidualNorm(_ksp[0], &final_resid); 	CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = KSPDestroy(&_ksp[0]);				CHKERRABORT(MPI_COMM_WORLD,ierr);

      UpdateTime = clock() - start_time;
      // ***************** END UPDATE ******************

      MatDestroy(&Pmat);
  }
  else if(_DirichletBCsHandlingMode==1) { // By elimination

    PetscVector* RESp=static_cast<PetscVector*> (_RES);
    Vec RES=RESp->vec();
    PetscMatrix* KKp=static_cast<PetscMatrix*>(_KK);
    Mat KK=KKp->mat();
              
    // ***************** ASSEMBLE *******************
    clock_t start_time=clock();
 
    PetscInt Asize=_indexai[1].size();

//     // generate IS
//     IS isA;
//     ierr=ISCreateGeneral(MPI_COMM_WORLD,Asize,&_indexai[1][0], PETSC_USE_POINTER ,&isA);	   CHKERRABORT(MPI_COMM_WORLD,ierr);

    IS isA=_isA[0];
    Vec Pr;
    ierr = VecGetSubVector(RES,isA,&Pr);				CHKERRABORT(MPI_COMM_WORLD,ierr);

    // Solution Vector Pw
    Vec Pw;
    ierr = VecDuplicate(Pr,&Pw);			        	CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecSet(Pw,0);				        	CHKERRABORT(MPI_COMM_WORLD,ierr);

    //Matrix A
    Mat A;
    ierr = MatGetSubMatrix(KK,isA,isA,MAT_INITIAL_MATRIX,&A);		CHKERRABORT(MPI_COMM_WORLD,ierr);

    AssemblyTime = clock()-start_time;
    // ***************** END ASSEMBLE ******************

    // ***************** SOLVE ******************
    start_time=clock();

    // Initialize the data structure with the matrix PA
    this->init(A,false,false);

    // Solve the linear system
    ierr = KSPSolve(_ksp[0],Pr,Pw);					CHKERRABORT(MPI_COMM_WORLD,ierr);

    // Get the number of iterations required for convergence
    ierr = KSPGetIterationNumber(_ksp[0], &its);			CHKERRABORT(MPI_COMM_WORLD,ierr);

    // Get the norm of the final residual to return to the user.
    ierr = KSPGetResidualNorm(_ksp[0], &final_resid); 			CHKERRABORT(MPI_COMM_WORLD,ierr);
 
    ierr = KSPDestroy(&_ksp[0]);				        CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = MatDestroy(&A);					        CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecRestoreSubVector(RES,isA,&Pr);				CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecDestroy(&Pr);					        CHKERRABORT(MPI_COMM_WORLD,ierr);
        
    SolveTime = clock()-start_time;
    // ***************** END SOLVE ******************

    // ***************** UPDATE ******************
    start_time=clock();

    _EPSC->zero();
    PetscVector* EPSCp=static_cast<PetscVector*> (_EPSC);
    Vec EPSC = EPSCp->vec();

    const PetscInt *ind;
    ierr = ISGetIndices(isA,&ind);				        CHKERRABORT(MPI_COMM_WORLD,ierr);
    PetscScalar *W[1];
    ierr = VecGetArray(Pw,W);						CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecSetValues(EPSC,Asize,ind,W[0],ADD_VALUES);	        CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecRestoreArray(Pw,W);				        CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = ISRestoreIndices(isA,&ind);					CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = VecDestroy(&Pw);					        CHKERRABORT(MPI_COMM_WORLD,ierr);
    //ierr = ISDestroy(&isA);						CHKERRABORT(MPI_COMM_WORLD,ierr);

    *_EPS += *_EPSC;

    _RESC->matrix_mult(*_EPSC,*_KK);
    *_RES -= *_RESC;

    _RES->close();
    _EPS->close();

    UpdateTime = clock() - start_time;
    // ***************** END UPDATE ******************
  }
  
  // *** Computational info ***
  cout << "Grid: " << _msh->GetGridNumber()<< "      SOLVER TIME:              " <<
    static_cast<double>( SearchTime + AssemblyTime + SolveTime + UpdateTime)/ CLOCKS_PER_SEC<<
    "  ITS: " << its << endl;

  return std::make_pair(its,final_resid);
    
    
    
}

void PetscLinearEquationSolver::clear() {
  
  for(unsigned i=0;i<_isA.size();i++){
    ISDestroy(&_isA[i]); 	
  }
  
  for(unsigned i=0;i<_isB.size();i++){
    ISDestroy(&_isB[i]); 	
  }
  
  
  if (this->initialized()) {
    this->_is_initialized = false;
    int ierr=0;
    ierr = KSPDestroy(&_ksp[0]);	CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Mimic PETSc default solver and preconditioner
    this->_solver_type  = GMRES;
    if (!this->_preconditioner)      {
      int i;
      MPI_Comm_size(MPI_COMM_WORLD,&i);
      if (i == 1) this->_preconditioner_type = ILU_PRECOND;
      else this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
    }
  }
}

// ==============================================================
void PetscLinearEquationSolver::init() {
  // Initialize the data structures if not done so already.
  if (!this->initialized()) {
    this->_is_initialized = true;
    int ierr=0;
    // Create the linear solver context
    ierr = KSPCreate(MPI_COMM_WORLD, &_ksp[0]);  			CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Create the preconditioner context
    ierr = KSPGetPC(_ksp[0], &_pc);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Have the Krylov subspace method use our good initial guess rather than 0
    ierr = KSPSetInitialGuessNonzero(_ksp[0], PETSC_TRUE);		CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Set user-specified  solver and preconditioner types
    this->set_petsc_solver_type();
    // Set the options from user-input
    // Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    //  These options will override those specified above as long as
    //  KSPSetFromOptions() is called _after_ any other customization  routines.
    ierr = KSPSetFromOptions(_ksp[0]);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
    // NOT NECESSARY!!!!
    //ierr = PCSetFromOptions (_pc);
    //CHKERRABORT(MPI_COMM_WORLD,ierr);
    // #endif

    // Notify PETSc of location to store residual history.
    // This needs to be called before any solves, since
    // it sets the residual history length to zero.  The default
    // behavior is for PETSc to allocate (internally) an array
    // of size 1000 to hold the residual norm history.
    ierr = KSPSetResidualHistory(_ksp[0],
                                 PETSC_NULL,   // pointer to the array which holds the history
                                 PETSC_DECIDE, // size of the array holding the history
                                 PETSC_TRUE);  // Whether or not to reset the history for each solve.
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);

    //If there is a preconditioner object we need to set the internal setup and apply routines
    if (this->_preconditioner) {
      PCShellSetContext(_pc,(void*)this->_preconditioner);
      PCShellSetSetUp(_pc,__libmesh_petsc_preconditioner_setup);
      PCShellSetApply(_pc,__libmesh_petsc_preconditioner_apply);
    }
  }
}


// ========================================================
void PetscLinearEquationSolver::init(Mat& matrix, const bool pc_flag, const bool Schur) {
  // Initialize the data structures if not done so already.
  //I have to initialize many times the data structure therefore the next line must be commented
  //   if (!this->initialized()) {
  this->_is_initialized = true;
  int ierr=0;
  int grid=_msh->GetGridNumber();

  // Since the matrix is sequential for us, we use as communicator MPI_COMM_WORLD (serial) 
  // instead of PETSC_COMM_WORLD (parallel)
  ierr = KSPCreate(MPI_COMM_WORLD,&_ksp[0]);					CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = KSPSetOperators(_ksp[0],matrix,matrix,DIFFERENT_NONZERO_PATTERN);	CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = KSPSetInitialGuessNonzero(_ksp[0], PETSC_TRUE);			CHKERRABORT(MPI_COMM_WORLD,ierr);

  if (!Schur || grid==0) {
    ierr = KSPSetTolerances(_ksp[0],_rtol[0],_abstol[0],_dtol[0],_maxits[0]);	CHKERRABORT(MPI_COMM_WORLD,ierr);
  } else {
    ierr = KSPSetTolerances(_ksp[0],_rtol[1],_abstol[1],_dtol[1],_maxits[1]);	CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  ierr = KSPSetFromOptions(_ksp[0]);						CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Set user-specified  solver and preconditioner types
  this->set_petsc_solver_type();

  if (!pc_flag) {
    ierr = KSPGetPC(_ksp[0],&_pc);						CHKERRABORT(MPI_COMM_WORLD,ierr);  
    PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);
    PetscReal zero = 1.e-16;
    ierr = PCFactorSetZeroPivot(_pc,zero); 					CHKERRABORT(MPI_COMM_WORLD,ierr); 
    ierr = PCFactorSetShiftType(_pc,MAT_SHIFT_NONZERO); 			CHKERRABORT(MPI_COMM_WORLD,ierr); 
  }
}

// ========================================================
void PetscLinearEquationSolver::init_schur(Mat& matrix) {
  // Initialize the data structures if not done so already.
  //I have to initialize many times the data structure therefore the next line must be commented
  //   if (!this->initialized()) {
  this->_is_initialized = true;
  int ierr=0;

  ierr = MatSchurComplementGetKSP(matrix,&_ksp[1]);				CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = KSPSetInitialGuessNonzero(_ksp[1], PETSC_TRUE);			CHKERRABORT(MPI_COMM_WORLD,ierr);
  

  ierr = KSPSetTolerances(_ksp[1],_rtol[0],_abstol[0],_dtol[0],_maxits[0]);	CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = KSPSetFromOptions(_ksp[1]);						CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Set user-specified  solver and preconditioner types
  this->set_petsc_solver_type2();

  ierr = KSPGetPC(_ksp[1], &_pc);						CHKERRABORT(MPI_COMM_WORLD,ierr);
  PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);
  PetscReal zero = 1.e-16;
  PCFactorSetZeroPivot(_pc,zero);
  PCFactorSetShiftType(_pc,MAT_SHIFT_NONZERO);

  //   }
}

// ========================================================
void PetscLinearEquationSolver::init(Mat& Amat, Mat& Pmat) {
  
  
  
  // Initialize the data structures if not done so already.
  //if (!this->initialized())    {
  this->_is_initialized = true;
  int ierr=0;
  // Create the linear solver context
  ierr = KSPCreate(MPI_COMM_WORLD, &_ksp[0]);					CHKERRABORT(MPI_COMM_WORLD,ierr);
  //ierr = PCCreate (MPI_COMM_WORLD, &_pc); CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Create the preconditioner context
  ierr = KSPGetPC(_ksp[0], &_pc);						CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Set operators. The input matrix works as the preconditioning matrix
  ierr = KSPSetOperators(_ksp[0], Amat, Pmat, SAME_NONZERO_PATTERN);	CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Have the Krylov subspace method use our good initial guess rather than 0
  ierr = KSPSetInitialGuessNonzero(_ksp[0], PETSC_TRUE);			CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Set user-specified  solver and preconditioner types
  this->set_petsc_solver_type();
    
  if (!this->same_preconditioner)  {
    ierr = KSPSetOperators(_ksp[0], Amat, Pmat, SAME_NONZERO_PATTERN);	CHKERRABORT(MPI_COMM_WORLD,ierr);
  } else {
    ierr = KSPSetOperators(_ksp[0], Amat, Pmat, SAME_PRECONDITIONER);	CHKERRABORT(MPI_COMM_WORLD,ierr);
  }

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances(_ksp[0],_rtol[0],_abstol[0],_dtol[0],_maxits[0]);	CHKERRABORT(MPI_COMM_WORLD,ierr);
  //   ierr = KSPSetTolerances(_ksp[0], , PETSC_DEFAULT,PETSC_DEFAULT, max_its);CHKERRABORT(MPI_COMM_WORLD,ierr);
  
    
  // Set the options from user-input
  // Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
  //  These options will override those specified above as long as
  //  KSPSetFromOptions() is called _after_ any other customization  routines.
  ierr = KSPSetFromOptions(_ksp[0]);						CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
  //ierr = PCSetFromOptions (_pc);	CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Notify PETSc of location to store residual history.
  // This needs to be called before any solves, since
  // it sets the residual history length to zero.  The default
  // behavior is for PETSc to allocate (internally) an array
  // of size 1000 to hold the residual norm history.
  ierr = KSPSetResidualHistory(_ksp[0],
			       PETSC_NULL,   // pointer to the array which holds the history
			       PETSC_DECIDE, // size of the array holding the history
			       PETSC_TRUE);  // Whether or not to reset the history for each solve.
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);
  PetscReal zero = 1.e-16;
  PCFactorSetZeroPivot(_pc,zero);
  PCFactorSetShiftType(_pc,MAT_SHIFT_NONZERO);
  //     if (this->_preconditioner) {
  //       this->_preconditioner->set_matrix(*matrix);
  //       PCShellSetContext(_pc,(void*)this->_preconditioner);
  //       PCShellSetSetUp(_pc,__libmesh_petsc_preconditioner_setup);
  //       PCShellSetApply(_pc,__libmesh_petsc_preconditioner_apply);
  //     }
  //}
}





// =========================================================================
void PetscLinearEquationSolver::get_residual_history(std::vector<double>& hist) {
  int ierr = 0;
  int its  = 0;
  // Fill the residual history vector with the residual norms
  // Note that GetResidualHistory() does not copy any values, it
  // simply sets the pointer p.  Note that for some Krylov subspace
  // methods, the number of residuals returned in the history
  // vector may be different from what you are expecting.  For
  // example, TFQMR returns two residual values per iteration step.
  PetscReal* p;
  ierr = KSPGetResidualHistory(_ksp[0], &p, &its);				  CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Check for early return
  if (its == 0) return;
  // Create space to store the result
  hist.resize(its);
  // Copy history into the vector provided by the user.
  for (int i=0; i<its; ++i) {
    hist[i] = *p;
    p++;
  }
}

// ======================================================
double PetscLinearEquationSolver::get_initial_residual() {
  int ierr = 0;
  int its  = 0;
  // Fill the residual history vector with the residual norms
  // Note that GetResidualHistory() does not copy any values, it
  // simply sets the pointer p.  Note that for some Krylov subspace
  // methods, the number of residuals returned in the history
  // vector may be different from what you are expecting.  For
  // example, TFQMR returns two residual values per iteration step.
  PetscReal* p;
  ierr = KSPGetResidualHistory(_ksp[0], &p, &its);					CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Check no residual history
  if (its == 0) {
    std::cerr << "No iterations have been performed, returning 0." << std::endl;
    return 0.;
  }
  // Otherwise, return the value pointed to by p.
  return *p;
}

// =================================================
void PetscLinearEquationSolver::set_petsc_solver_type() {
  int ierr = 0;
  switch (this->_solver_type) {
  case CG:
    ierr = KSPSetType(_ksp[0], (char*) KSPCG);						CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case CR:
    ierr = KSPSetType(_ksp[0], (char*) KSPCR);						CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case CGS:
    ierr = KSPSetType(_ksp[0], (char*) KSPCGS);						CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case BICG:
    ierr = KSPSetType(_ksp[0], (char*) KSPBICG);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case TCQMR:
    ierr = KSPSetType(_ksp[0], (char*) KSPTCQMR);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case TFQMR:
    ierr = KSPSetType(_ksp[0], (char*) KSPTFQMR);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case LSQR:
    ierr = KSPSetType(_ksp[0], (char*) KSPLSQR);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case BICGSTAB:
    ierr = KSPSetType(_ksp[0], (char*) KSPBCGS);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case MINRES:
    ierr = KSPSetType(_ksp[0], (char*) KSPMINRES);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case GMRES:
    ierr = KSPSetType(_ksp[0], (char*) KSPGMRES);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case RICHARDSON:
    ierr = KSPSetType(_ksp[0], (char*) KSPRICHARDSON);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case CHEBYSHEV:
    ierr = KSPSetType(_ksp[0], (char*) KSPCHEBYSHEV);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case PREONLY:
    ierr = KSPSetType(_ksp[0], (char*) KSPPREONLY);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  default:
    std::cerr << "ERROR:  Unsupported PETSC Solver: "
              << this->_solver_type               << std::endl
              << "Continuing with PETSC defaults" << std::endl;
  }
}


// =================================================
void PetscLinearEquationSolver::set_petsc_solver_type2() {
  int ierr = 0;
  switch (this->_solver_type) {
  case CG:
    ierr = KSPSetType(_ksp[1], (char*) KSPCG);						CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case CR:
    ierr = KSPSetType(_ksp[1], (char*) KSPCR);						CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case CGS:
    ierr = KSPSetType(_ksp[1], (char*) KSPCGS);						CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case BICG:
    ierr = KSPSetType(_ksp[1], (char*) KSPBICG);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case TCQMR:
    ierr = KSPSetType(_ksp[1], (char*) KSPTCQMR);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case TFQMR:
    ierr = KSPSetType(_ksp[1], (char*) KSPTFQMR);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case LSQR:
    ierr = KSPSetType(_ksp[1], (char*) KSPLSQR);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case BICGSTAB:
    ierr = KSPSetType(_ksp[1], (char*) KSPBCGS);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case MINRES:
    ierr = KSPSetType(_ksp[1], (char*) KSPMINRES);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case GMRES:
    ierr = KSPSetType(_ksp[1], (char*) KSPGMRES);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case RICHARDSON:
    ierr = KSPSetType(_ksp[1], (char*) KSPRICHARDSON);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case CHEBYSHEV:
    ierr = KSPSetType(_ksp[1], (char*) KSPCHEBYSHEV);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case PREONLY:
    ierr = KSPSetType(_ksp[1], (char*) KSPPREONLY);					CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  default:
    std::cerr << "ERROR:  Unsupported PETSC Solver: "
              << this->_solver_type               << std::endl
              << "Continuing with PETSC defaults" << std::endl;
  }
}

// =======================================================
void PetscLinearEquationSolver::print_converged_reason() {
  KSPConvergedReason reason;
  KSPGetConvergedReason(_ksp[0], &reason);
  std::cout << "Linear solver convergence/divergence reason: " << KSPConvergedReasons[reason] << std::endl;
}

#endif 

