//C++ include
#include <ctime>
#include <fstream>
#include <algorithm>
#include "LinSysPde.hpp"
#include "ElemType.hpp"
#include "ParalleltypeEnum.hpp"
#include "NumericVector.hpp"
#include "PetscVector.hpp"
#include "SparseRectangularMatrix.hpp"
#include "PetscRectangularMatrix.hpp"

using std::cout;
using std::endl;

//--------------------------------------------------------------------------------
lsysPde::lsysPde(mesh *other_msh){    
  _msh = other_msh;
  _is_symmetric = false;
  _stabilization = false;
  _compressibility = 0.;
  CC_flag=0;
}

//--------------------------------------------------------------------------------
lsysPde::~lsysPde() { }

//--------------------------------------------------------------------------------
void lsysPde::SetMatrixProperties(const bool property) {
  _is_symmetric = property;
}

//--------------------------------------------------------------------------------
bool lsysPde::GetMatrixProperties() {
  return _is_symmetric;
}

//--------------------------------------------------------------------------------
void lsysPde::AddStabilization(const bool stab, const double compressibility) {
  _stabilization = stab;
  _compressibility = compressibility;
}

//--------------------------------------------------------------------------------
double lsysPde::GetCompressibility() {
  return _compressibility;
};

//--------------------------------------------------------------------------------
bool lsysPde::GetStabilization() {
  return _stabilization;
};

//--------------------------------------------------------------------------------
unsigned lsysPde::GetIndex(const char name[]) {
  unsigned index=0;
  while (strcmp(_SolName[index],name)) {
    index++;
    if (index==_SolType.size()) {
      cout<<"error! invalid name entry GetIndex(...)"<<endl;
      exit(0);
    }
  }
  return index;
}

//--------------------------------------------------------------------------------
void lsysPde::SetBdcPointer(vector <NumericVector*> *Bdc_other){
    _Bdc=Bdc_other;
}

//--------------------------------------------------------------------------------
int lsysPde::InitMultigrid(const vector <unsigned> &_SolPdeIndex, const  vector <int> &SolType_other,  const vector <char*> &SolName_other ) {
   
  _SolType=SolType_other;
  _SolName=SolName_other;
  
  
  int ierr;
  KKIndex.resize(_SolPdeIndex.size()+1u);
  KKIndex[0]=0;
  for (unsigned i=1; i<KKIndex.size(); i++)
//     KKIndex[i]=KKIndex[i-1]+GetDofNumber(SolType[_SolPdeIndex[i-1]]);
  KKIndex[i]=KKIndex[i-1]+_msh->MetisOffset[_SolType[_SolPdeIndex[i-1]]][_msh->_nprocs];

  //-----------------------------------------------------------------------------------------------
  KKoffset.resize(_SolPdeIndex.size()+1);
  for(int i=0;i<_SolPdeIndex.size()+1;i++) {
    KKoffset[i].resize(_msh->nsubdom);
  }
  
   KKoffset[0][0]=0;
   for(int j=1; j<_SolPdeIndex.size()+1; j++) {
      unsigned indexSol=_SolPdeIndex[j-1];
      KKoffset[j][0] = KKoffset[j-1][0]+(_msh->MetisOffset[_SolType[indexSol]][1] - _msh->MetisOffset[_SolType[indexSol]][0]);
   }
  
  for(int i=1; i<_msh->nsubdom; i++) {
    KKoffset[0][i] = KKoffset[_SolPdeIndex.size()][i-1];
    for(int j=1; j<_SolPdeIndex.size()+1; j++) {
      unsigned indexSol=_SolPdeIndex[j-1];
      KKoffset[j][i] = KKoffset[j-1][i]+(_msh->MetisOffset[_SolType[indexSol]][i+1] - _msh->MetisOffset[_SolType[indexSol]][i]);
    }
  }
   
  //ghost size
  KKghostsize.resize(_msh->nsubdom,0);
  for(int i=0; i<_msh->nsubdom; i++) {
    for(int j=0; j<_SolPdeIndex.size(); j++) {
      unsigned indexSol=_SolPdeIndex[j];
      KKghostsize[i] += _msh->ghost_size[_SolType[indexSol]][i];
    }
  }
  
  //ghost nodes
  KKghost_nd.resize(_msh->nsubdom);
  KKghost_nd[0].resize(1);  KKghost_nd[0][0]=1;
  for(int i=1; i<_msh->nsubdom; i++) {
    KKghost_nd[i].resize(KKghostsize[i]);
  }
  
  
   for(int i=0; i<_msh->nsubdom; i++) {
     unsigned counter=0;
     for(int j=0; j<_SolPdeIndex.size(); j++) {
        unsigned indexSol=_SolPdeIndex[j];
	for(int k=0; k<_msh->ghost_size[_SolType[indexSol]][i];k++) {
	  //gambit ghost node
	  unsigned gmt_ghost_nd = _msh->ghost_nd[_SolType[indexSol]][i][k];
	  KKghost_nd[i][counter] =  GetKKDof(indexSol,j,gmt_ghost_nd);
	  counter++;
	}
     }
   }
   
  //-----------------------------------------------------------------------------------------------
  
  
//   PetscInt RESsize= KKIndex[KKIndex.size()-1];
// 
//   if(_msh->_nprocs==1) {  
//     ierr=VecCreateSeq(PETSC_COMM_SELF,RESsize,&RES);
//     CHKERRQ(ierr);
//     ierr=VecSetFromOptions(RES);
//     CHKERRQ(ierr);
//   } else {
//     PetscInt RES_local_size =KKoffset[KKIndex.size()-1][_msh->_iproc] - KKoffset[0][_msh->_iproc];
//     if(_msh->_iproc==0)
//       ierr=VecCreateGhost(PETSC_COMM_WORLD,RES_local_size,RESsize,0,PETSC_NULL,&RES); 
//     else
//       ierr=VecCreateGhost(PETSC_COMM_WORLD,RES_local_size,RESsize,KKghostsize[_msh->_iproc],&KKghost_nd[_msh->_iproc][0],&RES); 
//     CHKERRQ(ierr);
//   }
//   
// 
//   ierr=VecDuplicate(RES,&RESC);
//   CHKERRQ(ierr);
//   ierr=VecDuplicate(RES,&EPS);
//   CHKERRQ(ierr);
//   ierr=VecDuplicate(RES,&EPSC);
//   CHKERRQ(ierr);

  int EPSsize= KKIndex[KKIndex.size()-1];
  _EPS = NumericVector::build().release();
  if(_msh->_nprocs==1) { // IF SERIAL
    _EPS->init(EPSsize,EPSsize,false,SERIAL);
  } 
  else { // IF PARALLEL
    int EPS_local_size =KKoffset[KKIndex.size()-1][_msh->_iproc] - KKoffset[0][_msh->_iproc];
    _EPS->init(EPSsize,EPS_local_size, KKghost_nd[_msh->_iproc], false,GHOSTED);
  }
  
  _RES = NumericVector::build().release();
  _RES->init(*_EPS);
  
  _EPSC = NumericVector::build().release();
  _EPSC->init(*_EPS);
  
  _RESC = NumericVector::build().release();
  _RESC->init(*_EPS);
  
  PetscVector* EPSp=static_cast<PetscVector*> (_EPS);  //TODO
  PetscVector* EPSCp=static_cast<PetscVector*> (_EPSC);  //TODO
  PetscVector* RESp=static_cast<PetscVector*> (_RES);  //TODO
  PetscVector* RESCp=static_cast<PetscVector*> (_RESC); //TODO
  
  EPS=EPSp->vec(); //TODO
  EPSC=EPSCp->vec(); //TODO
  RES=RESp->vec(); //TODO
  RESC=RESCp->vec(); //TODO
  
  //--------------------------------------------------------------------------------------
  DrchKKdofs.resize(KKoffset[KKIndex.size()-1][_msh->_iproc] - KKoffset[0][_msh->_iproc]);
  unsigned counter=0;
  for(int k=0; k<_SolPdeIndex.size(); k++) {
    unsigned indexSol=_SolPdeIndex[k];
    unsigned soltype=_SolType[indexSol];
    if(soltype<3) {
      for(unsigned inode_mts=_msh->MetisOffset[soltype][_msh->_iproc]; 
	 inode_mts<_msh->MetisOffset[soltype][_msh->_iproc+1]; inode_mts++) {
	 if((*(*_Bdc)[indexSol])(inode_mts)<1.9) {
	   int local_mts = inode_mts-_msh->MetisOffset[soltype][_msh->_iproc];
	   int idof_kk = KKoffset[k][_msh->_iproc] +local_mts; 
	   DrchKKdofs[counter]=idof_kk;
	   counter++;
	 }
      }
    }
  } 
  DrchKKdofs.resize(counter);
 
  
  return 1;
}

//--------------------------------------------------------------------------------
int lsysPde::SetResZero() {
//   int ierr;
//   
//   if(_msh->_nprocs==1) {
//     ierr=VecSet(RES,0.);
//     CHKERRQ(ierr);
//   }  else {
//      /* Vectors that include ghost values require a special
//     handling.  */
//     Vec loc_vec; ierr = VecGhostGetLocalForm (RES,&loc_vec);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     ierr = VecSet (loc_vec, 0.);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     ierr = VecGhostRestoreLocalForm (RES,&loc_vec);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
// }
// return ierr;

  _RES->zero();
  return 1;
}

//--------------------------------------------------------------------------------
int lsysPde::SetEpsZero() {
//   int ierr;
//   if(_msh->_nprocs==1) {
//     ierr=VecSet(EPS,0.);
//     CHKERRQ(ierr);
//     ierr=VecSet(EPSC,0.);
//     CHKERRQ(ierr);
//   } else {
//      /* Vectors that include ghost values require a special
//     handling.  */
//     Vec loc_vec; ierr = VecGhostGetLocalForm (EPS,&loc_vec);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     ierr = VecSet (loc_vec, 0.);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     ierr = VecGhostRestoreLocalForm (EPS,&loc_vec);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     
//     ierr = VecGhostGetLocalForm (EPSC,&loc_vec);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     ierr = VecSet (loc_vec, 0.);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     ierr = VecGhostRestoreLocalForm (EPSC,&loc_vec);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     
//   }
//   return ierr;

  _EPS->zero();
  _EPSC->zero();
  return 1;
}

//--------------------------------------------------------------------------------
int lsysPde::SumEpsCToEps() {
//   int ierr;
//   ierr=VecAXPBY(EPS,1,1,EPSC);
//   CHKERRQ(ierr);
//   return 1;

  *_EPS += *_EPSC;

}

//--------------------------------------------------------------------------------
int lsysPde::UpdateResidual() {
//   int ierr;
//   ierr = MatMult(KK,EPSC,RESC);
//   CHKERRQ(ierr);
//   ierr = VecAXPBY(RES,-1.,1.,RESC);
//   CHKERRQ(ierr);
  
  
  _RESC->matrix_mult(*_EPSC,*_KK);
  *_RES -= *_RESC;
  
  return 1;
}

//-------------------------------------------------------------------------------------------
int lsysPde::DeallocateMatrix() {
  int ierr;
//   ierr=MatDestroy(&KK);
//   CHKERRQ(ierr);
//   
  delete _KK;
  
  if (_msh->GetGridNumber()>0) {
//     ierr=MatDestroy(&PP);
//     CHKERRQ(ierr);
    delete _PP;
  }
 /* 
  ierr=VecDestroy(&RES);
  CHKERRQ(ierr);
  ierr=VecDestroy(&RESC);
  CHKERRQ(ierr);

  ierr=VecDestroy(&EPSC);
  CHKERRQ(ierr);
  ierr=VecDestroy(&EPS);
  CHKERRQ(ierr);*/
  
  delete _EPS;
  delete _EPSC;
  delete _RES;
  delete _RESC;
  
  return ierr;
}

//-------------------------------------------------------------------------------------------
int lsysPde::AllocateMatrix() {

//   PetscInt KKsize=KKIndex[KKIndex.size()-1u];
//   PetscErrorCode ierr;
//   const unsigned dim = _msh->GetDimension();
//   PetscInt KK_UNIT_SIZE_ = pow(5,dim);
// 
//   if(_msh->_nprocs==1) {
//     //TODO optimize memory
//     ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,KKsize,KKsize,KK_UNIT_SIZE_*KKIndex.size(),PETSC_NULL,&KK);
//     CHKERRQ(ierr);
//     
//     ierr = MatSetFromOptions(KK);
//     CHKERRQ(ierr);
//     
//     ierr=MatSetOption(KK,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE); 
//     CHKERRQ(ierr);
//     
//   }
//   else {
// 
//     PetscInt KK_local_size =KKoffset[KKIndex.size()-1][_msh->_iproc] - KKoffset[0][_msh->_iproc];
//     
//     ierr = MatCreate(MPI_COMM_WORLD, &KK);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//     ierr = MatSetSizes(KK, KK_local_size, KK_local_size, KKsize,KKsize);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//     ierr = MatSetType(KK, MATMPIAIJ);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//     ierr = MatMPIAIJSetPreallocation(KK, KK_UNIT_SIZE_*KKIndex.size(), PETSC_NULL, KK_UNIT_SIZE_*KKIndex.size(), PETSC_NULL);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//   
//     ierr = MatSetFromOptions(KK);
//     CHKERRQ(ierr);
//     
//     ierr=MatSetOption(KK,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE); 
//     CHKERRQ(ierr);
//     
//   }

  
  const unsigned dim = _msh->GetDimension();
  int KK_UNIT_SIZE_ = pow(5,dim);
  int KK_size=KKIndex[KKIndex.size()-1u];
  int KK_local_size =KKoffset[KKIndex.size()-1][_msh->_iproc] - KKoffset[0][_msh->_iproc];
    
  _KK = SparseRectangularMatrix::build().release();
  _KK->init(KK_size,KK_size,KK_local_size,KK_local_size,KK_UNIT_SIZE_*KKIndex.size(),KK_UNIT_SIZE_*KKIndex.size());
  
  PetscRectangularMatrix* KKp=static_cast<PetscRectangularMatrix*>(_KK); //TODO
  KK=KKp->mat(); //TODO
  
  return 1;
}

//-------------------------------------------------------------------------------------------


