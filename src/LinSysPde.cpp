//C++ include
#include <ctime>
#include <fstream>
#include <algorithm>
#include "LinSysPde.hpp"
#include "ElemType.hpp"
#include "ParalleltypeEnum.hpp"
#include "NumericVector.hpp"
#include "SparseRectangularMatrix.hpp"
using std::cout;
using std::endl;

const unsigned lsysPDE::END_IND[5]= {0,1,3,4,5};

//--------------------------------------------------------------------------------
//lsysPDE::lsysPDE(const char infile[], vector < vector < double> > &vt,const double Lref):
//  mesh(infile, vt,Lref) {
    
lsysPDE::lsysPDE(const char infile[], vector < vector < double> > &vt,const double Lref){
  _msh = new mesh(infile, vt,Lref);
  _is_symmetric = false;
  _stabilization = false;
  _compressibility = 0.;
  for(int i=0;i<5;i++){
    Proj_mat_flag[i]=0;
  }
  CC_flag=0;
};

//--------------------------------------------------------------------------------

//lsysPDE::lsysPDE(const unsigned &igrid,elem *elc):
//  mesh(igrid,elc) {
    

lsysPDE::lsysPDE(const unsigned &igrid,elem *elc){    
  _msh = new mesh(igrid,elc); 
  _is_symmetric = false;
  _stabilization = false;
  _compressibility = 0.;
  for(int i=0;i<5;i++){
    Proj_mat_flag[i]=0;
  }
  CC_flag=0;
}

lsysPDE::lsysPDE(mesh *other_msh){    
  _msh = other_msh;
  _is_symmetric = false;
  _stabilization = false;
  _compressibility = 0.;
  for(int i=0;i<5;i++){
    Proj_mat_flag[i]=0;
  }
  CC_flag=0;
}

//--------------------------------------------------------------------------------
lsysPDE::~lsysPDE() {
  for (unsigned i=0; i<SolName.size(); i++) {
    delete [] SolName[i];
  }
  //delete _msh;
}

/**
 * This function flags the elements that will be refined
 **/
//--------------------------------------------------------------------------------
void lsysPDE::set_elr(const unsigned &test) {
 
  _msh->el->InitRefinedToZero();
  
  unsigned grid = _msh->GetGridNumber();
  unsigned nel = _msh->GetElementNumber();
   
  if (test==1) { //refine all next grid elements
    for (unsigned iel=0; iel<nel; iel++) {
      _msh->el->SetRefinedElementIndex(iel,1);
      _msh->el->AddToRefinedElementNumber(1);
      short unsigned elt=_msh->el->GetElementType(iel);
      _msh->el->AddToRefinedElementNumber(1,elt);
    }
  } else if (test==2) { //refine based of the position of the mid point;
    //serial loop
    for (unsigned iel=0; iel<nel; iel+=1) {
      unsigned nve=_msh->el->GetElementDofNumber(iel,0);
      double vtx=0.,vty=0.,vtz=0.;
      for ( unsigned i=0; i<nve; i++) {
        unsigned inode=_msh->el->GetElementVertexIndex(iel,i)-1u;
	unsigned inode_Metis=_msh->GetMetisDof(inode,2);
// 	vtx+=(*Sol_[0])(inode_Metis);  
// 	vty+=(*Sol_[1])(inode_Metis);
// 	vtz+=(*Sol_[2])(inode_Metis);
      }
      vtx/=nve;
      vty/=nve;
      vtz/=nve;
      if (_msh->_SetRefinementFlag(vtx,vty,vtz,_msh->el->GetElementGroup(iel),grid)) {
        _msh->el->SetRefinedElementIndex(iel,1);
        _msh->el->AddToRefinedElementNumber(1);
        short unsigned elt=_msh->el->GetElementType(iel);
        _msh->el->AddToRefinedElementNumber(1,elt);
      }
    }
  
//     //parallel loop
//     for(int isdom=_iproc; isdom<_iproc+1; isdom++) {
//       for(int iel=IS_Mts2Gmt_elem_offset[isdom]; iel < IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
//         unsigned kel = IS_Mts2Gmt_elem[iel];
//         unsigned nve=el->GetElementDofNumber(kel,0);
//         double vtx=0.,vty=0.,vtz=0.;
//         for ( unsigned i=0; i<nve; i++) {
//           unsigned inode=el->GetElementVertexIndex(kel,i)-1u;
// 	  unsigned inode_Metis=GetMetisDof(inode,2);
// 	  vtx+=(*Sol_[0])(inode_Metis);  
// 	  vty+=(*Sol_[1])(inode_Metis);
// 	  vtz+=(*Sol_[2])(inode_Metis);
//         }
//         vtx/=nve;
//         vty/=nve;
//         vtz/=nve;
//         if (_SetRefinementFlag(vtx,vty,vtz,el->GetElementGroup(kel),grid)) {
//           el->SetRefinedElementIndex(kel,1);
//           el->AddToRefinedElementNumber(1);
//           short unsigned elt=el->GetElementType(kel);
//           el->AddToRefinedElementNumber(1,elt);
//        }
//      } 
//    }
   
  } else if (test==3) { //refine all next grid even elements
    for (unsigned iel=0; iel<nel; iel+=2) {
      _msh->el->SetRefinedElementIndex(iel,1);
      _msh->el->AddToRefinedElementNumber(1);
      short unsigned elt=_msh->el->GetElementType(iel);
      _msh->el->AddToRefinedElementNumber(1,elt);
    }
  }
  _msh->el->AllocateChildrenElement(_msh->_ref_index);
}

//--------------------------------------------------------------------------------
void lsysPDE::GetCoordinates( vector < vector < double> > &vt){
  unsigned indexSol;
 
  indexSol=GetIndex("X");
  *Sol_[indexSol]=vt[0];
//   Sol_[indexSol]->close();
  
  indexSol=GetIndex("Y");
  *Sol_[indexSol]=vt[1];
//   Sol_[indexSol]->close();
  
  indexSol=GetIndex("Z");
  *Sol_[indexSol]=vt[2];
//   Sol_[indexSol]->close();
  
}

//--------------------------------------------------------------------------------
void lsysPDE::AddSolutionVector( const char name[], const char order[],
                                const unsigned& tmorder, const bool &PDE_type) {
  unsigned n=Sol_.size();

  //Bdc.resize(n+1u);
  SolType.resize(n+1u);
  SolName.resize(n+1u);
  SolTmOrder.resize(n+1u);
  SolTmOrder[n]=tmorder;

  //new
  Sol_.resize(n+1u);
  Res_.resize(n+1u);
  Eps_.resize(n+1u);
  Bdc_.resize(n+1u);
  ResEpsBdc_flag_.resize(n+1u);
  if(PDE_type) ResEpsBdc_flag_[n]=1;

  Sol_old_.resize(n+1u);

  if (!strcmp(order,"linear")) {
    SolType[n]=0;
  } else if (!strcmp(order,"quadratic")) {
    SolType[n]=1;
  } else if (!strcmp(order,"biquadratic")) {
    SolType[n]=2;
  } else if (!strcmp(order,"constant")) {
    SolType[n]=3;
  } else if (!strcmp(order,"disc_linear")) {
    SolType[n]=4;
  } else {
    cout<<"error! invalid order entry in AddSolutionVector(...)"<<endl;
    exit(0);
  }
  SolName[n]=new char [8];
  strcpy(SolName[n],name);

}



//--------------------------------------------------------------------------------
void lsysPDE::SetMatrixProperties(const bool property) {
  _is_symmetric = property;
}

//--------------------------------------------------------------------------------
bool lsysPDE::GetMatrixProperties() {
  return _is_symmetric;
}

//--------------------------------------------------------------------------------
void lsysPDE::AddStabilization(const bool stab, const double compressibility) {
  _stabilization = stab;
  _compressibility = compressibility;
}

//--------------------------------------------------------------------------------
double lsysPDE::GetCompressibility() {
  return _compressibility;
};

//--------------------------------------------------------------------------------
bool lsysPDE::GetStabilization() {
  return _stabilization;
};

//--------------------------------------------------------------------------------
void lsysPDE::ResizeSolutionVector(const char name[]) {

  unsigned i=GetIndex(name);
  
  Sol_[i] = NumericVector::build().release();
  if(_msh->_nprocs==1) {
    Sol_[i]->init(_msh->MetisOffset[SolType[i]][_msh->_nprocs],_msh->own_size[SolType[i]][_msh->_iproc],false,SERIAL);
  } else {
    if(SolType[i]<3) {
     if(_msh->ghost_size[SolType[i]][_msh->_iproc]!=0) { 
      Sol_[i]->init(_msh->MetisOffset[SolType[i]][_msh->_nprocs],_msh->own_size[SolType[i]][_msh->_iproc],_msh->ghost_nd_mts[SolType[i]][_msh->_iproc],
		    false,GHOSTED);
      } else {
      std::vector <int> fake_ghost(1,_msh->own_size[SolType[i]][_msh->_iproc]);
      Sol_[i]->init(_msh->MetisOffset[SolType[i]][_msh->_nprocs],_msh->own_size[SolType[i]][_msh->_iproc],fake_ghost,false,GHOSTED);
      }
    }
    else { //discontinuous pressure has not ghost nodes
      Sol_[i]->init(_msh->MetisOffset[SolType[i]][_msh->_nprocs],_msh->own_size[SolType[i]][_msh->_iproc],false,PARALLEL); 
      
    }
  }
  
  if (SolTmOrder[i]==2) {
    Sol_old_[i] = NumericVector::build().release();
    Sol_old_[i]->init(*Sol_[i]);
  }
    
  if(ResEpsBdc_flag_[i]) {
    
    Res_[i] = NumericVector::build().release();
    Res_[i]->init(*Sol_[i]);

    Eps_[i] = NumericVector::build().release();
    Eps_[i]->init(*Sol_[i]);
    
    Bdc_[i] = NumericVector::build().release();
    Bdc_[i]->init(*Sol_[i]);
    
    //Bdc[i].resize(MetisOffset[SolType[i]][_nprocs]);
  }
}

//--------------------------------------------------------------------------------
unsigned lsysPDE::GetIndex(const char name[]) {
  unsigned index=0;
  while (strcmp(SolName[index],name)) {
    index++;
    if (index==Res_.size()) {
      cout<<"error! invalid name entry GetIndex(...)"<<endl;
      exit(0);
    }
  }
  return index;
}

//--------------------------------------------------------------------------------
int lsysPDE::InitMultigrid(const vector <unsigned> &MGIndex) {

  int ierr;
  KKIndex.resize(MGIndex.size()+1u);
  KKIndex[0]=0;
  for (unsigned i=1; i<KKIndex.size(); i++)
//     KKIndex[i]=KKIndex[i-1]+GetDofNumber(SolType[MGIndex[i-1]]);
  KKIndex[i]=KKIndex[i-1]+_msh->MetisOffset[SolType[MGIndex[i-1]]][_msh->_nprocs];

  //-----------------------------------------------------------------------------------------------
  KKoffset.resize(MGIndex.size()+1);
  for(int i=0;i<MGIndex.size()+1;i++) {
    KKoffset[i].resize(_msh->nsubdom);
  }
  
   KKoffset[0][0]=0;
   for(int j=1; j<MGIndex.size()+1; j++) {
      unsigned indexSol=MGIndex[j-1];
      KKoffset[j][0] = KKoffset[j-1][0]+(_msh->MetisOffset[SolType[indexSol]][1] - _msh->MetisOffset[SolType[indexSol]][0]);
   }
  
  for(int i=1; i<_msh->nsubdom; i++) {
    KKoffset[0][i] = KKoffset[MGIndex.size()][i-1];
    for(int j=1; j<MGIndex.size()+1; j++) {
      unsigned indexSol=MGIndex[j-1];
      KKoffset[j][i] = KKoffset[j-1][i]+(_msh->MetisOffset[SolType[indexSol]][i+1] - _msh->MetisOffset[SolType[indexSol]][i]);
    }
  }
   
  //ghost size
  KKghostsize.resize(_msh->nsubdom,0);
  for(int i=0; i<_msh->nsubdom; i++) {
    for(int j=0; j<MGIndex.size(); j++) {
      unsigned indexSol=MGIndex[j];
      KKghostsize[i] += _msh->ghost_size[SolType[indexSol]][i];
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
     for(int j=0; j<MGIndex.size(); j++) {
        unsigned indexSol=MGIndex[j];
	for(int k=0; k<_msh->ghost_size[SolType[indexSol]][i];k++) {
	  //gambit ghost node
	  unsigned gmt_ghost_nd = _msh->ghost_nd[SolType[indexSol]][i][k];
	  KKghost_nd[i][counter] =  GetKKDof(indexSol,j,gmt_ghost_nd);
	  counter++;
	}
     }
   }
   
//    for(int i=0; i<nsubdom; i++) {
//      cout<<i<<" "<<KKghostsize[i]<<endl;
//      for(int j=0;j<KKghostsize[i];j++){
//        cout<<KKghost_nd[i][j]<<" ";
//      }
//      cout<<endl;
//    }
  
 
//   printing for debugging purpose
//   cout << "Grid: " << GetGridNumber() << endl;
//   for(int i=0; i<nsubdom; i++) {
//   cout << "nsubdom: " << i << endl;
//   for(int j=0; j<MGIndex.size()+1; j++) {
//    cout<< j << "  " << KKoffset[j][i]<<endl;
//     }
//   }
  
    //-----------------------------------------------------------------------------------------------
  
  
  PetscInt RESsize= KKIndex[KKIndex.size()-1];

  if(_msh->_nprocs==1) {  
    ierr=VecCreateSeq(PETSC_COMM_SELF,RESsize,&RES);
    CHKERRQ(ierr);
    ierr=VecSetFromOptions(RES);
    CHKERRQ(ierr);
  } else {
    PetscInt RES_local_size =KKoffset[KKIndex.size()-1][_msh->_iproc] - KKoffset[0][_msh->_iproc];
    if(_msh->_iproc==0)
      ierr=VecCreateGhost(PETSC_COMM_WORLD,RES_local_size,RESsize,0,PETSC_NULL,&RES); 
    else
      ierr=VecCreateGhost(PETSC_COMM_WORLD,RES_local_size,RESsize,KKghostsize[_msh->_iproc],&KKghost_nd[_msh->_iproc][0],&RES); 
    CHKERRQ(ierr);
  }
  

  ierr=VecDuplicate(RES,&RESC);
  CHKERRQ(ierr);
  ierr=VecDuplicate(RES,&EPS);
  CHKERRQ(ierr);
  ierr=VecDuplicate(RES,&EPSC);
  CHKERRQ(ierr);

  
  //--------------------------------------------------------------------------------------
  DrchKKdofs.resize(KKoffset[KKIndex.size()-1][_msh->_iproc] - KKoffset[0][_msh->_iproc]);
  unsigned counter=0;
  for(int k=0; k<MGIndex.size(); k++) {
    unsigned indexSol=MGIndex[k];
    unsigned soltype=SolType[indexSol];
    if(soltype<3) {
      for(unsigned inode_mts=_msh->MetisOffset[soltype][_msh->_iproc]; 
	 inode_mts<_msh->MetisOffset[soltype][_msh->_iproc+1]; inode_mts++) {
	 if((*Bdc_[indexSol])(inode_mts)<1.9) {
	   int local_mts = inode_mts-_msh->MetisOffset[soltype][_msh->_iproc];
	   int idof_kk = KKoffset[k][_msh->_iproc] +local_mts; 
	   DrchKKdofs[counter]=idof_kk;
	   counter++;
	 }
      }
    }
  } 
  DrchKKdofs.resize(counter);
  
//   for(int i=0;i< DrchKKdofs.size();i++){
//     cout<< DrchKKdofs[i]<<" ";
//   }
//   cout<<endl<<DrchKKdofs.size()<<endl;
  
  return 1;
}

//--------------------------------------------------------------------------------
int lsysPDE::SetResZero(const vector <unsigned> &MGIndex) {
  int ierr;
  
  if(_msh->_nprocs==1) {
    ierr=VecSet(RES,0.);
    CHKERRQ(ierr);
  }  else {
     /* Vectors that include ghost values require a special
    handling.  */
    Vec loc_vec; ierr = VecGhostGetLocalForm (RES,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecSet (loc_vec, 0.);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostRestoreLocalForm (RES,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  
  return ierr;
}

//--------------------------------------------------------------------------------
int lsysPDE::SetEpsZero(const vector <unsigned> &MGIndex) {
  int ierr;
  if(_msh->_nprocs==1) {
    ierr=VecSet(EPS,0.);
    CHKERRQ(ierr);
    ierr=VecSet(EPSC,0.);
    CHKERRQ(ierr);
  } else {
     /* Vectors that include ghost values require a special
    handling.  */
    Vec loc_vec; ierr = VecGhostGetLocalForm (EPS,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecSet (loc_vec, 0.);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostRestoreLocalForm (EPS,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    
    ierr = VecGhostGetLocalForm (EPSC,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecSet (loc_vec, 0.);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostRestoreLocalForm (EPSC,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    
  }
  return ierr;
}

//--------------------------------------------------------------------------------
int lsysPDE::SumEpsCToEps(const vector <unsigned> &MGIndex) {
  int ierr;
  ierr=VecAXPBY(EPS,1,1,EPSC);
  CHKERRQ(ierr);
  return 1;
}

//--------------------------------------------------------------------------------
int lsysPDE::UpdateResidual() {
  int ierr;
  ierr = MatMult(KK,EPSC,RESC);
  CHKERRQ(ierr);
  ierr = VecAXPBY(RES,-1.,1.,RESC);
  CHKERRQ(ierr);

  return 1;
}

//--------------------------------------------------------------------------------
int lsysPDE::SumEpsToSol(const vector <unsigned> &MGIndex) {

  PetscScalar* R;
  PetscScalar* E;
  PetscScalar zero=0.;
  int ierr;
  PetscScalar value;
  
  Vec RESloc;
  Vec EPSloc;

  
  if(_msh->_nprocs==1) {
    ierr = VecGetArray(RES,&R);
    CHKERRQ(ierr);
    ierr = VecGetArray(EPS,&E);
    CHKERRQ(ierr);
  } 
  else {
    ierr=VecGhostGetLocalForm(RES,&RESloc);
    CHKERRQ(ierr);
    ierr=VecGhostGetLocalForm(EPS,&EPSloc);
    CHKERRQ(ierr);
    ierr = VecGetArray(RESloc,&R);
    CHKERRQ(ierr);
    ierr = VecGetArray(EPSloc,&E);
    CHKERRQ(ierr);
  }
  
  for (unsigned k=0; k<MGIndex.size(); k++) {
    unsigned indexSol=MGIndex[k];
    unsigned soltype =  SolType[indexSol];
//     cout << "indexSol" << indexSol << endl;
//     
     int loc_size   = Eps_[indexSol]->local_size();
     int loc_offset_EPS = KKoffset[k][_msh->_iproc] - KKoffset[0][_msh->_iproc];
//      int offset_eps = Eps_[indexSol]->first_local_index();
     int glob_offset_eps = _msh->MetisOffset[soltype][_msh->_iproc];
     
//    if(_iproc==1) printf("iproc %d  local_size: %d    offset_EPS: %d    offset_eps: %d \n",_iproc,loc_size, offset_EPS, offset_eps);
//      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank %d   %d \n",_iproc,own_size[soltype][_iproc]);
//      PetscSynchronizedFlush(PETSC_COMM_WORLD);
     //ottimizzare utilizzando il puntatore 
     //porre il residuo uguale a zero direttamente nell'assembly
     for(int i=0; i<_msh->own_size[soltype][_msh->_iproc]; i++) {
       Eps_[indexSol]->set(i+glob_offset_eps,E[loc_offset_EPS+i]);
       if ((*Bdc_[indexSol])(i+glob_offset_eps)>1.1) Res_[indexSol]->set(i+glob_offset_eps,R[loc_offset_EPS+i]);
       else Res_[indexSol]->set(i+glob_offset_eps,zero);
    }
    
//     for (PetscInt i=KKIndex[k]; i<KKIndex[k+1]; i++) {
//       PetscInt inode=i-KKIndex[k];
//       PetscInt kkdof=GetKKDof(indexSol,k,inode);
//       PetscInt inode_Metis=GetMetisDof(inode,SolType[indexSol]);
//       if (Bdc[indexSol][inode]>1) {
//         Res_[indexSol]->set(inode_Metis,R[0][kkdof]);
//       } 
//       else {
//         Res_[indexSol]->set(inode_Metis,zero);
//       }
//       Eps_[indexSol]->set(inode_Metis,E[0][kkdof]);
//     }
// 
//     //close TODO
    Res_[indexSol]->close();
    Eps_[indexSol]->close();
  }

  if(_msh->_nprocs==1) {
    ierr = VecRestoreArray(RES,&R);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(EPS,&E);
    CHKERRQ(ierr);
  } else {
    ierr = VecRestoreArray(RESloc,&R);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(EPSloc,&E);
    CHKERRQ(ierr);
    ierr=VecGhostRestoreLocalForm(RES,&RESloc);
    CHKERRQ(ierr);
    ierr=VecGhostRestoreLocalForm(EPS,&EPSloc);
    CHKERRQ(ierr);
  }


  for (unsigned k=0; k<MGIndex.size(); k++) {
    unsigned indexSol=MGIndex[k];
    //Adding Eps to Solu for every variable
    Sol_[indexSol]->add(*Eps_[indexSol]);
//     cout << "pluto " << endl;
    Sol_[indexSol]->close();
  }
  
 

  return ierr;
}

//-------------------------------------------------------------------------------------------
int lsysPDE::DeallocateMatrix(const vector <unsigned> &MGIndex) {

  int ierr;
  ierr=MatDestroy(&KK);
  CHKERRQ(ierr);

  if (_msh->GetGridNumber()>0) {
    ierr=MatDestroy(&PP);
    CHKERRQ(ierr);
  }

  ierr=VecDestroy(&RES);
  CHKERRQ(ierr);
  ierr=VecDestroy(&RESC);
  CHKERRQ(ierr);

  ierr=VecDestroy(&EPSC);
  CHKERRQ(ierr);
  ierr=VecDestroy(&EPS);
  CHKERRQ(ierr);

  return ierr;
}

//-------------------------------------------------------------------------------------------
int lsysPDE::AllocateMatrix(const vector <unsigned> &MGIndex) {

  PetscInt KKsize=KKIndex[KKIndex.size()-1u];
  PetscErrorCode ierr;
  const unsigned dim = _msh->GetDimension();
  PetscInt KK_UNIT_SIZE_ = pow(5,dim);

  if(_msh->_nprocs==1) {
    //TODO optimize memory
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,KKsize,KKsize,KK_UNIT_SIZE_*KKIndex.size(),PETSC_NULL,&KK);
    CHKERRQ(ierr);
    
    ierr = MatSetFromOptions(KK);
    CHKERRQ(ierr);
    
    ierr=MatSetOption(KK,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE); 
    CHKERRQ(ierr);
    
  }
  else {

    PetscInt KK_local_size =KKoffset[KKIndex.size()-1][_msh->_iproc] - KKoffset[0][_msh->_iproc];
    
//     printf(" a %d  %d  %d ", _iproc, KK_local_size, KKsize );
    
    ierr = MatCreate(MPI_COMM_WORLD, &KK);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = MatSetSizes(KK, KK_local_size, KK_local_size, KKsize,KKsize);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = MatSetType(KK, MATMPIAIJ);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = MatMPIAIJSetPreallocation(KK, KK_UNIT_SIZE_*KKIndex.size(), PETSC_NULL, KK_UNIT_SIZE_*KKIndex.size(), PETSC_NULL);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  
    ierr = MatSetFromOptions(KK);
    CHKERRQ(ierr);
    
    ierr=MatSetOption(KK,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE); 
    CHKERRQ(ierr);
    
  }

  return ierr;
}

//-------------------------------------------------------------------------------------------
void lsysPDE::FreeSolutionVectors() {
  for (unsigned i=0; i<Sol_.size(); i++) {
    //Delete structures
    if(ResEpsBdc_flag_[i]){
      delete Res_[i];
      delete Eps_[i];
      delete Bdc_[i];
    }
    delete Sol_[i];

    //
    if (SolTmOrder[i]==2) {
      delete Sol_old_[i];
    }
  }
  
  for (unsigned i=0; i<5; i++) 
    if (Proj_mat_flag[i]) {
      delete Proj_mat[i];
  }
}

//-------------------------------------------------------------------------------------------
void lsysPDE::UpdateSolution() {

  for (unsigned i=0; i<Sol_.size(); i++) {

// cout << soldiscr[i]<<endl;
// Copy the old vector
    if (SolTmOrder[i]==2) {
      *(Sol_old_[i]) = *(Sol_[i]);
    }
  }

}
