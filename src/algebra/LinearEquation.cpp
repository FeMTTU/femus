/*=========================================================================

 Program: FEMUS
 Module: LinearEquation
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

#include <ctime>
#include <fstream>
#include <algorithm>
#include "LinearEquation.hpp"
#include "ElemType.hpp"
#include "ParalleltypeEnum.hpp"
#include "NumericVector.hpp"
#include "SparseMatrix.hpp"


namespace femus {



using std::cout;
using std::endl;

//--------------------------------------------------------------------------------
LinearEquation::LinearEquation(mesh *other_msh){    
  _msh = other_msh;
  _CC_flag=0;
  _EPS = NULL;
  _EPSC = NULL;
  _RES = NULL;
  _RESC = NULL;
  _PP = NULL;
  _RR = NULL;
  _KK = NULL;
  _CC = NULL;
}

//--------------------------------------------------------------------------------
LinearEquation::~LinearEquation() { }


//--------------------------------------------------------------------------------
unsigned LinearEquation::GetIndex(const char name[]) {
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


unsigned LinearEquation::GetKKDof(const unsigned &index_sol, const unsigned &kkindex_sol, 
				  const unsigned &idof_gmt) const {
  
   unsigned soltype =  _SolType[index_sol]; 
   unsigned isubdom = (soltype<3)?_msh->npart[idof_gmt]:(_msh->epart[idof_gmt % _msh->GetElementNumber()]);
   unsigned idof_metis = _msh->GetMetisDof(idof_gmt,soltype);   
   return KKoffset[kkindex_sol][isubdom] + idof_metis - _msh->MetisOffset[soltype][isubdom];
}

//--------------------------------------------------------------------------------
void LinearEquation::InitPde(const vector <unsigned> &SolPdeIndex_other, const  vector <int> &SolType_other,  
		     const vector <char*> &SolName_other, vector <NumericVector*> *Bdc_other, 
		     const unsigned &other_gridr, const unsigned &other_gridn) {
  _SolPdeIndex=SolPdeIndex_other;
  _gridr=other_gridr;
  _gridn=other_gridn;
  
  _SolType=SolType_other;
  _SolName=SolName_other;
  _Bdc=Bdc_other;
  
  int ierr;
  KKIndex.resize(_SolPdeIndex.size()+1u);
  KKIndex[0]=0;
  for (unsigned i=1; i<KKIndex.size(); i++)
  KKIndex[i]=KKIndex[i-1]+_msh->MetisOffset[_SolType[_SolPdeIndex[i-1]]][n_processors()];

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
  //KKghost_nd[0].resize(1);  KKghost_nd[0][0]=1;
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
  int EPSsize= KKIndex[KKIndex.size()-1];
  _EPS = NumericVector::build().release();
  if(n_processors()==1) { // IF SERIAL
    _EPS->init(EPSsize,EPSsize,false,SERIAL);
  } 
  else { // IF PARALLEL
    int EPS_local_size =KKoffset[KKIndex.size()-1][processor_id()] - KKoffset[0][processor_id()];
    _EPS->init(EPSsize,EPS_local_size, KKghost_nd[processor_id()], false,GHOSTED);
  }
    
  _RES = NumericVector::build().release();
  _RES->init(*_EPS);
  
  _EPSC = NumericVector::build().release();
  _EPSC->init(*_EPS);
  
  _RESC = NumericVector::build().release();
  _RESC->init(*_EPS);
  
  
  GetSparsityPatternSize();
  
  const unsigned dim = _msh->GetDimension();
  int KK_UNIT_SIZE_ = pow(5,dim);
  int KK_size=KKIndex[KKIndex.size()-1u];
  int KK_local_size =KKoffset[KKIndex.size()-1][processor_id()] - KKoffset[0][processor_id()];
    
 _KK = SparseMatrix::build().release();
 _KK->init(KK_size,KK_size,KK_local_size,KK_local_size,KK_UNIT_SIZE_*KKIndex.size(),KK_UNIT_SIZE_*KKIndex.size());
   
  unsigned igrid=_msh->GetGridNumber()+1;
  if(igrid>=_gridr && igrid<_gridn){
    _CC = SparseMatrix::build().release();
    _CC->init(KK_size,KK_size,KK_local_size,KK_local_size,KK_UNIT_SIZE_*KKIndex.size(),KK_UNIT_SIZE_*KKIndex.size());
  }

}

//--------------------------------------------------------------------------------
void LinearEquation::AddLevel(){
  _gridn++;
  if(!_CC){
    const unsigned dim = _msh->GetDimension();
    int KK_UNIT_SIZE_ = pow(5,dim);
    int KK_size=KKIndex[KKIndex.size()-1u];
    int KK_local_size =KKoffset[KKIndex.size()-1][processor_id()] - KKoffset[0][processor_id()];
   
    unsigned igrid=_msh->GetGridNumber()+1;
    if(igrid>=_gridr && igrid<_gridn){
      _CC = SparseMatrix::build().release();
      _CC->init(KK_size,KK_size,KK_local_size,KK_local_size,KK_UNIT_SIZE_*KKIndex.size(),KK_UNIT_SIZE_*KKIndex.size());
    }
  }
}

//--------------------------------------------------------------------------------
void LinearEquation::SetResZero() {
  _RES->zero();
}

//--------------------------------------------------------------------------------
void LinearEquation::SetEpsZero() {
  _EPS->zero();
  _EPSC->zero();
}

//--------------------------------------------------------------------------------
void LinearEquation::SumEpsCToEps() {
  *_EPS += *_EPSC;
}

//--------------------------------------------------------------------------------
void LinearEquation::UpdateResidual() {
  _RESC->matrix_mult(*_EPSC,*_KK);
  *_RES -= *_RESC;
}

//-------------------------------------------------------------------------------------------
void LinearEquation::DeletePde() {
  
  if(_KK)
    delete _KK;
  
  unsigned igrid=_msh->GetGridNumber()+1;
  if(igrid>=_gridr && igrid<_gridn){
    if(_CC)
      delete _CC;
  }
  
  
  if (_msh->GetGridNumber()>0) {
     if(_PP) 
       delete _PP;
     if(_RR) 
       delete _RR;
  }
  
  if(_EPS)
    delete _EPS;
  
  if(_EPSC)
    delete _EPSC;
  
  if(_RES)
    delete _RES;
  
  if(_RESC)
    delete _RESC;
  
}
#include "vector"
# include "map"
  void LinearEquation::GetSparsityPatternSize() {
    
    const unsigned dim = _msh->GetDimension();
   
    const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
   
    vector < vector < unsigned > > dofsVAR(_SolPdeIndex.size());
    vector < unsigned > end_ind(_SolPdeIndex.size());
       
    for(int i=0;i<_SolPdeIndex.size();i++){
      dofsVAR[i].reserve(max_size);
      end_ind[i] = _msh->GetEndIndex(_SolType[_SolPdeIndex[i]]);
    }
    
    // mesh and procs
    unsigned nel    = _msh->GetElementNumber();
    unsigned igrid  = _msh->GetGridNumber();
    unsigned iproc  = _msh->processor_id();
    int nprocs=	      _msh->n_processors();
     
    // *** element loop ***
    
    unsigned IndexStart= KKoffset[0][iproc];
    unsigned IndexEnd  = KKoffset[KKIndex.size()-1][iproc];
    unsigned Offset    = IndexEnd-IndexStart;
    
    
    
    int KK_UNIT_SIZE_ = pow(5,dim);
    
    vector < std::map < int, bool > > BlgToMe_BlgToMe(Offset);
    std::map < int, std::map <int, bool > > NnBlgToMe_BlgToMe;
        
    cout<<IndexStart<<" "<<IndexEnd<<" "<<Offset<<" "<<BlgToMe_BlgToMe[0].max_size()<<endl;
    
    for(int iel=_msh->IS_Mts2Gmt_elem_offset[iproc]; iel < _msh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

      unsigned kel        = _msh->IS_Mts2Gmt_elem[iel]; 
      short unsigned kelt = _msh->el->GetElementType(kel);
      vector < unsigned > nve(_SolPdeIndex.size());
      for(int i=0;i<_SolPdeIndex.size();i++){
	nve[i] = _msh->el->GetElementDofNumber(kel,end_ind[i]);
      }
            
      for(int i=0; i<_SolPdeIndex.size(); i++) {
	dofsVAR[i].resize(nve[i]);
      }
      
      for(int i=0; i<_SolPdeIndex.size(); i++) {
	unsigned ThisSolType=_SolType[_SolPdeIndex[i]];
	for (unsigned j=0;j<nve[i];j++) {
	  unsigned inode=(ThisSolType<3)?(_msh->el->GetElementVertexIndex(kel,j)-1u):(kel+j*nel);
	  unsigned inode_Metis=_msh->GetMetisDof(inode,ThisSolType);
	  dofsVAR[i][j]= GetKKDof(i,i,inode); 
	  //std::cout<<dofsVAR[i][j]<<" ";
	}
      }
      //std::cout<<std::endl;
      for(int i=0;i<_SolPdeIndex.size();i++){
	for(int inode=0;inode<nve[i];inode++){
	  unsigned idof_local=dofsVAR[i][inode] - IndexStart;
	  for(int j=0;j<_SolPdeIndex.size();j++){
	    for(int jnode=0;jnode<nve[j];jnode++){
	      unsigned jdof_local=dofsVAR[j][jnode] - IndexStart;
	      if(jdof_local >= 0 && jdof_local < Offset){ //jnode belongs to iproc 
		if(idof_local >= 0 && idof_local < Offset){ // inode belogs to iporc		  
 		  BlgToMe_BlgToMe[ idof_local ][ jdof_local] = 1;
		}
		else{
		  NnBlgToMe_BlgToMe[dofsVAR[i][inode]][jdof_local]=1.;
		}
	      }
	    } 
	  }
	}
      }     
    }
    
    
    
   NumericVector  *_sizeNBM = NumericVector::build().release();
   _sizeNBM->init(*_EPS);
   _sizeNBM->zero();
    
    for (std::map < int, std::map <int, bool > >::iterator it=NnBlgToMe_BlgToMe.begin(); it!=NnBlgToMe_BlgToMe.end(); ++it){
      _sizeNBM->set(it->first,it->second.size());
    }
    _sizeNBM->close();
    
    
    
    vector < int > d_nnz(Offset);
    vector < int > o_nnz(Offset);
    
        
    
    for(int i=0; i<Offset;i++){
     d_nnz[i]=BlgToMe_BlgToMe[i].size();
     o_nnz[i]=(*_sizeNBM)(IndexStart+i);
    }
    
    delete _sizeNBM;
  }
} 




