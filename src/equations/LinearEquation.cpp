//C++ include
#include <ctime>
#include <fstream>
#include <algorithm>
#include "LinearEquation.hpp"
#include "ElemType.hpp"
#include "ParalleltypeEnum.hpp"
#include "NumericVector.hpp"
#include "SparseMatrix.hpp"

using std::cout;
using std::endl;

//--------------------------------------------------------------------------------
lsysPde::lsysPde(mesh *other_msh){    
  _msh = other_msh;
  _is_symmetric = false;
  _stabilization = false;
  _compressibility = 0.;
  _CC_flag=0;
  _DirichletBCsHandlingMode = 0;
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
lsysPde::~lsysPde() { }

//--------------------------------------------------------------------------------
void lsysPde::set_dirichletBCsHandling(unsigned int DirichletBCsHandlingMode) {
  _DirichletBCsHandlingMode = DirichletBCsHandlingMode;
}


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
int lsysPde::InitPde(const vector <unsigned> &SolPdeIndex_other, const  vector <int> &SolType_other,  
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
  
  const unsigned dim = _msh->GetDimension();
  int KK_UNIT_SIZE_ = pow(5,dim);
  int KK_size=KKIndex[KKIndex.size()-1u];
  int KK_local_size =KKoffset[KKIndex.size()-1][_msh->_iproc] - KKoffset[0][_msh->_iproc];
    
 _KK = SparseMatrix::build().release();
 _KK->init(KK_size,KK_size,KK_local_size,KK_local_size,KK_UNIT_SIZE_*KKIndex.size(),KK_UNIT_SIZE_*KKIndex.size());
   
  unsigned igrid=_msh->GetGridNumber()+1;
  if(igrid>=_gridr && igrid<_gridn){
    _CC = SparseMatrix::build().release();
    _CC->init(KK_size,KK_size,KK_local_size,KK_local_size,KK_UNIT_SIZE_*KKIndex.size(),KK_UNIT_SIZE_*KKIndex.size());
  }
  return 1;
}

//--------------------------------------------------------------------------------
void lsysPde::SetResZero() {
  _RES->zero();
}

//--------------------------------------------------------------------------------
void lsysPde::SetEpsZero() {
  _EPS->zero();
  _EPSC->zero();
}

//--------------------------------------------------------------------------------
void lsysPde::SumEpsCToEps() {
  *_EPS += *_EPSC;
}

//--------------------------------------------------------------------------------
void lsysPde::UpdateResidual() {
  _RESC->matrix_mult(*_EPSC,*_KK);
  *_RES -= *_RESC;
}

//-------------------------------------------------------------------------------------------
void lsysPde::DeletePde() {
  
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
//-------------------------------------------------------------------------------------------


