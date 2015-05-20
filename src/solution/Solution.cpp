/*=========================================================================

 Program: FEMUS
 Module: Solution
 Authors: Eugenio Aulisa

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
#include "Solution.hpp"
#include "FemusDefault.hpp"
#include "ElemType.hpp"
#include "ParalleltypeEnum.hpp"
#include "NumericVector.hpp"


namespace femus {

using std::cout;
using std::endl;

/**
 *  Contructor
 **/
// ------------------------------------------------------------------
Solution::Solution(Mesh *other_msh){
 _msh = other_msh;
  for(int i=0;i<5;i++){
    _GradMat[i].resize(_msh->GetDimension());
    _AMR_flag=0;
  }
}

/**
 * Destructor
 **/
// ------------------------------------------------------------------
Solution::~Solution() {
  for (unsigned i=0; i<_SolName.size(); i++) {
    delete [] _SolName[i];
  }
}

/**
 * Add a new variable called 'name'
 */
void Solution::AddSolution( const char name[], const FEFamily fefamily, const FEOrder order,
			    const unsigned& tmorder, const bool &Pde_type) {

  unsigned n=_Sol.size();

  _SolType.resize(n+1u);
  _SolName.resize(n+1u);
  _SolTmOrder.resize(n+1u);
  _family.resize(n+1u);
  _order.resize(n+1u);

  _Sol.resize(n+1u);
  _Res.resize(n+1u);
  _Eps.resize(n+1u);

  _GradVec.resize(n+1u);
  _GradVec[n].resize(_msh->GetDimension());

  _Bdc.resize(n+1u);
  _ResEpsBdcFlag.resize(n+1u);

  _ResEpsBdcFlag[n]=Pde_type;
  _family[n] = fefamily;
  _order[n] = order;
  _SolType[n] = order - ((fefamily==LAGRANGE)?1:0) + fefamily*3;
  _SolTmOrder[n]=tmorder;
  _SolOld.resize(n+1u);
  _SolName[n]=new char [DEFAULT_SOL_NCHARS];
  strcpy(_SolName[n],name);

}

/**
 * Get the solution index for the variable called name
 **/
//-------------------------------------------------------------------
unsigned Solution::GetIndex(const char name[]) const {
  unsigned index=0;
  while (strcmp(_SolName[index],name)) {
    index++;
    if (index==_Res.size()) {
      cout<<"error! invalid name entry GetIndex(...)"<<endl;
      exit(0);
    }
  }
  return index;
}

/**
 * Allocate memory for the variable called name
 **/
// ------------------------------------------------------------------
void Solution::ResizeSolutionVector(const char name[]) {

  unsigned i=GetIndex(name);

  _Sol[i] = NumericVector::build().release();
  if(n_processors()==1) { // IF SERIAL
    _Sol[i]->init(_msh->MetisOffset[_SolType[i]][n_processors()],_msh->own_size[_SolType[i]][processor_id()],false,SERIAL);
  }
  else { // IF PARALLEL
    if(_SolType[i]<3) {
      if(_msh->ghost_size[_SolType[i]][processor_id()]!=0) {
	_Sol[i]->init(_msh->MetisOffset[_SolType[i]][n_processors()],_msh->own_size[_SolType[i]][processor_id()],
		      _msh->ghost_nd_mts[_SolType[i]][processor_id()],false,GHOSTED);
      }
      else {
	std::vector <int> fake_ghost(1,_msh->own_size[_SolType[i]][processor_id()]);
	_Sol[i]->init(_msh->MetisOffset[_SolType[i]][n_processors()],_msh->own_size[_SolType[i]][processor_id()],
		      fake_ghost,false,GHOSTED);
      }
    }
    else { //discontinuous pressure has no ghost nodes
      _Sol[i]->init(_msh->MetisOffset[_SolType[i]][n_processors()],_msh->own_size[_SolType[i]][processor_id()],false,PARALLEL);
    }
  }

  if (_SolTmOrder[i]==2) { // only if the variable is time dependent
    _SolOld[i] = NumericVector::build().release();
    _SolOld[i]->init(*_Sol[i]);
  }

  if(_ResEpsBdcFlag[i]) { //only if the variable is a Pde type

    _Res[i] = NumericVector::build().release();
    _Res[i]->init(*_Sol[i]);

    _Eps[i] = NumericVector::build().release();
    _Eps[i]->init(*_Sol[i]);

    _Bdc[i] = NumericVector::build().release();
    _Bdc[i]->init(*_Sol[i]);
  }
}


/** Init and set to zero The AMR Eps vector */
void Solution::InitAMREps(){
  _AMR_flag=1;
  _AMREps.resize(_Sol.size());
  for(int i=0;i<_Sol.size();i++){
    _AMREps[i] = NumericVector::build().release();
    _AMREps[i]->init(*_Sol[i]);
    _AMREps[i]->zero();
  }

}



/**
 * Deallocate memory for all variables
 **/
// ------------------------------------------------------------------
void Solution::FreeSolutionVectors() {
  for (unsigned i=0; i<_Sol.size(); i++) {
    delete _Sol[i];
    if(_ResEpsBdcFlag[i]){
      delete _Res[i];
      delete _Eps[i];
      delete _Bdc[i];
    }
    if (_SolTmOrder[i]==2) {
      delete _SolOld[i];
    }

    for(int j=0;j<_msh->GetDimension();j++){
      if(_GradVec[i][j]){
	delete _GradVec[i][j];
      }
    }

    if(_AMR_flag)
      delete _AMREps[i];
  }
  for (unsigned i=0; i<5; i++) {
    for(int j=0;j<_msh->GetDimension();j++){
      if(_GradMat[i][j]){
	delete _GradMat[i][j];
      }
    }
  }
}

/**
 * Update _Sol, _Res and _Eps based on EPS and RES
 **/
//--------------------------------------------------------------------------------
void Solution::UpdateSolAndRes(const vector <unsigned> &_SolPdeIndex,  NumericVector* _EPS,  NumericVector* _RES,
			       const vector <vector <unsigned> > &KKoffset){

  PetscScalar zero=0.;
  for (unsigned k=0; k<_SolPdeIndex.size(); k++) {
    unsigned indexSol=_SolPdeIndex[k];
    unsigned soltype =  _SolType[indexSol];

    int loc_offset_EPS = KKoffset[k][processor_id()];

    int glob_offset_eps = _msh->MetisOffset[soltype][processor_id()];

    vector <int> index(_msh->own_size[soltype][processor_id()]);
    for(int i=0; i<_msh->own_size[soltype][processor_id()]; i++) {
      index[i]=loc_offset_EPS+i;
    }
    vector <double> valueEPS(_msh->own_size[soltype][processor_id()]);
    _EPS->get(index,valueEPS);
    vector <double> valueRES(_msh->own_size[soltype][processor_id()]);
    _RES->get(index,valueRES);

    for(int i=0; i<_msh->own_size[soltype][processor_id()]; i++) {
      _Eps[indexSol]->set(i+glob_offset_eps,valueEPS[i]);
      if ((*_Bdc[indexSol])(i+glob_offset_eps)>1.1) _Res[indexSol]->set(i+glob_offset_eps,valueRES[i]);
      else _Res[indexSol]->set(i+glob_offset_eps,zero);
    }
    _Res[indexSol]->close();
    _Eps[indexSol]->close();
  }

  for (unsigned k=0; k<_SolPdeIndex.size(); k++) {
    unsigned indexSol=_SolPdeIndex[k];
    _Sol[indexSol]->add(*_Eps[indexSol]);
    _Sol[indexSol]->close();
    if(_AMR_flag){
      _AMREps[indexSol]->add(*_Eps[indexSol]);
      _AMREps[indexSol]->close();
    }

  }

}

/**
 * Update _Sol
 **/

void Solution::UpdateSol(const vector <unsigned> &_SolPdeIndex,  NumericVector* _EPS, const vector <vector <unsigned> > &KKoffset){

  PetscScalar zero=0.;
  for (unsigned k=0; k<_SolPdeIndex.size(); k++) {
    unsigned indexSol=_SolPdeIndex[k];
    unsigned soltype =  _SolType[indexSol];

    int loc_offset_EPS = KKoffset[k][processor_id()];

    int glob_offset_eps = _msh->MetisOffset[soltype][processor_id()];

    vector <int> index(_msh->own_size[soltype][processor_id()]);
    for(int i=0; i<_msh->own_size[soltype][processor_id()]; i++) {
      index[i]=loc_offset_EPS+i;
    }
    vector <double> valueEPS(_msh->own_size[soltype][processor_id()]);
    _EPS->get(index,valueEPS);
    //vector <double> valueRES(_msh->own_size[soltype][processor_id()]);
    //_RES->get(index,valueRES);

    for(int i=0; i<_msh->own_size[soltype][processor_id()]; i++) {
      _Eps[indexSol]->set(i+glob_offset_eps,valueEPS[i]);
      //if ((*_Bdc[indexSol])(i+glob_offset_eps)>1.1) _Res[indexSol]->set(i+glob_offset_eps,valueRES[i]);
      //else _Res[indexSol]->set(i+glob_offset_eps,zero);
    }
    //_Res[indexSol]->close();
    _Eps[indexSol]->close();
  }

  for (unsigned k=0; k<_SolPdeIndex.size(); k++) {
    unsigned indexSol=_SolPdeIndex[k];
    _Sol[indexSol]->add(*_Eps[indexSol]);
    _Sol[indexSol]->close();
    if(_AMR_flag){
      _AMREps[indexSol]->add(*_Eps[indexSol]);
      _AMREps[indexSol]->close();
    }

  }

}


/**
 * Update _Res
 **/
//--------------------------------------------------------------------------------
void Solution::UpdateRes(const vector <unsigned> &_SolPdeIndex, NumericVector* _RES, const vector <vector <unsigned> > &KKoffset) {

  PetscScalar zero=0.;
  for (unsigned k=0; k<_SolPdeIndex.size(); k++) {
    unsigned indexSol=_SolPdeIndex[k];
    unsigned soltype =  _SolType[indexSol];

    int loc_offset_RES = KKoffset[k][processor_id()];

    int glob_offset_res = _msh->MetisOffset[soltype][processor_id()];

    vector <int> index(_msh->own_size[soltype][processor_id()]);
    for(int i=0; i<_msh->own_size[soltype][processor_id()]; i++) {
      index[i]=loc_offset_RES+i;
    }

    vector <double> valueRES(_msh->own_size[soltype][processor_id()]);
    _RES->get(index,valueRES);

    for(int i=0; i<_msh->own_size[soltype][processor_id()]; i++) {
      if ((*_Bdc[indexSol])(i+glob_offset_res)>1.1) {
	_Res[indexSol]->set(i+glob_offset_res,valueRES[i]);
      }
      else {
	_Res[indexSol]->set(i+glob_offset_res,zero);
      }
    }
    _Res[indexSol]->close();
  }

}

bool Solution::FlagAMRRegionBasedOnl2(const vector <unsigned> &SolIndex,const double &AMRthreshold){

  vector <double> SolMax(SolIndex.size());
  vector <unsigned> SolType(SolIndex.size());
  vector <unsigned> SolEndInd(SolIndex.size());

  unsigned END_IND[5]= {0,1,1,4,5};

  for (unsigned k=0; k<SolIndex.size(); k++) {
    double EPSMAX = _AMREps[SolIndex[k]]->linfty_norm ();
    double SOLMAX = _Sol[SolIndex[k]]->linfty_norm ();
    cout << "Current maximum relative change = " <<EPSMAX/SOLMAX << endl << endl;
    SolMax[k] = AMRthreshold * SOLMAX;
    SolType[k] = _SolType[SolIndex[k]];
    SolEndInd[k]   = END_IND[SolType[k]];
  }

  Solution* AMR = _msh->_coordinate;
  unsigned  AMRIndex= AMR->GetIndex("AMR");
  AMR->_Sol[AMRIndex]->zero();

  unsigned nel= _msh->GetNumberOfElements();

  NumericVector *counter_vec;
  counter_vec = NumericVector::build().release();

  if(_nprocs==1) {
    counter_vec->init(_nprocs,1,false,SERIAL);
  }
  else {
    counter_vec->init(_nprocs,1,false,PARALLEL);
  }
  counter_vec->zero();

  for (int iel_metis=_msh->IS_Mts2Gmt_elem_offset[_iproc]; iel_metis < _msh->IS_Mts2Gmt_elem_offset[_iproc+1]; iel_metis++) {

    unsigned kel = _msh->IS_Mts2Gmt_elem[iel_metis];
    short unsigned kelt=_msh->el->GetElementType(kel);

    for (unsigned k=0; k<SolIndex.size(); k++) {

      if(SolType[k]<3){

      unsigned nve=_msh->el->GetElementDofNumber(kel,SolEndInd[k]);
      for(unsigned i=0; i<nve; i++) {
	unsigned inode=(SolType[k]<3)?(_msh->el->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
	unsigned inode_metis=_msh->GetMetisDof(inode,SolType[k]);
	double value = (*_AMREps[SolIndex[k]])(inode_metis);
	if(fabs(value)>SolMax[k]){
	  counter_vec->add(_iproc,1.);
	  AMR->_Sol[AMRIndex]->set(iel_metis,1.);
	  k=SolIndex.size();
	  i=nve;
	}
      }
    }
    }
  }
  AMR->_Sol[AMRIndex]->close();

  counter_vec->close();
  double counter=counter_vec->l1_norm();
  bool test=(counter<=_nprocs)?1:0;

  delete counter_vec;
  return test;

}



bool Solution::FlagAMRRegionBasedOnSemiNorm(const vector <unsigned> &SolIndex,const double &AMRthreshold){

  vector <double> GradSolMax(SolIndex.size());
  vector <unsigned> SolType(SolIndex.size());
  vector <unsigned> SolEndInd(SolIndex.size());
  unsigned dim=_msh->GetDimension();

  unsigned nel= _msh->GetNumberOfElements();

  for (unsigned k=0; k<SolIndex.size(); k++) {

    if(SolType[k]<3){

      SolType[k] = _SolType[SolIndex[k]];

      BuildGradMatrixStructure(SolType[k]);

      GradSolMax[k]=0.;
      for(int i=0;i<dim;i++){

	if(_GradVec[SolIndex[k]][i]==0){
	  _GradVec[SolIndex[k]][i] = NumericVector::build().release();
	  if(n_processors()==1) { // IF SERIAL
	    _GradVec[SolIndex[k]][i]->init(_msh->MetisOffset[3][n_processors()],_msh->own_size[3][processor_id()],false,SERIAL);
	  }
	  else { //discontinuous pressure has no ghost nodes
	    _GradVec[SolIndex[k]][i]->init(_msh->MetisOffset[3][n_processors()],_msh->own_size[3][processor_id()],false,PARALLEL);

	  }
	}

	_GradVec[SolIndex[k]][i]->matrix_mult(*_Sol[SolIndex[k]],*_GradMat[SolType[k]][i]);
	double GradSolMaxi=_GradVec[SolIndex[k]][i]->linfty_norm();
	GradSolMax[k] =GradSolMaxi*GradSolMaxi;
	_GradVec[SolIndex[k]][i]->close();
	_GradVec[SolIndex[k]][i]->matrix_mult(*_Eps[SolIndex[k]],*_GradMat[SolType[k]][i]);
      }
      GradSolMax[k]=AMRthreshold*(GradSolMax[k]);
    }
  }

  Solution* AMR = _msh->_coordinate;
  unsigned  AMRIndex= AMR->GetIndex("AMR");

  AMR->_Sol[AMRIndex]->zero();

  NumericVector *counter_vec;
  counter_vec = NumericVector::build().release();

  if(_nprocs==1) {
    counter_vec->init(_nprocs,1,false,SERIAL);
  }
  else {
    counter_vec->init(_nprocs,1,false,PARALLEL);
  }
  counter_vec->zero();

  for (int iel_metis=_msh->IS_Mts2Gmt_elem_offset[_iproc]; iel_metis < _msh->IS_Mts2Gmt_elem_offset[_iproc+1]; iel_metis++) {

    for (unsigned k=0; k<SolIndex.size(); k++) {

      if(SolType[k]<3){
	double value=0.;
	for(int i=0;i<dim;i++){
	  double valuei = (*_GradVec[SolIndex[k]][i])(iel_metis);
	  value+=valuei*valuei;
	}
	value=sqrt(value);
	if(fabs(value)>GradSolMax[k]){
	  counter_vec->add(_iproc,1.);
	  AMR->_Sol[AMRIndex]->set(iel_metis,1.);
	  k=SolIndex.size();
	}
      }
    }
  }
  AMR->_Sol[AMRIndex]->close();

  counter_vec->close();
  double counter=counter_vec->l1_norm();
  bool test=(counter<=_nprocs)?1:0;

  return test;
}


void Solution::BuildGradMatrixStructure(unsigned SolType) {

 if(SolType<3 && _GradMat[SolType][0]==0){

    unsigned dim=_msh->GetDimension();

    int nr     = _msh->MetisOffset[3][_nprocs];
    int nc     = _msh->MetisOffset[SolType][_nprocs];
    int nr_loc = _msh->own_size[3][_iproc];
    int nc_loc = _msh->own_size[SolType][_iproc];

    for(int i=0;i<dim;i++){
      _GradMat[SolType][i] = SparseMatrix::build().release();
      _GradMat[SolType][i]->init(nr,nc,nr_loc,nc_loc,27,27);
    }

//     // Begin build elem type structure
//     const elem_type *type_elem[6];
//     if(dim==3){
//       if(SolType==0){
// 	type_elem[0]=new const elem_type_3D("hex","linear","zero");
// 	type_elem[1]=new const elem_type_3D("tet","linear","zero");
// 	type_elem[2]=new const elem_type_3D("wedge","linear","zero");
//       }
//       else if(SolType==1){
// 	type_elem[0]=new const elem_type_3D("hex","quadratic","zero");
// 	type_elem[1]=new const elem_type_3D("tet","quadratic","zero");
// 	type_elem[2]=new const elem_type_3D("wedge","quadratic","zero");
//       }
//       else{
// 	type_elem[0]=new const elem_type_3D("hex","biquadratic","zero");
// 	type_elem[1]=new const elem_type_3D("tet","biquadratic","zero");
// 	type_elem[2]=new const elem_type_3D("wedge","biquadratic","zero");
//       }
//     }
//     else if(dim==2){
//       if(SolType==0){
// 	type_elem[3]=new const elem_type_2D("quad","linear","zero");
// 	type_elem[4]=new const elem_type_2D("tri","linear","zero");
//       }
//       else if(SolType==1){
// 	type_elem[3]=new const elem_type_2D("quad","quadratic","zero");
// 	type_elem[4]=new const elem_type_2D("tri","quadratic","zero");
//       }
//       else{
// 	type_elem[3]=new const elem_type_2D("quad","biquadratic","zero");
// 	type_elem[4]=new const elem_type_2D("tri","biquadratic","zero");
//       }
//     }
//     else if(dim==1){
//       if(SolType==0){
// 	type_elem[5]=new const elem_type_1D("line","linear","zero");
//       }
//       else if(SolType==1){
// 	type_elem[5]=new const elem_type_1D("line","quadratic","zero");
//       }
//       else{
// 	type_elem[5]=new const elem_type_1D("line","biquadratic","zero");
//       }
//     }

    vector< vector < double> > coordinates(dim);
    vector< int > column_dofs;
    vector< int > row_dof;
    vector <double> phi;
    vector <double> gradphi;
    vector <double> nablaphi;
    double weight;
    vector< vector< double> > B(dim);


    const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));

    for(int i=0; i<dim; i++)
      coordinates[i].reserve(max_size);

    row_dof.reserve(1);
    column_dofs.reserve(max_size);

    phi.reserve(max_size);
    gradphi.reserve(max_size*dim);
    nablaphi.reserve(max_size*(3*(dim-1)+!(dim-1)));
    for(int i=0;i<dim;i++){
      B[i].reserve(max_size);
    }
    // Set to zeto all the entries of the Global Matrix
    for(int i=0;i<dim;i++){
      _GradMat[SolType][i]->zero();
    }

    unsigned nel= _msh->GetNumberOfElements();

    for (int iel=_msh->IS_Mts2Gmt_elem_offset[_iproc]; iel < _msh->IS_Mts2Gmt_elem_offset[_iproc+1]; iel++) {

      row_dof.resize(1);
      row_dof[0]=iel;

      unsigned kel = _msh->IS_Mts2Gmt_elem[iel];
      short unsigned kelt=_msh->el->GetElementType(kel);

      unsigned nve=_msh->el->GetElementDofNumber(kel,SolType);

      // resize
      column_dofs.resize(nve);
      phi.resize(nve);
      gradphi.resize(nve*dim);
      nablaphi.resize(nve*(3*(dim-1)+!(dim-1)));
      for(int i=0; i<dim; i++) {
        coordinates[i].resize(nve);
      }

      // set to zero all the entries of the FE matrices
      for(int i=0;i<dim;i++){
	B[i].resize(nve);
	memset(&B[i][0],0,nve*sizeof(double));
      }

      for( unsigned i=0; i<nve; i++) {
	unsigned inode=_msh->el->GetElementVertexIndex(kel,i)-1u;
	unsigned inode_coord_metis=_msh->GetMetisDof(inode,2);
	column_dofs[i]=_msh->GetMetisDof(inode,SolType);
        for(unsigned ivar=0; ivar<dim; ivar++) {
          coordinates[ivar][i]=(*_msh->_coordinate->_Sol[ivar])(inode_coord_metis);
        }
      }
      _msh->_finiteElement[kelt][SolType]->Jacobian(coordinates,0,weight,phi,gradphi, nablaphi );


      for(int i=0;i<nve;i++){
	for(int j=0;j<dim;j++){
	  B[j][i]=gradphi[i*dim+j];
	}
      }

      for(int i=0;i<dim;i++){
        _GradMat[SolType][i]->add_matrix_blocked(B[i],row_dof,column_dofs);
      }
    }


    // End build elem type structure
    for(int i=0;i<dim;i++){
      _GradMat[SolType][i]->close();
    }

//     // Begin free elem type structures
//     if(dim==3){
//       delete type_elem[0];
//       delete type_elem[1];
//       delete type_elem[2];
//     }
//     else if(dim==2){
//       delete type_elem[3];
//       delete type_elem[4];
//     }
//     else if(dim==1){
//       delete type_elem[5];
//     }
//     // End free elem type structures


  }



}












// ------------------------------------------------------------------
void Solution::UpdateSolution() {
  for (unsigned i=0; i<_Sol.size(); i++) {
    // Copy the old vector
    if (_SolTmOrder[i]==2) {
      *(_SolOld[i]) = *(_Sol[i]);
    }
  }
}


} //end namespace femus


