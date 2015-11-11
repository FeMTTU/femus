/*=========================================================================

 Program: FEMUS
 Module: Mesh
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
#include "Mesh.hpp"
#include "MeshGeneration.hpp"
#include "MeshMetisPartitioning.hpp"
#include "GambitIO.hpp"
#include "SalomeIO.hpp"
#include "NumericVector.hpp"

// C++ includes
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <algorithm>


namespace femus {

using std::cout;
using std::endl;
using std::min;
using std::sort;
using std::map;

bool Mesh::_IsUserRefinementFunctionDefined = false;

unsigned Mesh::_dimension=2;
unsigned Mesh::_ref_index=4;  // 8*DIM[2]+4*DIM[1]+2*DIM[0];
unsigned Mesh::_face_index=2; // 4*DIM[2]+2*DIM[1]+1*DIM[0];

//------------------------------------------------------------------------------------------------------
  Mesh::Mesh(){

    _coarseMsh = NULL;

    for(int i=0;i<5;i++){
      _ProjCoarseToFine[i]=NULL;
    }

    for (int itype=0; itype<3; itype++) {
      for (int jtype=0; jtype<3; jtype++) {
	_ProjQitoQj[itype][jtype] = NULL;
      }
    }
  }


  Mesh::~Mesh(){
    delete el;
    _coordinate->FreeSolutionVectors();
    delete _coordinate;

    for (int itype=0; itype<3; itype++) {
      for (int jtype=0; jtype<3; jtype++) {
	if(_ProjQitoQj[itype][jtype]){
	  delete _ProjQitoQj[itype][jtype];
	  _ProjQitoQj[itype][jtype] = NULL;
	}
      }
    }

    for (unsigned i=0; i<5; i++) {
      if (_ProjCoarseToFine[i]) {
	delete _ProjCoarseToFine[i];
	_ProjCoarseToFine[i]=NULL;
      }
    }
  }

/// print Mesh info
void Mesh::PrintInfo() {

 std::cout << " Mesh Level        : " << _level  << std::endl;
 std::cout << " Number of elements: " << _nelem  << std::endl;
 std::cout << " Number of nodes   : " << _nnodes << std::endl;

}

/**
 *  This function generates the coarse Mesh level, $l_0$, from an input Mesh file (Now only the Gambit Neutral File)
 **/
void Mesh::ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool> &type_elem_flag) {

  _coords.resize(3);

  _level = 0;

  if(name.rfind(".neu") < name.size())
  {
    GambitIO(*this).read(name,_coords,Lref,type_elem_flag);
  }
  else if(name.rfind(".med") < name.size()) {
    SalomeIO(*this).read(name,_coords,Lref,type_elem_flag);
  }
  else
  {
    std::cerr << " ERROR: Unrecognized file extension: " << name
	      << "\n   I understand the following:\n\n"
	      << "     *.neu -- Gambit Neutral File\n"
              << std::endl;
  }

  el->SetNodeNumber(_nnodes);

  std::vector < int > epart(GetNumberOfElements());
  MeshMetisPartitioning meshmetispartitioning(*this);
  meshmetispartitioning.DoPartition(epart, false);
  FillISvector(epart);
  epart.resize(0);


  BuildAdjVtx();
  Buildkel();

  _coordinate = new Solution(this);

  _coordinate->AddSolution("X",LAGRANGE,SECOND,1,0);
  _coordinate->AddSolution("Y",LAGRANGE,SECOND,1,0);
  _coordinate->AddSolution("Z",LAGRANGE,SECOND,1,0);

  _coordinate->ResizeSolutionVector("X");
  _coordinate->ResizeSolutionVector("Y");
  _coordinate->ResizeSolutionVector("Z");

  _coordinate->GetSolutionName("X") = _coords[0];
  _coordinate->GetSolutionName("Y") = _coords[1];
  _coordinate->GetSolutionName("Z") = _coords[2];

  _coordinate->AddSolution("AMR",DISCONTINOUS_POLYNOMIAL,ZERO,1,0);

  _coordinate->ResizeSolutionVector("AMR");

};

/**
 *  This function generates the coarse Box Mesh level using the built-in generator
 **/
void Mesh::GenerateCoarseBoxMesh(
        const unsigned int nx, const unsigned int ny, const unsigned int nz,
        const double xmin, const double xmax,
        const double ymin, const double ymax,
        const double zmin, const double zmax,
        const ElemType type, std::vector<bool> &type_elem_flag) {

  _coords.resize(3);

  _level = 0;

  MeshTools::Generation::BuildBox(*this,_coords,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,type,type_elem_flag);

  el->SetNodeNumber(_nnodes);

  std::vector < int > epart(GetNumberOfElements());
  MeshMetisPartitioning meshmetispartitioning(*this);
  meshmetispartitioning.DoPartition(epart, false);
  FillISvector(epart);
  epart.resize(0);

  BuildAdjVtx();

  Buildkel();

  _coordinate = new Solution(this);

  _coordinate->AddSolution("X",LAGRANGE,SECOND,1,0);
  _coordinate->AddSolution("Y",LAGRANGE,SECOND,1,0);
  _coordinate->AddSolution("Z",LAGRANGE,SECOND,1,0);

  _coordinate->ResizeSolutionVector("X");
  _coordinate->ResizeSolutionVector("Y");
  _coordinate->ResizeSolutionVector("Z");

  _coordinate->GetSolutionName("X") = _coords[0];
  _coordinate->GetSolutionName("Y") = _coords[1];
  _coordinate->GetSolutionName("Z") = _coords[2];

  _coordinate->AddSolution("AMR",DISCONTINOUS_POLYNOMIAL,ZERO,1,0);

  _coordinate->ResizeSolutionVector("AMR");

}


/**
 * This function searches all the elements around all the vertices
 **/
void Mesh::BuildAdjVtx() {
  el->AllocateVertexElementMemory();
  for (unsigned iel=0; iel<_nelem; iel++) {
    for (unsigned inode=0; inode < el->GetElementDofNumber(iel,0); inode++) {
      unsigned irow=el->GetElementVertexIndex(iel,inode)-1u;
      unsigned jcol=0;
      while ( 0 != el->GetVertexElementIndex(irow,jcol) ) jcol++;
      el->SetVertexElementIndex(irow,jcol,iel+1u);
    }
  }
}

/**
 * This function stores the element adiacent to the element face (iel,iface)
 * and stores it in kel[iel][iface]
 **/
void Mesh::Buildkel() {
  for (unsigned iel=0; iel<el->GetElementNumber(); iel++) {
    for (unsigned iface=0; iface<el->GetElementFaceNumber(iel); iface++) {
      if ( el->GetFaceElementIndex(iel,iface) <= 0) {//TODO probably just == -1
        unsigned i1=el->GetFaceVertexIndex(iel,iface,0);
        unsigned i2=el->GetFaceVertexIndex(iel,iface,1);
        unsigned i3=el->GetFaceVertexIndex(iel,iface,2);
        for (unsigned j=0; j< el->GetVertexElementNumber(i1-1u); j++) {
          unsigned jel= el->GetVertexElementIndex(i1-1u,j)-1u;
          if (jel > iel) {
            for (unsigned jface=0; jface<el->GetElementFaceNumber(jel); jface++) {
              if ( el->GetFaceElementIndex(jel,jface) <= 0) {
                unsigned j1=el->GetFaceVertexIndex(jel,jface,0);
                unsigned j2=el->GetFaceVertexIndex(jel,jface,1);
                unsigned j3=el->GetFaceVertexIndex(jel,jface,2);
                unsigned j4=el->GetFaceVertexIndex(jel,jface,3);
                if ((Mesh::_dimension==3 &&
                     (i1==j1 || i1==j2 || i1==j3 ||  i1==j4 )&&
                     (i2==j1 || i2==j2 || i2==j3 ||  i2==j4 )&&
                     (i3==j1 || i3==j2 || i3==j3 ||  i3==j4 ))||
                    (Mesh::_dimension==2 &&
                     (i1==j1 || i1==j2 )&&
                     (i2==j1 || i2==j2 ))||
                    (Mesh::_dimension==1 &&
                     (i1==j1))
                   ) {
                  el->SetFaceElementIndex(iel,iface,jel+1u);
                  el->SetFaceElementIndex(jel,jface,iel+1u);
                }
              }
            }
          }
        }
      }
    }
  }
}


void Mesh::AllocateAndMarkStructureNode() {
  el->AllocateNodeRegion();
  for (unsigned iel=0; iel<_nelem; iel++) {

    int flag_mat = el->GetElementMaterial(iel);

    if (flag_mat==4) {
      unsigned nve=el->GetElementDofNumber(iel,2);
      for ( unsigned i=0; i<nve; i++) {
        unsigned inode=el->GetElementVertexIndex(iel,i)-1u;
        el->SetNodeRegion(inode, 1);
      }
    }
  }
}


void Mesh::SetFiniteElementPtr(const elem_type * OtherFiniteElement[6][5]){
  for(int i=0;i<6;i++)
    for(int j=0;j<5;j++)
      _finiteElement[i][j] = OtherFiniteElement[i][j];
}


  // *******************************************************

//dof map: piecewise liner 0, quadratic 1, bi-quadratic 2, piecewise constant 3, piecewise linear discontinuous 4

void Mesh::FillISvector(vector < int > &epart) {

  //BEGIN Initialization for k = 0,1,2,3,4
  //resize the vector IS_Gmt2Mts_dof and dof
  for(int k=0;k<5;k++) {
    IS_Gmt2Mts_dof_offset[k].resize(_nprocs+1);
  }
  IS_Mts2Gmt_elem.resize(GetNumberOfElements());
  IS_Mts2Gmt_elem_offset.resize(_nprocs+1);
  IS_Mts2Gmt_elem_offset[0] = 0;

  //END Initialization for k = 0,1,2,3,4

  //BEGIN building the  metis2Gambit_elem and  k = 3,4
  unsigned counter = 0;
  for(int isdom = 0; isdom < _nprocs; isdom++) { // isdom = iprocess
    for(unsigned iel = 0; iel < GetNumberOfElements(); iel++){
      if( epart[iel] == isdom ){
	//filling the Metis to Mesh element mapping
	IS_Mts2Gmt_elem[ counter ] = iel;
        counter++;
	IS_Mts2Gmt_elem_offset[isdom + 1] = counter;
      }
    }
  }
  epart.resize(0);

  if( GetLevel() == 0 ){
    el->ReorderMeshElements(IS_Mts2Gmt_elem, NULL);
  }
  else{
    el->ReorderMeshElements(IS_Mts2Gmt_elem, _coarseMsh->el);
  }

  IS_Gmt2Mts_dof[3].assign(GetNumberOfElements(), 0);
  IS_Gmt2Mts_dof[4].assign(GetNumberOfElements() * (_dimension + 1), 0);

  for(int isdom = 0; isdom < _nprocs; isdom++){
    unsigned localSize = IS_Mts2Gmt_elem_offset[isdom+1] - IS_Mts2Gmt_elem_offset[isdom];
    unsigned offsetPWLD = IS_Mts2Gmt_elem_offset[isdom] * (_dimension + 1);
    for(unsigned iel = IS_Mts2Gmt_elem_offset[isdom]; iel < IS_Mts2Gmt_elem_offset[isdom+1]; iel++){
      IS_Mts2Gmt_elem[iel] = iel;
      IS_Gmt2Mts_dof[3][iel] = iel;
      //piecewise linear discontinuous
      unsigned locIel = iel - IS_Mts2Gmt_elem_offset[isdom];
      for(unsigned k = 0; k < _dimension + 1; k++){
        unsigned locKel = ( k * localSize ) + locIel;
        unsigned kel = offsetPWLD + locKel;
        IS_Gmt2Mts_dof[4][iel + k * GetNumberOfElements()] =  kel;
      }
    }
  }

  // ghost vs owned nodes: 3 and 4 have no ghost nodes
  for(unsigned k = 3; k < 5; k++){
    own_size[k].assign(_nprocs,0);
  }

  for(int isdom = 0; isdom < _nprocs; isdom++){
    own_size[3][isdom] = IS_Mts2Gmt_elem_offset[isdom+1] - IS_Mts2Gmt_elem_offset[isdom];
    own_size[4][isdom] = (IS_Mts2Gmt_elem_offset[isdom+1] - IS_Mts2Gmt_elem_offset[isdom])*(_dimension+1);
  }

  for(int k = 3; k < 5; k++) {
    ghost_nd[k].resize(_nprocs);
    ghost_nd_mts[k].resize(_nprocs);
    for(int isdom = 0; isdom < _nprocs; isdom++) {
      ghost_nd[k][isdom].resize( 0 );
      ghost_nd_mts[k][isdom].resize( 0 );
    }
  }
  //END building the  metis2Gambit_elem and  k = 3,4

  //BEGIN building for k = 0,1,2
  vector < unsigned > npart;
  npart.reserve(GetNumberOfNodes());

  MetisOffset.resize(5);

  // Initialization for k = 0,1,2
  npart.assign(GetNumberOfNodes(),_nprocs);
  IS_Gmt2Mts_dof[2].assign(GetNumberOfNodes(), 0);

  for( unsigned k = 0; k < 3; k++){
    own_size[k].assign(_nprocs,0);
  }
  counter = 0;
  for(int isdom = 0; isdom < _nprocs; isdom++){
    for( unsigned k = 0; k < 3; k++){
      for( unsigned iel = IS_Mts2Gmt_elem_offset[isdom]; iel < IS_Mts2Gmt_elem_offset[isdom+1]; iel++){
	unsigned nodeStart = (k == 0) ? 0 : el->GetElementDofNumber(iel,k-1);
	unsigned nodeEnd = el->GetElementDofNumber(iel,k);
	for ( unsigned inode = nodeStart; inode < nodeEnd; inode++) {
	  unsigned ii = el->GetElementVertexIndex(iel,inode) - 1;
	  if(npart[ii] > isdom) {
	    npart[ii] = isdom;
	    IS_Gmt2Mts_dof[2][ii] = counter;
	    counter++;
	    for( int j = k; j < 3; j++){
	      own_size[j][isdom]++;
	    }
	  }
	}
      }
    }
  }

  npart.resize(0);

  MetisOffset[2].resize(_nprocs+1);
  MetisOffset[2][0]=0;
  for(int i = 1 ;i <= _nprocs; i++){
    MetisOffset[2][i]= MetisOffset[2][i-1] + own_size[2][i-1];
  }

  el->ReorderMeshNodes( IS_Gmt2Mts_dof[2]);

  if( GetLevel() == 0 ){
    vector <double> coords_temp;
    for(int i = 0;i < 3; i++){
      coords_temp = _coords[i];
        for(unsigned j = 0; j < GetNumberOfNodes(); j++) {
	  _coords[i][IS_Gmt2Mts_dof[2][j]] = coords_temp[j];
      }
    }
  }

  for(unsigned j=0;j<GetNumberOfNodes();j++) {
    IS_Gmt2Mts_dof[2][j]=j;
  }
  //END building for k = 2, but incomplete for k = 0, 1

  //BEGIN ghost nodes search k = 0, 1, 2
  for(int k = 0; k < 3; k++){
    ghost_nd[k].resize(_nprocs);
    ghost_nd_mts[k].resize(_nprocs);
    for(int isdom = 0; isdom < _nprocs; isdom++){
      std::map < unsigned, bool > ghostMap;
      for(unsigned iel = IS_Mts2Gmt_elem_offset[isdom]; iel < IS_Mts2Gmt_elem_offset[isdom+1]; iel++){
	for (unsigned inode = 0; inode < el->GetElementDofNumber(iel,k); inode++) {
	  unsigned ii = el->GetElementVertexIndex(iel,inode)-1;
	  if(ii < MetisOffset[2][isdom]){
	    ghostMap[ii] = true;
	  }
	}
      }
      ghost_nd[k][isdom].resize( ghostMap.size() );
      ghost_nd_mts[k][isdom].resize( ghostMap.size() );
      unsigned counter = 0;
      for( std::map < unsigned, bool >::iterator it = ghostMap.begin(); it != ghostMap.end(); it++ ){
	ghost_nd[k][isdom][counter] = it->first;
	ghost_nd_mts[k][isdom][counter] = it->first;
	counter++;
      }
    }
  }
  //END ghost nodes search k = 0, 1, 2


  //BEGIN completing k = 0, 1
  
  for(unsigned k = 0; k < 2; k++){

    IS_Gmt2Mts_dof[k].assign(GetNumberOfNodes(), GetNumberOfNodes());
    std::vector < unsigned > ownedGhostCounter( _nprocs , 0);
    unsigned counter = 0;

    
    for(int isdom = 0; isdom < _nprocs; isdom++){

      //owned nodes
      for(unsigned inode = MetisOffset[2][isdom]; inode < own_size[k][isdom] + MetisOffset[2][isdom]; inode++) {
	IS_Gmt2Mts_dof[k][inode] = counter;
	counter++;
      }

      //ghost nodes
      ghost_nd_mts[k][isdom].reserve(ghost_nd[k][isdom].size());
      ghost_nd_mts[k][isdom].resize(0);
      
      for (unsigned inode = 0; inode < ghost_nd[k][isdom].size(); inode++){
	unsigned ghostNode = ghost_nd[k][isdom][inode];

// 	unsigned isdom0 = 0;
// 	unsigned isdom1 = isdom ;
// 	unsigned ksdom  = isdom /2;
// 	while( ghostNode < MetisOffset[2][ksdom] || ghostNode >= MetisOffset[2][ksdom + 1] ){
// 	  if( ghostNode < MetisOffset[2][ksdom] ) isdom1 = ksdom;
// 	  else isdom0 = ksdom + 1;
// 	  ksdom = ( isdom0 + isdom1 ) / 2;
// 	}

	unsigned ksdom = IsdomBisectionSearch(ghostNode, 2);
	
	int upperBound = MetisOffset[2][ksdom] + own_size[k][ksdom];
	if( ghostNode < upperBound || IS_Gmt2Mts_dof[k][ ghostNode ] != GetNumberOfNodes()){ // real ghost nodes
	  unsigned ghostSize = ghost_nd_mts[k][isdom].size();
	  ghost_nd_mts[k][isdom].resize(ghostSize + 1);
	  ghost_nd_mts[k][isdom][ghostSize] = IS_Gmt2Mts_dof[k][ ghostNode ];
	}
	else { // owned ghost nodes
	  IS_Gmt2Mts_dof[k][ ghostNode ] = counter;
	  _ownedGhostMap[k][ ghostNode ] = counter;
	  counter++;
	  ownedGhostCounter[isdom]++;

          for(unsigned jnode = inode; jnode < ghost_nd[k][isdom].size()-1; jnode++ ){
	    ghost_nd[k][isdom][jnode] = ghost_nd[k][isdom][jnode + 1];
	  }

          ghost_nd[k][isdom].resize(ghost_nd[k][isdom].size()-1);
	  inode--;
	}
      }
    }

    _originalOwnSize[k].resize(_nprocs);
    
    for(int isdom = 0; isdom < _nprocs; isdom++){
      _originalOwnSize[k][isdom] = own_size[k][isdom];
      own_size[k][isdom] += ownedGhostCounter[isdom];
    }
  }
  //END completing for k = 0, 1


  //BEGIN Initilize and set all the Offsets
  MetisOffset.resize(5);
  for(int k=0;k<5;k++) {
    MetisOffset[k].resize(_nprocs+1);
    MetisOffset[k][0]=0;
    for(int i = 1 ;i <= _nprocs; i++){
      MetisOffset[k][i]= MetisOffset[k][i-1] + own_size[k][i-1];
    }
  }
  //END Initilize and set all the Offsets

}


  // *******************************************************
  unsigned Mesh::IsdomBisectionSearch(const unsigned &dof, const short unsigned &solType) const{

    unsigned isdom0 = 0;
    unsigned isdom1 = _nprocs ;
    unsigned isdom = _iproc;

    while( dof < MetisOffset[solType][isdom] || dof >= MetisOffset[solType][isdom + 1] ){
      if( dof < MetisOffset[solType][isdom] ) isdom1 = isdom;
      else isdom0 = isdom + 1;
      isdom = ( isdom0 + isdom1 ) / 2;
    }
  
    return isdom;
  }
  // *******************************************************

  unsigned Mesh::GetMetisDof(const unsigned &i, const unsigned &iel, const short unsigned &solType) const {

    unsigned dof;
    switch(solType){
      case 0: // linear Lagrange
      {
        unsigned iNode = el->GetMeshDof(iel, i, solType);
	unsigned isdom = IsdomBisectionSearch(iNode, 2);
	if(iNode < MetisOffset[2][isdom]+_originalOwnSize[0][isdom]){
 	  dof = (iNode - MetisOffset[2][isdom]) + MetisOffset[0][isdom];
	}
	else{
	  dof = _ownedGhostMap[0].find(iNode)->second;
	}
      }
	break;
      case 1: // quadratic Lagrange
       	{
        unsigned iNode = el->GetMeshDof(iel, i, solType);
	unsigned isdom = IsdomBisectionSearch(iNode, 2);
	if(iNode < MetisOffset[2][isdom]+_originalOwnSize[1][isdom]){
 	  dof = (iNode - MetisOffset[2][isdom]) + MetisOffset[1][isdom];
	}
	else{
	  dof = _ownedGhostMap[1].find(iNode)->second;
	}
      }
	break;
	
      case 2: // bi-quadratic Lagrange
        dof = el->GetMeshDof(iel, i, solType);
        break;
      case 3: // piecewise constant
        dof = iel;
        break;
      case 4: // piecewise linear discontinuous
	unsigned isdom = IsdomBisectionSearch(iel, 3);
	unsigned offset = IS_Mts2Gmt_elem_offset[isdom];
        unsigned offsetp1 = IS_Mts2Gmt_elem_offset[isdom + 1];
        unsigned ownSize = offsetp1 - offset;
        unsigned offsetPWLD = offset * (_dimension + 1);
        unsigned locIel = iel - offset;
        dof = offsetPWLD + ( i * ownSize ) + locIel;
        break;
      }
    return dof;
  }

  // *******************************************************

SparseMatrix* Mesh::GetQitoQjProjection(const unsigned& itype, const unsigned& jtype) {
  if(itype < 3 && jtype < 3){
    if(!_ProjQitoQj[itype][jtype]){
      BuildQitoQjProjection(itype, jtype);
    }
  }
  else{
    std::cout<<"Wrong argument range in function"
	     <<"Mesh::GetLagrangeProjectionMatrix(const unsigned& itype, const unsigned& jtype)"<<std::cout;
    abort();
  }
  return _ProjQitoQj[itype][jtype];
}

void Mesh::BuildQitoQjProjection(const unsigned& itype, const unsigned& jtype){

  unsigned ni = MetisOffset[itype][_nprocs];
  unsigned ni_loc = own_size[itype][_iproc];

  unsigned nj = MetisOffset[jtype][_nprocs];
  unsigned nj_loc = own_size[itype][_iproc];

  NumericVector *NNZ_d = NumericVector::build().release();
  if(1 == _nprocs) { // IF SERIAL
    NNZ_d->init(ni, ni_loc, false, SERIAL);
  }
  else{
    NNZ_d->init(ni, ni_loc, ghost_nd_mts[itype][processor_id()], false, GHOSTED);
  }
  NNZ_d->zero();

  NumericVector *NNZ_o = NumericVector::build().release();
  NNZ_o->init(*NNZ_d);
  NNZ_o->zero();

  for(unsigned isdom = _iproc; isdom < _iproc+1; isdom++) {
    for (unsigned iel_mts = IS_Mts2Gmt_elem_offset[isdom]; iel_mts < IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++){
      unsigned iel = IS_Mts2Gmt_elem[iel_mts];
      short unsigned ielt = el->GetElementType(iel);
      _finiteElement[ielt][jtype]->GetSparsityPatternSize(*this, iel, NNZ_d, NNZ_o, itype);
    }
  }

  NNZ_d->close();
  NNZ_o->close();

  unsigned offset = MetisOffset[itype][_iproc];

  vector < int > nnz_d(ni_loc);
  vector < int > nnz_o(ni_loc);
  for(unsigned i = 0; i < ni_loc; i++){
    nnz_d[i] = static_cast < int > ((*NNZ_d)(offset+i));
    nnz_o[i] = static_cast < int > ((*NNZ_o)(offset+i));
  }

  _ProjQitoQj[itype][jtype] = SparseMatrix::build().release();
  _ProjQitoQj[itype][jtype]->init(ni, nj, own_size[itype][_iproc], own_size[jtype][_iproc], nnz_d, nnz_o);
  for(unsigned isdom = _iproc; isdom < _iproc+1; isdom++) {
    for (unsigned iel_mts = IS_Mts2Gmt_elem_offset[isdom]; iel_mts < IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++){
      unsigned iel = IS_Mts2Gmt_elem[iel_mts];
      short unsigned ielt = el->GetElementType(iel);
      _finiteElement[ielt][jtype]->BuildProlongation(*this, iel, _ProjQitoQj[itype][jtype], NNZ_d, NNZ_o,itype);
    }
  }
  _ProjQitoQj[itype][jtype]->close();

  delete NNZ_d;
  delete NNZ_o;
}



SparseMatrix* Mesh::GetCoarseToFineProjection(const unsigned& solType){

  if( solType >= 5 ){
    std::cout<<"Wrong argument range in function \"GetCoarseToFineProjection\": "
	     <<"solType is greater then SolTypeMax"<<std::endl;
    abort();
  }

  if(!_ProjCoarseToFine[solType])
    BuildCoarseToFineProjection(solType);

  return _ProjCoarseToFine[solType];
}



void Mesh::BuildCoarseToFineProjection(const unsigned& solType){

  if (!_coarseMsh) {
    std::cout<<"Error! In function \"BuildCoarseToFineProjection\": the coarse mesh has not been set"<<std::endl;
    abort();
  }

  if( !_ProjCoarseToFine[solType] ){

    int nf     = MetisOffset[solType][_nprocs];
    int nc     = _coarseMsh->MetisOffset[solType][_nprocs];
    int nf_loc = own_size[solType][_iproc];
    int nc_loc = _coarseMsh->own_size[solType][_iproc];

    //build matrix sparsity pattern size
    NumericVector *NNZ_d = NumericVector::build().release();
    if(n_processors()==1) { // IF SERIAL
      NNZ_d->init(nf, nf_loc, false, SERIAL);
    }
    else { // IF PARALLEL
      if(solType<3) { // GHOST nodes only for Lagrange FE families
	NNZ_d->init(nf, nf_loc, ghost_nd_mts[solType][processor_id()], false, GHOSTED);
      }
      else { //piecewise discontinuous variables have no ghost nodes
	NNZ_d->init(nf, nf_loc, false, PARALLEL);
      }
    }
    NNZ_d->zero();

    NumericVector *NNZ_o = NumericVector::build().release();
    NNZ_o->init(*NNZ_d);
    NNZ_o->zero();

    for(int isdom=_iproc; isdom<_iproc+1; isdom++) {
      for (int iel_mts=_coarseMsh->IS_Mts2Gmt_elem_offset[isdom];iel_mts < _coarseMsh->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
	unsigned iel = _coarseMsh->IS_Mts2Gmt_elem[iel_mts];
	short unsigned ielt=_coarseMsh->el->GetElementType(iel);
	_finiteElement[ielt][solType]->GetSparsityPatternSize( *this, *_coarseMsh, iel, NNZ_d, NNZ_o);
      }
    }
    NNZ_d->close();
    NNZ_o->close();

    unsigned offset = MetisOffset[solType][_iproc];
    vector <int> nnz_d(nf_loc);
    vector <int> nnz_o(nf_loc);
    for(int i=0; i<nf_loc;i++){
      nnz_d[i]=static_cast <int> (floor((*NNZ_d)(offset+i)+0.5));
      nnz_o[i]=static_cast <int> (floor((*NNZ_o)(offset+i)+0.5));
    }
    delete NNZ_d;
    delete NNZ_o;

    //build matrix
    _ProjCoarseToFine[solType] = SparseMatrix::build().release();
    _ProjCoarseToFine[solType]->init(nf,nc,nf_loc,nc_loc,nnz_d,nnz_o);

    // loop on the coarse grid
    for(int isdom=_iproc; isdom<_iproc+1; isdom++) {
      for (int iel_mts=_coarseMsh->IS_Mts2Gmt_elem_offset[isdom];
	   iel_mts < _coarseMsh->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
	unsigned iel = _coarseMsh->IS_Mts2Gmt_elem[iel_mts];
	short unsigned ielt=_coarseMsh->el->GetElementType(iel);
	_finiteElement[ielt][solType]->BuildProlongation(*this, *_coarseMsh,iel, _ProjCoarseToFine[solType]);
      }
    }
    _ProjCoarseToFine[solType]->close();
  }
}








} //end namespace femus


