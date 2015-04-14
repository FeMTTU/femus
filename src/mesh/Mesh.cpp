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

bool Mesh::_TestSetRefinementFlag=0;

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
  
 std::cout << " Mesh Level        : " << _level << std::endl; 
 std::cout << "   Number of elements: " << _nelem << std::endl; 
 std::cout << "   Number of nodes   : " << _nnodes << std::endl;
  
}

/**
 *  This function generates the coarse Mesh level, $l_0$, from an input Mesh file (Now only the Gambit Neutral File)
 **/
void Mesh::ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool> &type_elem_flag) {
    
  vector <vector <double> > coords(3);  
    
  _level = 0;

  if(name.rfind(".neu") < name.size())
  {
    GambitIO(*this).read(name,coords,Lref,type_elem_flag);
  }
  else if(name.rfind(".med") < name.size()) {
    SalomeIO(*this).read(name,coords,Lref,type_elem_flag);
  }
  else
  {
    std::cerr << " ERROR: Unrecognized file extension: " << name
	      << "\n   I understand the following:\n\n"
	      << "     *.neu -- Gambit Neutral File\n"
              << std::endl;
  }
  
  ReorderMeshDofs(coords);
  
  BuildAdjVtx();
  
  Buildkel();
  
  MeshMetisPartitioning meshmetispartitioning(*this);
  meshmetispartitioning.DoPartition();
  
  FillISvector();
  
  vector <double> coords_temp;
  for(int i=0;i<3;i++){
    coords_temp=coords[i];
    for(unsigned j=0;j<_nnodes;j++) {
      coords[i][GetMetisDof(j,2)]=coords_temp[j];
    }
  }
  
  _coordinate = new Solution(this);
 
  _coordinate->AddSolution("X",LAGRANGE,SECOND,1,0); 
  _coordinate->AddSolution("Y",LAGRANGE,SECOND,1,0); 
  _coordinate->AddSolution("Z",LAGRANGE,SECOND,1,0); 
  
  _coordinate->ResizeSolutionVector("X");
  _coordinate->ResizeSolutionVector("Y");
  _coordinate->ResizeSolutionVector("Z");
    
  _coordinate->GetSolutionName("X") = coords[0];
  _coordinate->GetSolutionName("Y") = coords[1];
  _coordinate->GetSolutionName("Z") = coords[2];
  
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
  
  vector <vector < double > > coords(3);  
    
  _level = 0;
    
  MeshTools::Generation::BuildBox(*this,coords,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,type,type_elem_flag);
  
  ReorderMeshDofs(coords);
  
  BuildAdjVtx();
  
  Buildkel();
  
  MeshMetisPartitioning meshmetispartitioning(*this);
  meshmetispartitioning.DoPartition();
  
  FillISvector();
  
  vector <double> coords_temp;
  for(int i=0;i<3;i++){
    coords_temp=coords[i];
    for(unsigned j=0;j<_nnodes;j++) {
      coords[i][GetMetisDof(j,2)]=coords_temp[j];
    }
  }
  
  _coordinate = new Solution(this);
 
  _coordinate->AddSolution("X",LAGRANGE,SECOND,1,0); 
  _coordinate->AddSolution("Y",LAGRANGE,SECOND,1,0); 
  _coordinate->AddSolution("Z",LAGRANGE,SECOND,1,0); 
  
  _coordinate->ResizeSolutionVector("X");
  _coordinate->ResizeSolutionVector("Y");
  _coordinate->ResizeSolutionVector("Z");
    
  _coordinate->GetSolutionName("X") = coords[0];
  _coordinate->GetSolutionName("Y") = coords[1];
  _coordinate->GetSolutionName("Z") = coords[2];
  
  _coordinate->AddSolution("AMR",DISCONTINOUS_POLYNOMIAL,ZERO,1,0); 
  
  _coordinate->ResizeSolutionVector("AMR");
  
}  
  

//------------------------------------------------------------------------------------------------------
void Mesh::ReorderMeshDofs(vector < vector < double> > &coords) {
  
  vector <unsigned> dof_index;
  dof_index.resize(_nnodes);
  for(unsigned i=0;i<_nnodes;i++){
    dof_index[i]=i+1;
  }
  //reorder vertices and mid-points vs central points
  for (unsigned iel=0; iel<_nelem; iel++) {
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,1); inode++) {
      for (unsigned jel=0; jel<_nelem; jel++) {
	for (unsigned jnode=el->GetElementDofNumber(jel,1); jnode<el->GetElementDofNumber(jel,2); jnode++) { 
	  unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	  unsigned jj=el->GetElementVertexIndex(jel,jnode)-1;
	  unsigned i0=dof_index[ii];
          unsigned i1=dof_index[jj];
	  if(i0>i1){
	    dof_index[ii]=i1;
	    dof_index[jj]=i0; 
	  }
	}
      }
    }
  }
  //reorder vertices vs mid-points
  for (unsigned iel=0; iel<_nelem; iel++) {
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
      for (unsigned jel=0; jel<_nelem; jel++) {
        for (unsigned jnode=el->GetElementDofNumber(jel,0); jnode<el->GetElementDofNumber(jel,1); jnode++) {
          unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	  unsigned jj=el->GetElementVertexIndex(jel,jnode)-1;
	  unsigned i0=dof_index[ii];
          unsigned i1=dof_index[jj];
	  if(i0>i1){
	    dof_index[ii]=i1;
	    dof_index[jj]=i0; 
	  }
	}
      }
    }
  }
  
  // update all
  for (unsigned iel=0; iel<_nelem; iel++) {
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,2); inode++) {
      unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
      el->SetElementVertexIndex(iel,inode,dof_index[ii]);
    }
  }
  vector <double> coords_temp;
  for(int i=0;i<3;i++) {
    coords_temp=coords[i];
    for(unsigned j=0;j<_nnodes;j++) {
      coords[i][dof_index[j]-1]=coords_temp[j];
    }
  }
  // **************  end reoreder mesh dofs **************
 
  el->SetNodeNumber(_nnodes);

  unsigned nv0=0;
  for (unsigned iel=0; iel<_nelem; iel++)
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
      unsigned i0=el->GetElementVertexIndex(iel,inode);
      if (nv0<i0) nv0=i0;
  }
  el->SetVertexNodeNumber(nv0);

  unsigned nv1=0;
  for (unsigned iel=0; iel<_nelem; iel++)
    for (unsigned inode=el->GetElementDofNumber(iel,0); inode<el->GetElementDofNumber(iel,1); inode++) {
      unsigned i1=el->GetElementVertexIndex(iel,inode);
      if (nv1<i1) nv1=i1;
  }
  el->SetMidpointNodeNumber(nv1-nv0);

  el->SetCentralNodeNumber(_nnodes-nv1);
  
}  

/**
 * This function searches all the elements around all the vertices
 **/
void Mesh::BuildAdjVtx() {
  el->AllocateVertexElementMemory();
  for (unsigned iel=0; iel<_nelem; iel++) {
    for (unsigned inode=0; inode < el->GetElementDofNumber(iel,0); inode++) {
      unsigned ii=el->GetElementVertexIndex(iel,inode)-1u;
      unsigned jj=0;
      while ( 0 != el->GetVertexElementIndex(ii,jj) ) jj++;
      el->SetVertexElementIndex(ii,jj,iel+1u);
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
      if ( el->GetFaceElementIndex(iel,iface) <= 0) {
        unsigned i1=el->GetFaceVertexIndex(iel,iface,0);
        unsigned i2=el->GetFaceVertexIndex(iel,iface,1);
        unsigned i3=el->GetFaceVertexIndex(iel,iface,2);
        for (unsigned j=0; j< el->GetVertexElementNumber(i1-1u); j++) {
          unsigned jel= el->GetVertexElementIndex(i1-1u,j)-1u;
          if (jel>iel) {
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


/**
 * This function returns the number of Mesh nodes for different type of elemets
 **/
unsigned Mesh::GetDofNumber(const unsigned type) const {

  switch (type) {
  case 0:
    return el->GetVertexNodeNumber();
    break;
  case 1:
    return el->GetVertexNodeNumber()+el->GetMidpointNodeNumber();
    break;
  case 2:
    return _nnodes;
    break;
  case 3:
    return _nelem;
    break;
  case 4:
    return _nelem*(2+Mesh::_dimension-1);
    break;
  }
  return 0;
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



void Mesh::FillISvector() {
  
   //dof map: piecewise liner 0, quadratic 1, biquadratic 2, piecewise constant 3, picewise discontinous linear 4 
  
  //resize the vector IS_Gmt2Mts_dof and dof
  for(int k=0;k<5;k++) {
    IS_Gmt2Mts_dof[k].resize(GetDofNumber(k)); 
    IS_Gmt2Mts_dof_offset[k].resize(nsubdom+1);
  }
  IS_Mts2Gmt_elem.resize(_nelem);
  IS_Mts2Gmt_elem_offset.resize(nsubdom+1);
  
  // I 
  for(unsigned i=0;i<_nnodes;i++) {
    npart[i] = nsubdom;
  }
  
  IS_Mts2Gmt_elem_offset[0] = 0;
  vector <unsigned> IS_Gmt2Mts_dof_counter(5,0);
  
   for(int k=0;k<5;k++) {
     IS_Gmt2Mts_dof[k].assign(GetDofNumber(k),GetDofNumber(k)-1); 
     //TODO for domain decomposition pourposes! the non existing dofs point to the last dof!!!!!!
   }
  
  IS_Gmt2Mts_dof_counter[3]=0;
  IS_Gmt2Mts_dof_counter[4]=0;
  
  for(int isdom=0;isdom<nsubdom;isdom++) {
    for(unsigned iel=0;iel<_nelem;iel++){
      if(epart[iel]==isdom){
	//filling the piecewise IS_Mts2Gmt_elem metis->gambit
	IS_Mts2Gmt_elem[ IS_Gmt2Mts_dof_counter[3] ] = iel;
	IS_Gmt2Mts_dof[3][iel] = IS_Gmt2Mts_dof_counter[3];
	IS_Gmt2Mts_dof_counter[3]++;
	IS_Mts2Gmt_elem_offset[isdom+1]=IS_Gmt2Mts_dof_counter[3];
	// linear+quadratic+biquadratic
	for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
	  unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	  if(npart[ii]>isdom) {
	    npart[ii]=isdom;
	    IS_Gmt2Mts_dof[0][ii]=IS_Gmt2Mts_dof_counter[0];
	    IS_Gmt2Mts_dof_counter[0]++;
	    IS_Gmt2Mts_dof[1][ii]=IS_Gmt2Mts_dof_counter[1];
	    IS_Gmt2Mts_dof_counter[1]++;
	    IS_Gmt2Mts_dof[2][ii]=IS_Gmt2Mts_dof_counter[2];
	    IS_Gmt2Mts_dof_counter[2]++;
	  }
	}
	// quadratic+biquadratic
	for (unsigned inode=el->GetElementDofNumber(iel,0); inode<el->GetElementDofNumber(iel,1); inode++) {
	  unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	  if(npart[ii]>isdom){
	    npart[ii]=isdom;
	    IS_Gmt2Mts_dof[1][ii]=IS_Gmt2Mts_dof_counter[1];
	    IS_Gmt2Mts_dof_counter[1]++;
	    IS_Gmt2Mts_dof[2][ii]=IS_Gmt2Mts_dof_counter[2];
	    IS_Gmt2Mts_dof_counter[2]++;
	  }
	}
	// biquadratic
	for (unsigned inode=el->GetElementDofNumber(iel,1); inode<el->GetElementDofNumber(iel,2); inode++) {
	  unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	  if(npart[ii]>isdom){
	    npart[ii]=isdom;
	    IS_Gmt2Mts_dof[2][ii]=IS_Gmt2Mts_dof_counter[2];
	    IS_Gmt2Mts_dof_counter[2]++;
	  }
	}	
      }
    }
    for(unsigned k_dim=0;k_dim<_dimension+1;k_dim++){
      for(unsigned iel=0;iel<_nelem;iel++){
     	if(epart[iel]==isdom){
	  IS_Gmt2Mts_dof[4][iel+k_dim*_nelem]=IS_Gmt2Mts_dof_counter[4];
	  IS_Gmt2Mts_dof_counter[4]++;
	}
      }
    }
  }  
   
  // ghost vs own nodes
  vector <unsigned short> node_count(_nnodes,0);
  
  for(unsigned k=0;k<5;k++){
    ghost_size[k].assign(nsubdom,0);
    own_size[k].assign(nsubdom,0);
  }
  
  for(int isdom=0;isdom<nsubdom;isdom++){
    
    own_size[3][isdom] = IS_Mts2Gmt_elem_offset[isdom+1]-IS_Mts2Gmt_elem_offset[isdom];
    own_size[4][isdom] = (IS_Mts2Gmt_elem_offset[isdom+1]-IS_Mts2Gmt_elem_offset[isdom])*(_dimension+1);
    
    for(unsigned i=IS_Mts2Gmt_elem_offset[isdom];i<IS_Mts2Gmt_elem_offset[isdom+1];i++){
      unsigned iel=IS_Mts2Gmt_elem[i]; 
      
       for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	if(node_count[ii]<isdom+1){
	  node_count[ii]=isdom+1;
	  if( npart[ii] != isdom){
	    ghost_size[0][isdom]++;
	    ghost_size[1][isdom]++;
	    ghost_size[2][isdom]++;
	  }else {
	    own_size[0][isdom]++;
	    own_size[1][isdom]++;
	    own_size[2][isdom]++;
	  }
	}
      }
      
      for (unsigned inode=el->GetElementDofNumber(iel,0); inode<el->GetElementDofNumber(iel,1); inode++) {
	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	if(node_count[ii]<isdom+1){
	  node_count[ii]=isdom+1;
	  if( npart[ii] != isdom){
	    ghost_size[1][isdom]++;
	    ghost_size[2][isdom]++;
	  }else {
	    own_size[1][isdom]++;
	    own_size[2][isdom]++;
	  }
	}
      }
      
      
      for (unsigned inode=el->GetElementDofNumber(iel,1); inode<el->GetElementDofNumber(iel,2); inode++) {
	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	if(node_count[ii]<isdom+1){
	  node_count[ii]=isdom+1;
	  if( npart[ii] != isdom){
	    ghost_size[2][isdom]++;
	  }else {
	    own_size[2][isdom]++;
	  }
	}
      }
    }
  }
  
  for(int k=0; k<5; k++) {
    ghost_nd[k].resize(nsubdom);
    ghost_nd_mts[k].resize(nsubdom);
    for(int isdom=0; isdom<nsubdom; isdom++) {
      ghost_nd[k][isdom].resize(ghost_size[k][isdom]);
      ghost_nd_mts[k][isdom].resize(ghost_size[k][isdom]);    
    }
  } 
  
  node_count.assign (_nnodes,0);  
  for(int k=0; k<5; k++) {
    ghost_size[k].assign(nsubdom,0);
  }
  
  for(int isdom=0;isdom<nsubdom;isdom++) {
    for(unsigned i=IS_Mts2Gmt_elem_offset[isdom];i<IS_Mts2Gmt_elem_offset[isdom+1];i++){
      unsigned iel=IS_Mts2Gmt_elem[i]; 
      
      
      for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	if(node_count[ii]<isdom+1){
	  node_count[ii]=isdom+1;
	  if(npart[ii] != isdom){
 	    ghost_nd_mts[0][isdom][ghost_size[0][isdom]]=IS_Gmt2Mts_dof[0][ii];
            ghost_nd[0][isdom][ghost_size[0][isdom]]=ii;
	    ghost_size[0][isdom]++;
 	    ghost_nd_mts[1][isdom][ghost_size[1][isdom]]=IS_Gmt2Mts_dof[1][ii];
	    ghost_nd[1][isdom][ghost_size[1][isdom]]=ii;
	    ghost_size[1][isdom]++;
 	    ghost_nd_mts[2][isdom][ghost_size[2][isdom]]=IS_Gmt2Mts_dof[2][ii];
	    ghost_nd[2][isdom][ghost_size[2][isdom]]=ii;
	    ghost_size[2][isdom]++;
	  }
	}
      }
      
      for (unsigned inode=el->GetElementDofNumber(iel,0); inode<el->GetElementDofNumber(iel,1); inode++) {
	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	if(node_count[ii]<isdom+1){
	  node_count[ii]=isdom+1;
	  if( npart[ii] != isdom){
 	    ghost_nd_mts[1][isdom][ghost_size[1][isdom]]=IS_Gmt2Mts_dof[1][ii];
	    ghost_nd[1][isdom][ghost_size[1][isdom]]=ii;
	    ghost_size[1][isdom]++;
	    ghost_nd_mts[2][isdom][ghost_size[2][isdom]]=IS_Gmt2Mts_dof[2][ii];
	    ghost_nd[2][isdom][ghost_size[2][isdom]]=ii;
	    ghost_size[2][isdom]++;
	  }
	}
      }
      
      for(unsigned inode=el->GetElementDofNumber(iel,1); inode<el->GetElementDofNumber(iel,2); inode++) {
	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	if(node_count[ii]<isdom+1){
	  node_count[ii]=isdom+1;
	  if( npart[ii] != isdom){
 	    ghost_nd_mts[2][isdom][ghost_size[2][isdom]]=IS_Gmt2Mts_dof[2][ii];
            ghost_nd[2][isdom][ghost_size[2][isdom]]=ii;
	    ghost_size[2][isdom]++;
	  }
	}
      } 
    }
  }
    
  MetisOffset.resize(5);
  for(int i=0;i<5;i++) 
    MetisOffset[i].resize(nsubdom+1);
  
  MetisOffset[0][0]=0;
  MetisOffset[1][0]=0;
  MetisOffset[2][0]=0;
  MetisOffset[3][0]=0;
  MetisOffset[4][0]=0;
  
  for(int i = 1 ;i <= nsubdom; i++){
    MetisOffset[0][i]= MetisOffset[0][i-1] + own_size[0][i-1];
    MetisOffset[1][i]= MetisOffset[1][i-1] + own_size[1][i-1];
    MetisOffset[2][i]= MetisOffset[2][i-1] + own_size[2][i-1];
    MetisOffset[3][i]= IS_Mts2Gmt_elem_offset[i];
    MetisOffset[4][i]= IS_Mts2Gmt_elem_offset[i]*(_dimension+1);
  } 
  
}

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
  NumericVector *NNZ_o = NumericVector::build().release();
  if(1 == _nprocs) { // IF SERIAL
    NNZ_d->init(ni, ni_loc, false, SERIAL);
    NNZ_o->init(ni, ni_loc, false, SERIAL);
  } 
  else{
    NNZ_d->init(ni, ni_loc, false, PARALLEL); 
    NNZ_o->init(ni, ni_loc, false, PARALLEL); 
  }
  NNZ_d->zero();
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
            
  delete NNZ_d;
  delete NNZ_o;
	
  _ProjQitoQj[itype][jtype] = SparseMatrix::build().release();
  _ProjQitoQj[itype][jtype]->init(ni, nj, own_size[itype][_iproc], own_size[jtype][_iproc], nnz_d, nnz_o);
  for(unsigned isdom = _iproc; isdom < _iproc+1; isdom++) {
    for (unsigned iel_mts = IS_Mts2Gmt_elem_offset[isdom]; iel_mts < IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++){
      unsigned iel = IS_Mts2Gmt_elem[iel_mts];
      short unsigned ielt = el->GetElementType(iel);
      _finiteElement[ielt][jtype]->BuildProlongation(*this, iel, _ProjQitoQj[itype][jtype], itype);	  
    }
  }
  _ProjQitoQj[itype][jtype]->close();
}



SparseMatrix* Mesh::GetCoarseToFineProjection(const unsigned& solType){
  
  if( solType > 4 ){
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
      if(solType<3) {
	if(ghost_size[solType][processor_id()]!=0) { 
	  NNZ_d->init(nf, nf_loc, ghost_nd_mts[solType][processor_id()], false, GHOSTED);
	} 
	else { 
	  std::vector < int > fake_ghost(1,nf_loc);
	  NNZ_d->init(nf, nf_loc, fake_ghost, false, GHOSTED);
	}
      }
      else { //discontinuous pressure has no ghost nodes
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
	if(_coarseMsh->el->GetRefinedElementIndex(iel)){ //only if the coarse element has been refined
	  short unsigned ielt=_coarseMsh->el->GetElementType(iel);
	  _finiteElement[ielt][solType]->GetSparsityPatternSize( *this, *_coarseMsh, iel, NNZ_d, NNZ_o); 
	}
      }
    }
    NNZ_d->close();
    NNZ_o->close();
    
    unsigned offset = MetisOffset[solType][_iproc];
    vector <int> nnz_d(nf_loc);
    vector <int> nnz_o(nf_loc);
    for(int i=0; i<nf_loc;i++){
      nnz_d[i]=static_cast <int> ((*NNZ_d)(offset+i));
      nnz_o[i]=static_cast <int> ((*NNZ_o)(offset+i));
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
	if(_coarseMsh->el->GetRefinedElementIndex(iel)){ //only if the coarse element has been refined
	  short unsigned ielt=_coarseMsh->el->GetElementType(iel);
	  _finiteElement[ielt][solType]->BuildProlongation(*this, *_coarseMsh,iel, _ProjCoarseToFine[solType]); 
	}
      }
    }
    _ProjCoarseToFine[solType]->close();
  }
}








} //end namespace femus


