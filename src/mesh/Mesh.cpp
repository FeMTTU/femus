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

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "mpi.h"
#include <algorithm>
#include "SparseMatrix.hpp"
#include "NumericVector.hpp"
#include "map"
#include "Mesh.hpp"
#include "metis.h"
#include "GambitIO.hpp"

namespace femus {
  
using std::cout;
using std::endl;
using std::min;
using std::sort;
using std::map;

bool Mesh::_TestSetRefinementFlag=0;

const unsigned Mesh::_END_IND[5]= {0,1,3,4,5};

unsigned Mesh::_dimension=2;
unsigned Mesh::_ref_index=4;  // 8*DIM[2]+4*DIM[1]+2*DIM[0];
unsigned Mesh::_face_index=2; // 4*DIM[2]+2*DIM[1]+1*DIM[0];


/**
 *  This function generates the coarse Mesh level, $l_0$, from an input Mesh file (Now only the Gambit Neutral File)
 **/
void Mesh::ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool> &type_elem_flag) {
    
  MPI_Comm_rank(MPI_COMM_WORLD, &_iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &_nprocs);
  
  vector <vector <double> > vt;  
  vt.resize(3);
    
  _grid=0;

  if(name.rfind(".neu") < name.size())
  {
    GambitIO(*this).read(name,vt,Lref,type_elem_flag);
  }
  else
  {
    std::cerr << " ERROR: Unrecognized file extension: " << name
	      << "\n   I understand the following:\n\n"
	      << "     *.neu -- Gambit Neutral File\n"
              << std::endl;
  }
  
  // connectivity: find all the element near the vertices
  BuildAdjVtx();
  Buildkel();
  
 if (_nprocs>=1) GenerateMetisMeshPartition();
  vector <double> vt_temp;
  for(int i=0;i<3;i++){
    vt_temp=vt[i];
    for(unsigned j=0;j<nvt;j++) {
      vt[i][GetMetisDof(j,2)]=vt_temp[j];
    }
  }
  
  _coordinate = new Solution(this);
  _coordinate->AddSolution("X",LAGRANGE,SECOND,1,0); 
  _coordinate->AddSolution("Y",LAGRANGE,SECOND,1,0); 
  _coordinate->AddSolution("Z",LAGRANGE,SECOND,1,0); 
  
  _coordinate->ResizeSolutionVector("X");
  _coordinate->ResizeSolutionVector("Y");
  _coordinate->ResizeSolutionVector("Z");
    
  _coordinate->SetCoarseCoordinates(vt);
  
  _coordinate->AddSolution("AMR",DISCONTINOUS_POLYNOMIAL,ZERO,1,0); 
  _coordinate->ResizeSolutionVector("AMR");
    
};


//-------------------------------------------------------------------
void Mesh::FlagAllElementsToBeRefined() {
  
   el->InitRefinedToZero();
   
   //refine all next grid elements
   for (unsigned iel=0; iel<nel; iel++) {
     el->SetRefinedElementIndex(iel,1);
     el->AddToRefinedElementNumber(1);
     short unsigned elt=el->GetElementType(iel);
     el->AddToRefinedElementNumber(1,elt);
   }

   el->AllocateChildrenElement(_ref_index);
}

//-------------------------------------------------------------------
void Mesh::FlagElementsToBeRefinedByUserDefinedFunction() {
     el->InitRefinedToZero();
   
    //refine based on the function SetRefinementFlag defined in the main;
    // the Mesh is serial, we cannot in parallel use the coordinates to selectively refine
    std::vector<double> X_local;
    std::vector<double> Y_local;
    std::vector<double> Z_local;
    _coordinate->_Sol[0]->localize_to_all(X_local);
    _coordinate->_Sol[1]->localize_to_all(Y_local);
    _coordinate->_Sol[2]->localize_to_all(Z_local);
  
    for (unsigned iel=0; iel<nel; iel+=1) {
      unsigned nve=el->GetElementDofNumber(iel,0);
      double vtx=0.,vty=0.,vtz=0.;
      for ( unsigned i=0; i<nve; i++) {
	unsigned inode=el->GetElementVertexIndex(iel,i)-1u;
	unsigned inode_Metis=GetMetisDof(inode,2);
	vtx+=X_local[inode_Metis];  
	vty+=Y_local[inode_Metis]; 
	vtz+=Z_local[inode_Metis]; 
      }
      vtx/=nve;
      vty/=nve;
      vtz/=nve;
      if(!el->GetRefinedElementIndex(iel)){
	if (_SetRefinementFlag(vtx,vty,vtz,el->GetElementGroup(iel),_grid)) {
	  el->SetRefinedElementIndex(iel,1);
	  el->AddToRefinedElementNumber(1);
	  short unsigned elt=el->GetElementType(iel);
	  el->AddToRefinedElementNumber(1,elt);
	}
      }
    }
    el->AllocateChildrenElement(_ref_index);
}



//-------------------------------------------------------------------
void Mesh::FlagElementsToBeRefinedByAMR() {
    
       
    if(_TestSetRefinementFlag){
      for (int iel_metis=IS_Mts2Gmt_elem_offset[_iproc]; iel_metis < IS_Mts2Gmt_elem_offset[_iproc+1]; iel_metis++) {
	unsigned kel = IS_Mts2Gmt_elem[iel_metis];
	short unsigned kelt=el->GetElementType(kel);
        unsigned nve=el->GetElementDofNumber(kel,0);
	double vtx=0.,vty=0.,vtz=0.;
	for(unsigned i=0; i<nve; i++) {
	  unsigned inode=el->GetElementVertexIndex(kel,i)-1u;
	  unsigned inode_metis=GetMetisDof(inode,2);
      	  vtx+= (*_coordinate->_Sol[0])(inode_metis);
	  vty+= (*_coordinate->_Sol[1])(inode_metis);
	  vtz+= (*_coordinate->_Sol[2])(inode_metis);
	}
	vtx/=nve;
	vty/=nve;
	vtz/=nve;
	if( (*_coordinate->_Sol[3])(iel_metis) < 0.5 &&
	    _SetRefinementFlag(vtx,vty,vtz,el->GetElementGroup(kel),_grid) ) {
	    _coordinate->_Sol[3]->set(iel_metis,1.);
	}
      }
      _coordinate->_Sol[3]->close();
    }
       
    std::vector<double> AMR_local;
    _coordinate->_Sol[3]->localize_to_all(AMR_local);
  
    el->InitRefinedToZero();
    
    for (unsigned iel_metis=0; iel_metis<nel; iel_metis++) {
      if(AMR_local[iel_metis]>0.5){
	unsigned iel=IS_Mts2Gmt_elem[iel_metis];
	el->SetRefinedElementIndex(iel,1);
	el->AddToRefinedElementNumber(1);
	short unsigned elt=el->GetElementType(iel);
	el->AddToRefinedElementNumber(1,elt);
      }   
    }
    el->AllocateChildrenElement(_ref_index);
}







//-------------------------------------------------------------------
void Mesh::FlagOnlyEvenElementsToBeRefined() {
  
   el->InitRefinedToZero();

   //refine all next grid even elements
   for (unsigned iel=0; iel<nel; iel+=2) {
     el->SetRefinedElementIndex(iel,1);
     el->AddToRefinedElementNumber(1);
     short unsigned elt=el->GetElementType(iel);
     el->AddToRefinedElementNumber(1,elt);
   }
   el->AllocateChildrenElement(_ref_index);
}


/**
 *  This function generates a finer Mesh level, $l_i$, from a coarser Mesh level $l_{i-1}$, $i>0$
 **/

void Mesh::SetFiniteElementPtr(const elem_type * OtherFiniteElement[6][5]){
  
  for(int i=0;i<6;i++)
    for(int j=0;j<5;j++)
      _finiteElement[i][j] = OtherFiniteElement[i][j];
}

//------------------------------------------------------------------------------------------------------
void Mesh::RefineMesh(const unsigned & igrid, Mesh *mshc, const elem_type *otherFiniteElement[6][5]) {
  
  SetFiniteElementPtr(otherFiniteElement);
    
  elem *elc=mshc->el;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &_iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &_nprocs);
  
  const unsigned fine2CoarseVertexMapping[6][8][8]= { // coarse Mesh dof = f2CVM[element type][fine element][fine vertex]
    { {1,9,25,12,17,21,27,24},
      {9,2,10,25,21,18,22,27},
      {25,10,3,11,27,22,19,23},
      {12,25,11,4,24,27,23,20},
      {17,21,27,24,5,13,26,16},
      {21,18,22,27,13,6,14,26},
      {27,22,19,23,26,14,7,15},
      {24,27,23,20,16,26,15,8} },
    { {1,5,7,8},
      {2,6,5,9},
      {3,7,6,10},
      {8,9,10,4},
      {5,6,8,9},
      {6,10,8,9},
      {5,8,6,7},
      {6,8,10,7} },
    { {1,7,9,13,16,18},
      {2,8,7,14,17,16},
      {3,9,8,15,18,17},
      {7,8,9,16,17,18},
      {13,16,18,4,10,12},
      {14,17,16,5,11,10},
      {15,18,17,6,12,11},
      {16,17,18,10,11,12} },
    { {1,5,9,8},
      {5,2,6,9},
      {9,6,3,7},
      {8,9,7,4} },
    { {1,4,6},
      {4,2,5},
      {6,5,3},
      {5,6,4} },
    { {1,3},
      {3,2} }
  };

  const unsigned coarse2FineFaceMapping[6][6][4][2]= { // fine element,fine face=c2FFM[element type][coarse face][face split index][0,1]
    {
      { {0,0},{1,0},{4,0},{5,0} },
      { {1,1},{2,1},{5,1},{6,1} },
      { {2,2},{3,2},{6,2},{7,2} },
      { {3,3},{0,3},{7,3},{4,3} },
      { {0,4},{1,4},{2,4},{3,4} },
      { {4,5},{5,5},{6,5},{7,5} }
    },
    { 
      { {0,0},{1,0},{2,0},{6,3} },
      { {0,1},{1,3},{3,1},{4,3} },
      { {1,1},{2,3},{3,2},{5,1} },
      { {2,1},{0,3},{3,3},{7,2} }
    },
    { 
      { {0,0},{1,2},{4,0},{5,2} },
      { {1,0},{2,2},{5,0},{6,2} },
      { {2,0},{0,2},{6,0},{4,2} },
      { {0,3},{1,3},{2,3},{3,3} },
      { {4,4},{5,4},{6,4},{7,4} }
    },
    { 
      { {0,0},{1,0} },
      { {1,1},{2,1} },
      { {2,2},{3,2} },
      { {3,3},{0,3} }
    },
    { 
      { {0,0},{1,0} },
      { {1,1},{2,1} },
      { {2,2},{0,2} }
    },
    { 
      { {0,0} },
      { {1,1} }
    }
  };

  const unsigned edge2VerticesMapping[6][12][2]= { // vertex1,vertex2=e2VM[element type][edge][0,1]
    {
      {0,1},{1,2},{2,3},{3,0},
      {4,5},{5,6},{6,7},{7,4},
      {0,4},{1,5},{2,6},{3,7}
    },
    { 
      {0,1},{1,2},{2,0},
      {0,3},{1,3},{2,3}
    },
    { 
      {0,1},{1,2},{2,0},
      {3,4},{4,5},{5,3},
      {0,3},{1,4},{2,5}
    },
    {
      {0,1},{1,2},{2,3},{3,0}
    },
    {
      {0,1},{1,2},{2,0}
    },
    { 
      {0,1}      
    }
  };

  const unsigned vertices2EdgeMapping[6][7][8]= { // edge =v2EM[element type][vertex1][vertex2] with vertex1<vertex2
    { {0,8,0,11,16},
      {0,0,9,0,0,17},
      {0,0,0,10,0,0,18},
      {0,0,0,0,0,0,0,19},
      {0,0,0,0,0,12,0,15},
      {0,0,0,0,0,0,13},
      {0,0,0,0,0,0,0,14}
    },
    { 
      {0,4,6,7},
      {0,0,5,8},
      {0,0,0,9}
    },
    { 
      {0,6,8,12},
      {0,0,7,0,13},
      {0,0,0,0,0,14},
      {0,0,0,0,9,11},
      {0,0,0,0,0,10}
    },
    { 
      {0,4,0,7},
      {0,0,5},
      {0,0,0,6}
    },
    { 
      {0,3,5},
      {0,0,4}
    },
    {
      {0,2}
    }
  };
  _grid=igrid;

//   nel=elc->GetRefinedElementNumber()*REF_INDEX;
  nel=elc->GetRefinedElementNumber()*_ref_index;
  el=new elem(elc,Mesh::_ref_index);



  unsigned jel=0;
  //divide each coarse element in 8(3D), 4(2D) or 2(1D) fine elemets and find all the vertices

  el->SetElementGroupNumber(elc->GetElementGroupNumber());
  el->SetNumberElementFather(elc->GetElementNumber());

  for (unsigned iel=0; iel<elc->GetElementNumber(); iel++) {
    if ( elc->GetRefinedElementIndex(iel) ) {
      elc->SetRefinedElementIndex(iel,jel+1u);
      unsigned elt=elc->GetElementType(iel);
      // project element type
//       for(unsigned j=0;j<REF_INDEX;j++){
      for (unsigned j=0; j<_ref_index; j++) {
        el->SetElementType(jel+j,elt);
        el->SetElementFather(jel+j,iel);
	elc->SetChildElement(iel,j,jel+j);
      }

      unsigned elg=elc->GetElementGroup(iel);
      unsigned elmat=elc->GetElementMaterial(iel);
      // project element group
//       for(unsigned j=0;j<REF_INDEX;j++){
      for (unsigned j=0; j<_ref_index; j++) {
        el->SetElementGroup(jel+j,elg);
	el->SetElementMaterial(jel+j,elmat);
      }

      // project vertex indeces
//       for(unsigned j=0;j<REF_INDEX;j++)
      for (unsigned j=0; j<_ref_index; j++)
        for (unsigned inode=0; inode<elc->GetElementDofNumber(iel,0); inode++)
          el->SetElementVertexIndex(jel+j,inode,elc->GetElementVertexIndex(iel,fine2CoarseVertexMapping[elt][j][inode]-1u));
      // project face indeces
      for (unsigned iface=0; iface<elc->GetElementFaceNumber(iel); iface++) {
        int value=elc->GetFaceElementIndex(iel,iface);
        if (0>value)
// 	  for(unsigned jface=0;jface<FACE_INDEX;jface++)
          for (unsigned jface=0; jface<_face_index; jface++)
            el->SetFaceElementIndex(jel+coarse2FineFaceMapping[elt][iface][jface][0],coarse2FineFaceMapping[elt][iface][jface][1], value);
      }
      // update element numbers
//       jel+=REF_INDEX;
      jel+=_ref_index;
//       el->AddToElementNumber(REF_INDEX,elt);
      el->AddToElementNumber(_ref_index,elt);
    }
  }

  nvt=elc->GetNodeNumber();
  el->SetVertexNodeNumber(nvt);

  
  //find all the elements near each vertex
  BuildAdjVtx();
  //initialize to zero all the middle edge points
  for (unsigned iel=0; iel<nel; iel++)
    for (unsigned inode=el->GetElementDofNumber(iel,0); inode<el->GetElementDofNumber(iel,1); inode++)
      el->SetElementVertexIndex(iel,inode,0);
  //find all the middle edge points
  for (unsigned iel=0; iel<nel; iel++) {
    unsigned ielt=el->GetElementType(iel);
    unsigned istart=el->GetElementDofNumber(iel,0);
    unsigned iend=el->GetElementDofNumber(iel,1);
    for (unsigned inode=istart; inode<iend; inode++)
      if (0==el->GetElementVertexIndex(iel,inode)) {
        nvt++;
        el->SetElementVertexIndex(iel,inode,nvt);
        unsigned im=el->GetElementVertexIndex(iel,edge2VerticesMapping[ielt][inode-istart][0]);
        unsigned ip=el->GetElementVertexIndex(iel,edge2VerticesMapping[ielt][inode-istart][1]);
        //find all the near elements which share the same middle edge point
        for (unsigned j=0; j<el->GetVertexElementNumber(im-1u); j++) {
          unsigned jel=el->GetVertexElementIndex(im-1u,j)-1u;
          if (jel>iel) {
            unsigned jm=0,jp=0;
            unsigned jelt=el->GetElementType(jel);
            for (unsigned jnode=0; jnode<el->GetElementDofNumber(jel,0); jnode++) {
              if (el->GetElementVertexIndex(jel,jnode)==im) {
                jm=jnode+1u;
                break;
              }
            }
            if (jm!=0) {
              for (unsigned jnode=0; jnode<el->GetElementDofNumber(jel,0); jnode++) {
                if (el->GetElementVertexIndex(jel,jnode)==ip) {
                  jp=jnode+1u;
                  break;
                }
              }
              if (jp!=0) {
                if (jp<jm) {
                  unsigned tp=jp;
                  jp=jm;
                  jm=tp;
                }
                el->SetElementVertexIndex(jel,vertices2EdgeMapping[jelt][--jm][--jp],nvt);
              }
            }
          }
        }
      }
  }
  el->SetMidpointNodeNumber(nvt - el->GetVertexNodeNumber());
  // Now build kmid
  Buildkmid();
  Buildkel();
  
  //for parallel computations
  if (_nprocs>=1) GenerateMetisMeshPartition();
    
  // build Mesh coordinates by projecting the coarse coordinats
  _coordinate = new Solution(this);
  _coordinate->AddSolution("X",LAGRANGE,SECOND,1,0); 
  _coordinate->AddSolution("Y",LAGRANGE,SECOND,1,0); 
  _coordinate->AddSolution("Z",LAGRANGE,SECOND,1,0); 
  
  _coordinate->ResizeSolutionVector("X");
  _coordinate->ResizeSolutionVector("Y");
  _coordinate->ResizeSolutionVector("Z");    
  
  _coordinate->AddSolution("AMR",DISCONTINOUS_POLYNOMIAL,ZERO,1,0); 
  _coordinate->ResizeSolutionVector("AMR");
     
  //build projection Matrix
  unsigned thisSolType=2;
  if(_coordinate->_ProjMatFlag[thisSolType]==0){ 
    _coordinate->_ProjMatFlag[thisSolType]=1;

    int nf     = MetisOffset[thisSolType][_nprocs];
    int nc     = mshc->MetisOffset[thisSolType][_nprocs];
    int nf_loc = own_size[thisSolType][_iproc];
    int nc_loc = mshc->own_size[thisSolType][_iproc]; 

    //build matrix sparsity pattern size for efficient storage
    
    NumericVector *NNZ_d = NumericVector::build().release();
    NNZ_d->init(*_coordinate->_Sol[0]);
    NNZ_d->zero();
    
    NumericVector *NNZ_o = NumericVector::build().release();
    NNZ_o->init(*_coordinate->_Sol[0]);
    NNZ_o->zero();
        
    for(int isdom=_iproc; isdom<_iproc+1; isdom++) {
      for (int iel_mts=mshc->IS_Mts2Gmt_elem_offset[isdom];iel_mts < mshc->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
	unsigned iel = mshc->IS_Mts2Gmt_elem[iel_mts];
	if(mshc->el->GetRefinedElementIndex(iel)){ //only if the coarse element has been refined
	  short unsigned ielt=mshc->el->GetElementType(iel);
	  _finiteElement[ielt][thisSolType]->GetSparsityPatternSize(*this, *mshc, iel, NNZ_d, NNZ_o); 
	}
      }
    }
    NNZ_d->close();
    NNZ_o->close();
    
    unsigned offset = MetisOffset[thisSolType][_iproc];
    vector <int> nnz_d(nf_loc);
    vector <int> nnz_o(nf_loc);
    for(int i=0; i<nf_loc;i++){
      nnz_d[i]=static_cast <int> ((*NNZ_d)(offset+i));
      nnz_o[i]=static_cast <int> ((*NNZ_o)(offset+i));
    }
            
    delete NNZ_d;
    delete NNZ_o;
    
    //build matrix
    _coordinate->_ProjMat[thisSolType] = SparseMatrix::build().release();
    _coordinate->_ProjMat[thisSolType]->init(nf,nc,nf_loc,nc_loc,nnz_d,nnz_o);
    // loop on the coarse grid 
    for(int isdom=_iproc; isdom<_iproc+1; isdom++) {
      for (int iel_mts=mshc->IS_Mts2Gmt_elem_offset[isdom]; 
	   iel_mts < mshc->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
	unsigned iel = mshc->IS_Mts2Gmt_elem[iel_mts];
	if( mshc->el->GetRefinedElementIndex(iel)){ //only if the coarse element has been refined
	  short unsigned ielt= mshc->el->GetElementType(iel);
	  _finiteElement[ielt][thisSolType]->BuildProlongation(*this,*mshc,iel,_coordinate->_ProjMat[thisSolType]); 
	}
      }
    }
    _coordinate->_ProjMat[thisSolType]->close();
  }
      
  _coordinate->_Sol[0]->matrix_mult(*mshc->_coordinate->_Sol[0],*_coordinate->_ProjMat[thisSolType]);
  _coordinate->_Sol[1]->matrix_mult(*mshc->_coordinate->_Sol[1],*_coordinate->_ProjMat[thisSolType]);
  _coordinate->_Sol[2]->matrix_mult(*mshc->_coordinate->_Sol[2],*_coordinate->_ProjMat[thisSolType]);
  _coordinate->_Sol[0]->close();
  _coordinate->_Sol[1]->close();
  _coordinate->_Sol[2]->close();     
         
}

//------------------------------------------------------------------------------------------------------

 Mesh::~Mesh(){
    delete el;
    _coordinate->FreeSolutionVectors(); 
    delete _coordinate;
    if(_nprocs>=1){
      delete [] epart;
      delete [] npart;
    }
  }

/// print Mesh info
void Mesh::PrintInfo() {
  
 std::cout << " Mesh Level        : " << _grid << std::endl; 
 std::cout << "   Number of elements: " << nel << std::endl; 
 std::cout << "   Number of nodes   : " << nvt << std::endl;
  
}


//------------------------------------------------------------------------------------------------------
void Mesh::GenerateMetisMeshPartition(){
   
  unsigned eind_size = el->GetElementNumber("Hex")*NVE[0][3] + el->GetElementNumber("Tet")*NVE[1][3] 
                     + el->GetElementNumber("Wedge")*NVE[2][3] + el->GetElementNumber("Quad")*NVE[3][3] 
                     + el->GetElementNumber("Triangle")*NVE[4][3] + el->GetElementNumber("Line")*NVE[5][3];

  idx_t *eptr=new idx_t [nel+1];
  idx_t *eind=new idx_t [eind_size];
  
  idx_t objval;
  idx_t options[METIS_NOPTIONS]; 
  
  METIS_SetDefaultOptions(options);
    
  options[METIS_OPTION_NUMBERING]= 0;
  options[METIS_OPTION_DBGLVL]   = 0;
  options[METIS_OPTION_CTYPE]    = METIS_CTYPE_SHEM; 
  options[METIS_OPTION_PTYPE]	 = METIS_PTYPE_KWAY;
  options[METIS_OPTION_IPTYPE]   = METIS_IPTYPE_RANDOM;
  options[METIS_OPTION_CONTIG]   = 0;//params->contig;
  options[METIS_OPTION_MINCONN]  = 1;
  options[METIS_OPTION_NITER]    = 10;
  options[METIS_OPTION_UFACTOR]  = 100;
  
  eptr[0]=0;
  unsigned counter=0;
  for (unsigned iel=0; iel<nel; iel++) {
    unsigned ielt=el->GetElementType(iel);
    eptr[iel+1]=eptr[iel]+NVE[ielt][3];
    
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,3); inode++){
      eind[counter]=el->GetElementVertexIndex(iel,inode)-1;
    
      counter++;
    }
    
  }
  
  idx_t mnel = nel; 
  idx_t mnvt = nvt;
  idx_t ncommon = _dimension+1;
  nsubdom = _nprocs;
    
  epart=new idx_t [nel];
  npart=new idx_t [nvt];
  
  if(nsubdom!=1) {
  //I call the Mesh partioning function of Metis library (output is epart(own elem) and npart (own nodes))
  int err = METIS_PartMeshDual(&mnel, &mnvt, eptr, eind, NULL, NULL, &ncommon, &nsubdom, NULL, options, &objval, epart, npart);
  
  if(err==METIS_OK) {
    cout << " METIS PARTITIONING IS OK " << endl;
  }
  else if(err==METIS_ERROR_INPUT) {
    cout << " METIS_ERROR_INPUT " << endl;
    exit(1);
  }
  else if (err==METIS_ERROR_MEMORY) {
    cout << " METIS_ERROR_MEMORY " << endl;
    exit(2);
  }
  else {
    cout << " METIS_GENERIC_ERROR " << endl;
    exit(3);
   }
  }
  else {
    //serial computation
    for(unsigned i=0;i<nel;i++) {
      epart[i]=0;
    } 
  }
  
  delete [] eptr;
  delete [] eind;

   //dof map: piecewise liner 0, quadratic 1, biquadratic 2, piecewise constant 3, picewise discontinous linear 4 
  
  //resize the vector IS_Gmt2Mts_dof and dof
  for(int k=0;k<5;k++) {
    IS_Gmt2Mts_dof[k].resize(GetDofNumber(k)); 
    IS_Gmt2Mts_dof_offset[k].resize(nsubdom+1);
  }
  IS_Mts2Gmt_elem.resize(nel);
  IS_Mts2Gmt_elem_offset.resize(nsubdom+1);
  
  // I 
  for(unsigned i=0;i<nvt;i++) {
    npart[i]=nsubdom;
  }
  
  IS_Mts2Gmt_elem_offset[0]=0;
  vector <unsigned> IS_Gmt2Mts_dof_counter(5,0);
  
   for(int k=0;k<5;k++) {
     IS_Gmt2Mts_dof[k].assign(GetDofNumber(k),GetDofNumber(k)-1); 
     //TODO for domain decomposition pourposes! the non existing dofs point to the last dof!!!!!!
   }
  
  
  IS_Gmt2Mts_dof_counter[3]=0;
  IS_Gmt2Mts_dof_counter[4]=0;
  
  for(int isdom=0;isdom<nsubdom;isdom++) {
    for(unsigned iel=0;iel<nel;iel++){
      if(epart[iel]==isdom){
	//filling the piecewise IS_Mts2Gmt_elem metis->gambit
	IS_Mts2Gmt_elem[ IS_Gmt2Mts_dof_counter[3] ]=iel;
	IS_Gmt2Mts_dof[3][iel]=IS_Gmt2Mts_dof_counter[3];
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
	for (unsigned inode=el->GetElementDofNumber(iel,1); inode<el->GetElementDofNumber(iel,3); inode++) {
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
      for(unsigned iel=0;iel<nel;iel++){
     	if(epart[iel]==isdom){
	  IS_Gmt2Mts_dof[4][iel+k_dim*nel]=IS_Gmt2Mts_dof_counter[4];
	  IS_Gmt2Mts_dof_counter[4]++;
	}
      }
    }
  }  
   
  // ghost vs own nodes
  vector <unsigned short> node_count(nvt,0);
  
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
      
      
      for (unsigned inode=el->GetElementDofNumber(iel,1); inode<el->GetElementDofNumber(iel,3); inode++) {
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
  
  node_count.assign (nvt,0);  
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
      
      for(unsigned inode=el->GetElementDofNumber(iel,1); inode<el->GetElementDofNumber(iel,3); inode++) {
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
  
  for(int i=1;i<=nsubdom;i++){
    MetisOffset[0][i]= MetisOffset[0][i-1]+own_size[0][i-1];
    MetisOffset[1][i]= MetisOffset[1][i-1]+own_size[1][i-1];
    MetisOffset[2][i]= MetisOffset[2][i-1]+own_size[2][i-1];
    MetisOffset[3][i]= IS_Mts2Gmt_elem_offset[i];
    MetisOffset[4][i]= IS_Mts2Gmt_elem_offset[i]*(_dimension+1);
    
  }
  
   
  return; 
  
}

/**
 * This function searches all the elements around all the vertices
 **/
void Mesh::BuildAdjVtx() {
  el->AllocateVertexElementMemory();
  for (unsigned iel=0; iel<nel; iel++) {
    for (unsigned inode=0; inode < el->GetElementDofNumber(iel,0); inode++) {
      unsigned ii=el->GetElementVertexIndex(iel,inode)-1u;
      unsigned jj=0;
      while ( 0 != el->GetVertexElementIndex(ii,jj) ) jj++;
      el->SetVertexElementIndex(ii,jj,iel+1u);
    }
  }
}


/**
 * This function generates kmid for hex and wedge elements
 **/
void Mesh::Buildkmid() {
  for (unsigned iel=0; iel<el->GetElementNumber(); iel++)
    for (unsigned inode=el->GetElementDofNumber(iel,1); inode<el->GetElementDofNumber(iel,2); inode++)
      el->SetElementVertexIndex(iel,inode,0);

  for (unsigned iel=0; iel<el->GetElementNumber(); iel++) {
    for (unsigned iface=0; iface<el->GetElementFaceNumber(iel,0); iface++) {
      unsigned inode=el->GetElementDofNumber(iel,1)+iface;
      if ( 0==el->GetElementVertexIndex(iel,inode) ) {
        el->SetElementVertexIndex(iel,inode,++nvt);
        unsigned i1=el->GetFaceVertexIndex(iel,iface,0);
        unsigned i2=el->GetFaceVertexIndex(iel,iface,1);
        unsigned i3=el->GetFaceVertexIndex(iel,iface,2);
        for (unsigned j=0; j< el->GetVertexElementNumber(i1-1u); j++) {
          unsigned jel= el->GetVertexElementIndex(i1-1u,j)-1u;
          if (jel>iel) {
            for (unsigned jface=0; jface<el->GetElementFaceNumber(jel,0); jface++) {
              unsigned jnode=el->GetElementDofNumber(jel,1)+jface;
              if ( 0==el->GetElementVertexIndex(jel,jnode) ) {
                unsigned j1=el->GetFaceVertexIndex(jel,jface,0);
                unsigned j2=el->GetFaceVertexIndex(jel,jface,1);
                unsigned j3=el->GetFaceVertexIndex(jel,jface,2);
                unsigned j4=el->GetFaceVertexIndex(jel,jface,3);
                if ((i1==j1 || i1==j2 || i1==j3 ||  i1==j4 )&&
                    (i2==j1 || i2==j2 || i2==j3 ||  i2==j4 )&&
                    (i3==j1 || i3==j2 || i3==j3 ||  i3==j4 )) {
                  el->SetElementVertexIndex(jel,jnode,nvt);
                }
              }
            }
          }
        }
      }
    }
  }

  for (unsigned iel=0; iel<el->GetElementNumber(); iel++) {
    if (0==el->GetElementType(iel)) {
      el->SetElementVertexIndex(iel,26,++nvt);
    }
    if (3==el->GetElementType(iel)) {
      el->SetElementVertexIndex(iel,8,++nvt);
    }
  }
  el->SetNodeNumber(nvt);

  unsigned nv0= el->GetVertexNodeNumber();
  unsigned nv1= el->GetMidpointNodeNumber();
  el->SetCentralNodeNumber(nvt-nv0-nv1);

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
// 		if((DIM[2]==1 &&
                if ((Mesh::_dimension==3 &&
                     (i1==j1 || i1==j2 || i1==j3 ||  i1==j4 )&&
                     (i2==j1 || i2==j2 || i2==j3 ||  i2==j4 )&&
                     (i3==j1 || i3==j2 || i3==j3 ||  i3==j4 ))||
// 		   (DIM[1]==1 &&
                    (Mesh::_dimension==2 &&
                     (i1==j1 || i1==j2 )&&
                     (i2==j1 || i2==j2 ))||
// 		   (DIM[0]==1 &&
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

// if(Mesh::_dimension != 2) exit(1);
  switch (type) {
  case 0:
    return el->GetVertexNodeNumber();
    break;
  case 1:
    return el->GetVertexNodeNumber()+el->GetMidpointNodeNumber();
    break;
  case 2:
    return nvt;
    break;
  case 3:
    return nel;
    break;
  case 4:
//     return nel*(2+DIMENSION);
    return nel*(2+Mesh::_dimension-1);
    break;
  }
  return 0;
}


/**
 * This function copies the refined element index vector in other_vector
 **/
void Mesh::copy_elr(vector <unsigned> &other_vec) const {
  for (unsigned i=0; i<nel; i++)
    other_vec[i]=el->GetRefinedElementIndex(i);
}


void Mesh::AllocateAndMarkStructureNode() {
  el->AllocateNodeRegion();
  for (unsigned iel=0; iel<nel; iel++) {

    int flag_mat = el->GetElementMaterial(iel);

    if (flag_mat==4) {
      unsigned nve=el->GetElementDofNumber(iel);
      for ( unsigned i=0; i<nve; i++) {
        unsigned inode=el->GetElementVertexIndex(iel,i)-1u;
        el->SetNodeRegion(inode, 1);
      }
    }
  }
  return;
}


void Mesh::GenerateVankaPartitions_FAST( const unsigned &block_size, vector < vector< unsigned > > &block_elements,
					 vector <unsigned> &block_type_range){
  unsigned iproc=processor_id();
  unsigned ElemOffset    = IS_Mts2Gmt_elem_offset[iproc];
  unsigned ElemOffsetp1  = IS_Mts2Gmt_elem_offset[iproc+1];
  unsigned OwnedElements = ElemOffsetp1 - ElemOffset;
  unsigned reminder = OwnedElements % block_size;
  
  unsigned nbloks = (0 == reminder)? OwnedElements/block_size : OwnedElements/block_size+1 ;
    
  //vector < unsigned > size(nbloks,0);
  block_elements.resize(nbloks);
  
  for(int i = 0; i < block_elements.size(); i++) 
    block_elements[i].resize(block_size);
  if  (0 != reminder ){
    block_elements[ block_elements.size()-1].resize(reminder); 
  }
  
  for (unsigned iel=0; iel<OwnedElements; iel++) {
    block_elements[iel/block_size][iel%block_size]=iel+ElemOffset;
  }
  block_type_range.resize(2);
  block_type_range[0]=block_elements.size();
  block_type_range[1]=block_elements.size();
} 

void Mesh::GenerateVankaPartitions_FSI( const unsigned &block_size, vector < vector< unsigned > > &block_elements,
					vector <unsigned> &block_type_range){

  unsigned iproc=processor_id();
  unsigned ElemOffset    = IS_Mts2Gmt_elem_offset[iproc];
  unsigned ElemOffsetp1  = IS_Mts2Gmt_elem_offset[iproc+1];
  unsigned OwnedElements = ElemOffsetp1 - ElemOffset;

  block_elements.resize(2); 
  block_elements[0].resize(OwnedElements); 
  block_elements[1].resize(OwnedElements); 
  
  block_type_range.resize(2);
  block_type_range[0]=1;
  block_type_range[1]=2;
  
  
  unsigned counter_f=0;
  unsigned counter_s=0;
  for (unsigned iel_mts = ElemOffset; iel_mts < ElemOffsetp1; iel_mts++) {
    unsigned kel        = IS_Mts2Gmt_elem[iel_mts]; 
    unsigned flag_mat   = el->GetElementMaterial(kel);
    if(2 == flag_mat){
      block_elements[0][counter_f]=iel_mts;
      counter_f++;
    } 
    else{
      block_elements[1][counter_s]=iel_mts;
      counter_s++;
    }
  }  
  
  if(counter_s==0){
    block_elements.erase (block_elements.begin()+1);
    block_type_range[1]=1;
  }
  else{
    block_elements[1].resize(counter_s);
  }
  
  if(counter_f==0){
    block_elements.erase (block_elements.begin());
    block_type_range[0]=0;
    block_type_range[1]=1;
  }
  else{
    block_elements[0].resize(counter_f); 
  }
  
} 

const unsigned Mesh::GetElementMaterial(unsigned &kel) const {
  unsigned flag_mat = el->GetElementMaterial(kel);
  if(flag_mat==2){
    unsigned nve = el->GetElementDofNumber(kel,2);
    for(int i=0;i<nve;i++){
      unsigned inode=el->GetElementVertexIndex(kel,i)-1u;
      if(1 == el->GetNodeRegion(inode)) flag_mat=4;
    }
  }
  return flag_mat;
}



void Mesh::GenerateVankaPartitions_FSI1( const unsigned *block_size, vector < vector< unsigned > > &block_elements,
					 vector <unsigned> &block_type_range){

  unsigned iproc=processor_id();
  unsigned ElemOffset    = IS_Mts2Gmt_elem_offset[iproc];
  unsigned ElemOffsetp1  = IS_Mts2Gmt_elem_offset[iproc+1];
  unsigned OwnedElements = ElemOffsetp1 - ElemOffset;

  unsigned counter[2]={0,0};
  for (unsigned iel_mts = ElemOffset; iel_mts < ElemOffsetp1; iel_mts++) {
    unsigned kel        = IS_Mts2Gmt_elem[iel_mts]; 
    unsigned flag_mat   = el->GetElementMaterial(kel);
//    unsigned flag_mat   = GetElementMaterial(kel);
    if(2 == flag_mat){
      counter[1]++;
    } 
  }
  counter[0]=OwnedElements - counter[1];
  
  block_type_range.resize(2);
 
  unsigned flag_block[2]={4,2};
  
  unsigned block_start=0;
  unsigned iblock=0;
  while(iblock < 2){
    if(counter[iblock] !=0 ){
      unsigned reminder = counter[iblock] % block_size[iblock];
      unsigned blocks = (0 == reminder)? counter[iblock]/block_size[iblock] : counter[iblock]/block_size[iblock] + 1 ;
      block_elements.resize(block_start+blocks);
  
      for(int i = 0; i < blocks; i++) 
	block_elements[block_start + i].resize(block_size[iblock]);
      if  (0 != reminder ){
	block_elements[block_start + (blocks-1u)].resize(reminder); 
      }
      
      unsigned counter=0;
      for (unsigned iel_mts = ElemOffset; iel_mts < ElemOffsetp1; iel_mts++) {
	unsigned kel        = IS_Mts2Gmt_elem[iel_mts]; 
 	unsigned flag_mat   = el->GetElementMaterial(kel);
//	unsigned flag_mat   = GetElementMaterial(kel);
	if( flag_block[iblock] == flag_mat ){
	  block_elements[ block_start + (counter / block_size[iblock]) ][ counter % block_size[iblock] ]=iel_mts;
	  counter++;
	}	 
      }
      block_type_range[iblock]=block_start+blocks;
      block_start += blocks;
    }
    else{
      block_type_range[iblock]=block_start;
    }
    iblock++;
  }

} 

void Mesh::GenerateVankaPartitions_METIS( const unsigned &vnk_blck, vector < vector< unsigned > > &block_elements){
   
 
  
  unsigned iproc=processor_id();
   
  unsigned ElemOffset    = IS_Mts2Gmt_elem_offset[iproc];
  unsigned ElemOffsetp1  = IS_Mts2Gmt_elem_offset[iproc+1];
  unsigned OwnedElements = ElemOffsetp1 - ElemOffset;
    
  vector < idx_t > connectivity_index_size(OwnedElements+1);
  vector < idx_t > connectivity(0); 
  
  connectivity.reserve(static_cast < int > (pow(3,_dimension)) * OwnedElements + 1);
     
  map <unsigned, unsigned> node_map;
    
  connectivity_index_size[0] = 0;
  unsigned counter = 0;
  unsigned node_counter = 0;
  for (unsigned iel_metis = ElemOffset; iel_metis < ElemOffsetp1; iel_metis++) {
    unsigned iel        = IS_Mts2Gmt_elem[iel_metis]; 
    short unsigned ielt = el->GetElementType(iel);
    unsigned nve        = el->GetElementDofNumber(iel,3);
       
    connectivity_index_size[ iel_metis - ElemOffset + 1 ]=connectivity_index_size[iel_metis - ElemOffset]+nve;
    connectivity.resize( connectivity.size() + nve );
    for (unsigned j = 0; j < nve; j++){
      unsigned jnode = el->GetElementVertexIndex(iel,j)-1;
      map <unsigned, unsigned>::iterator it = node_map.find(jnode);
      if(it != node_map.end() ){
	connectivity[counter] = it->second;
      }
      else{
	node_map[jnode]       = node_counter;
	connectivity[counter] = node_counter;
	node_counter++;
      }
    }
  }		     
 
  //cout<<1<<endl;
 
  idx_t objval;
  idx_t options[METIS_NOPTIONS]; 
  
  METIS_SetDefaultOptions(options);
    
  options[METIS_OPTION_NUMBERING]= 0;
  options[METIS_OPTION_DBGLVL]   = 0;
  options[METIS_OPTION_CTYPE]    = METIS_CTYPE_SHEM; 
  options[METIS_OPTION_PTYPE]	 = METIS_PTYPE_KWAY;
  options[METIS_OPTION_IPTYPE]   = METIS_IPTYPE_RANDOM;
  options[METIS_OPTION_CONTIG]   = 0;//params->contig;
  options[METIS_OPTION_MINCONN]  = 1;
  options[METIS_OPTION_NITER]    = 10;
  options[METIS_OPTION_UFACTOR]  = 100;
  
  idx_t mnel = OwnedElements; 
  idx_t mnvt = node_map.size();
  idx_t ncommon = _dimension+1;
  idx_t nsubdom = OwnedElements/vnk_blck;
    
  vector < idx_t > epart(OwnedElements);
  vector < idx_t > npart(node_map.size());
  
  if(nsubdom!=1) {
  //I call the Mesh partioning function of Metis library (output is epart(own elem) and npart (own nodes))
  int err = METIS_PartMeshDual(&mnel, &mnvt, &connectivity_index_size[0], &connectivity[0], 
			       NULL, NULL, &ncommon, &nsubdom, NULL, options, &objval, &epart[0], &npart[0]);
  
    if(err==METIS_OK) {
      cout << " METIS PARTITIONING IS OK " << endl;
    }
    else if(err==METIS_ERROR_INPUT) {
      cout << " METIS_ERROR_INPUT " << endl;
      exit(1);
    }
    else if (err==METIS_ERROR_MEMORY) {
      cout << " METIS_ERROR_MEMORY " << endl;
      exit(2);
    }
    else {
      cout << " METIS_GENERIC_ERROR " << endl;
      exit(3);
    }
  }
  else {
    //serial computation
    for(unsigned i=0;i<OwnedElements;i++) {
      epart[i]=0;
    } 
  }
  
  //cout<<2<<endl;
  
  vector < unsigned > size(nsubdom,0);
  block_elements.resize(nsubdom);
  
  //cout<<4<<endl;
  for(int i=0;i<block_elements.size();i++) 
    block_elements[i].resize(2 * vnk_blck );
   //for(int i=0; i<block_elements.size(); i++) {
    //cout<<i<<" "<<block_elements[i].size() << " "<< vnk_blck<<" "<< OwnedElements<<" "<<2 * OwnedElements / vnk_blck +1 << endl;
   //}
  for(int iel=0; iel<OwnedElements; iel++){
    block_elements[ epart[iel] ] [size[epart[iel]]] = iel + ElemOffset;
    (size[epart[iel]])++;
  }
  
  //cout<<6<<endl;
  
  for(int i=0; i<block_elements.size(); i++) {
    //cout<<size[i]<< " "<<block_elements[i].size() << endl;
    block_elements[i].resize( size[i]);
  }
  //cout<<7<<endl;
  
}




} //end namespace femus


