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

#include "Mesh.hpp"
#include "metis.h"


namespace femus {
  
using std::cout;
using std::endl;
using std::min;
using std::sort;


bool mesh::_TestSetRefinementFlag=0;

const unsigned mesh::_END_IND[5]= {0,1,3,4,5};

unsigned mesh::_dimension=2;
unsigned mesh::_ref_index=4;  // 8*DIM[2]+4*DIM[1]+2*DIM[0];
unsigned mesh::_face_index=2; // 4*DIM[2]+2*DIM[1]+1*DIM[0];


mesh::mesh() 
{

}


/**
 *  This function generates the coarse mesh level, $l_0$, from an input mesh file (Now only the Gambit Neutral File)
 **/
void mesh::ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool> &type_elem_flag) {
  
  MPI_Comm_rank(MPI_COMM_WORLD, &_iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &_nprocs);
  
  vector <vector <double> > vt;  
  vt.resize(3);
    
  grid=0;

  if(name.rfind(".neu") < name.size())
  {
    ReadGambit(name,vt,Lref,type_elem_flag);
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
  
 if (_nprocs>=1) generate_metis_mesh_partition();
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
void mesh::FlagAllElementsToBeRefined() {
  
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
void mesh::FlagElementsToBeRefinedByUserDefinedFunction() {
     el->InitRefinedToZero();
   
    //refine based on the function SetRefinementFlag defined in the main;
    // the mesh is serial, we cannot in parallel use the coordinates to selectively refine
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
	if (_SetRefinementFlag(vtx,vty,vtz,el->GetElementGroup(iel),grid)) {
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
void mesh::FlagElementsToBeRefinedByAMR() {
    
       
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
	    _SetRefinementFlag(vtx,vty,vtz,el->GetElementGroup(kel),grid) ) {
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
void mesh::FlagOnlyEvenElementsToBeRefined() {
  
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
 *  This function generates a finer mesh level, $l_i$, from a coarser mesh level $l_{i-1}$, $i>0$
 **/
//------------------------------------------------------------------------------------------------------
void mesh::RefineMesh(const unsigned & igrid, mesh *mshc, const elem_type* type_elem[6][5]) {
  
  elem *elc=mshc->el;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &_iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &_nprocs);
  
  const unsigned vertex_index[6][8][8]= {{{1,9,25,12,17,21,27,24},{9,2,10,25,21,18,22,27},
      {25,10,3,11,27,22,19,23},{12,25,11,4,24,27,23,20},
      {17,21,27,24,5,13,26,16},{21,18,22,27,13,6,14,26},
      {27,22,19,23,26,14,7,15},{24,27,23,20,16,26,15,8}
    },
    { {1,5,7,8},{2,6,5,9},{3,7,6,10},{8,9,10,4},
      {5,6,8,9},{6,10,8,9},{5,8,6,7},{6,8,10,7}
    },
    { {1,7,9,13,16,18},{2,8,7,14,17,16},
      {3,9,8,15,18,17},{7,8,9,16,17,18},
      {13,16,18,4,10,12},{14,17,16,5,11,10},
      {15,18,17,6,12,11},{16,17,18,10,11,12}
    },
    {{1,5,9,8},{5,2,6,9},{9,6,3,7},{8,9,7,4}},
    {{1,4,6},{2,5,4},{3,6,5},{4,5,6}},
    {{1,3},{3,2}}
  };

  const unsigned face_index[6][6][4][2]= {{{{0,0},{1,0},{4,0},{5,0}},
      {{1,1},{2,1},{5,1},{6,1}},
      {{2,2},{3,2},{6,2},{7,2}},
      {{3,3},{0,3},{7,3},{4,3}},
      {{0,4},{1,4},{2,4},{3,4}},
      {{4,5},{5,5},{6,5},{7,5}}
    },
    { {{0,0},{1,0},{2,0},{6,3}},
      {{0,1},{1,3},{3,1},{4,3}},
      {{1,1},{2,3},{3,2},{5,1}},
      {{2,1},{0,3},{3,3},{7,2}}
    },
    { {{0,0},{1,2},{4,0},{5,2}},
      {{1,0},{2,2},{5,0},{6,2}},
      {{2,0},{0,2},{6,0},{4,2}},
      {{0,3},{1,3},{2,3},{3,3}},
      {{4,4},{5,4},{6,4},{7,4}}
    },
    { {{0,0},{1,0}},
      {{1,1},{2,1}},
      {{2,2},{3,2}},
      {{3,3},{0,3}}
    },
    { {{0,0},{1,2}},
      {{1,0},{2,2}},
      {{2,0},{0,2}}
    },
    { {{0,0}},
      {{1,1}}
    }
  };

  const unsigned midpoint_index[6][12][2]= {{{0,1},{1,2},{2,3},{3,0},
      {4,5},{5,6},{6,7},{7,4},
      {0,4},{1,5},{2,6},{3,7}
    },
    { {0,1},{1,2},{2,0},
      {0,3},{1,3},{2,3}
    },
    { {0,1},{1,2},{2,0},
      {3,4},{4,5},{5,3},
      {0,3},{1,4},{2,5}
    },
    {{0,1},{1,2},{2,3},{3,0}},
    {{0,1},{1,2},{2,0}},
    {{0,1}}
  };

  const unsigned midpoint_index2[6][7][8]= {{{0,8,0,11,16},
      {0,0,9,0,0,17},
      {0,0,0,10,0,0,18},
      {0,0,0,0,0,0,0,19},
      {0,0,0,0,0,12,0,15},
      {0,0,0,0,0,0,13},
      {0,0,0,0,0,0,0,14}
    },
    { {0,4,6,7},
      {0,0,5,8},
      {0,0,0,9}
    },
    { {0,6,8,12},
      {0,0,7,0,13},
      {0,0,0,0,0,14},
      {0,0,0,0,9,11},
      {0,0,0,0,0,10}
    },
    { {0,4,0,7},
      {0,0,5},
      {0,0,0,6}
    },
    { {0,3,5},
      {0,0,4}
    },
    {{0,2}}
  };
  grid=igrid;

//   nel=elc->GetRefinedElementNumber()*REF_INDEX;
  nel=elc->GetRefinedElementNumber()*_ref_index;
  el=new elem(elc,mesh::_ref_index);



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
          el->SetElementVertexIndex(jel+j,inode,elc->GetElementVertexIndex(iel,vertex_index[elt][j][inode]-1u));
      // project face indeces
      for (unsigned iface=0; iface<elc->GetElementFaceNumber(iel); iface++) {
        int value=elc->GetFaceElementIndex(iel,iface);
        if (0>value)
// 	  for(unsigned jface=0;jface<FACE_INDEX;jface++)
          for (unsigned jface=0; jface<_face_index; jface++)
            el->SetFaceElementIndex(jel+face_index[elt][iface][jface][0],face_index[elt][iface][jface][1], value);
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
        unsigned im=el->GetElementVertexIndex(iel,midpoint_index[ielt][inode-istart][0]);
        unsigned ip=el->GetElementVertexIndex(iel,midpoint_index[ielt][inode-istart][1]);
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
                el->SetElementVertexIndex(jel,midpoint_index2[jelt][--jm][--jp],nvt);
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
  if (_nprocs>=1) generate_metis_mesh_partition();
    
  // build mesh coordinates by projecting the coarse coordinats
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
	  type_elem[ielt][thisSolType]->GetSparsityPatternSize(*this, *mshc, iel, NNZ_d, NNZ_o); 
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
	  type_elem[ielt][thisSolType]->BuildProlongation(*this,*mshc,iel,_coordinate->_ProjMat[thisSolType]); 
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

 mesh::~mesh(){
    delete el;
    _coordinate->FreeSolutionVectors(); 
    delete _coordinate;
    if(_nprocs>=1){
      delete [] epart;
      delete [] npart;
    }
  }

/// print mesh info
void mesh::print_info() {
  
 std::cout << " Mesh Level        : " << grid << std::endl; 
 std::cout << "   Number of elements: " << nel << std::endl; 
 std::cout << "   Number of nodes   : " << nvt << std::endl;
  
}


//------------------------------------------------------------------------------------------------------
void mesh::generate_metis_mesh_partition(){
   
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
  //I call the mesh partioning function of Metis library (output is epart(own elem) and npart (own nodes))
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
 * This function read the data form the gambit file.
 * It is used in the constructor of the coarse mesh.
 **/
void mesh::ReadGambit(const std::string& name, vector < vector < double> > &vt, const double Lref, std::vector<bool> &type_elem_flag) {
  // set the file
  const unsigned GambitVertexIndex[6][27]= {{
      4,16,0,15,23,11,7,19,3,
      12,20,8,25,26,24,14,22,10,
      5,17,1,13,21,9,6,18,2
    },
    {
      0,4,1,6,5,
      2,7,8,9,3
    },
    {
      3, 11,5, 9, 10,4,
      12,17,14,15,16,13,
      0, 8, 2, 6, 7, 1
    },
    {0,4,1,5,2,6,3,7,8},
    {0,3,1,4,2,5},
    {0,2,1}
  };

  const unsigned GambitFaceIndex[6][6]= {{0,4,2,5,3,1},
    {0,1,2,3},
    {2,1,0,4,3},
    {0,1,2,3},
    {0,1,2},
    {0,1}
  };

  std::ifstream inf;
  std::string str2;
  unsigned ngroup;
  unsigned nbcd;
  unsigned dim;
  double x,y,z;
//   double _lref = 1.;

  grid=0;
  // read control data ******************** A
  inf.open(name.c_str());
  if (!inf) {
    cout<<"Generic-mesh file "<< name << " can not read parameters\n";
    exit(0);
  }
  str2="0";
  while (str2.compare("NDFVL") != 0) inf >> str2;
  inf >> nvt >> nel >>  ngroup >> nbcd >> dim >> str2 ;
  mesh::_dimension = dim;
  mesh::_ref_index = pow(2,mesh::_dimension);  // 8*DIM[2]+4*DIM[1]+2*DIM[0];
  mesh::_face_index = pow(2,mesh::_dimension-1u); // 4*DIM[2]+2*DIM[1]+1*DIM[0];
//   cout << " ref index " << _ref_index << "    " << REF_INDEX << endl;
//   cout << " face index " << _face_index << "    " << FACE_INDEX << endl;
//   cout << " mesh dim " << mesh::_dimension <<  endl;
  inf >> str2;
  if (str2.compare("ENDOFSECTION") != 0) {
    cout<<"error control data mesh"<<endl;
    exit(0);
  }
//   cout << "***************" << _dimension << endl;
  inf.close();
  // end read control data **************** A
  // read ELEMENT/cell ******************** B
  inf.open(name.c_str());
  if (!inf) {
    cout<<"Generic-mesh file "<< name << " cannot read elements\n";
    exit(0);
  }
  el= new elem(nel);
  while (str2.compare("ELEMENTS/CELLS") != 0) inf >> str2;
  inf >> str2;
  for (unsigned iel=0; iel<nel; iel++) {
    el->SetElementGroup(iel,1);
    unsigned nve;
    inf >> str2 >> str2 >> nve;
    if (nve==27) {
      type_elem_flag[0]=type_elem_flag[3]=true;
      el->AddToElementNumber(1,"Hex");
      el->SetElementType(iel,0);
    } else if (nve==10) {
      type_elem_flag[1]=type_elem_flag[4]=true;
      el->AddToElementNumber(1,"Tet");
      el->SetElementType(iel,1);
    } else if (nve==18) {
      type_elem_flag[2]=type_elem_flag[3]=type_elem_flag[4]=true;
      el->AddToElementNumber(1,"Wedge");
      el->SetElementType(iel,2);
    } else if (nve==9) {
      type_elem_flag[3]=true;
      el->AddToElementNumber(1,"Quad");
      el->SetElementType(iel,3);
    }
    else if (nve==6 && mesh::_dimension==2) {
      type_elem_flag[4]=true;
      el->AddToElementNumber(1,"Triangle");
      el->SetElementType(iel,4);
    }
    else if (nve==3 && mesh::_dimension==1) {
      el->AddToElementNumber(1,"Line");
      el->SetElementType(iel,5);
    } else {
      cout<<"Error! Invalid element type in reading Gambit File!"<<endl;
      cout<<"Error! Use a second order discretization"<<endl;
      exit(0);
    }
    for (unsigned i=0; i<nve; i++) {
      unsigned inode=GambitVertexIndex[el->GetElementType(iel)][i];
      double value;
      inf>>value;
      el->SetElementVertexIndex(iel,inode,value);
    }
  }
  inf >> str2;
  if (str2.compare("ENDOFSECTION") != 0) {
    cout<<"error element data mesh"<<endl;
    exit(0);
  }
  inf.close();

  // end read  ELEMENT/CELL **************** B

  // read NODAL COORDINATES **************** C
  inf.open(name.c_str());
  if (!inf) {
    cout<<"Generic-mesh file "<< name << " cannot read nodes\n";
    exit(0);
  }
  while (str2.compare("COORDINATES") != 0) inf >> str2;
  inf >> str2;  // 2.0.4
  vt[0].resize(nvt);
  vt[1].resize(nvt);
  vt[2].resize(nvt);

  if (mesh::_dimension==3) {
    for (unsigned j=0; j<nvt; j++) {
      inf >> str2 >> x >> y >> z;
      vt[0][j] = x/Lref;
      vt[1][j] = y/Lref;
      vt[2][j] = z/Lref;
    }
  }

  else if (mesh::_dimension==2) {
    for (unsigned j=0; j<nvt; j++) {
      inf >> str2 >> x >> y;
      vt[0][j] = x/Lref;
      vt[1][j] = y/Lref;
      vt[2][j] = 0.;
    }
  }

  else if (mesh::_dimension==1) {
    for (unsigned j=0; j<nvt; j++) {
      inf >> str2 >> x;
      vt[0][j] = x/Lref;
      vt[1][j]=0.;
      vt[2][j]=0.;
    }
  }
  inf >> str2; // "ENDOFSECTION"
  if (str2.compare("ENDOFSECTION") != 0) {
    cout<<"error node data mesh 1"<<endl;
    exit(0);
  }
  inf.close();
  // end read NODAL COORDINATES ************* C

  // read GROUP **************** E
  inf.open(name.c_str());
  if (!inf) {
    cout<<"Generic-mesh file "<< name << " cannot read group\n";
    exit(0);
  }
  el->SetElementGroupNumber(ngroup);
  for (unsigned k=0; k<ngroup; k++) {
    int ngel;
    int name;
    int mat;
    while (str2.compare("GROUP:") != 0) inf >> str2;
    inf >> str2 >> str2 >> ngel >> str2 >> mat >>str2 >> str2 >>name>> str2;
//     cout<<ngel<<endl;
//     cout<<name<<endl;
//     cout << mat << endl;
    for (int i=0; i<ngel; i++) {
      int iel;
      inf >> iel;
      el->SetElementGroup(iel-1,name);
      el->SetElementMaterial(iel-1,mat);
    }
    inf >> str2;
    if (str2.compare("ENDOFSECTION") != 0) {
      cout<<"error group data mesh"<<endl;
      exit(0);
    }
  }
  inf.close();
  // end read boundary **************** E

  // read boundary **************** D
  inf.open(name.c_str());
  if (!inf) {
    cout<<"Generic-mesh file "<< name << " cannot read boudary\n";
    exit(0);
  }
  for (unsigned k=0; k<nbcd; k++) {
    while (str2.compare("CONDITIONS") != 0) inf >> str2;
    inf >> str2;
    int value;
    unsigned nface;
    inf >>value>>str2>>nface>>str2>>str2;
    value=-value-1;
    for (unsigned i=0; i<nface; i++) {
      unsigned iel,iface;
      inf>>iel>>str2>>iface;
      iel--;
      iface=GambitFaceIndex[el->GetElementType(iel)][iface-1u];
      el->SetFaceElementIndex(iel,iface,value);
    }
    inf >> str2;
    if (str2.compare("ENDOFSECTION") != 0) {
      cout<<"error boundary data mesh"<<endl;
      exit(0);
    }
  }
  inf.close();
  // end read boundary **************** D

  
  
  //*************** start reorder mesh dofs **************
  //(1)linear (2)quadratic (3)biquaratic
  
  vector <unsigned> dof_index;
  dof_index.resize(nvt);
  for(unsigned i=0;i<nvt;i++){
    dof_index[i]=i+1;
  }
  //reorder vertices and mid-points vs central points
  for (unsigned iel=0; iel<nel; iel++) {
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,1); inode++) {
      for (unsigned jel=0; jel<nel; jel++) {
	for (unsigned jnode=el->GetElementDofNumber(jel,1); jnode<el->GetElementDofNumber(jel,3); jnode++) { 
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
  for (unsigned iel=0; iel<nel; iel++) {
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
      for (unsigned jel=0; jel<nel; jel++) {
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
  for (unsigned iel=0; iel<nel; iel++) {
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,3); inode++) {
      unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
      el->SetElementVertexIndex(iel,inode,dof_index[ii]);
    }
  }
  vector <double> vt_temp;
  for(int i=0;i<3;i++){
    vt_temp=vt[i];
    for(unsigned j=0;j<nvt;j++){
      vt[i][dof_index[j]-1]=vt_temp[j];
    }
  }
  // **************  end reoreder mesh dofs **************
 
  el->SetNodeNumber(nvt);

  unsigned nv0=0;
  for (unsigned iel=0; iel<nel; iel++)
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
      unsigned i0=el->GetElementVertexIndex(iel,inode);
      if (nv0<i0) nv0=i0;
  }
  el->SetVertexNodeNumber(nv0);

  unsigned nv1=0;
  for (unsigned iel=0; iel<nel; iel++)
    for (unsigned inode=el->GetElementDofNumber(iel,0); inode<el->GetElementDofNumber(iel,1); inode++) {
      unsigned i1=el->GetElementVertexIndex(iel,inode);
      if (nv1<i1) nv1=i1;
  }
  el->SetMidpointNodeNumber(nv1-nv0);

  el->SetCentralNodeNumber(nvt-nv1);
  
};



     /**
       * A useful inline function which replaces the #defines
       * used previously.  Not private since this is a namespace,
       * but would be if this were a class.  The first one returns
       * the proper node number for 2D elements while the second
       * one returns the node number for 3D elements.
       */
unsigned int mesh::idx(const ElemType type, const unsigned int nx, const unsigned int i, const unsigned int j) {
	switch(type)
	  {
// 	  case INVALID_ELEM:
// 	  case QUAD4:
// 	  case TRI3:
// 	    {
// 	      return i + j*(nx+1);
// 	      break;
// 	    }
// 
// 	  case QUAD8:
	  case QUAD9:
 	  case TRI6:
	    {
	      return i + j*(2*nx+1);
	      break;
	    }

	  default:
	    {
	      std::cout << "ERROR: Unrecognized or Not Supported 2D element type." << std::endl;
	      exit(1);
	    }
	  }

	return -1; // invalid_uint
}
      

      // Same as the function above, but for 3D elements
      unsigned int mesh::idx(const ElemType type,
		       const unsigned int nx,
		       const unsigned int ny,
		       const unsigned int i,
		       const unsigned int j,
		       const unsigned int k)
      {
	switch(type)
	  {
// 	  case INVALID_ELEM:
// 	  case HEX8:
// 	  case PRISM6:
// 	    {
// 	      return i + (nx+1)*(j + k*(ny+1));
// 	      break;
// 	    }

// 	  case HEX20:
	  case HEX27:
// 	  case TET4:  // TET4's are created from an initial HEX27 discretization
// 	  case TET10: // TET10's are created from an initial HEX27 discretization
// 	  case PYRAMID5: // PYRAMID5's are created from an initial HEX27 discretization
// 	  case PRISM15:
// 	  case PRISM18:
	    {
	      return i + (2*nx+1)*(j + k*(2*ny+1));
	      break;
	    }

	  default:
	    {
	      std::cout << "ERROR: Unrecognized element type." << std::endl;
	      exit(1);
	    }
	  }

	return -1;
 }

// ------------------------------------------------------------
// MeshTools::Generation function for mesh generation
void mesh::BuildBrick(const unsigned int nx,
	              const unsigned int ny,
	              const unsigned int nz,
		      const double xmin, const double xmax,
		      const double ymin, const double ymax,
		      const double zmin, const double zmax,
		      const ElemType type,
		      std::vector<bool> &type_elem_flag ) 
{
  
  MPI_Comm_rank(MPI_COMM_WORLD, &_iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &_nprocs);
  
  vector <vector <double> > vt;  
  vt.resize(3);
    
  grid=0;
  
//   // Clear the mesh and start from scratch
//   mesh.clear();

  if (nz != 0)
    _dimension = 3;
  else if (ny != 0)
    _dimension = 2;
  else if (nx != 0)
    _dimension = 1;
  else
    _dimension = 0;
  
  unsigned ngroup;
  unsigned nbcd;

  switch (_dimension)
    {

      //---------------------------------------------------------------------
       // Build a 1D line
     case 1:
       {
        assert(nx!=0);
	assert(ny==0);
	assert(nz==0);
	assert(xmin<xmax);

         // Reserve elements
         switch (type)
         {
//           case INVALID_ELEM:
//           case EDGE2:
           case EDGE3:
//           case EDGE4:
             {
	       nel    =  nx;
	      
	       ngroup = 1;
	       nbcd   = 2;
	       mesh::_ref_index = pow(2,mesh::_dimension);
	       mesh::_face_index = pow(2,mesh::_dimension-1u);
 	       break;
             }
// 
            default:
            {
 	      std::cerr << "ERROR: Unrecognized 1D element type." << std::endl;
 	    }
 	  }
 
         // Reserve nodes
         switch (type)
           {
//           case INVALID_ELEM:
//           case EDGE2:
//             {
//               mesh.reserve_nodes(nx+1);
//               break;
//             }
// 
           case EDGE3:
             {
	       nvt = (2*nx+1);
 	       break;
             }
// 
//           case EDGE4:
//             {
//               mesh.reserve_nodes(3*nx+1);
//               break;
//             }
// 
           default:
             {
               std::cerr << "ERROR: Unrecognized 1D element type." << std::endl;
             }
           }


        // Build the nodes, depends on whether we're using linears,
        // quadratics or cubics and whether using uniform grid or Gauss-Lobatto
        unsigned int node_id = 0;
         switch(type)
         {
//           case INVALID_ELEM:
//           case EDGE2:
//             {
//               for (unsigned int i=0; i<=nx; i++)
//               {
//                 if (gauss_lobatto_grid)
//                   mesh.add_point (Point(0.5*(std::cos(libMesh::pi*static_cast<double>(nx-i)/static_cast<double>(nx))+1.0),
//                         0,
//                         0), node_id++);
//                 else
//                   mesh.add_point (Point(static_cast<double>(i)/static_cast<double>(nx),
//                         0,
//                         0), node_id++);
//               }
//               break;
//             }
// 
           case EDGE3:
             {
	      vt[0].resize(nvt);
              vt[1].resize(nvt);
              vt[2].resize(nvt);
	       
               for (unsigned int i=0; i<=2*nx; i++)
               {
//                 if (gauss_lobatto_grid)
//                 {
//                   // The x location of the point.
//                   double x=0.;
// 
//                   // Shortcut quantities (do not depend on i)
//                   const double c = std::cos( libMesh::pi*i / static_cast<double>(2*nx) );
// 
//                   // If i is even, compute a normal Gauss-Lobatto point
//                   if (i%2 == 0)
//                     x = 0.5*(1.0 - c);
// 
//                   // Otherwise, it is the average of the previous and next points
//                   else
//                   {
//                     double cmin = std::cos( libMesh::pi*(i-1) / static_cast<double>(2*nx) );
//                     double cmax = std::cos( libMesh::pi*(i+1) / static_cast<double>(2*nx) );
// 
//                     double gl_xmin = 0.5*(1.0 - cmin);
//                     double gl_xmax = 0.5*(1.0 - cmax);
//                     x = 0.5*(gl_xmin + gl_xmax);
//                   }
// 
//                   mesh.add_point (Point(x,0.,0.), node_id++);
//                 }
//                 else
		   vt[0][node_id] = (static_cast<double>(i)/static_cast<double>(2*nx))*(xmax-xmin) + xmin;
                   vt[1][node_id] = 0.;
                   vt[2][node_id] = 0.;
  
                   node_id++;

               }
               break;
             }
 
           default:
             {
               std::cerr << "ERROR: Unrecognized 1D element type." << std::endl;
             }
 
         }
 
         // Build the elements of the mesh
       unsigned iel = 0;
	 el= new elem(nel);
	 el->SetElementGroupNumber(1);
         // Build the elements.  Each one is a bit different.
         switch(type)
           {
//             case INVALID_ELEM:
//             case EDGE2:
//               {
//                 for (unsigned int i=0; i<nx; i++)
//                 {
//                   Elem* elem = mesh.add_elem (new Edge2);
//                   elem->set_node(0) = mesh.node_ptr(i);
//                   elem->set_node(1) = mesh.node_ptr(i+1);
// 
//                   if (i == 0)
//                     mesh.boundary_info->add_side(elem, 0, 0);
// 
//                   if (i == (nx-1))
//                     mesh.boundary_info->add_side(elem, 1, 1);
//                 }
//               break;
//               }
// 
             case EDGE3:
               {
		 
		 
		unsigned LocalToGlobalNodePerElement[3];
		
		for (unsigned int i=0; i<nx; i++)
 		  {
		  
                  el->SetElementGroup(iel,1);
		  el->SetElementMaterial(iel, 2); 
 		  type_elem_flag[5]=true;
                  el->AddToElementNumber(1,"Line");
                  el->SetElementType(iel,5);
  
		  LocalToGlobalNodePerElement[0] = 2*i   + 1; 
		  LocalToGlobalNodePerElement[1] = 2*i+2 + 1;
		  LocalToGlobalNodePerElement[2] = 2*i+1 + 1;
		  
                  // connectivity
                  for (unsigned iloc=0; iloc<3; iloc++) {
                    el->SetElementVertexIndex(iel,iloc,LocalToGlobalNodePerElement[iloc]);
                  }
                                   
                  if (i == 0) {
		    el->SetFaceElementIndex(iel,0,-2);
		    _boundaryinfo.insert( std::pair<unsigned int, std::string>(0,"left"));
		  }

		  if (i == (nx-1)) {
		    el->SetFaceElementIndex(iel,1,-3);
		    _boundaryinfo.insert( std::pair<unsigned int, std::string>(1,"right")); 
		  }
		  
                    
                  iel++;
		  }

               break;
               }
 
             default:
               {
                 std::cerr << "ERROR: Unrecognized 1D element type." << std::endl;
               }
           }

// 	// Scale the nodal positions
// 	for (unsigned int p=0; p<mesh.n_nodes(); p++)
// 	  mesh.node(p)(0) = (mesh.node(p)(0))*(xmax-xmin) + xmin;
// 
//         // Add sideset names to boundary info
//         mesh.boundary_info->sideset_name(0) = "left";
//         mesh.boundary_info->sideset_name(1) = "right";
// 
//         // Add nodeset names to boundary info
//         mesh.boundary_info->nodeset_name(0) = "left";
//         mesh.boundary_info->nodeset_name(1) = "right";

	std::cout << "Error: NotImplemented " << std::endl; 
	break;
      }

      //---------------------------------------------------------------------
      // Build a 2D quadrilateral
     case 2:
       {
 	assert (nx != 0);
 	assert (ny != 0);
 	assert (nz == 0);
 	assert (xmin < xmax);
 	assert (ymin < ymax);
 
 	switch (type)
 	  {
// 	  case INVALID_ELEM:
// 	  case QUAD4:
// 	  case QUAD8:
 	  case QUAD9:
 	    {
 	      nel    =  nx*ny;
	      
	      ngroup = 1;
	      nbcd   = 4;
	      mesh::_ref_index = pow(2,mesh::_dimension);
	      mesh::_face_index = pow(2,mesh::_dimension-1u);
 	      break;
 	    }
// 
// 	  case TRI3:
 	  case TRI6:
 	    {
	      nel    =  2*nx*ny;
	      
	      ngroup = 1;
	      nbcd   = 4;
	      mesh::_ref_index = pow(2,mesh::_dimension);
	      mesh::_face_index = pow(2,mesh::_dimension-1u);
 	      break;
 	    }
// 
 	  default:
 	    {
 	      std::cout << "ERROR: Unrecognized or NotSupported 2D element type." << std::endl;
 	      exit(1);
 	    }
 	  }

// 	// Reserve nodes.  The quadratic element types
// 	// need to reserve more nodes than the linear types.
 	switch (type)
 	  {
// 	  case INVALID_ELEM:
// 	  case QUAD4:
// 	  case TRI3:
// 	    {
// 	      mesh.reserve_nodes( (nx+1)*(ny+1) );
// 	      break;
// 	    }
// 
// 	  case QUAD8:
 	  case QUAD9:
	  case TRI6:
 	    {
 	      nvt = (2*nx+1)*(2*ny+1) ;
 	      break;
 	    }
 
 
 	  default:
 	    {
 	      std::cout << "ERROR: Unrecognized or not Supported 2D element type." << std::endl;
 	      exit(1);
 	    }
 	  }
 
 
 
	// Build the nodes. Depends on whether you are using a linear
	// or quadratic element, and whether you are using a uniform
	// grid [ or the Gauss-Lobatto grid points - NOTIMPLENTED ].
        unsigned int node_id = 0;
	switch (type)
	  {
// 	  case INVALID_ELEM:
// 	  case QUAD4:
// 	  case TRI3:
// 	    {
// 	      for (unsigned int j=0; j<=ny; j++)
// 		for (unsigned int i=0; i<=nx; i++)
// 		  {
// 		    if (gauss_lobatto_grid)
// 		      {
// 			mesh.add_point (Point(0.5*(1.0 - std::cos(libMesh::pi*static_cast<double>(i)/static_cast<double>(nx))),
// 					      0.5*(1.0 - std::cos(libMesh::pi*static_cast<double>(j)/static_cast<double>(ny))),
// 					      0.), node_id++);
// 		      }
// 
// 		    else
// 		      mesh.add_point (Point(static_cast<double>(i)/static_cast<double>(nx),
// 					    static_cast<double>(j)/static_cast<double>(ny),
// 					    0.), node_id++);
// 		  }
// 
// 	      break;
// 	    }
// 
// 	  case QUAD8:
 	  case QUAD9:
 	  case TRI6:
 	    {
	      vt[0].resize(nvt);
              vt[1].resize(nvt);
              vt[2].resize(nvt);
	      
	      
	      for (unsigned int j=0; j<=(2*ny); j++)
		for (unsigned int i=0; i<=(2*nx); i++)
		  {
// 		    if (gauss_lobatto_grid) // NOT YET IMPLENTED
// 		      {
// 			// The x,y locations of the point.
// 			double x=0., y=0.;
// 
// 			// Shortcut quantities (do not depend on i,j)
// 			const double a = std::cos( libMesh::pi / static_cast<double>(2*nx) );
// 			const double b = std::cos( libMesh::pi / static_cast<double>(2*ny) );
// 
// 			// Shortcut quantities (depend on i,j)
// 			const double c = std::cos( libMesh::pi*i / static_cast<double>(2*nx) );
// 			const double d = std::cos( libMesh::pi*j / static_cast<double>(2*ny) );
// 
// 			// If i is even, compute a normal Gauss-Lobatto point
// 			if (i%2 == 0)
// 			  x = 0.5*(1.0 - c);
// 
// 			// Otherwise, it is the average of the previous and next points
// 			else
// 			  x = 0.5*(1.0 - a*c);
// 
// 			// If j is even, compute a normal Gauss-Lobatto point
// 			if (j%2 == 0)
// 			  y = 0.5*(1.0 - d);
// 
// 			// Otherwise, it is the average of the previous and next points
// 			else
// 			  y = 0.5*(1.0 - b*d);
// 
// 
// 			mesh.add_point (Point(x,y,0.), node_id++);
// 		      }


// 		    else
		          vt[0][node_id] = (static_cast<double>(i)/static_cast<double>(2*nx))*(xmax-xmin) + xmin;
                          vt[1][node_id] = (static_cast<double>(j)/static_cast<double>(2*ny))*(ymax-ymin) + ymin;
                          vt[2][node_id] = 0.;
			  
//			  std::cout << "inode: " << node_id <<  " x: " << vt[0][node_id] << "   y: " << vt[1][node_id] << std::endl;
			  
			  node_id++;
 
		}

 	      break;
 	    }

 	  default:
 	    {
 	      std::cout << "ERROR: Unrecognized or NotSupported 2D element type." << std::endl;
 	      exit(1);
 	    }
 	  }

 	  	    
	unsigned iel = 0;
	el= new elem(nel);
	el->SetElementGroupNumber(1);
 	// Build the elements.  Each one is a bit different.
 	switch (type)
 	  {
// 
// 	  case INVALID_ELEM:
// 	  case QUAD4:
// 	    {
// 	      for (unsigned int j=0; j<ny; j++)
// 		for (unsigned int i=0; i<nx; i++)
// 		  {
// 		    Elem* elem = mesh.add_elem(new Quad4);
// 
// 		    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
// 		    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+1,j)  );
// 		    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i+1,j+1));
// 		    elem->set_node(3) = mesh.node_ptr(idx(type,nx,i,j+1)  );
// 
// 		    if (j == 0)
// 		      mesh.boundary_info->add_side(elem, 0, 0);
// 
// 		    if (j == (ny-1))
// 		      mesh.boundary_info->add_side(elem, 2, 2);
// 
// 		    if (i == 0)
// 		      mesh.boundary_info->add_side(elem, 3, 3);
// 
// 		    if (i == (nx-1))
// 		      mesh.boundary_info->add_side(elem, 1, 1);
// 		  }
// 	      break;
// 	    }
// 
// 
// 	  case TRI3:
// 	    {
// 	      for (unsigned int j=0; j<ny; j++)
// 		for (unsigned int i=0; i<nx; i++)
// 		  {
// 		    Elem* elem = NULL;
// 
// 		    // Add first Tri3
// 		    elem = mesh.add_elem(new Tri3);
// 
// 		    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
// 		    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+1,j)  );
// 		    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i+1,j+1));
// 
// 		    if (j == 0)
// 		      mesh.boundary_info->add_side(elem, 0, 0);
// 
// 		    if (i == (nx-1))
// 		      mesh.boundary_info->add_side(elem, 1, 1);
// 
// 		    // Add second Tri3
// 		    elem = mesh.add_elem(new Tri3);
// 
// 		    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
// 		    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+1,j+1));
// 		    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i,j+1)  );
// 
// 		    if (j == (ny-1))
// 		      mesh.boundary_info->add_side(elem, 1, 2);
// 
// 		    if (i == 0)
// 		      mesh.boundary_info->add_side(elem, 2, 3);
// 		  }
// 	      break;
// 	    }
// 
// 
// 
// 	  case QUAD8:
	  case QUAD9:
	    {
	      
		unsigned LocalToGlobalNodePerElement[9];
		
		for (unsigned int j=0; j<(2*ny); j += 2)
 		  for (unsigned int i=0; i<(2*nx); i += 2)
 		  {
		  
                  el->SetElementGroup(iel,1);
		  el->SetElementMaterial(iel, 2); 
 		  type_elem_flag[3]=true;
                  el->AddToElementNumber(1,"Quad");
                  el->SetElementType(iel,3);
  
		  LocalToGlobalNodePerElement[0] = idx(type,nx,i,j) + 1; 
		  LocalToGlobalNodePerElement[1] = idx(type,nx,i+2,j) + 1;
		  LocalToGlobalNodePerElement[2] = idx(type,nx,i+2,j+2) + 1;
		  LocalToGlobalNodePerElement[3] = idx(type,nx,i,j+2) + 1;
		  LocalToGlobalNodePerElement[4] = idx(type,nx,i+1,j) + 1;
		  LocalToGlobalNodePerElement[5] = idx(type,nx,i+2,j+1) + 1;
		  LocalToGlobalNodePerElement[6] = idx(type,nx,i+1,j+2) + 1;
		  LocalToGlobalNodePerElement[7] = idx(type,nx,i,j+1) + 1;
		  LocalToGlobalNodePerElement[8] = idx(type,nx,i+1,j+1) + 1;
		  
                  // connectivity
                  for (unsigned iloc=0; iloc<9; iloc++) {
                    el->SetElementVertexIndex(iel,iloc,LocalToGlobalNodePerElement[iloc]);
                  }
                                   
                  if (j == 0) {
		    el->SetFaceElementIndex(iel,0,-2);
		    _boundaryinfo.insert( std::pair<unsigned int, std::string>(0,"bottom"));
		  }

                  if (j == 2*(ny-1)) {
		    el->SetFaceElementIndex(iel,2,-4);
		    _boundaryinfo.insert( std::pair<unsigned int, std::string>(2,"top"));
		  }
 
                  if (i == 0) {
		    el->SetFaceElementIndex(iel,3,-5);
		    _boundaryinfo.insert( std::pair<unsigned int, std::string>(3,"left"));
		  }

                  if (i == 2*(nx-1)) {
		    el->SetFaceElementIndex(iel,1,-3);
		    _boundaryinfo.insert(std::pair<unsigned int, std::string>(1,"right"));
		  }
                    
                  iel++;
               }
                  
	      break;
	    }

	  case TRI6:
	    {
	      
	      	unsigned LocalToGlobalNodePerElement[6];
		
		for (unsigned int j=0; j<(2*ny); j += 2)
 		  for (unsigned int i=0; i<(2*nx); i += 2)
 		  {
		  
		  // Add first Tri6  
                  el->SetElementGroup(iel,1);
		  el->SetElementMaterial(iel, 2); // 2 == fluid
 		  type_elem_flag[4]=true;
                  el->AddToElementNumber(1,"Triangle");
                  el->SetElementType(iel,4);

		
		  LocalToGlobalNodePerElement[0] = idx(type,nx,i,j) + 1;
		  LocalToGlobalNodePerElement[1] = idx(type,nx,i+2,j) + 1; 
		  LocalToGlobalNodePerElement[2] = idx(type,nx,i+2,j+2) + 1;
		  LocalToGlobalNodePerElement[3] = idx(type,nx,i+1,j) + 1; 
		  LocalToGlobalNodePerElement[4] = idx(type,nx,i+2,j+1) + 1;
		  LocalToGlobalNodePerElement[5] = idx(type,nx,i+1,j+1) + 1;

                  // connectivity
                  for (unsigned iloc=0; iloc<6; iloc++) {
                    el->SetElementVertexIndex(iel,iloc,LocalToGlobalNodePerElement[iloc]);
                  }
                    
                  if (j == 0) {
		    el->SetFaceElementIndex(iel,0,-2);
		    _boundaryinfo.insert( std::pair<unsigned int, std::string>(0,"bottom")); 
		  }
 
 		  if (i == 2*(nx-1)) {
		    el->SetFaceElementIndex(iel,1,-3);
		    _boundaryinfo.insert( std::pair<unsigned int, std::string>(1,"right"));
		  }

                  iel++;
		    
		    // Add second Tri6  
		    
		    el->SetElementGroup(iel,1);
		    el->SetElementMaterial(iel, 2); // 2 == fluid
 		    type_elem_flag[4]=true;
                    el->AddToElementNumber(1,"Triangle");
                    el->SetElementType(iel,4);
  
		    LocalToGlobalNodePerElement[0] = idx(type,nx,i,j) + 1;
		    LocalToGlobalNodePerElement[1] = idx(type,nx,i+2,j+2) + 1;
		    LocalToGlobalNodePerElement[2] = idx(type,nx,i,j+2) + 1;
		    LocalToGlobalNodePerElement[3] = idx(type,nx,i+1,j+1) + 1;
		    LocalToGlobalNodePerElement[4] = idx(type,nx,i+1,j+2) + 1;
		    LocalToGlobalNodePerElement[5] = idx(type,nx,i,j+1) + 1;
		    
		    // connectivity
                    for (unsigned iloc=0; iloc<6; iloc++) {
                      el->SetElementVertexIndex(iel,iloc,LocalToGlobalNodePerElement[iloc]);
                    }
	    
		    if (j == 2*(ny-1)) {
		      el->SetFaceElementIndex(iel,1,-4);
		      _boundaryinfo.insert( std::pair<unsigned int, std::string>(2,"top"));
		    }
		    
                    if (i == 0) {
		      el->SetFaceElementIndex(iel,2,-5);
		      _boundaryinfo.insert( std::pair<unsigned int, std::string>(3,"left"));
		    }
		    
		    iel++;
		    
		  }
	      
	      break;
	    }


	  default:
	    {
	      std::cout << "ERROR: Unrecognized or Not Supported 2D element type." << std::endl;
	      exit(1);
	    }
	  }
 
 
//         // Add sideset names to boundary info
//         mesh.boundary_info->sideset_name(0) = "bottom";
//         mesh.boundary_info->sideset_name(1) = "right";
//         mesh.boundary_info->sideset_name(2) = "top";
//         mesh.boundary_info->sideset_name(3) = "left";
// 
//         // Add nodeset names to boundary info
//         mesh.boundary_info->nodeset_name(0) = "bottom";
// 	mesh.boundary_info->nodeset_name(1) = "right";
// 	mesh.boundary_info->nodeset_name(2) = "top";
// 	mesh.boundary_info->nodeset_name(3) = "left";
// 
 	break;
       }

     //---------------------------------------------------------------------
     // Build a 3D mesh using hexahedral or prismatic elements.
     case 3:
       {
// 	libmesh_assert_not_equal_to (nx, 0);
// 	libmesh_assert_not_equal_to (ny, 0);
// 	libmesh_assert_not_equal_to (nz, 0);
// 	libmesh_assert_less (xmin, xmax);
// 	libmesh_assert_less (ymin, ymax);
// 	libmesh_assert_less (zmin, zmax);

	assert(nx!=0);
	assert(ny!=0);
	assert(nz!=0);
	assert(xmin<xmax);
	assert(ymin<ymax);
	assert(zmin<zmax);
 
 	// Reserve elements.  Meshes with prismatic elements require
 	// twice as many elements.
 	switch (type)
 	  {
// 	  case INVALID_ELEM:
// 	  case HEX8:
// 	  case HEX20:
 	  case HEX27:
// 	  case TET4:  // TET4's are created from an initial HEX27 discretization
// 	  case TET10: // TET10's are created from an initial HEX27 discretization
// 	  case PYRAMID5: // PYRAMID5's are created from an initial HEX27 discretization
	    {
	      nel    =  nx*ny*nz;
	      ngroup = 1;
	      nbcd   = 6;
	      mesh::_ref_index = pow(2,mesh::_dimension);
	      mesh::_face_index = pow(2,mesh::_dimension-1u);
              break;
            }
// 
// 	  case PRISM6:
// 	  case PRISM15:
// 	  case PRISM18:
// 	    {
// 	      mesh.reserve_elem(2*nx*ny*nz);
// 	      break;
// 	    }
// 
 	  default:
 	    {
 	      std::cerr << "ERROR: Unrecognized 3D element type." << std::endl;
 	      exit(1);
 	    }
 	  }


 	// Reserve nodes.  Quadratic elements need twice as many nodes as linear elements.
 	switch (type)
 	  {
// 	  case INVALID_ELEM:
// 	  case HEX8:
// 	  case PRISM6:
// 	    {
// 	      mesh.reserve_nodes( (nx+1)*(ny+1)*(nz+1) );
// 	      break;
// 	    }
// 
// 	  case HEX20:
 	  case HEX27:
// 	  case TET4: // TET4's are created from an initial HEX27 discretization
// 	  case TET10: // TET10's are created from an initial HEX27 discretization
// 	  case PYRAMID5: // PYRAMID5's are created from an initial HEX27 discretization
// 	  case PRISM15:
// 	  case PRISM18:
 	    {
// 	      // FYI: The resulting TET4 mesh will have exactly
// 	      // 5*(nx*ny*nz) + 2*(nx*ny + nx*nz + ny*nz) + (nx+ny+nz) + 1
// 	      // nodes once the additional mid-edge nodes for the HEX27 discretization
// 	      // have been deleted.
 	      
	      //mesh.reserve_nodes( (2*nx+1)*(2*ny+1)*(2*nz+1) );
	      nvt = (2*nx+1)*(2*ny+1)*(2*nz+1); 
 	      break;
 	    }

 	  default:
 	    {
 	      std::cerr << "ERROR: Unrecognized 3D element type." << std::endl;
 	      exit(1);
 	    }
 	  }
 
 
 	// Build the nodes.
        unsigned int node_id = 0;
	
	vt[0].resize(nvt);
        vt[1].resize(nvt);
        vt[2].resize(nvt);
	
 	switch (type)
 	  {
// 	  case INVALID_ELEM:
// 	  case HEX8:
// 	  case PRISM6:
// 	    {
// 	      for (unsigned int k=0; k<=nz; k++)
// 		for (unsigned int j=0; j<=ny; j++)
// 		  for (unsigned int i=0; i<=nx; i++)
// 		    {
// 		      if (gauss_lobatto_grid)
// 			{
// 			  mesh.add_point (Point(0.5*(1.0 - std::cos(libMesh::pi*static_cast<double>(i)/static_cast<double>(nx))),
// 						0.5*(1.0 - std::cos(libMesh::pi*static_cast<double>(j)/static_cast<double>(ny))),
// 						0.5*(1.0 - std::cos(libMesh::pi*static_cast<double>(k)/static_cast<double>(nz)))), node_id++);
// 			}
// 
// 		      else
// 			mesh.add_point(Point(static_cast<double>(i)/static_cast<double>(nx),
// 					     static_cast<double>(j)/static_cast<double>(ny),
// 					     static_cast<double>(k)/static_cast<double>(nz)), node_id++);
// 		    }
// 	      break;
// 	    }
// 
// 	  case HEX20:
 	  case HEX27:
// 	  case TET4: // TET4's are created from an initial HEX27 discretization
// 	  case TET10: // TET10's are created from an initial HEX27 discretization
// 	  case PYRAMID5: // PYRAMID5's are created from an initial HEX27 discretization
// 	  case PRISM15:
// 	  case PRISM18:
 	    {
	      for (unsigned int k=0; k<=(2*nz); k++)
 		for (unsigned int j=0; j<=(2*ny); j++)
 		  for (unsigned int i=0; i<=(2*nx); i++)
 		    {
// 		      if (gauss_lobatto_grid)
// 			{
// 			  // The x,y locations of the point.
// 			  double x=0., y=0., z=0.;
// 
// 			  // Shortcut quantities (do not depend on i,j)
// 			  const double a = std::cos( libMesh::pi / static_cast<double>(2*nx) );
// 			  const double b = std::cos( libMesh::pi / static_cast<double>(2*ny) );
// 
// 			  // Shortcut quantities (depend on i,j)
// 			  const double c = std::cos( libMesh::pi*i / static_cast<double>(2*nx) );
// 			  const double d = std::cos( libMesh::pi*j / static_cast<double>(2*ny) );
// 
// 			  // Additional shortcut quantities (for 3D)
// 			  const double e = std::cos( libMesh::pi / static_cast<double>(2*nz) );
// 			  const double f = std::cos( libMesh::pi*k / static_cast<double>(2*nz) );
// 
// 			  // If i is even, compute a normal Gauss-Lobatto point
// 			  if (i%2 == 0)
// 			    x = 0.5*(1.0 - c);
// 
// 			  // Otherwise, it is the average of the previous and next points
// 			  else
// 			    x = 0.5*(1.0 - a*c);
// 
// 			  // If j is even, compute a normal Gauss-Lobatto point
// 			  if (j%2 == 0)
// 			    y = 0.5*(1.0 - d);
// 
// 			  // Otherwise, it is the average of the previous and next points
// 			  else
// 			    y = 0.5*(1.0 - b*d);
// 
// 			  // If k is even, compute a normal Gauss-Lobatto point
// 			  if (k%2 == 0)
// 			    z = 0.5*(1.0 - f);
// 
// 			  // Otherwise, it is the average of the previous and next points
// 			  else
// 			    z = 0.5*(1.0 - e*f);
// 
// 
// 			  mesh.add_point (Point(x,y,z), node_id++);
// 			}
// 
// 		      else
// 			mesh.add_point(Point(static_cast<double>(i)/static_cast<double>(2*nx),
// 					     static_cast<double>(j)/static_cast<double>(2*ny),
// 					     static_cast<double>(k)/static_cast<double>(2*nz)), node_id++);
		          vt[0][node_id] = (static_cast<double>(i)/static_cast<double>(2*nx))*(xmax-xmin) + xmin;
                          vt[1][node_id] = (static_cast<double>(j)/static_cast<double>(2*ny))*(ymax-ymin) + ymin;
                          vt[2][node_id] = (static_cast<double>(k)/static_cast<double>(2*nz))*(zmax-zmin) + zmin;;
			  
//			  std::cout << "inode: " << node_id <<  " x: " << vt[0][node_id] << "   y: " << vt[1][node_id] << std::endl;
			  
			  node_id++;
 		    }
	      
 	      break;
 	    }
 
 
 	  default:
 	    {
 	      std::cerr << "ERROR: Unrecognized 3D element type." << std::endl;
 	      exit(1);
 	    }
 	  }
 
 
 
 	// Build the elements.
        unsigned iel = 0;
	el= new elem(nel);
	el->SetElementGroupNumber(1);
 	switch (type)
 	  {
// 	  case INVALID_ELEM:
// 	  case HEX8:
// 	    {
// 	      for (unsigned int k=0; k<nz; k++)
// 		for (unsigned int j=0; j<ny; j++)
// 		  for (unsigned int i=0; i<nx; i++)
// 		    {
// 		      Elem* elem = mesh.add_elem(new Hex8);
// 
// 		      elem->set_node(0) = mesh.node_ptr(idx(type,nx,ny,i,j,k)      );
// 		      elem->set_node(1) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k)    );
// 		      elem->set_node(2) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
// 		      elem->set_node(3) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k)    );
// 		      elem->set_node(4) = mesh.node_ptr(idx(type,nx,ny,i,j,k+1)    );
// 		      elem->set_node(5) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k+1)  );
// 		      elem->set_node(6) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
// 		      elem->set_node(7) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k+1)  );
// 
// 		      if (k == 0)
// 			mesh.boundary_info->add_side(elem, 0, 0);
// 
// 		      if (k == (nz-1))
// 			mesh.boundary_info->add_side(elem, 5, 5);
// 
// 		      if (j == 0)
// 			mesh.boundary_info->add_side(elem, 1, 1);
// 
// 		      if (j == (ny-1))
// 			mesh.boundary_info->add_side(elem, 3, 3);
// 
// 		      if (i == 0)
// 			mesh.boundary_info->add_side(elem, 4, 4);
// 
// 		      if (i == (nx-1))
// 			mesh.boundary_info->add_side(elem, 2, 2);
// 		    }
// 	      break;
// 	    }
// 
// 
// 
// 
// 	  case PRISM6:
// 	    {
// 	      for (unsigned int k=0; k<nz; k++)
// 		for (unsigned int j=0; j<ny; j++)
// 		  for (unsigned int i=0; i<nx; i++)
// 		    {
// 		      // First Prism
// 		      Elem* elem = NULL;
// 		      elem = mesh.add_elem(new Prism6);
// 
// 		      elem->set_node(0) = mesh.node_ptr(idx(type,nx,ny,i,j,k)      );
// 		      elem->set_node(1) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k)    );
// 		      elem->set_node(2) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k)    );
// 		      elem->set_node(3) = mesh.node_ptr(idx(type,nx,ny,i,j,k+1)    );
// 		      elem->set_node(4) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k+1)  );
// 		      elem->set_node(5) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k+1)  );
// 
// 		      // Add sides for first prism to boundary info object
// 		      if (i==0)
// 			mesh.boundary_info->add_side(elem, 3, 4);
// 
// 		      if (j==0)
// 			mesh.boundary_info->add_side(elem, 1, 1);
// 
// 		      if (k==0)
// 			mesh.boundary_info->add_side(elem, 0, 0);
// 
// 		      if (k == (nz-1))
// 			mesh.boundary_info->add_side(elem, 4, 5);
// 
// 		      // Second Prism
// 		      elem = mesh.add_elem(new Prism6);
// 
// 		      elem->set_node(0) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k)    );
// 		      elem->set_node(1) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
// 		      elem->set_node(2) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k)    );
// 		      elem->set_node(3) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k+1)  );
// 		      elem->set_node(4) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
// 		      elem->set_node(5) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k+1)  );
// 
// 		      // Add sides for second prism to boundary info object
// 		      if (i == (nx-1))
// 			mesh.boundary_info->add_side(elem, 1, 2);
// 
// 		      if (j == (ny-1))
// 			mesh.boundary_info->add_side(elem, 2, 3);
// 
// 		      if (k==0)
// 			mesh.boundary_info->add_side(elem, 0, 0);
// 
// 		      if (k == (nz-1))
// 			mesh.boundary_info->add_side(elem, 4, 5);
// 		    }
// 	      break;
// 	    }
// 
// 
// 
// 
// 
// 
// 	  case HEX20:
 	  case HEX27:
// 	  case TET4: // TET4's are created from an initial HEX27 discretization
// 	  case TET10: // TET10's are created from an initial HEX27 discretization
// 	  case PYRAMID5: // PYRAMID5's are created from an initial HEX27 discretization
 	    {
	      double LocalToGlobalNodePerElement[27];
	      
 	      for (unsigned int k=0; k<(2*nz); k += 2)
 		for (unsigned int j=0; j<(2*ny); j += 2)
 		  for (unsigned int i=0; i<(2*nx); i += 2)
 		    {
// 		      Elem* elem = (type == HEX20) ?
// 			mesh.add_elem(new Hex20) :
// 			mesh.add_elem(new Hex27);
		      
		      el->SetElementGroup(iel,1);
		      el->SetElementMaterial(iel, 2); 
 		      type_elem_flag[0]=true; // hex
		      type_elem_flag[3]=true; // quad face
                      el->AddToElementNumber(1,"Hex");
                      el->SetElementType(iel,0);
  
		      LocalToGlobalNodePerElement[0] = idx(type,nx,ny,i,  j,  k) + 1; 
		      LocalToGlobalNodePerElement[1] = idx(type,nx,ny,i+2,j,  k) + 1;
		      LocalToGlobalNodePerElement[2] = idx(type,nx,ny,i+2,j+2,k) + 1;
		      LocalToGlobalNodePerElement[3] = idx(type,nx,ny,i,  j+2,k) + 1;
		      LocalToGlobalNodePerElement[4] = idx(type,nx,ny,i,  j,  k+2) + 1;
		      LocalToGlobalNodePerElement[5] = idx(type,nx,ny,i+2,j,  k+2) + 1;
		      LocalToGlobalNodePerElement[6] = idx(type,nx,ny,i+2,j+2,k+2) + 1;
		      LocalToGlobalNodePerElement[7] = idx(type,nx,ny,i,  j+2,k+2) + 1;
		      LocalToGlobalNodePerElement[8] = idx(type,nx,ny,i+1,j,  k) + 1;
		      LocalToGlobalNodePerElement[9] = idx(type,nx,ny,i+2,j+1,k) + 1; 
		      LocalToGlobalNodePerElement[10] = idx(type,nx,ny,i+1,j+2,k) + 1;
		      LocalToGlobalNodePerElement[11] = idx(type,nx,ny,i,  j+1,k) + 1;
		      
		      // mid-point - up
		      LocalToGlobalNodePerElement[12] = idx(type,nx,ny,i+1,j,  k+2) + 1;
		      LocalToGlobalNodePerElement[13] = idx(type,nx,ny,i+2,j+1,k+2) + 1;
		      LocalToGlobalNodePerElement[14] = idx(type,nx,ny,i+1,j+2,k+2) + 1;
		      LocalToGlobalNodePerElement[15] = idx(type,nx,ny,i,  j+1,k+2) + 1 ;
		      
		      // vertices - middle
		      LocalToGlobalNodePerElement[16] = idx(type,nx,ny,i,  j,  k+1) + 1;
		      LocalToGlobalNodePerElement[17] = idx(type,nx,ny,i+2,j,  k+1) + 1;
		      LocalToGlobalNodePerElement[18] = idx(type,nx,ny,i+2,j+2,k+1) + 1;
		      LocalToGlobalNodePerElement[19] = idx(type,nx,ny,i,  j+2,k+1) + 1;
		      
		      // mid-point - middle
		      LocalToGlobalNodePerElement[20] = idx(type,nx,ny,i+1,j,  k+1) + 1;
		      LocalToGlobalNodePerElement[21] = idx(type,nx,ny,i+2,j+1,k+1) + 1;
		      LocalToGlobalNodePerElement[22] = idx(type,nx,ny,i+1,j+2,k+1) + 1;
		      LocalToGlobalNodePerElement[23] = idx(type,nx,ny,i,  j+1,k+1) + 1;
		      
		      // center - bottom
		      LocalToGlobalNodePerElement[24] = idx(type,nx,ny,i+1,j+1,k) + 1;
		      
		      // center - top
		      LocalToGlobalNodePerElement[25] = idx(type,nx,ny,i+1,j+1,k+2) + 1;
		      
		      // center - middle
		      LocalToGlobalNodePerElement[26] = idx(type,nx,ny,i+1,j+1,k+1) + 1;
		      
                      // connectivity
                      for (unsigned iloc=0; iloc<27; iloc++) {
                        el->SetElementVertexIndex(iel,iloc,LocalToGlobalNodePerElement[iloc]);
                      }
                      
                      if (k == 0) {
			el->SetFaceElementIndex(iel,0,-1);
			_boundaryinfo.insert( std::pair<unsigned int, std::string>(0,"bottom"));
		      }

 
 		      if (k == 2*(nz-1)) {
			el->SetFaceElementIndex(iel,5,-6);
			_boundaryinfo.insert( std::pair<unsigned int, std::string>(5,"top"));
		      }


 		      if (j == 0) {
			el->SetFaceElementIndex(iel,1,-2);
			_boundaryinfo.insert( std::pair<unsigned int, std::string>(1,"front"));
		      }

 
 		      if (j == 2*(ny-1)) {
			el->SetFaceElementIndex(iel,3,-4);
			_boundaryinfo.insert( std::pair<unsigned int, std::string>(3,"behind"));
		      }


 		      if (i == 0) {
			el->SetFaceElementIndex(iel,4,-5);
			_boundaryinfo.insert( std::pair<unsigned int, std::string>(4,"left"));
		      }


 		      if (i == 2*(nx-1)) {
			el->SetFaceElementIndex(iel,2,-3);
			_boundaryinfo.insert( std::pair<unsigned int, std::string>(2,"right"));
		      }

                      
		      iel++;
 
// 		      elem->set_node(0)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k)  );
// 		      elem->set_node(1)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k)  );
// 		      elem->set_node(2)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k)  );
// 		      elem->set_node(3)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k)  );
		      
// 		      elem->set_node(4)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+2));
// 		      elem->set_node(5)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+2));
// 		      elem->set_node(6)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+2));
// 		      elem->set_node(7)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+2));
// 		      elem->set_node(8)  = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k)  );
		      
// 		      elem->set_node(9)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k)  );
// 		      elem->set_node(10) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k)  );
// 		      elem->set_node(11) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k)  );
// 		      elem->set_node(12) = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+1));
		      
// 		      elem->set_node(13) = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+1));
// 		      elem->set_node(14) = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+1));
// 		      elem->set_node(15) = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+1));
// 		      elem->set_node(16) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+2));
// 		      elem->set_node(17) = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+2));
// 		      elem->set_node(18) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+2));
// 		      elem->set_node(19) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+2));
// 		      if ((type == HEX27) || (type == TET4) || (type == TET10) || (type == PYRAMID5))
// 			{
// 			  elem->set_node(20) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
// 			  elem->set_node(21) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+1));
// 			  elem->set_node(22) = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+1));
// 			  elem->set_node(23) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+1));
// 			  elem->set_node(24) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+1));
// 			  elem->set_node(25) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+2));
// 			  elem->set_node(26) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
// 			}
// 
// 
// 		      if (k == 0)
// 			mesh.boundary_info->add_side(elem, 0, 0);
// 
// 		      if (k == 2*(nz-1))
// 			mesh.boundary_info->add_side(elem, 5, 5);
// 
// 		      if (j == 0)
// 			mesh.boundary_info->add_side(elem, 1, 1);
// 
// 		      if (j == 2*(ny-1))
// 			mesh.boundary_info->add_side(elem, 3, 3);
// 
// 		      if (i == 0)
// 			mesh.boundary_info->add_side(elem, 4, 4);
// 
// 		      if (i == 2*(nx-1))
// 			mesh.boundary_info->add_side(elem, 2, 2);
  		    }
 	      break;
 	    }
// 
// 
// 
// 
// 	  case PRISM15:
// 	  case PRISM18:
// 	    {
// 	      for (unsigned int k=0; k<(2*nz); k += 2)
// 		for (unsigned int j=0; j<(2*ny); j += 2)
// 		  for (unsigned int i=0; i<(2*nx); i += 2)
// 		    {
// 		      // First Prism
// 		      Elem* elem = NULL;
// 		      elem = ((type == PRISM15) ?
// 			      mesh.add_elem(new Prism15) :
// 			      mesh.add_elem(new Prism18));
// 
// 		      elem->set_node(0)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k)  );
// 		      elem->set_node(1)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k)  );
// 		      elem->set_node(2)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k)  );
// 		      elem->set_node(3)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+2));
// 		      elem->set_node(4)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+2));
// 		      elem->set_node(5)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+2));
// 		      elem->set_node(6)  = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k)  );
// 		      elem->set_node(7)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
// 		      elem->set_node(8)  = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k)  );
// 		      elem->set_node(9)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+1));
// 		      elem->set_node(10) = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+1));
// 		      elem->set_node(11) = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+1));
// 		      elem->set_node(12) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+2));
// 		      elem->set_node(13) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+2));
// 		      elem->set_node(14) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+2));
// 		      if (type == PRISM18)
// 			{
// 			  elem->set_node(15) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+1));
// 			  elem->set_node(16) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
// 			  elem->set_node(17) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+1));
// 			}
// 
// 		      // Add sides for first prism to boundary info object
// 		      if (i==0)
// 			mesh.boundary_info->add_side(elem, 3, 4);
// 
// 		      if (j==0)
// 			mesh.boundary_info->add_side(elem, 1, 1);
// 
// 		      if (k==0)
// 			mesh.boundary_info->add_side(elem, 0, 0);
// 
// 		      if (k == 2*(nz-1))
// 			mesh.boundary_info->add_side(elem, 4, 5);
// 
// 
// 		      // Second Prism
// 		      elem = ((type == PRISM15) ?
// 			      mesh.add_elem(new Prism15) :
// 			      mesh.add_elem(new Prism18));
// 
// 		      elem->set_node(0)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,k)     );
// 		      elem->set_node(1)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k)   );
// 		      elem->set_node(2)  = mesh.node_ptr(idx(type,nx,ny,i,j+2,k)     );
// 		      elem->set_node(3)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,k+2)   );
// 		      elem->set_node(4)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+2) );
// 		      elem->set_node(5)  = mesh.node_ptr(idx(type,nx,ny,i,j+2,k+2)   );
// 		      elem->set_node(6)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k)  );
// 		      elem->set_node(7)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k)  );
// 		      elem->set_node(8)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
// 		      elem->set_node(9)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,k+1)  );
// 		      elem->set_node(10) = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+1));
// 		      elem->set_node(11) = mesh.node_ptr(idx(type,nx,ny,i,j+2,k+1)  );
// 		      elem->set_node(12) = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+2));
// 		      elem->set_node(13) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+2));
// 		      elem->set_node(14) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+2));
// 		      if (type == PRISM18)
// 			{
// 			  elem->set_node(15)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+1));
// 			  elem->set_node(16)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+1));
// 			  elem->set_node(17)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
// 			}
// 
// 		      // Add sides for second prism to boundary info object
// 		      if (i == 2*(nx-1))
// 			mesh.boundary_info->add_side(elem, 1, 2);
// 
// 		      if (j == 2*(ny-1))
// 			mesh.boundary_info->add_side(elem, 2, 3);
// 
// 		      if (k==0)
// 			mesh.boundary_info->add_side(elem, 0, 0);
// 
// 		      if (k == 2*(nz-1))
// 			mesh.boundary_info->add_side(elem, 4, 5);
// 
//  		    }
//  	      break;
//  	    }
// 

 	  default:
 	    {
 	      std::cerr << "ERROR: Unrecognized 3D element type." << std::endl;
 	      exit(1);
 	    }
 	  }

// 	//.......................................
// 	// Scale the nodal positions
// 	for (unsigned int p=0; p<mesh.n_nodes(); p++)
// 	  {
// 	    mesh.node(p)(0) = (mesh.node(p)(0))*(xmax-xmin) + xmin;
// 	    mesh.node(p)(1) = (mesh.node(p)(1))*(ymax-ymin) + ymin;
// 	    mesh.node(p)(2) = (mesh.node(p)(2))*(zmax-zmin) + zmin;
// 	  }
// 
// 
// 
// 
// 	// Additional work for tets and pyramids: we take the existing
// 	// HEX27 discretization and split each element into 24
// 	// sub-tets or 6 sub-pyramids.
// 	//
// 	// 24 isn't the minimum-possible number of tets, but it
// 	// obviates any concerns about the edge orientations between
// 	// the various elements.
// 	if ((type == TET4) ||
// 	    (type == TET10) ||
// 	    (type == PYRAMID5))
// 	  {
// 	    // Temporary storage for new elements. (24 tets per hex, 6 pyramids)
// 	    std::vector<Elem*> new_elements;
// 
// 	    if ((type == TET4) || (type == TET10))
// 	      new_elements.reserve(24*mesh.n_elem());
// 	    else
// 	      new_elements.reserve(6*mesh.n_elem());
// 
// 	    // Create tetrahedra or pyramids
// 	    {
// 	      MeshBase::element_iterator       el     = mesh.elements_begin();
// 	      const MeshBase::element_iterator end_el = mesh.elements_end();
// 
// 	      for ( ; el != end_el;  ++el)
// 		{
// 		  // Get a pointer to the HEX27 element.
// 		  Elem* base_hex = *el;
// 
// 		  // Get a pointer to the node located at the HEX27 centroid
// 		  Node* apex_node = base_hex->get_node(26);
// 
// 		  for (unsigned int s=0; s<base_hex->n_sides(); ++s)
// 		    {
// 		      // Get the boundary ID for this side
// 		      boundary_id_type b_id = mesh.boundary_info->boundary_id(*el, s);
// 
// 		      // Need to build the full-ordered side!
// 		      AutoPtr<Elem> side = base_hex->build_side(s);
// 
// 		      if ((type == TET4) || (type == TET10))
// 			{
// 			  // Build 4 sub-tets per side
// 			  for (unsigned int sub_tet=0; sub_tet<4; ++sub_tet)
// 			    {
// 			      new_elements.push_back( new Tet4 );
// 			      Elem* sub_elem = new_elements.back();
// 			      sub_elem->set_node(0) = side->get_node(sub_tet);
// 			      sub_elem->set_node(1) = side->get_node(8);                           // centroid of the face
// 			      sub_elem->set_node(2) = side->get_node(sub_tet==3 ? 0 : sub_tet+1 ); // wrap-around
// 			      sub_elem->set_node(3) = apex_node;                                   // apex node always used!
// 
// 			      // If the original hex was a boundary hex, add the new sub_tet's side
// 			      // 0 with the same b_id.  Note: the tets are all aligned so that their
// 			      // side 0 is on the boundary.
// 			      if (b_id != BoundaryInfo::invalid_id)
// 				mesh.boundary_info->add_side(sub_elem, 0, b_id);
// 			    }
// 			} // end if ((type == TET4) || (type == TET10))
// 
// 		      else // type==PYRAMID5
// 			{
// 			  // Build 1 sub-pyramid per side.
// 			  new_elements.push_back(new Pyramid5);
// 			  Elem* sub_elem = new_elements.back();
// 
// 			  // Set the base.  Note that since the apex is *inside* the base_hex,
// 			  // and the pyramid uses a counter-clockwise base numbering, we need to
// 			  // reverse the [1] and [3] node indices.
// 			  sub_elem->set_node(0) = side->get_node(0);
// 			  sub_elem->set_node(1) = side->get_node(3);
// 			  sub_elem->set_node(2) = side->get_node(2);
// 			  sub_elem->set_node(3) = side->get_node(1);
// 
// 			  // Set the apex
// 			  sub_elem->set_node(4) = apex_node;
// 
// 			  // If the original hex was a boundary hex, add the new sub_pyr's side
// 			  // 4 (the square base) with the same b_id.
// 			  if (b_id != BoundaryInfo::invalid_id)
// 			    mesh.boundary_info->add_side(sub_elem, 4, b_id);
// 			} // end else type==PYRAMID5
// 		    }
// 		}
// 	    }
// 
// 
// 	    // Delete the original HEX27 elements from the mesh, and the boundary info structure.
// 	    {
// 	      MeshBase::element_iterator       el     = mesh.elements_begin();
// 	      const MeshBase::element_iterator end_el = mesh.elements_end();
// 
// 	      for ( ; el != end_el;  ++el)
// 		{
// 		  mesh.boundary_info->remove(*el); // Safe even if *el has no boundary info.
// 		  mesh.delete_elem(*el);
// 		}
// 	    }
// 
// 	    // Add the new elements
// 	    for (unsigned int i=0; i<new_elements.size(); ++i)
// 	      mesh.add_elem(new_elements[i]);
// 
// 	  } // end if (type == TET4,TET10,PYRAMID5
// 
// 
// 	// Use all_second_order to convert the TET4's to TET10's
// 	if (type == TET10)
// 	  {
// 	    mesh.all_second_order();
// 	  }
// 
//         // Add sideset names to boundary info (Z axis out of the screen)
//         mesh.boundary_info->sideset_name(0) = "back";
//         mesh.boundary_info->sideset_name(1) = "bottom";
//         mesh.boundary_info->sideset_name(2) = "right";
//         mesh.boundary_info->sideset_name(3) = "top";
//         mesh.boundary_info->sideset_name(4) = "left";
//         mesh.boundary_info->sideset_name(5) = "front";
// 
//         // Add nodeset names to boundary info
//         mesh.boundary_info->nodeset_name(0) = "back";
// 	mesh.boundary_info->nodeset_name(1) = "bottom";
// 	mesh.boundary_info->nodeset_name(2) = "right";
// 	mesh.boundary_info->nodeset_name(3) = "top";
// 	mesh.boundary_info->nodeset_name(4) = "left";
// 	mesh.boundary_info->nodeset_name(5) = "front";

        break;
      } // end case dim==3

    default:
      {
	std::cout << " Error! " << std::endl;
      }
    }
    
    
    //*************** start reorder mesh dofs **************
  //(1)linear (2)quadratic (3)biquaratic
  
  vector <unsigned> dof_index;
  dof_index.resize(nvt);
  for(unsigned i=0;i<nvt;i++){
    dof_index[i]=i+1;
  }
  //reorder vertices and mid-points vs central points
  for (unsigned iel=0; iel<nel; iel++) {
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,1); inode++) {
      for (unsigned jel=0; jel<nel; jel++) {
	for (unsigned jnode=el->GetElementDofNumber(jel,1); jnode<el->GetElementDofNumber(jel,3); jnode++) { 
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
  for (unsigned iel=0; iel<nel; iel++) {
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
      for (unsigned jel=0; jel<nel; jel++) {
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
  for (unsigned iel=0; iel<nel; iel++) {
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,3); inode++) {
      unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
      el->SetElementVertexIndex(iel,inode,dof_index[ii]);
    }
  }
  vector <double> vt_temp;
  for(int i=0;i<3;i++){
    vt_temp=vt[i];
    for(unsigned j=0;j<nvt;j++){
      vt[i][dof_index[j]-1]=vt_temp[j];
    }
  }
  // **************  end reoreder mesh dofs **************
 
  el->SetNodeNumber(nvt);
 
  unsigned nv0=0;
  for (unsigned iel=0; iel<nel; iel++)
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
      unsigned i0=el->GetElementVertexIndex(iel,inode);
      if (nv0<i0) nv0=i0;
  }
  el->SetVertexNodeNumber(nv0);

  unsigned nv1=0;
  for (unsigned iel=0; iel<nel; iel++)
    for (unsigned inode=el->GetElementDofNumber(iel,0); inode<el->GetElementDofNumber(iel,1); inode++) {
      unsigned i1=el->GetElementVertexIndex(iel,inode);
      if (nv1<i1) nv1=i1;
  }
  el->SetMidpointNodeNumber(nv1-nv0);

  el->SetCentralNodeNumber(nvt-nv1);

  
  // connectivity: find all the element near the vertices
  BuildAdjVtx();
  Buildkel();
  
 if (_nprocs>=1) generate_metis_mesh_partition();
  vector <double> vt_temp2;
  for(int i=0;i<3;i++){
    vt_temp2=vt[i];
    for(unsigned j=0;j<nvt;j++) {
      vt[i][GetMetisDof(j,2)]=vt_temp2[j];
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
  
 
}


/**
 * This function searches all the elements around all the vertices
 **/
void mesh::BuildAdjVtx() {
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
void mesh::Buildkmid() {
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
void mesh::Buildkel() {
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
                if ((mesh::_dimension==3 &&
                     (i1==j1 || i1==j2 || i1==j3 ||  i1==j4 )&&
                     (i2==j1 || i2==j2 || i2==j3 ||  i2==j4 )&&
                     (i3==j1 || i3==j2 || i3==j3 ||  i3==j4 ))||
// 		   (DIM[1]==1 &&
                    (mesh::_dimension==2 &&
                     (i1==j1 || i1==j2 )&&
                     (i2==j1 || i2==j2 ))||
// 		   (DIM[0]==1 &&
                    (mesh::_dimension==1 &&
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


unsigned mesh::GetDimension() {
  return mesh::_dimension;
}

unsigned mesh::GetRefIndex() {
  return mesh::_ref_index;
}


/**
 * This function returns the number of mesh nodes for different type of elemets
 **/
unsigned mesh::GetDofNumber(const unsigned type) const {

// if(mesh::_dimension != 2) exit(1);
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
    return nel*(2+mesh::_dimension-1);
    break;
  }
  return 0;
}

unsigned mesh::GetElementNumber() const {
  return nel;
}
unsigned mesh::GetGridNumber() const {
  return grid;
}


/**
 * This function copies the refined element index vector in other_vector
 **/
void mesh::copy_elr(vector <unsigned> &other_vec) const {
  for (unsigned i=0; i<nel; i++)
    other_vec[i]=el->GetRefinedElementIndex(i);
}


void mesh::AllocateAndMarkStructureNode() {
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

unsigned mesh::GetEndIndex(const unsigned i) const{
  return _END_IND[i];
};

#include "map"
using std::map;



void mesh::GenerateVankaPartitions_FAST( const unsigned &block_size, vector < vector< unsigned > > &block_elements,
					 vector <unsigned> &block_element_type){
  unsigned iproc=processor_id();
  unsigned ElemOffset    = IS_Mts2Gmt_elem_offset[iproc];
  unsigned ElemOffsetp1  = IS_Mts2Gmt_elem_offset[iproc+1];
  unsigned OwnedElements = ElemOffsetp1 - ElemOffset;
  unsigned reminder = OwnedElements % block_size;
  
  unsigned nbloks = (0 == reminder)? OwnedElements/block_size : OwnedElements/block_size+1 ;
    
  vector < unsigned > size(nbloks,0);
  block_elements.resize(nbloks);
  
  for(int i = 0; i < block_elements.size(); i++) 
    block_elements[i].resize(block_size);
  if  (0 != reminder ){
    block_elements[ block_elements.size()-1].resize(reminder); 
  }
  
  for (unsigned iel=0; iel<OwnedElements; iel++) {
    block_elements[iel/block_size][iel%block_size]=iel+ElemOffset;
  }
  block_element_type.resize(2);
  block_element_type[0]=block_elements.size();
  block_element_type[1]=block_elements.size();
} 

void mesh::GenerateVankaPartitions_FSI( const unsigned &block_size, vector < vector< unsigned > > &block_elements,
					vector <unsigned> &block_element_type){

  unsigned iproc=processor_id();
  unsigned ElemOffset    = IS_Mts2Gmt_elem_offset[iproc];
  unsigned ElemOffsetp1  = IS_Mts2Gmt_elem_offset[iproc+1];
  unsigned OwnedElements = ElemOffsetp1 - ElemOffset;

  block_elements.resize(2); 
  block_elements[0].resize(OwnedElements); 
  block_elements[1].resize(OwnedElements); 
  
  block_element_type.resize(2);
  block_element_type[0]=1;
  block_element_type[1]=2;
  
  
  unsigned counter_f=0;
  unsigned counter_s=0;
  for (unsigned iel_mts = ElemOffset; iel_mts < ElemOffsetp1; iel_mts++) {
    unsigned kel        = IS_Mts2Gmt_elem[iel_mts]; 
    unsigned flag_mat   = el->GetElementMaterial(kel);
    if(2 == flag_mat){
      block_elements[0][counter_f]=iel_mts;
      counter_f++;
      
      
//       bool test_brd=0;
//       unsigned nve = el->GetElementDofNumber(kel,3);
//       for (unsigned i=0;i<nve;i++) {
// 	unsigned inode=el->GetElementVertexIndex(kel,i)-1u;
// 	if( el->GetNodeRegion(inode) ){
// 	  test_brd=1;
// 	  break;
// 	}
//       }
//       if(0==test_brd){
// 	block_elements[0][counter_f]=iel_mts;
// 	counter_f++;
//       }
//       else{
// 	block_elements[1][counter_s]=iel_mts;
// 	counter_s++;
//       }
    } 
    else{
      block_elements[1][counter_s]=iel_mts;
      counter_s++;
    }
  }  
  
  if(counter_s==0){
    block_elements.erase (block_elements.begin()+1);
    block_element_type[1]=1;
  }
  else{
    block_elements[1].resize(counter_s);
  }
  
  if(counter_f==0){
    block_elements.erase (block_elements.begin());
    block_element_type[0]=0;
    block_element_type[1]=1;
  }
  else{
    block_elements[0].resize(counter_f); 
  }
  
} 


void mesh::GenerateVankaPartitions_METIS( const unsigned &vnk_blck, vector < vector< unsigned > > &block_elements){
   
 
  
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
  //I call the mesh partioning function of Metis library (output is epart(own elem) and npart (own nodes))
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


