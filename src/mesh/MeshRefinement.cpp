/*=========================================================================

 Program: FEMUS
 Module: MeshRefinement
 Authors: Simone Bnà, Eugenio Aulisa
 
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
#include "MeshMetisPartitioning.hpp"
#include "MeshRefinement.hpp"
#include "NumericVector.hpp"


namespace femus {
  
//-------------------------------------------------------------------
MeshRefinement::MeshRefinement(Mesh& mesh): _mesh(mesh) {
  
}

//-------------------------------------------------------------------
MeshRefinement::~MeshRefinement() {
  
}
  
//-------------------------------------------------------------------
void MeshRefinement::FlagAllElementsToBeRefined() {
  
   _mesh.el->InitRefinedToZero();
   
   //refine all next grid elements
   for (unsigned iel=0; iel<_mesh.GetNumberOfElements(); iel++) {
     _mesh.el->SetRefinedElementIndex(iel,1);
     _mesh.el->AddToRefinedElementNumber(1);
     short unsigned elt=_mesh.el->GetElementType(iel);
     _mesh.el->AddToRefinedElementNumber(1,elt);
   }

   _mesh.el->AllocateChildrenElement(_mesh.GetRefIndex());
}

//-------------------------------------------------------------------
void MeshRefinement::FlagElementsToBeRefinedByUserDefinedFunction() {
  
    _mesh.el->InitRefinedToZero();
   
    //refine based on the function SetRefinementFlag defined in the main;
    // the Mesh is serial, we cannot in parallel use the coordinates to selectively refine
    std::vector<double> X_local;
    std::vector<double> Y_local;
    std::vector<double> Z_local;
    _mesh._coordinate->_Sol[0]->localize_to_all(X_local);
    _mesh._coordinate->_Sol[1]->localize_to_all(Y_local);
    _mesh._coordinate->_Sol[2]->localize_to_all(Z_local);
  
    for (unsigned iel=0; iel<_mesh.GetNumberOfElements(); iel+=1) {
      unsigned nve=_mesh.el->GetElementDofNumber(iel,0);
      std::vector < double > vtx(3,0.);
      for ( unsigned i=0; i<nve; i++) {
	unsigned inode=_mesh.el->GetElementVertexIndex(iel,i)-1u;
	unsigned inode_Metis=_mesh.GetMetisDof(inode,2);
	vtx[0]+= (*_mesh._coordinate->_Sol[0])(inode_Metis);
	vtx[1]+= (*_mesh._coordinate->_Sol[1])(inode_Metis);
	vtx[2]+= (*_mesh._coordinate->_Sol[2])(inode_Metis);
      }
      vtx[0]/=nve;
      vtx[1]/=nve;
      vtx[2]/=nve;
      if(!_mesh.el->GetRefinedElementIndex(iel)){
	if (_mesh._SetRefinementFlag(vtx,_mesh.el->GetElementGroup(iel),_mesh.GetLevel())) {
	  _mesh.el->SetRefinedElementIndex(iel,1);
	  _mesh.el->AddToRefinedElementNumber(1);
	  short unsigned elt=_mesh.el->GetElementType(iel);
	  _mesh.el->AddToRefinedElementNumber(1,elt);
	}
      }
    }
    _mesh.el->AllocateChildrenElement(_mesh.GetRefIndex());
}



//-------------------------------------------------------------------
void MeshRefinement::FlagElementsToBeRefinedByAMR() {
    
       
    if(_mesh._TestSetRefinementFlag){
      for (int iel_metis=_mesh.IS_Mts2Gmt_elem_offset[_iproc]; iel_metis < _mesh.IS_Mts2Gmt_elem_offset[_iproc+1]; iel_metis++) {
	unsigned kel = _mesh.IS_Mts2Gmt_elem[iel_metis];
	short unsigned kelt=_mesh.el->GetElementType(kel);
        unsigned nve=_mesh.el->GetElementDofNumber(kel,0);
	std::vector < double > vtx(3,0.);
	for(unsigned i=0; i<nve; i++) {
	  unsigned inode=_mesh.el->GetElementVertexIndex(kel,i)-1u;
	  unsigned inode_metis=_mesh.GetMetisDof(inode,2);
	  vtx[0]+= (*_mesh._coordinate->_Sol[0])(inode_metis);
	  vtx[1]+= (*_mesh._coordinate->_Sol[1])(inode_metis);
	  vtx[2]+= (*_mesh._coordinate->_Sol[2])(inode_metis);
	}
	vtx[0]/=nve;
	vtx[1]/=nve;
	vtx[2]/=nve;
	if( (*_mesh._coordinate->_Sol[3])(iel_metis) < 0.5 &&
	    _mesh._SetRefinementFlag(vtx,_mesh.el->GetElementGroup(kel),_mesh.GetLevel()) ) {
	    _mesh._coordinate->_Sol[3]->set(iel_metis,1.);
	}
      }
      _mesh._coordinate->_Sol[3]->close();
    }
       
    std::vector<double> AMR_local;
    _mesh._coordinate->_Sol[3]->localize_to_all(AMR_local);
  
    _mesh.el->InitRefinedToZero();
    
    for (unsigned iel_metis=0; iel_metis<_mesh.GetNumberOfElements(); iel_metis++) {
      if(AMR_local[iel_metis]>0.5){
	unsigned iel=_mesh.IS_Mts2Gmt_elem[iel_metis];
	_mesh.el->SetRefinedElementIndex(iel,1);
	_mesh.el->AddToRefinedElementNumber(1);
	short unsigned elt=_mesh.el->GetElementType(iel);
	_mesh.el->AddToRefinedElementNumber(1,elt);
      }   
    }
    _mesh.el->AllocateChildrenElement(_mesh.GetRefIndex());
}


//-------------------------------------------------------------------
void MeshRefinement::FlagOnlyEvenElementsToBeRefined() {
  
   _mesh.el->InitRefinedToZero();

   //refine all next grid even elements
   for (unsigned iel=0; iel<_mesh.GetNumberOfElements(); iel+=2) {
     _mesh.el->SetRefinedElementIndex(iel,1);
     _mesh.el->AddToRefinedElementNumber(1);
     short unsigned elt=_mesh.el->GetElementType(iel);
     _mesh.el->AddToRefinedElementNumber(1,elt);
   }
   _mesh.el->AllocateChildrenElement(_mesh.GetRefIndex());
}


//---------------------------------------------------------------------------------------------------------------
void MeshRefinement::RefineMesh(const unsigned & igrid, Mesh *mshc, const elem_type *otherFiniteElement[6][5]) {
  
  _mesh.SetCoarseMesh(mshc);
  
  _mesh.SetFiniteElementPtr(otherFiniteElement);
    
  elem *elc=mshc->el;
    
  _mesh.SetLevel(igrid);
  //_grid=igrid;


  int nelem=elc->GetRefinedElementNumber()*_mesh.GetRefIndex();
  _mesh.SetNumberOfElements(nelem);
  _mesh.el=new elem(elc,_mesh.GetRefIndex());

  unsigned jel=0;
  //divide each coarse element in 8(3D), 4(2D) or 2(1D) fine elemets and find all the vertices

  _mesh.el->SetElementGroupNumber(elc->GetElementGroupNumber());
  _mesh.el->SetNumberElementFather(elc->GetElementNumber());

  bool AMR = false;
  for (unsigned iel=0; iel<elc->GetElementNumber(); iel++) {
    if ( elc->GetRefinedElementIndex(iel) ) {
      elc->SetRefinedElementIndex(iel,jel+1u);
      unsigned elt=elc->GetElementType(iel);
      // project element type
      for (unsigned j=0; j<_mesh.GetRefIndex(); j++) {
        _mesh.el->SetElementType(jel+j,elt);
        _mesh.el->SetElementFather(jel+j,iel);
	elc->SetChildElement(iel,j,jel+j);
      }

      unsigned elg=elc->GetElementGroup(iel);
      unsigned elmat=elc->GetElementMaterial(iel);
      // project element group
//       for(unsigned j=0;j<REF_INDEX;j++){
      for (unsigned j=0; j<_mesh.GetRefIndex(); j++) {
        _mesh.el->SetElementGroup(jel+j,elg);
	_mesh.el->SetElementMaterial(jel+j,elmat);
      }

      // project vertex indeces
//       for(unsigned j=0;j<REF_INDEX;j++)
      for (unsigned j=0; j<_mesh.GetRefIndex(); j++)
        for (unsigned inode=0; inode<elc->GetElementDofNumber(iel,0); inode++)
          _mesh.el->SetElementVertexIndex(jel+j,inode,elc->GetElementVertexIndex(iel,fine2CoarseVertexMapping[elt][j][inode]-1u));
      // project face indeces
      for (unsigned iface=0; iface<elc->GetElementFaceNumber(iel); iface++) {
        int value=elc->GetFaceElementIndex(iel,iface);
        if (0>value)
// 	  for(unsigned jface=0;jface<FACE_INDEX;jface++)
          for (unsigned jface=0; jface<_mesh.GetFaceIndex(); jface++)
            _mesh.el->SetFaceElementIndex(jel+coarse2FineFaceMapping[elt][iface][jface][0],coarse2FineFaceMapping[elt][iface][jface][1], value);
      }
      // update element numbers
//       jel+=REF_INDEX;
      jel+=_mesh.GetRefIndex();
//       el->AddToElementNumber(REF_INDEX,elt);
      _mesh.el->AddToElementNumber(_mesh.GetRefIndex(),elt);
    }
    else {
      AMR=true;
    }
  }


  int ncoarsenodes=elc->GetNodeNumber();
  _mesh.SetNumberOfNodes(ncoarsenodes);
  _mesh.el->SetVertexNodeNumber(ncoarsenodes);
  int nnodes = _mesh.GetNumberOfNodes();
  
  //find all the elements near each vertex
  _mesh.BuildAdjVtx();
  //initialize to zero all the middle edge points
  for (unsigned iel=0; iel<_mesh.GetNumberOfElements(); iel++)
    for (unsigned inode=_mesh.el->GetElementDofNumber(iel,0); inode<_mesh.el->GetElementDofNumber(iel,1); inode++)
      _mesh.el->SetElementVertexIndex(iel,inode,0);
  //find all the middle edge points
  for (unsigned iel=0; iel<_mesh.GetNumberOfElements(); iel++) {
    unsigned ielt=_mesh.el->GetElementType(iel);
    unsigned istart=_mesh.el->GetElementDofNumber(iel,0);
    unsigned iend=_mesh.el->GetElementDofNumber(iel,1);
    for (unsigned inode=istart; inode<iend; inode++)
      if (0==_mesh.el->GetElementVertexIndex(iel,inode)) {
	nnodes++;
	_mesh.SetNumberOfNodes(nnodes);
        _mesh.el->SetElementVertexIndex(iel,inode,nnodes);
        unsigned im=_mesh.el->GetElementVertexIndex(iel,edge2VerticesMapping[ielt][inode-istart][0]);
        unsigned ip=_mesh.el->GetElementVertexIndex(iel,edge2VerticesMapping[ielt][inode-istart][1]);
        //find all the near elements which share the same middle edge point
        for (unsigned j=0; j<_mesh.el->GetVertexElementNumber(im-1u); j++) {
          unsigned jel=_mesh.el->GetVertexElementIndex(im-1u,j)-1u;
          if (jel>iel) {
            unsigned jm=0,jp=0;
            unsigned jelt=_mesh.el->GetElementType(jel);
            for (unsigned jnode=0; jnode<_mesh.el->GetElementDofNumber(jel,0); jnode++) {
              if (_mesh.el->GetElementVertexIndex(jel,jnode)==im) {
                jm=jnode+1u;
                break;
              }
            }
            if (jm!=0) {
              for (unsigned jnode=0; jnode<_mesh.el->GetElementDofNumber(jel,0); jnode++) {
                if (_mesh.el->GetElementVertexIndex(jel,jnode)==ip) {
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
                _mesh.el->SetElementVertexIndex(jel,vertices2EdgeMapping[jelt][--jm][--jp],nnodes);
              }
            }
          }
        }
      }
  }
  _mesh.el->SetMidpointNodeNumber(_mesh.GetNumberOfNodes() - _mesh.el->GetVertexNodeNumber());
  

  Buildkmid();
  
  _mesh.Buildkel();
  
  MeshMetisPartitioning meshmetispartitioning(_mesh);
  if( AMR == true ){
    meshmetispartitioning.DoPartition();
  }
  else{
    meshmetispartitioning.DoPartition(*mshc);
  }
  
  _mesh.FillISvector();
    
  // build Mesh coordinates by projecting the coarse coordinats
  _mesh._coordinate = new Solution(&_mesh);
  _mesh._coordinate->AddSolution("X",LAGRANGE,SECOND,1,0); 
  _mesh._coordinate->AddSolution("Y",LAGRANGE,SECOND,1,0); 
  _mesh._coordinate->AddSolution("Z",LAGRANGE,SECOND,1,0); 
  
  _mesh._coordinate->ResizeSolutionVector("X");
  _mesh._coordinate->ResizeSolutionVector("Y");
  _mesh._coordinate->ResizeSolutionVector("Z");    
  
  _mesh._coordinate->AddSolution("AMR",DISCONTINOUS_POLYNOMIAL,ZERO,1,0); 
  _mesh._coordinate->ResizeSolutionVector("AMR");
 
  unsigned solType=2;
     
  _mesh._coordinate->_Sol[0]->matrix_mult(*mshc->_coordinate->_Sol[0],*_mesh.GetCoarseToFineProjection(solType));
  _mesh._coordinate->_Sol[1]->matrix_mult(*mshc->_coordinate->_Sol[1],*_mesh.GetCoarseToFineProjection(solType));
  _mesh._coordinate->_Sol[2]->matrix_mult(*mshc->_coordinate->_Sol[2],*_mesh.GetCoarseToFineProjection(solType));
  _mesh._coordinate->_Sol[0]->close();
  _mesh._coordinate->_Sol[1]->close();
  _mesh._coordinate->_Sol[2]->close();     
         
}


/**
 * This function generates kmid for hex and wedge elements
 **/
void MeshRefinement::Buildkmid() {
  
  unsigned int nnodes = _mesh.GetNumberOfNodes();
  
  for (unsigned iel=0; iel<_mesh.el->GetElementNumber(); iel++)
    for (unsigned inode=_mesh.el->GetElementDofNumber(iel,1); inode<_mesh.el->GetElementDofNumber(iel,2); inode++)
      _mesh.el->SetElementVertexIndex(iel,inode,0);

  for (unsigned iel=0; iel<_mesh.el->GetElementNumber(); iel++) {
    for (unsigned iface=0; iface<_mesh.el->GetElementFaceNumber(iel,0); iface++) {
      unsigned inode=_mesh.el->GetElementDofNumber(iel,1)+iface;
      if ( 0==_mesh.el->GetElementVertexIndex(iel,inode) ) {
        _mesh.el->SetElementVertexIndex(iel,inode,++nnodes);
        unsigned i1=_mesh.el->GetFaceVertexIndex(iel,iface,0);
        unsigned i2=_mesh.el->GetFaceVertexIndex(iel,iface,1);
        unsigned i3=_mesh.el->GetFaceVertexIndex(iel,iface,2);
        for (unsigned j=0; j< _mesh.el->GetVertexElementNumber(i1-1u); j++) {
          unsigned jel= _mesh.el->GetVertexElementIndex(i1-1u,j)-1u;
          if (jel>iel) {
            for (unsigned jface=0; jface<_mesh.el->GetElementFaceNumber(jel,0); jface++) {
              unsigned jnode=_mesh.el->GetElementDofNumber(jel,1)+jface;
              if ( 0==_mesh.el->GetElementVertexIndex(jel,jnode) ) {
                unsigned j1=_mesh.el->GetFaceVertexIndex(jel,jface,0);
                unsigned j2=_mesh.el->GetFaceVertexIndex(jel,jface,1);
                unsigned j3=_mesh.el->GetFaceVertexIndex(jel,jface,2);
                unsigned j4=_mesh.el->GetFaceVertexIndex(jel,jface,3);
                if ((i1==j1 || i1==j2 || i1==j3 ||  i1==j4 )&&
                    (i2==j1 || i2==j2 || i2==j3 ||  i2==j4 )&&
                    (i3==j1 || i3==j2 || i3==j3 ||  i3==j4 )) {
                  _mesh.el->SetElementVertexIndex(jel,jnode,nnodes);
                }
              }
            }
          }
        }
      }
    }
  }

  for (unsigned iel=0; iel<_mesh.el->GetElementNumber(); iel++) {
    if (0==_mesh.el->GetElementType(iel)) {
      _mesh.el->SetElementVertexIndex(iel,26,++nnodes);
    }
    if (3==_mesh.el->GetElementType(iel)) {
      _mesh.el->SetElementVertexIndex(iel,8,++nnodes);
    }
  }
  _mesh.el->SetNodeNumber(nnodes);

  unsigned nv0= _mesh.el->GetVertexNodeNumber();
  unsigned nv1= _mesh.el->GetMidpointNodeNumber();
  _mesh.el->SetCentralNodeNumber(nnodes-nv0-nv1);
  
  _mesh.SetNumberOfNodes(nnodes);

}


}