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
void MeshRefinement::FlagElementsToBeRefined() {

    if(_mesh._IsUserRefinementFunctionDefined){
      for (int kel=_mesh._elementOffset[_iproc]; kel < _mesh._elementOffset[_iproc+1]; kel++) {
	if( _mesh.GetLevel() == 0 || _mesh.el->IsFatherRefined(kel)  ){
	  short unsigned kelt=_mesh.el->GetElementType(kel);
	  unsigned nve=_mesh.el->GetElementDofNumber(kel,0);
	  std::vector < double > vtx(3,0.);
	  for(unsigned i=0; i<nve; i++) {
	    //unsigned inode=_mesh.el->GetElementVertexIndex(kel,i)-1u;
	    unsigned inode_metis=_mesh.GetSolutionDof(i,kel,2);
	    vtx[0]+= (*_mesh._topology->_Sol[0])(inode_metis);
	    vtx[1]+= (*_mesh._topology->_Sol[1])(inode_metis);
	    vtx[2]+= (*_mesh._topology->_Sol[2])(inode_metis);
	  }
	  vtx[0]/=nve;
	  vtx[1]/=nve;
	  vtx[2]/=nve;
	  if( (*_mesh._topology->_Sol[3])(kel) < 0.5 &&
	      _mesh._SetRefinementFlag(vtx,_mesh.GetElementGroup(kel),_mesh.GetLevel()) ) {
	      _mesh._topology->_Sol[3]->set(kel,1.);
	  }
	}
      }
      _mesh._topology->_Sol[3]->close();
    }

    std::vector<double> AMR_local;
    _mesh._topology->_Sol[3]->localize_to_all(AMR_local);

    _mesh.el->InitRefinedToZero();

    for (unsigned iel=0; iel<_mesh.GetNumberOfElements(); iel++) {
      if(AMR_local[iel]>0.5){
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

  for (int iel=_mesh._elementOffset[_iproc]; iel < _mesh._elementOffset[_iproc+1]; iel++) {
    if( _mesh.GetLevel() == 0 || _mesh.el->IsFatherRefined(iel)){
      if( (*_mesh._topology->_Sol[3])(iel) < 0.5 && iel%2 == 0) {
	_mesh._topology->_Sol[3]->set(iel,1.);
      }
    }
  }

  std::vector<double> AMR_local;
  _mesh._topology->_Sol[3]->localize_to_all(AMR_local);

  _mesh.el->InitRefinedToZero();

  for (unsigned iel = 0; iel < _mesh.GetNumberOfElements(); iel++) {
    if(AMR_local[iel]>0.5){
      _mesh.el->SetRefinedElementIndex(iel,1);
      _mesh.el->AddToRefinedElementNumber(1);
      short unsigned elt=_mesh.el->GetElementType(iel);
      _mesh.el->AddToRefinedElementNumber(1,elt);
    }
  }
  _mesh.el->AllocateChildrenElement(_mesh.GetRefIndex());
}


//---------------------------------------------------------------------------------------------------------------
void MeshRefinement::RefineMesh(const unsigned & igrid, Mesh *mshc, const elem_type *otherFiniteElement[6][5]) {

  _mesh.SetCoarseMesh(mshc);

  _mesh.SetFiniteElementPtr(otherFiniteElement);

  elem *elc = mshc->el;

  _mesh.SetLevel(igrid);

  // total number of elements on the fine level
  int nelem = elc->GetRefinedElementNumber() * _mesh.GetRefIndex(); // refined
  nelem += elc->GetElementNumber() - elc->GetRefinedElementNumber(); // not-refined

  _mesh.SetNumberOfElements(nelem);
  _mesh.el = new elem( elc, _mesh.GetRefIndex() );

  unsigned jel=0;
  //divide each coarse element in 8(3D), 4(2D) or 2(1D) fine elements and find all the vertices

  _mesh.el->SetElementGroupNumber(elc->GetElementGroupNumber());
  _mesh.el->SetNumberElementFather(elc->GetElementNumber()); // setta il num di elementi padre per il mesh fine

  bool AMR = false;
  for (unsigned iel=0; iel<elc->GetElementNumber(); iel++) {
    if ( elc->GetRefinedElementIndex(iel) ) {
      elc->SetRefinedElementIndex(iel,jel+1u);
      unsigned elt=elc->GetElementType(iel);
      // project element type
      for (unsigned j=0; j<_mesh.GetRefIndex(); j++) {
        _mesh.el->SetElementType(jel+j,elt);
        _mesh.el-> SetIfFatherIsRefined(jel+j, true);
	elc->SetChildElement(iel,j,jel+j);
      }

      //unsigned elg=elc->GetElementGroup(iel);
      //unsigned elmat=elc->GetElementMaterial(iel);
      // project element group
      //for (unsigned j=0; j<_mesh.GetRefIndex(); j++) {
        //_mesh.el->SetElementGroup(jel+j,elg);
	//_mesh.el->SetElementMaterial(jel+j,elmat);
      //}

      // project vertex indeces
      for (unsigned j=0; j<_mesh.GetRefIndex(); j++)
        for (unsigned inode=0; inode<elc->GetElementDofNumber(iel,0); inode++)
          _mesh.el->SetElementVertexIndex(jel+j,inode,elc->GetElementVertexIndex(iel,fine2CoarseVertexMapping[elt][j][inode]-1u));
      // project face indeces
      for (unsigned iface=0; iface<elc->GetElementFaceNumber(iel); iface++) {
        int value=elc->GetFaceElementIndex(iel,iface);
        if (0>value)
          for (unsigned jface=0; jface<_mesh.GetFaceIndex(); jface++)
            _mesh.el->SetFaceElementIndex(jel+coarse2FineFaceMapping[elt][iface][jface][0],coarse2FineFaceMapping[elt][iface][jface][1], value);
      }
      // update element numbers
      jel+=_mesh.GetRefIndex();
      _mesh.el->AddToElementNumber(_mesh.GetRefIndex(),elt);
    }
    else {
      AMR = true;
      //elc->SetRefinedElementIndex(iel,jel+1u); // to understand
      unsigned elt=elc->GetElementType(iel);

      // project element type
      _mesh.el->SetElementType(jel,elt);
      _mesh.el-> SetIfFatherIsRefined(jel, false);
      elc->SetChildElement(iel,0,jel);

      //unsigned elg = elc->GetElementGroup(iel);
      //unsigned elmat = elc->GetElementMaterial(iel);

      // project element group
      //_mesh.el->SetElementGroup(jel,elg);
      //_mesh.el->SetElementMaterial(jel,elmat);

      // project nodes indeces
      for (unsigned inode=0; inode<elc->GetElementDofNumber(iel,2); inode++)
        _mesh.el->SetElementVertexIndex(jel,inode,elc->GetElementVertexIndex(iel,inode));

      // project face indeces
      for (unsigned iface=0; iface<elc->GetElementFaceNumber(iel); iface++) {
        int value = elc->GetFaceElementIndex(iel,iface);
        if ( 0 > value){
          _mesh.el->SetFaceElementIndex(jel, iface, value);
	}
      }
      // update element numbers
      jel++;
      _mesh.el->AddToElementNumber(1,elt);
    }
  }

  int nnodes = elc->GetNodeNumber();
  _mesh.SetNumberOfNodes(nnodes);
  _mesh.el->SetNodeNumber(nnodes);

  //find all the elements near each vertex
  _mesh.BuildAdjVtx(); //TODO
  //initialize to zero all the middle edge points
  for (unsigned iel=0; iel<_mesh.GetNumberOfElements(); iel++){
    if ( _mesh.el->IsFatherRefined( iel ) ) {
      for (unsigned inode=_mesh.el->GetElementDofNumber(iel,0); inode<_mesh.el->GetElementDofNumber(iel,1); inode++){
	_mesh.el->SetElementVertexIndex(iel,inode,0);
      }
    }
  }
  //find all the middle edge points
  for (unsigned iel=0; iel<_mesh.GetNumberOfElements(); iel++) {
    if ( _mesh.el->IsFatherRefined(iel) ) {
      unsigned ielt=_mesh.el->GetElementType(iel);
      unsigned istart=_mesh.el->GetElementDofNumber(iel,0);
      unsigned iend=_mesh.el->GetElementDofNumber(iel,1);
      for (unsigned inode=istart; inode<iend; inode++) {
	if (0 == _mesh.el->GetElementVertexIndex(iel, inode)) {
	  nnodes++;
	  _mesh.el->SetElementVertexIndex(iel, inode, nnodes);
	  unsigned im=_mesh.el->GetElementVertexIndex(iel,edge2VerticesMapping[ielt][inode-istart][0]);
	  unsigned ip=_mesh.el->GetElementVertexIndex(iel,edge2VerticesMapping[ielt][inode-istart][1]);
	  //find all the near elements which share the same middle edge point
	  for (unsigned j=0; j<_mesh.el->GetVertexElementNumber(im-1u); j++) {
	    unsigned jel=_mesh.el->GetVertexElementIndex(im-1u,j)-1u;
	    if ( _mesh.el->IsFatherRefined(jel) && jel > iel ) {  // to skip coarse elements
	      unsigned jm=0,jp=0;
	      unsigned jelt=_mesh.el->GetElementType(jel);
	      for (unsigned jnode=0; jnode<_mesh.el->GetElementDofNumber(jel,0); jnode++) {
		if (_mesh.el->GetElementVertexIndex(jel,jnode)==im) {
		  jm=jnode+1u;
		  break;
		}
	      }
	      if (jm!=0) { //TODO this can be changed and put inside (by Sara)
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
    }
  }
  _mesh.SetNumberOfNodes(nnodes);
  _mesh.el->SetNodeNumber(nnodes);

  Buildkmid();

  std::vector < int > partition;
  partition.reserve(_mesh.GetNumberOfNodes());
  partition.resize(_mesh.GetNumberOfElements());
   
  MeshMetisPartitioning meshMetisPartitioning(_mesh);
  if( AMR == true ){
    meshMetisPartitioning.DoPartition(partition, AMR);
  }
  else{
    meshMetisPartitioning.DoPartition(partition, *mshc);
  }
  _mesh.FillISvector(partition);
  partition.resize(0);

  _mesh.BuildAdjVtx(); //TODO

  _mesh.Buildkel();



  // build Mesh coordinates by projecting the coarse coordinats
  _mesh._topology = new Solution(&_mesh);
  _mesh._topology->AddSolution("X",LAGRANGE,SECOND,1,0);
  _mesh._topology->AddSolution("Y",LAGRANGE,SECOND,1,0);
  _mesh._topology->AddSolution("Z",LAGRANGE,SECOND,1,0);

  _mesh._topology->ResizeSolutionVector("X");
  _mesh._topology->ResizeSolutionVector("Y");
  _mesh._topology->ResizeSolutionVector("Z");

  _mesh._topology->AddSolution("AMR",DISCONTINOUS_POLYNOMIAL,ZERO,1,0);
  _mesh._topology->ResizeSolutionVector("AMR");

  unsigned solType=2;

  _mesh._topology->_Sol[0]->matrix_mult(*mshc->_topology->_Sol[0],*_mesh.GetCoarseToFineProjection(solType));
  _mesh._topology->_Sol[1]->matrix_mult(*mshc->_topology->_Sol[1],*_mesh.GetCoarseToFineProjection(solType));
  _mesh._topology->_Sol[2]->matrix_mult(*mshc->_topology->_Sol[2],*_mesh.GetCoarseToFineProjection(solType));
  _mesh._topology->_Sol[0]->close();
  _mesh._topology->_Sol[1]->close();
  _mesh._topology->_Sol[2]->close();
  
  _mesh._topology->AddSolution("Material", DISCONTINOUS_POLYNOMIAL, ZERO, 1 , 0);
  _mesh._topology->ResizeSolutionVector("Material");
  NumericVector &materialf =  _mesh._topology->GetSolutionName("Material");
  NumericVector &materialc =   mshc->_topology->GetSolutionName("Material");
  materialf.matrix_mult(materialc, *_mesh.GetCoarseToFineProjection(3));
  materialf.close(); 
  
  _mesh._topology->AddSolution("Group", DISCONTINOUS_POLYNOMIAL, ZERO, 1 , 0);
  _mesh._topology->ResizeSolutionVector("Group");
  NumericVector &groupf =  _mesh._topology->GetSolutionName("Group");
  NumericVector &groupc =   mshc->_topology->GetSolutionName("Group");
  groupf.matrix_mult(groupc, *_mesh.GetCoarseToFineProjection(3));
  groupf.close();
  
  _mesh._topology->AddSolution("Type", DISCONTINOUS_POLYNOMIAL, ZERO, 1 , 0);
  _mesh._topology->ResizeSolutionVector("Type");
  NumericVector &typef =  _mesh._topology->GetSolutionName("Type");
  NumericVector &typec =   mshc->_topology->GetSolutionName("Type");
  typef.matrix_mult(typec, *_mesh.GetCoarseToFineProjection(3));
  typef.close();
  
}


/**
 * This function generates face (for hex and wedge elements) and element (for hex and quad) dofs
 **/
void MeshRefinement::Buildkmid() {

  unsigned int nnodes = _mesh.GetNumberOfNodes();

  //intialize to zero
  for (unsigned iel=0; iel<_mesh.el->GetElementNumber(); iel++){
    if ( _mesh.el->IsFatherRefined(iel) ) {
      for (unsigned inode=_mesh.el->GetElementDofNumber(iel,1); inode<_mesh.el->GetElementDofNumber(iel,2); inode++){
	_mesh.el->SetElementVertexIndex(iel,inode,0);
      }
    }
  }

  // generate face dofs for hex and wedge elements
  for (unsigned iel=0; iel<_mesh.el->GetElementNumber(); iel++) {
    if ( _mesh.el->IsFatherRefined(iel) ) {
      for (unsigned iface=0; iface<_mesh.el->GetElementFaceNumber(iel,0); iface++) { // I think is on all the faces that are quads
	unsigned inode=_mesh.el->GetElementDofNumber(iel,1)+iface;
	if ( 0 == _mesh.el->GetElementVertexIndex(iel,inode) ) {
	  _mesh.el->SetElementVertexIndex(iel,inode,++nnodes);
	  unsigned i1=_mesh.el->GetFaceVertexIndex(iel,iface,0);
	  unsigned i2=_mesh.el->GetFaceVertexIndex(iel,iface,1);
	  unsigned i3=_mesh.el->GetFaceVertexIndex(iel,iface,2);
	  for (unsigned j=0; j< _mesh.el->GetVertexElementNumber(i1-1u); j++) {
	    unsigned jel= _mesh.el->GetVertexElementIndex(i1-1u,j)-1u;
	    if ( _mesh.el->IsFatherRefined(jel) && jel > iel ) {
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
  }

  // generates element dofs for hex and quad elements
  for (unsigned iel=0; iel<_mesh.el->GetElementNumber(); iel++) {
    if ( _mesh.el->IsFatherRefined(iel) ) {
      if (0==_mesh.el->GetElementType(iel)) { //hex
	_mesh.el->SetElementVertexIndex(iel,26,++nnodes);
      }
      if (3==_mesh.el->GetElementType(iel)) {//quad
	_mesh.el->SetElementVertexIndex(iel,8,++nnodes);
      }
    }
  }
  _mesh.el->SetNodeNumber(nnodes);
  _mesh.SetNumberOfNodes(nnodes);

}


}