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

#include "MeshRefinement.hpp"
#include "Mesh.hpp"
#include "MeshMetisPartitioning.hpp"

#include <climits>


namespace femus {

//-------------------------------------------------------------------
  MeshRefinement::MeshRefinement(Mesh& mesh): _mesh(mesh) {

  }

//-------------------------------------------------------------------
  MeshRefinement::~MeshRefinement() {

  }

//-------------------------------------------------------------------
  void MeshRefinement::FlagAllElementsToBeRefined() {
    FlagElementsToRefine(0);
  }

//-------------------------------------------------------------------
  bool MeshRefinement::FlagElementsToBeRefined(const double & treshold, NumericVector& error) {
    return FlagElementsToRefineBasedOnError(treshold, error);
  }
//-------------------------------------------------------------------
  bool MeshRefinement::FlagElementsToBeRefined() {
    FlagElementsToRefine(1);
    return true;
  }
//-------------------------------------------------------------------
  void MeshRefinement::FlagOnlyEvenElementsToBeRefined() {
    FlagElementsToRefine(2);
  }
//-------------------------------------------------------------------
  void MeshRefinement::FlagElementsToRefine(const unsigned& type) {

    //BEGIN temporary parallel vector initialization
    NumericVector* numberOfRefinedElement;
    numberOfRefinedElement = NumericVector::build().release();

    if(_nprocs == 1) numberOfRefinedElement->init(_nprocs, 1, false, SERIAL);
    else numberOfRefinedElement->init(_nprocs, 1, false, PARALLEL);

    numberOfRefinedElement->zero();
    //END temporary parallel vector initialization

    //BEGIN flag element to be refined
    if(type == 0) {   // Flag all element
      for(int iel = _mesh._elementOffset[_iproc]; iel < _mesh._elementOffset[_iproc + 1]; iel++) {
        if(_mesh.el->GetIfElementCanBeRefined(iel)) {
          _mesh.GetTopology()->_Sol[_mesh.GetAmrIndex()]->set(iel, 1.);
          numberOfRefinedElement->add(_iproc, 1.);
        }
      }
    }
    else if(type == 1) {   // Flag AMR elements
      for(int iel = _mesh._elementOffset[_iproc]; iel < _mesh._elementOffset[_iproc + 1]; iel++) {
        if(_mesh.el->GetIfElementCanBeRefined(iel)) {
          if((*_mesh.GetTopology()->_Sol[ _mesh.GetAmrIndex() ])(iel) > 0.5) {
            numberOfRefinedElement->add(_iproc, 1.);
          }
          else if( _mesh.get_is_refinement_function_defined() ) {
            short unsigned ielt = _mesh.GetElementType(iel);
            unsigned nve = _mesh.GetElementDofNumber(iel, 0);
            std::vector < double > x(3, 0.);

            for(unsigned i = 0; i < nve; i++) {
              unsigned inode_metis = _mesh.GetSolutionDof(i, iel, 2);
              x[0] += (*_mesh.GetTopology()->_Sol[0])(inode_metis);
              x[1] += (*_mesh.GetTopology()->_Sol[1])(inode_metis);
              x[2] += (*_mesh.GetTopology()->_Sol[2])(inode_metis);
            }

            x[0] /= nve;
            x[1] /= nve;
            x[2] /= nve;

            if(_mesh._SetRefinementFlag(x, _mesh.GetElementGroup(iel), _mesh.GetLevel())) {
              _mesh.GetTopology()->_Sol[ _mesh.GetAmrIndex() ]->set(iel, 1.);
              numberOfRefinedElement->add(_iproc, 1.);
            }
          }
        }
        else {
          _mesh.GetTopology()->_Sol[ _mesh.GetAmrIndex() ]->set(iel, 0.);
        }
      }
    }
    else if(type == 2) {   // Flag only even elements (for debugging purposes)
      for(int iel = _mesh._elementOffset[_iproc]; iel < _mesh._elementOffset[_iproc + 1]; iel++) {
        if(_mesh.el->GetIfElementCanBeRefined(iel)) {
          if((*_mesh.GetTopology()->_Sol[_mesh.GetAmrIndex()])(iel) < 0.5 && iel % 2 == 0) {
            _mesh.GetTopology()->_Sol[_mesh.GetAmrIndex()]->set(iel, 1.);
            numberOfRefinedElement->add(_iproc, 1.);
          }
        }
      }
    }

    _mesh.GetTopology()->_Sol[_mesh.GetAmrIndex()]->close();
    //END flag element to be refined

    //BEGIN update elem
    numberOfRefinedElement->close();
    double totalNumber = numberOfRefinedElement->l1_norm();
    _mesh.el->SetRefinedElementNumber(static_cast < unsigned >(totalNumber + 0.25));
    delete numberOfRefinedElement;

    //END update elem
  }


  bool MeshRefinement::FlagElementsToRefineBasedOnError(const double& treshold, NumericVector& error) {

    unsigned type = 2;

    //BEGIN temporary parallel vector initialization
    NumericVector* numberOfRefinedElement;
    numberOfRefinedElement = NumericVector::build().release();

    if(_nprocs == 1) numberOfRefinedElement->init(_nprocs, 1, false, SERIAL);
    else numberOfRefinedElement->init(_nprocs, 1, false, PARALLEL);

    numberOfRefinedElement->zero();
    //END temporary parallel vector initialization


    //BEGIN flag element to be refined
    for(int iel = _mesh._elementOffset[_iproc]; iel < _mesh._elementOffset[_iproc + 1]; iel++) {
      if(_mesh.el->GetIfElementCanBeRefined(iel)) {
        if((*_mesh.GetTopology()->_Sol[_mesh.GetAmrIndex()])(iel) < 0.5 && error(iel) > treshold) {
          _mesh.GetTopology()->_Sol[_mesh.GetAmrIndex()]->set(iel, 1.);
          numberOfRefinedElement->add(_iproc, 1.);
        }
      }
    }

    _mesh.GetTopology()->_Sol[_mesh.GetAmrIndex()]->close();
    //END flag element to be refined

    //BEGIN update elem
    numberOfRefinedElement->close();
    double totalNumber = numberOfRefinedElement->l1_norm();
    _mesh.el->SetRefinedElementNumber(static_cast < unsigned >(totalNumber + 0.25));
    delete numberOfRefinedElement;

    bool elementsHaveBeenRefined = true;
    if(totalNumber < _nprocs) {
      elementsHaveBeenRefined = false;
    }
    return elementsHaveBeenRefined;
    //END update elem
  }









//---------------------------------------------------------------------------------------------------------------
/// 
void MeshRefinement::RefineMesh(const unsigned& igrid, Mesh* mshc, /*const*/ elem_type* otherFiniteElement[N_GEOM_ELS][NFE_FAMS]) {


//==== Equivalent of: ReadCoarseMeshBeforePartitioning - BEGIN ======== 

      
//==== Level - BEGIN ==============================
    _mesh.SetLevel(igrid);
//==== Level - END ==============================
    
    
//==== AMR - BEGIN ==============================
    _mesh.SetIfHomogeneous(true);
//==== AMR - END ==============================

    
//==== info from the coarse mesh - BEGIN ==============================
    
    const unsigned dim = mshc->GetDimension();
    
    // geom el, refinement - BEGIN ********************
    _mesh.SetRefinementCellAndFaceIndices(dim);
    // geom el, refinement - END ******************** 

    _mesh.get_prol_matrices().SetCoarseMesh(mshc);

    elem* elc = mshc->el;

    _mesh.SetFiniteElementPtr(otherFiniteElement);
    
//==== info from the coarse mesh - END ==============================

    
//====== BEGIN ELEMENTS  ==============================

    // total number of elements on the fine level
     int nelem = elem::InitializeNumberOfElementsFromCoarseList(elc, _mesh.GetRefIndex() );

    _mesh.SetNumberOfElements(nelem);

    unsigned elementOffsetCoarse   = mshc->_elementOffset[_iproc];
    unsigned elementOffsetCoarseP1 = mshc->_elementOffset[_iproc + 1];

    std::vector < double > coarseLocalizedAmrVector;
    mshc->GetTopology()->_Sol[mshc->GetAmrIndex()]->localize_to_all(coarseLocalizedAmrVector);

    mshc->el->AllocateChildrenElement   (_mesh.GetRefIndex(), mshc);
    mshc->el->AllocateChildrenElementDof(_mesh.GetRefIndex(), mshc);

    _mesh.el = new elem(elc, mshc->GetDimension(), _mesh.GetRefIndex(), coarseLocalizedAmrVector);


    unsigned jel = 0;
    //divide each coarse element in 8(3D), 4(2D) or 2(1D) fine elements and find all the vertices

    _mesh.el->SetElementGroupNumber(elc->GetElementGroupNumber());

    bool AMR = false;

    
    
    for(unsigned isdom = 0; isdom < _nprocs; isdom++) {
        
      elc->LocalizeElementDof(isdom);
      elc->LocalizeElementNearFace(isdom);
      elc->LocalizeElement_Level_Type_Group_Material(isdom);
      
      for(unsigned iel = mshc->_elementOffset[isdom]; iel < mshc->_elementOffset[isdom + 1]; iel++) {
          
        if(static_cast < unsigned short >(coarseLocalizedAmrVector[iel] + 0.25) == 1) {
            
          unsigned elt = elc->GetElementType(iel);
          // project element type, group, material; child element -----------------
          for(unsigned j = 0; j < _mesh.GetRefIndex(); j++) {
            _mesh.el->SetElementType(jel + j, elc->GetElementType(iel));
            _mesh.el->SetElementGroup(jel + j, elc->GetElementGroup(iel));
            
	    unsigned gr_mat = elc->GetElementMaterial(iel);
	    _mesh.el->SetElementMaterial(jel + j, gr_mat);
        
	        
            _mesh.el->SetElementLevel(jel + j, elc->GetElementLevel(iel) + 1);
            if(iel >= elementOffsetCoarse && iel < elementOffsetCoarseP1) {
              elc->SetChildElement(iel, j, jel + j);
            }
          }

          // project vertex indices -----------------
          for(unsigned j = 0; j < _mesh.GetRefIndex(); j++)
            for(unsigned inode = 0; inode < elc->GetNVE(elt, CONTINUOUS_LINEAR); inode++) {
              unsigned jDof =  otherFiniteElement[elt][ CONTINUOUS_LINEAR ]->GetBasis()->GetFine2CoarseVertexMapping(j, inode);
              _mesh.el->SetElementDofIndex(jel + j, inode,  elc->GetElementDofIndex(iel, jDof));
            }

          // project face indices -----------------
          for(unsigned iface = 0; iface <  elc->GetNFC(elt, 1); iface++) {
            int value = elc->GetFaceElementIndex(iel, iface);

            if(value < -1)
              for(unsigned jface = 0; jface < _mesh.GetRefFaceIndex(); jface++)
                _mesh.el->SetFaceElementIndex(jel + coarse2FineFaceMapping[elt][iface][jface][0], coarse2FineFaceMapping[elt][iface][jface][1], value);
          }

          // update element numbers -----------------
          jel += _mesh.GetRefIndex();
          _mesh.el->AddToElementNumber(_mesh.GetRefIndex(), elt);
        }
        else {
            
          _mesh.SetIfHomogeneous(false);
          AMR = true;
          
          unsigned elt = elc->GetElementType(iel);
          
          // project element type, group, material; child element  -----------------
          _mesh.el->SetElementType(jel, elc->GetElementType(iel));
          _mesh.el->SetElementGroup(jel , elc->GetElementGroup(iel));
          _mesh.el->SetElementMaterial(jel, elc->GetElementMaterial(iel));
	  
	  unsigned gr_mat = elc->GetElementMaterial(iel);
	  _mesh.el->SetElementMaterial(jel, gr_mat);
	  
          _mesh.el->SetElementLevel(jel, elc->GetElementLevel(iel));
          if(iel >= elementOffsetCoarse && iel < elementOffsetCoarseP1) {
            elc->SetChildElement(iel, 0, jel);
          }

          // project nodes indices -----------------
          for(unsigned inode = 0; inode < elc->GetNVE(elt, CONTINUOUS_BIQUADRATIC); inode++)
            _mesh.el->SetElementDofIndex(jel, inode, elc->GetElementDofIndex(iel, inode));

          // project face indices -----------------
          for(unsigned iface = 0; iface <  elc->GetNFC(elt, 1); iface++) {
            int value = elc->GetFaceElementIndex(iel, iface);

            if(value < -1) {
              _mesh.el->SetFaceElementIndex(jel, iface, value);
            }
          }

          // update element numbers -----------------
          jel++;
          _mesh.el->AddToElementNumber(1, elt);
        }
        
      }
      
      elc->FreeLocalizedElementDof();
      elc->FreeLocalizedElementNearFace();
      elc->FreeLocalizedElement_Level_Type_Group_Material();
      
    }
    
   
    std::vector< double > ().swap(coarseLocalizedAmrVector);

//====== END ELEMENTS ==============================
    
    
//====== BEGIN NODES  ==============================
    
    int nnodes = elc->GetNodeNumber();
    _mesh.SetNumberOfNodes(nnodes);
    _mesh.el->SetNodeNumber(nnodes);

    //find all the elements near each vertex
    _mesh.el->BuildElementNearVertex();

    //initialize to UINT_MAX all the middle edge points
    for(unsigned iel = 0; iel < _mesh.GetNumberOfElements(); iel++) {
      if(_mesh.el->GetIfFatherHasBeenRefined(iel)) {
        for(unsigned inode = _mesh.el->GetElementDofNumber(iel, 0); inode < _mesh.el->GetElementDofNumber(iel, 1); inode++) {
          _mesh.el->SetElementDofIndex(iel, inode, UINT_MAX);
        }
      }
    }

    //find all the middle edge points
    for(unsigned iel = 0; iel < _mesh.GetNumberOfElements(); iel++) {
      if(_mesh.el->GetIfFatherHasBeenRefined(iel)) {
        unsigned ielt = _mesh.el->GetElementType(iel);
        unsigned istart = _mesh.el->GetElementDofNumber(iel, 0);
        unsigned iend = _mesh.el->GetElementDofNumber(iel, 1);

        for(unsigned inode = istart; inode < iend; inode++) {
          if(UINT_MAX == _mesh.el->GetElementDofIndex(iel, inode)) {
            nnodes++;
            _mesh.el->SetElementDofIndex(iel, inode, nnodes  - 1u);
            unsigned im = _mesh.el->GetElementDofIndex(iel, edge2VerticesMapping[ielt][inode - istart][0]);
            unsigned ip = _mesh.el->GetElementDofIndex(iel, edge2VerticesMapping[ielt][inode - istart][1]);

            //find all the near elements which share the same middle edge point
            for(unsigned j = 0; j < _mesh.el->GetElementNearVertexNumber(im); j++) {
              unsigned jel = _mesh.el->GetElementNearVertex(im, j);

              if(_mesh.el->GetIfFatherHasBeenRefined(jel) && jel > iel) {     // to skip coarse elements
                unsigned jm = 0, jp = 0;
                unsigned jelt = _mesh.el->GetElementType(jel);

                for(unsigned jnode = 0; jnode < _mesh.el->GetElementDofNumber(jel, 0); jnode++) {
                  if(_mesh.el->GetElementDofIndex(jel, jnode) == im) {
                    jm = jnode + 1u;
                    break;
                  }
                }

                if(jm != 0) {   //TODO this can be changed and put inside (by Sara)
                  for(unsigned jnode = 0; jnode < _mesh.el->GetElementDofNumber(jel, 0); jnode++) {
                    if(_mesh.el->GetElementDofIndex(jel, jnode) == ip) {
                      jp = jnode + 1u;
                      break;
                    }
                  }

                  if(jp != 0) {
                    if(jp < jm) {
                      unsigned tp = jp;
                      jp = jm;
                      jm = tp;
                    }

                    _mesh.el->SetElementDofIndex(jel, vertices2EdgeMapping[jelt][--jm][--jp], nnodes  - 1u);
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

    
    AddFaceDofAndElementDof();
//====== END NODES  ==============================



//==== Equivalent of: ReadCoarseMeshBeforePartitioning - END ======== 




    
//==== Partition: BEGIN ======== 

    std::vector < unsigned > partition = _mesh.PartitionForElements_refinement(AMR, mshc);
    
    _mesh.FillISvectorDofMapAllFEFamilies(partition);
    
    std::vector<unsigned> ().swap(partition);

//==== Partition: END ======== 
    
    
    
//====================================
//==== SetChildElementDof ======== 
//====================================
    elc->SetChildElementDof(_mesh.el); ///@todo you can only do this after the reordering that took place in FillISVector

    
    
//==== BuildElem_NearFace_NearElem_using_NearVertex - BEGIN ======== 
    _mesh.GetMeshElements()->BuildElem_NearFace_NearElem_using_NearVertex();
    
    _mesh.GetMeshElements()->ScatterElement_Level_Type_Group_Material___NearFace();
    
    _mesh.GetMeshElements()->ScatterElementDof();
//==== BuildElem_NearFace_NearElem_using_NearVertex - END ======== 
  
    
//====================================
    
//====  Topology, Coordinates - BEGIN ======== 
    _mesh.Topology_InitializeCoordinates();

    // build Mesh coordinates by projecting the coarse coordinates
    unsigned solType = CONTINUOUS_BIQUADRATIC;

    _mesh.GetTopology()->_Sol[0]->matrix_mult(*mshc->GetTopology()->_Sol[0], *_mesh.get_prol_matrices().GetCoarseToFineProjection(solType, _mesh));
    _mesh.GetTopology()->_Sol[1]->matrix_mult(*mshc->GetTopology()->_Sol[1], *_mesh.get_prol_matrices().GetCoarseToFineProjection(solType, _mesh));
    _mesh.GetTopology()->_Sol[2]->matrix_mult(*mshc->GetTopology()->_Sol[2], *_mesh.get_prol_matrices().GetCoarseToFineProjection(solType, _mesh));
    _mesh.GetTopology()->_Sol[0]->close();
    _mesh.GetTopology()->_Sol[1]->close();
    _mesh.GetTopology()->_Sol[2]->close();
//====  Topology, Coordinates - END ======== 

    

//====  Topology, AMR - BEGIN ======== 
    _mesh.Topology_InitializeAMR();
//====  Topology, AMR - END ======== 

    
//====  Topology, Solidnodeflag - BEGIN ======== 
    _mesh.Topology_InitializeSolidNodeFlag();
    _mesh.Topology_FillSolidNodeFlag();
//====  Topology, Solidnodeflag - END ======== 

    
    
//==== AMR (only uses Topology in one point, does not modify it) - BEGIN ========
    _mesh.InitializeAndPossiblyFillAmrRestriction(AMR); 
//==== AMR (only uses Topology in one point, does not modify it) - END ======== 
    
    
//====================================
//====  CharacteristicLength ======== 
//====================================
    _mesh.SetCharacteristicLength( mshc->GetCharacteristicLength() );


//====  Print Info ======== 
    _mesh.PrintInfo();


  }


  /**
   * This function generates face dof (for hex and wedge elements) and element dof (for hex and quad)
   **/
  void MeshRefinement::AddFaceDofAndElementDof() {

    unsigned int nnodes = _mesh.GetNumberOfNodes();

    //intialize to UINT_MAX
    for(unsigned iel = 0; iel < _mesh.el->GetElementNumber(); iel++) {
      if(_mesh.el->GetIfFatherHasBeenRefined(iel)) {
        for(unsigned inode = _mesh.el->GetElementDofNumber(iel, 1); inode < _mesh.el->GetElementDofNumber(iel, 2); inode++) {
          _mesh.el->SetElementDofIndex(iel, inode, UINT_MAX);
        }
      }
    }

    // generate face dofs on quad faces for hex and wedge elements
    for(unsigned iel = 0; iel < _mesh.el->GetElementNumber(); iel++) {
      if(_mesh.el->GetIfFatherHasBeenRefined(iel)) {
        for(unsigned iface = 0; iface < _mesh.el->GetElementFaceNumber(iel, 0); iface++) {     // I think is on all the faces that are quads
          unsigned inode = _mesh.el->GetElementDofNumber(iel, 1) + iface;

          if(UINT_MAX == _mesh.el->GetElementDofIndex(iel, inode)) {
            _mesh.el->SetElementDofIndex(iel, inode, ++nnodes  - 1u);
            unsigned i1 = _mesh.el->GetFaceVertexIndex(iel, iface, 0);
            unsigned i2 = _mesh.el->GetFaceVertexIndex(iel, iface, 1);
            unsigned i3 = _mesh.el->GetFaceVertexIndex(iel, iface, 2);

            for(unsigned j = 0; j < _mesh.el->GetElementNearVertexNumber(i1); j++) {
              unsigned jel = _mesh.el->GetElementNearVertex(i1, j);

              if(_mesh.el->GetIfFatherHasBeenRefined(jel) && jel > iel) {
                for(unsigned jface = 0; jface < _mesh.el->GetElementFaceNumber(jel, 0); jface++) {
                  unsigned jnode = _mesh.el->GetElementDofNumber(jel, 1) + jface;

                  if(UINT_MAX == _mesh.el->GetElementDofIndex(jel, jnode)) {
                    unsigned j1 = _mesh.el->GetFaceVertexIndex(jel, jface, 0);
                    unsigned j2 = _mesh.el->GetFaceVertexIndex(jel, jface, 1);
                    unsigned j3 = _mesh.el->GetFaceVertexIndex(jel, jface, 2);
                    unsigned j4 = _mesh.el->GetFaceVertexIndex(jel, jface, 3);

                    if((i1 == j1 || i1 == j2 || i1 == j3 ||  i1 == j4) &&
                        (i2 == j1 || i2 == j2 || i2 == j3 ||  i2 == j4) &&
                        (i3 == j1 || i3 == j2 || i3 == j3 ||  i3 == j4)) {
                      _mesh.el->SetElementDofIndex(jel, jnode, nnodes  - 1u);
                    }
                  }
                }
              }
            }
          }
        }
        // generate face dofs on tri faces for tet and wedge elements
        unsigned elementType = _mesh.el->GetElementType(iel);
        if(elementType == 1 || elementType == 2) {
          for(unsigned iface = _mesh.el->GetElementFaceNumber(iel, 0); iface < _mesh.el->GetElementFaceNumber(iel, 1); iface++) {       //on all the faces that are triangles
            unsigned inode = _mesh.el->GetElementDofNumber(iel, 1) + iface;
            if(UINT_MAX == _mesh.el->GetElementDofIndex(iel, inode)) {
              _mesh.el->SetElementDofIndex(iel, inode, ++nnodes  - 1u);
              unsigned i1 = _mesh.el->GetFaceVertexIndex(iel, iface, 0);
              unsigned i2 = _mesh.el->GetFaceVertexIndex(iel, iface, 1);
              unsigned i3 = _mesh.el->GetFaceVertexIndex(iel, iface, 2);
              for(unsigned j = 0; j < _mesh.el->GetElementNearVertexNumber(i1); j++) {
                unsigned jel = _mesh.el->GetElementNearVertex(i1, j);
                if(_mesh.el->GetIfFatherHasBeenRefined(jel) && jel > iel) {
                  for(unsigned jface = _mesh.el->GetElementFaceNumber(jel, 0); jface < _mesh.el->GetElementFaceNumber(jel, 1); jface++) {
                    unsigned jnode = _mesh.el->GetElementDofNumber(jel, 1) + jface;
                    if(UINT_MAX == _mesh.el->GetElementDofIndex(jel, jnode)) {
                      unsigned j1 = _mesh.el->GetFaceVertexIndex(jel, jface, 0);
                      unsigned j2 = _mesh.el->GetFaceVertexIndex(jel, jface, 1);
                      unsigned j3 = _mesh.el->GetFaceVertexIndex(jel, jface, 2);
                      if((i1 == j1 || i1 == j2 || i1 == j3) &&
                          (i2 == j1 || i2 == j2 || i2 == j3) &&
                          (i3 == j1 || i3 == j2 || i3 == j3)) {
                        _mesh.el->SetElementDofIndex(jel, jnode, nnodes  - 1u);
                      }
                    }
                  }
                }
              }
            }
          }
        }

      }
    }

    // generates element dofs for hex, tet, wedge, quad and triangle elements
    for(unsigned iel = 0; iel < _mesh.el->GetElementNumber(); iel++) {
      if(_mesh.el->GetIfFatherHasBeenRefined(iel)) {
        if(0 == _mesh.el->GetElementType(iel)) {     //hex
          _mesh.el->SetElementDofIndex(iel, 26, ++nnodes - 1u);
        }
        if(1 == _mesh.el->GetElementType(iel)) {     //tet
          _mesh.el->SetElementDofIndex(iel, 14, ++nnodes - 1u);
        }
        if(2 == _mesh.el->GetElementType(iel)) {     //wedge
          _mesh.el->SetElementDofIndex(iel, 20, ++nnodes - 1u);
        }
        if(3 == _mesh.el->GetElementType(iel)) {     //quad
          _mesh.el->SetElementDofIndex(iel, 8, ++nnodes - 1u);
        }
        if(4 == _mesh.el->GetElementType(iel)) {     //triangle
          _mesh.el->SetElementDofIndex(iel, 6, ++nnodes - 1u);
        }
      }
    }

    _mesh.el->SetNodeNumber(nnodes);
    _mesh.SetNumberOfNodes(nnodes);

  }


}
