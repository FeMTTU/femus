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
#include "GeomElTypeEnum.hpp"

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
    return FlagElementsToRefineBaseOnError(treshold, error);
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
          _mesh._topology->_Sol[_mesh.GetAmrIndex()]->set(iel, 1.);
          numberOfRefinedElement->add(_iproc, 1.);
        }
      }
    }
    else if(type == 1) {   // Flag AMR elements
      for(int iel = _mesh._elementOffset[_iproc]; iel < _mesh._elementOffset[_iproc + 1]; iel++) {
        if(_mesh.el->GetIfElementCanBeRefined(iel)) {
          if((*_mesh._topology->_Sol[ _mesh.GetAmrIndex() ])(iel) > 0.5) {
            numberOfRefinedElement->add(_iproc, 1.);
          }
          else if(_mesh._IsUserRefinementFunctionDefined) {
            short unsigned ielt = _mesh.GetElementType(iel);
            unsigned nve = _mesh.GetElementDofNumber(iel, 0);
            std::vector < double > x(3, 0.);

            for(unsigned i = 0; i < nve; i++) {
              unsigned inode_metis = _mesh.GetSolutionDof(i, iel, 2);
              x[0] += (*_mesh._topology->_Sol[0])(inode_metis);
              x[1] += (*_mesh._topology->_Sol[1])(inode_metis);
              x[2] += (*_mesh._topology->_Sol[2])(inode_metis);
            }

            x[0] /= nve;
            x[1] /= nve;
            x[2] /= nve;

            if(_mesh._SetRefinementFlag(x, _mesh.GetElementGroup(iel), _mesh.GetLevel())) {
              _mesh._topology->_Sol[ _mesh.GetAmrIndex() ]->set(iel, 1.);
              numberOfRefinedElement->add(_iproc, 1.);
            }
          }
        }
        else {
          _mesh._topology->_Sol[ _mesh.GetAmrIndex() ]->set(iel, 0.);
        }
      }
    }
    else if(type == 2) {   // Flag only even elements (for debugging purposes)
      for(int iel = _mesh._elementOffset[_iproc]; iel < _mesh._elementOffset[_iproc + 1]; iel++) {
        if(_mesh.el->GetIfElementCanBeRefined(iel)) {
          if((*_mesh._topology->_Sol[_mesh.GetAmrIndex()])(iel) < 0.5 && iel % 2 == 0) {
            _mesh._topology->_Sol[_mesh.GetAmrIndex()]->set(iel, 1.);
            numberOfRefinedElement->add(_iproc, 1.);
          }
        }
      }
    }

    _mesh._topology->_Sol[_mesh.GetAmrIndex()]->close();
    //END flag element to be refined

    //BEGIN update elem
    numberOfRefinedElement->close();
    double totalNumber = numberOfRefinedElement->l1_norm();
    _mesh.el->SetRefinedElementNumber(static_cast < unsigned >(totalNumber + 0.25));
    delete numberOfRefinedElement;

    //END update elem
  }


  bool MeshRefinement::FlagElementsToRefineBaseOnError(const double& treshold, NumericVector& error) {

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
        if((*_mesh._topology->_Sol[_mesh.GetAmrIndex()])(iel) < 0.5 && error(iel) > treshold) {
          _mesh._topology->_Sol[_mesh.GetAmrIndex()]->set(iel, 1.);
          numberOfRefinedElement->add(_iproc, 1.);
        }
      }
    }

    _mesh._topology->_Sol[_mesh.GetAmrIndex()]->close();
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
  void MeshRefinement::RefineMesh(const unsigned& igrid, Mesh* mshc, const elem_type* otherFiniteElement[6][5]) {

    _mesh.SetIfHomogeneous(true);

    _mesh.SetCoarseMesh(mshc);

    _mesh.SetFiniteElementPtr(otherFiniteElement);

    elem* elc = mshc->el;

    _mesh.SetLevel(igrid);

    // total number of elements on the fine level
    int nelem = elc->GetRefinedElementNumber() * _mesh.GetRefIndex(); // refined
    nelem += elc->GetElementNumber() - elc->GetRefinedElementNumber(); // not-refined

    unsigned elementOffsetCoarse   = mshc->_elementOffset[_iproc];
    unsigned elementOffsetCoarseP1 = mshc->_elementOffset[_iproc + 1];

    _mesh.SetNumberOfElements(nelem);

    vector < double > coarseLocalizedAmrVector;
    mshc->_topology->_Sol[mshc->GetAmrIndex()]->localize_to_all(coarseLocalizedAmrVector);

    mshc->el->AllocateChildrenElement(_mesh.GetRefIndex(), mshc);

    _mesh.el = new elem(elc, _mesh.GetRefIndex(), coarseLocalizedAmrVector);

    _mesh.SetCharacteristicLength( mshc->GetCharacteristicLength() );
    
    unsigned jel = 0;
    //divide each coarse element in 8(3D), 4(2D) or 2(1D) fine elements and find all the vertices

    _mesh.el->SetElementGroupNumber(elc->GetElementGroupNumber());

    bool AMR = false;

    std::vector < unsigned > materialElementCounter(3,0);
    
    for(unsigned isdom = 0; isdom < _nprocs; isdom++) {
      elc->LocalizeElementDof(isdom);
      elc->LocalizeElementNearFace(isdom);
      elc->LocalizeElementQuantities(isdom);
      for(unsigned iel = mshc->_elementOffset[isdom]; iel < mshc->_elementOffset[isdom + 1]; iel++) {
        if(static_cast < unsigned short >(coarseLocalizedAmrVector[iel] + 0.25) == 1) {
          unsigned elt = elc->GetElementType(iel);
          // project element type
          for(unsigned j = 0; j < _mesh.GetRefIndex(); j++) {
            _mesh.el->SetElementType(jel + j, elc->GetElementType(iel));
            _mesh.el->SetElementGroup(jel + j, elc->GetElementGroup(iel));
            
	    unsigned gr_mat = elc->GetElementMaterial(iel);
	    _mesh.el->SetElementMaterial(jel + j, gr_mat);
	    if( gr_mat == 2) materialElementCounter[0] += 1;
	    else if(gr_mat == 3 ) materialElementCounter[1] += 1;
	    else materialElementCounter[2] += 1;
	        
            _mesh.el->SetElementLevel(jel + j, elc->GetElementLevel(iel) + 1);
            if(iel >= elementOffsetCoarse && iel < elementOffsetCoarseP1) {
              elc->SetChildElement(iel, j, jel + j);
            }
          }

          // project vertex indeces
          for(unsigned j = 0; j < _mesh.GetRefIndex(); j++)
            for(unsigned inode = 0; inode < elc->GetNVE(elt, 0); inode++) {
              unsigned jDof =  otherFiniteElement[elt][0]->GetBasis()->GetFine2CoarseVertexMapping(j, inode);
              _mesh.el->SetElementDofIndex(jel + j, inode,  elc->GetElementDofIndex(iel, jDof));
            }

          // project face indeces
          for(unsigned iface = 0; iface <  elc->GetNFC(elt, 1); iface++) {
            int value = elc->GetFaceElementIndex(iel, iface);

            if(value < -1)
              for(unsigned jface = 0; jface < _mesh.GetFaceIndex(); jface++)
                _mesh.el->SetFaceElementIndex(jel + coarse2FineFaceMapping[elt][iface][jface][0], coarse2FineFaceMapping[elt][iface][jface][1], value);
          }

          // update element numbers
          jel += _mesh.GetRefIndex();
          _mesh.el->AddToElementNumber(_mesh.GetRefIndex(), elt);
        }
        else {
          _mesh.SetIfHomogeneous(false);
          AMR = true;
          unsigned elt = elc->GetElementType(iel);
          _mesh.el->SetElementType(jel, elc->GetElementType(iel));
          _mesh.el->SetElementGroup(jel , elc->GetElementGroup(iel));
          _mesh.el->SetElementMaterial(jel, elc->GetElementMaterial(iel));
	  
	  unsigned gr_mat = elc->GetElementMaterial(iel);
	  _mesh.el->SetElementMaterial(jel, gr_mat);
	  if( gr_mat == 2) materialElementCounter[0] += 1;
	  else if(gr_mat == 3 ) materialElementCounter[1] += 1;
	  else materialElementCounter[2] += 1;
	  
          _mesh.el->SetElementLevel(jel, elc->GetElementLevel(iel));
          if(iel >= elementOffsetCoarse && iel < elementOffsetCoarseP1) {
            elc->SetChildElement(iel, 0, jel);
          }

          // project nodes indeces
          for(unsigned inode = 0; inode < elc->GetNVE(elt, 2); inode++)
            _mesh.el->SetElementDofIndex(jel, inode, elc->GetElementDofIndex(iel, inode));

          // project face indeces
          for(unsigned iface = 0; iface <  elc->GetNFC(elt, 1); iface++) {
            int value = elc->GetFaceElementIndex(iel, iface);

            if(value < -1) {
              _mesh.el->SetFaceElementIndex(jel, iface, value);
            }
          }

          // update element numbers
          jel++;
          _mesh.el->AddToElementNumber(1, elt);
        }
      }
      elc->FreeLocalizedElementDof();
      elc->FreeLocalizedElementNearFace();
      elc->FreeLocalizedElementQuantities();
    }
    
    _mesh.el->SetMaterialElementCounter(materialElementCounter);
    
    
    std::vector<unsigned> MaterialElementCounter = _mesh.el->GetMaterialElementCounter();
    //std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAA\n";
    //std::cout << MaterialElementCounter[0]<<" "<< MaterialElementCounter[1]<<" "<< MaterialElementCounter[2]<<" \n";
   
    coarseLocalizedAmrVector.resize(0);
    //coarseLocalizedElementType.resize(0);

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

    Buildkmid();

    std::vector < unsigned > partition;
    partition.reserve(_mesh.GetNumberOfNodes());
    partition.resize(_mesh.GetNumberOfElements());

    MeshMetisPartitioning meshMetisPartitioning(_mesh);

    if(AMR == true) {
      meshMetisPartitioning.DoPartition(partition, AMR);
    }
    else {
      meshMetisPartitioning.DoPartition(partition, *mshc);
    }

    _mesh.FillISvector(partition);
    partition.resize(0);

    elc->SetChildElementDof(_mesh.el);

    _mesh.el->DeleteElementNearVertex();
    _mesh.el->BuildElementNearVertex();

    _mesh.Buildkel();

    // build Mesh coordinates by projecting the coarse coordinats
    _mesh._topology = new Solution(&_mesh);
    _mesh._topology->AddSolution("X", LAGRANGE, SECOND, 1, 0);
    _mesh._topology->AddSolution("Y", LAGRANGE, SECOND, 1, 0);
    _mesh._topology->AddSolution("Z", LAGRANGE, SECOND, 1, 0);

    _mesh._topology->ResizeSolutionVector("X");
    _mesh._topology->ResizeSolutionVector("Y");
    _mesh._topology->ResizeSolutionVector("Z");

    _mesh._topology->AddSolution("AMR", DISCONTINUOUS_POLYNOMIAL, ZERO, 1, 0);
    _mesh._topology->ResizeSolutionVector("AMR");

    unsigned solType = 2;

    _mesh._topology->_Sol[0]->matrix_mult(*mshc->_topology->_Sol[0], *_mesh.GetCoarseToFineProjection(solType));
    _mesh._topology->_Sol[1]->matrix_mult(*mshc->_topology->_Sol[1], *_mesh.GetCoarseToFineProjection(solType));
    _mesh._topology->_Sol[2]->matrix_mult(*mshc->_topology->_Sol[2], *_mesh.GetCoarseToFineProjection(solType));
    _mesh._topology->_Sol[0]->close();
    _mesh._topology->_Sol[1]->close();
    _mesh._topology->_Sol[2]->close();

    _mesh.el->BuildElementNearElement();
    _mesh.el->DeleteElementNearVertex();

    _mesh._topology->AddSolution("solidMrk", LAGRANGE, SECOND, 1, 0);
    _mesh.AllocateAndMarkStructureNode();
    
    _mesh.el->ScatterElementQuantities();
    _mesh.el->ScatterElementDof();
    _mesh.el->ScatterElementNearFace();

    std::vector < std::map < unsigned,  std::map < unsigned, double  > > >& restriction = _mesh.GetAmrRestrictionMap();
    if(AMR) {
      _mesh.el->GetAMRRestriction(&_mesh);
//       for(unsigned soltype = 0; soltype < 3; soltype++) {
//         std::cout << "solution type = " << soltype << std::endl;
//         for(std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restriction[soltype].begin(); it1 != restriction[soltype].end(); it1++) {
//           std::cout << it1->first << "\t";
//           for(std::map<unsigned, double> ::iterator it2 = restriction[soltype][it1->first].begin(); it2 != restriction[soltype][it1->first].end(); it2++) {
//             std::cout << it2->first << " (" << it2->second << ")  ";
// 
//           }
//           std::cout << std::endl;
//         }
//       }
    }
    else{
      restriction.resize(3);
    }
    
    _mesh.PrintInfo();
  }


  /**
   * This function generates face dof (for hex and wedge elements) and element dof (for hex and quad)
   **/
  void MeshRefinement::Buildkmid() {

    unsigned int nnodes = _mesh.GetNumberOfNodes();
    //std::cout << "nnodes before buildkmid= "  << nnodes << std::endl;

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
    //std::cout << "nnodes after buildkmid= "  << nnodes << std::endl;

  }


}
