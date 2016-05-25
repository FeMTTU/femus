/*=========================================================================

 Program: FEMuS
 Module: MeshRefinement
 Authors: Simone Bnà, Eugenio Aulisa

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_mesh_MeshRefinement_hpp__
#define __femus_mesh_MeshRefinement_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "ParallelObject.hpp"

namespace femus {


class Mesh;
class elem_type;

/**
 * This is the \p MeshRefinement class.  This class implements
 * adaptive, selective and global mesh refinement algorithms for a \p Mesh.
*/

class MeshRefinement : public ParallelObject {

public:

    /** Constructor */
    MeshRefinement(Mesh& mesh);

    /** destructor */
    ~MeshRefinement();

    /** Refinement functions */

    /** This function generates a finer mesh level, $l_i$, from a coarser mesh level $l_{i-1}$, $i>0$ */
    void RefineMesh(const unsigned &igrid, Mesh *mshc, const elem_type* otheFiniteElement[6][5]);

    /** Flag all the elements to be refined */
    void FlagAllElementsToBeRefined();
    bool FlagElementsToBeRefined(const double & treshold, NumericVector& error);

    /** Flag all the even elements to be refined */
    void FlagOnlyEvenElementsToBeRefined();

    /** Flag the elements to be refined in according to AMR criteria */
    bool FlagElementsToBeRefined();


private:

    void FlagElementsToRefine(const unsigned& type);
    bool FlagElementsToRefineBaseOnError(const double& treshold, NumericVector& error);

    /** To be added */
    void Buildkmid();


    Mesh& _mesh;                 //< reference to the mesh which is built by refinement

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
      { {0,0},{1,0},{2,0},{4,0} }, 
      { {0,1},{1,1},{3,1},{5,1} }, 
      { {1,2},{2,2},{3,2},{6,2} }, 
      { {2,3},{0,3},{3,3},{7,3} }  
    },
    {
      { {0,0},{1,0},{4,0},{5,0} },
      { {1,1},{2,1},{5,1},{6,1} },
      { {2,2},{0,2},{6,2},{4,2} },
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


}   //end namespace femus



#endif