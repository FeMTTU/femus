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

    /** Flag all the even elements to be refined */
    void FlagOnlyEvenElementsToBeRefined();

    /** Flag the elements to be refined in according to AMR criteria */
    void FlagElementsToBeRefined();


private:

    void FlagElementsToRefine(const unsigned& type);

    /** To be added */
    void Buildkmid();


    Mesh& _mesh;                 //< reference to the mesh which is built by refinement

};

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


}   //end namespace femus



#endif