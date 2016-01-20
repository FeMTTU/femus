/*=========================================================================

 Program: FEMuS
 Module: Elem
 Authors: Eugenio Aulisa

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_mesh_Elem_hpp__
#define __femus_mesh_Elem_hpp__

#include <vector>
#include <map>

#include "Mesh.hpp"

namespace femus {

  class Mesh;
  /**
   * The elem class
  */
  class elem {


    public:

      /** constructors */
      elem(const unsigned& other_nel);

      elem(elem* elc, const unsigned refindex, const std::vector < double >& coarseAmrLocal, const std::vector < double >& localizedElementType);

      void SharpMemoryAllocation();

      /** destructor */
      ~elem();

      void DeleteGroupAndMaterial();

      void DeleteElementType();

      void ScatterElementCanBeRefinedVector();

      void ScatterElementNearFace();
      void LocalizeElementNearFaceFromOneToAll( const unsigned& jproc );
      void FreeLocalizedElementNearFace();

      void ScatterElementDof();
      void LocalizeElementDofFromOneToAll( const unsigned &jproc);
      void LocalizeElementDofFromOneToOne( const unsigned &jproc, const unsigned &kproc );
      void FreeLocalizedElementDof();

      // reorder the element according to the new element mapping
      void ReorderMeshElements(const std::vector < unsigned >& elementMapping);

      // reorder the nodes according to the new node mapping
      void ReorderMeshNodes(const std::vector < unsigned >& nodeMapping);

      /** To be Added */
      unsigned GetElementDofNumber(const unsigned& iel, const unsigned& type) const;

      /** Return the local->global node number */
      unsigned GetElementDofIndex(const unsigned& iel, const unsigned& inode)const;

      /** To be Added */
      void SetElementDofIndex(const unsigned& iel, const unsigned& inode, const unsigned& value);

      /** To be Added */
      unsigned GetFaceVertexIndex(const unsigned& iel, const unsigned& iface, const unsigned& inode) const;

      /** To be Added */
      short unsigned GetElementType(const unsigned& iel) const;

      /** To be Added */
      void SetElementType(const unsigned& iel, const short unsigned& value);

      /** To be Added */
      short unsigned GetElementGroup(const unsigned& iel) const;

      /** To be Added */
      void SetElementGroup(const unsigned& iel, const short unsigned& value);

      /** To be Added */
      void SetElementMaterial(const unsigned& iel, const short unsigned& value);

      /** To be Added */
      short unsigned GetElementMaterial(const unsigned& iel) const;

      /** To be Added */
      unsigned GetElementGroupNumber() const;

      /** To be Added */
      void SetElementGroupNumber(const unsigned& value);

      /** To be Added */
      int GetFaceElementIndex(const unsigned& iel, const unsigned& iface) const;

      int GetBoundaryIndex(const unsigned& iel, const unsigned& iface) const;

      /** To be Added */
      void SetFaceElementIndex(const unsigned& iel, const unsigned& iface, const int& value);

      /** To be Added */
      unsigned GetIndex(const char name[]) const;

      /** To be Added */
      unsigned GetElementNumber(const char* name = "All") const;

      /** To be Added */
      void AddToElementNumber(const unsigned& value, const char name[]);

      /** To be Added */
      void AddToElementNumber(const unsigned& value, short unsigned ielt);

      /** To be Added */
      unsigned GetRefinedElementNumber() const {
        return _nelr;
      };

      /** To be Added */
      unsigned GetRefinedElementTypeNumber(const unsigned& ielt) const {
        return _nelrt[ielt];
      };

      /** To be Added */
      void SetRefinedElementNumber(const unsigned& value) {
        _nelr = value;
      };

      /** To be Added */
      void SetRefinedElemenTypeNumber(const unsigned& value, const unsigned& ielt) {
        _nelrt[ielt] = value;
      };

      /** To be Added */
      unsigned GetNodeNumber()const;

      /** To be Added */
      void SetNodeNumber(const unsigned& value);

      /** To be Added */
      unsigned GetElementFaceNumber(const unsigned& iel, const unsigned& type = 1)const;

      /** To be Added */
      void BuildElementNearVertex();

      /** To be Added */
      void SetChildElementDof(const unsigned& ref_index, Mesh* msh, const elem* elf);

      unsigned GetChildElementDof(const unsigned& iel, const unsigned& i0, const unsigned i1) const {
        return _childElemDof[iel - _elementOffset[_iproc] ][i0][i1];
      }

      void DeleteElementNearVertex();

      /** To be Added */
      unsigned GetElementNearVertexNumber(const unsigned& inode)const;

      /** To be Added */
      unsigned GetElementNearVertex(const unsigned& inode, const unsigned& jnode)const;

      void BuildLocalElementNearVertex();

      const std::vector<unsigned>& GetLocalElementNearVertex(const unsigned& inode)  {
        return _localElementNearVertexMap[inode];
      };

      /** To be Added */
      void SetIfElementCanBeRefined(const unsigned& iel, const bool& refined);

      /** To be Added */
      bool GetIfElementCanBeRefined(const unsigned& iel) const;

      /** To be Added */
      bool GetIfFatherHasBeenRefined(const unsigned& iel) {
        return GetIfElementCanBeRefined(iel);
      }

      /** To be Added */
      void AllocateChildrenElement(const unsigned& ref_index, Mesh* msh);

      /** To be Added */
      void SetChildElement(const unsigned& iel, const unsigned& json, const unsigned& value);

      /** To be Added */
      unsigned GetChildElement(const unsigned& iel, const unsigned& json) const;

      const unsigned GetNVE(const unsigned& elementType, const unsigned& doftype) const;

      const unsigned GetNFACENODES(const unsigned& elementType, const unsigned& jface, const unsigned& dof) const;

      const unsigned GetNFC(const unsigned& elementType, const unsigned& type) const;

      const unsigned GetIG(const unsigned& elementType, const unsigned& iface, const unsigned& jnode) const;

      void SetElementOffsets( const std::vector < unsigned > & elementOffset, const unsigned &iproc, const unsigned &nprocs) {
        _elementOffset = elementOffset;
        _elementOwned = elementOffset[iproc + 1] - elementOffset[iproc];
        _iproc = iproc;
        _nprocs = nprocs;
      }

    private:
      unsigned _iproc;
      unsigned _nprocs;


      // member data
      int** _elementNearFace;
      int* _elementNearFaceMemory;
      unsigned _elementNearFaceMemorySize;

      unsigned _elementNearFaceOffset;
      bool _elementNearFaceIsScattered;

      bool _elementNearFaceIsLocalizedFromJproc;
      unsigned _jprocElementNearFaceIsLocalizedFrom;

      int** _localElementNearFace; //element -> nodes
      int* _localElementNearFaceMemory;
      unsigned _localElementNearFaceMemorySize;


      std::map< unsigned, std::vector< unsigned > > _localElementNearVertexMap;
      unsigned** _elementNearVertex; //node->element
      unsigned* _elementNearVertexMemory;
      unsigned* _elementNearVertexNumber;

      unsigned** _elementDof; //element -> nodes
      unsigned* _elementDofMemory;
      unsigned _elementDofMemorySize;

      unsigned _elementDofOffset;
      bool _elementDofIsScattered;

      bool _elementDofIsLocalizedFromJproc;
      unsigned _jprocElementDofIsLocalizedFrom;

      unsigned** _localElementDof; //element -> nodes
      unsigned* _localElementDofMemory;
      unsigned _localElementDofMemorySize;

      unsigned** _childElem;
      unsigned* _childElemMemory;
      unsigned _childElemMemorySize;
      bool _childElemFlag;

      unsigned*** _childElemDof; //element -> nodes
      unsigned** _childElemDofMemoryPointer; //element -> nodes
      unsigned* _childElemDofMemory;
      unsigned _childElemDofMemorySize;

      std::vector < unsigned > _elementOffset;
      unsigned _elementOwned;

      short unsigned* _elementType, *_elementGroup, *_elementMaterial; //element

      unsigned _nvt;
      unsigned _nel, _nelt[6];
      unsigned _nelr, _nelrt[6];
      unsigned _ngroup;

      bool* _elementCanBeRefined; //element
      bool _elementCanBeRefinedIsScattered;

      unsigned _level;

      elem* _coarseElem;

  };

//linear, quadratic, biquadratic, picewise costant, picewise linear discontinuous
  const unsigned NVE[6][5] = {
    {8, 20, 27, 1, 4}, //hex
    {4, 10, 10, 1, 4}, //tet
    {6, 15, 18, 1, 4}, //wedge
    {4, 8, 9, 1, 3}, //quad
    {3, 6, 6, 1, 3}, //tri
    {2, 3, 3, 1, 2}  //line
  };

//number of dof objects, or "dof carriers" for every geometric element and every FE family
  const unsigned NDOFOBJS[6][5] = {
    {8, 20, 27, 1, 1}, //hex
    {4, 10, 10, 1, 1}, //tet
    {6, 15, 18, 1, 1}, //wedge
    {4, 8, 9, 1, 1}, //quad
    {3, 6, 6, 1, 1}, //tri
    {2, 3, 3, 1, 1}  //line
  };

  /**
   * Number of elements obtained with one refinement
  **/
  const unsigned NRE[6] = {8, 8, 8, 4, 4, 2};

  /**
   * Number of FACES(3D), edges(2D) or point-extrema(1D) for each considered element
   **/
  const unsigned NFC[6][2] = {
    {6, 6},
    {0, 4},
    {3, 5},
    {0, 4},
    {0, 3},
    {0, 2}
  };

  /**
   * Node ordering for each element face(3D), edge(2D) or point-extrema(1D) position for each considered element
   **/
  const unsigned ig[6][6][9] = {
    { {0, 1, 5, 4, 8, 17, 12, 16, 20},
      {1, 2, 6, 5, 9, 18, 13, 17, 21},
      {2, 3, 7, 6, 10, 19, 14, 18, 22},
      {3, 0, 4, 7, 11, 16, 15, 19, 23},
      {0, 3, 2, 1, 11, 10, 9, 8, 24},
      {4, 5, 6, 7, 12, 13, 14, 15, 25}
    },
    { {0, 2, 1, 6, 5, 4},
      {0, 1, 3, 4, 8, 7},
      {1, 2, 3, 5, 9, 8},
      {2, 0, 3, 6, 7, 9}
    },
    { {0, 1, 4, 3, 6, 13, 9, 12, 15},
      {1, 2, 5, 4, 7, 14, 10, 13, 16},
      {2, 0, 3, 5, 8, 12, 11, 14, 17},
      {0, 2, 1, 8, 7, 6},
      {3, 4, 5, 9, 10, 11}
    },
    { {0, 1, 4},
      {1, 2, 5},
      {2, 3, 6},
      {3, 0, 7}
    },
    { {0, 1, 3},
      {1, 2, 4},
      {2, 0, 5}
    },
    { {0},
      {1}
    }
  };


  const unsigned NFACENODES[6][6][3] =
  {
    { {4, 8, 9}, // Hex
      {4, 8, 9},
      {4, 8, 9},
      {4, 8, 9},
      {4, 8, 9},
      {4, 8, 9}
    },
    { {3, 6, 6}, // Tet
      {3, 6, 6},
      {3, 6, 6},
      {3, 6, 6}
    },
    { {4, 8, 9}, // Wedge
      {4, 8, 9},
      {4, 8, 9},
      {3, 6, 6},
      {3, 6, 6}
    },
    { {2, 3, 3},
      {2, 3, 3}, // Quad
      {2, 3, 3},
      {2, 3, 3}
    },
    { {2, 3, 3}, // Tri
      {2, 3, 3},
      {2, 3, 3}
    },
    { {1, 1, 1}, // Line
      {1, 1, 1}
    }
  };


} //end namespace femus


const unsigned referenceElementDirection[6][3][2] = { //Endpoint1, Endpoint2 =rEED[elemem type][direction][0,1]
  {
    {23, 21}, {20, 22}, {24, 25}
  },
  {
    {0, 1}, {0, 2}, {0, 3}
  },
  {
    {12, 13}, {12, 14}, {0, 3}
  },
  {
    {7, 5}, {4, 6}
  },
  {
    {0, 1}, {0, 2}
  },
  {
    {0, 1}
  }
};

// const unsigned referenceElementPoint[6]={26,0,12,8,0,2};

#endif


// ********************  class solver**************************

/*

//         7------14-------6
//        /|              /|
//       / |             / |
//     15  |   25      13  |
//     /  19      22   /  18
//    /    |          /    |
//   4------12-------5     |
//   | 23  |   26    | 21  |
//   |     3------10-|-----2
//   |    /          |    /
//  16   /  20      17   /
//   | 11      24    |  9
//   | /             | /
//   |/              |/
//   0-------8-------1




//            3
//           /|\
//          / | \
//         /  |  \
//        9   |   8
//       /    |    \
//      /     |     \
//     /      7      \
//    2-------|5------1
//     \      |      /
//      \     |     /
//       \    |    /
//        6   |   4
//         \  |  /
//          \ | /
//           \|/
//            0

//           5
//          /|\
//         / | \
//        /  |  \
//      11   |  10
//      /   14    \
//     /     |     \
//    /      |      \
//   3--------9------4
//   |  17   |  16   |
//   |       2       |
//   |      / \      |
//   |     /   \     |
//  12    / 15  \   13
//   |   8       7   |
//   |  /         \  |
//   | /           \ |
//   |/             \|
//   0-------6-------1

//      3-----6-----2
//      |           |
//      |           |
//      7     8     5
//      |           |
//      |           |
//      0-----4-----1


//      2
//      | \
//      |   \
//      5     4
//      |       \
//      |         \
//      0-----3----1


//
//	0-----2-----1
//
*/
