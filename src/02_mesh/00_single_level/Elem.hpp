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


#include "GeomElTypeEnum.hpp"
#include "FElemTypeEnum_list.hpp"
#include "MyVector.hpp"
#include "MyMatrix.hpp"
#include "PolynomialBases.hpp"


#include <vector>
#include <map>


namespace femus {

  //Forward declarations  
  class Mesh;
  /**
   * The elem class: it contains the list of all Mesh Geometric Elements, along with several Element-based and also Node-based properties
   * @todo I believe it would even be more linear if this class did not have any function at all involving the Mesh pointer, there are very few in any case
  */
  class elem {


// === CONSTR-DESTR - BEGIN =================
    public:

      /** constructors */
      elem(const unsigned& other_nel, const unsigned dim_in);

      //elem(elem* elc, const unsigned refindex, const std::vector < double >& coarseAmrLocal, const std::vector < double >& localizedElementType);
      elem(elem* elc, const unsigned dim_in, const unsigned refindex, const std::vector < double >& coarseAmrLocal);

      /** destructor */
      ~elem();
// === CONSTR-DESTR - END =================
      
      
    public:


      // reorder the element according to the new element mapping
      void ReorderMeshElements(const std::vector < unsigned >& elementMapping);

      // reorder the nodes according to the new node mapping
      void ReorderMeshNodes(const std::vector < unsigned >& nodeMapping);
      
      

      
      


      void BuildMeshElemStructures();
      

      //BEGIN _ElementLevel functions
      void ResizeElementQuantities(const unsigned nel_in, const unsigned level_in) {
       _elementLevel.resize(nel_in, level_in);
       _elementType.resize(nel_in);
       _elementGroup.resize(nel_in);
       _elementMaterial.resize(nel_in);
      }
      
      void ScatterElementQuantities() {
        _elementLevel.scatter(_elementOffset);
        _elementType.scatter(_elementOffset);
        _elementGroup.scatter(_elementOffset);
        _elementMaterial.scatter(_elementOffset);
      }
      
      void LocalizeElementQuantities(const unsigned &lproc) {
        _elementLevel.broadcast(lproc);
        _elementType.broadcast(lproc);
        _elementGroup.broadcast(lproc);
        _elementMaterial.broadcast(lproc);
      }
      
      void FreeLocalizedElementQuantities() {
        _elementLevel.clearBroadcast();
        _elementType.clearBroadcast();
        _elementGroup.clearBroadcast();
        _elementMaterial.clearBroadcast();
      }

      bool GetIfElementCanBeRefined(const unsigned& iel) {
        return (_elementLevel[iel] == _level) ? true : false;
      }
      bool GetIfFatherHasBeenRefined(const unsigned& iel) {
        return GetIfElementCanBeRefined(iel);
      }
      //END _ElementLevel functions


      /** To be Added */
      void SetChildElementDof(elem* elf);

      unsigned GetChildElementDof(const unsigned& iel, const unsigned& i0, const unsigned i1);
      
      /** To be Added */
      void AllocateChildrenElement(const unsigned int& refindex, const Mesh* msh);

      /** To be Added */
      void SetChildElement(const unsigned& iel, const unsigned& json, const unsigned& value);

      /** To be Added */
      unsigned GetChildElement(const unsigned& iel, const unsigned& json);
      

// === Geometric Element, Single - BEGIN =================
  public:
    
      /** To be Added */
      const unsigned GetElementFaceNumber(const unsigned& iel, const unsigned& type = 1) const;
      
      const unsigned GetNFC(const unsigned& elementType, const unsigned& type) const;

      const unsigned GetIG(const unsigned& elementType, const unsigned& iface, const unsigned& jnode) const;
      
      const unsigned GetNRE(const unsigned& elementType) const { return NRE[elementType]; }
      
      const unsigned GetReferenceElementDirection(const unsigned& elementType, const unsigned dir, const unsigned node) const { 
        return referenceElementDirection[elementType][dir][node]; }
      
  private:

  /**
   * Number of FACES(3D), edges(2D) or point-extrema(1D) for each considered element
   * The 1st number is the quadrilaterals
   * The 2nd number is the total number of faces, such that the different "2nd - 1st" is the number of triangular faces
   **/
  const unsigned NFC[N_GEOM_ELS][2] = {
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
  const unsigned ig[ N_GEOM_ELS ][ MAXIMUM_NUMBER_OF_FACES_PER_GEOM_EL ][ MAXIMUM_NUMBER_OF_NODES_PER_FACE ] = {
    { {0, 1, 5, 4, 8, 17, 12, 16, 20},
      {1, 2, 6, 5, 9, 18, 13, 17, 21},
      {2, 3, 7, 6, 10, 19, 14, 18, 22},
      {3, 0, 4, 7, 11, 16, 15, 19, 23},
      {0, 3, 2, 1, 11, 10, 9, 8, 24},
      {4, 5, 6, 7, 12, 13, 14, 15, 25}
    },
    { {0, 2, 1, 6, 5, 4, 10},
      {0, 1, 3, 4, 8, 7, 11},
      {1, 2, 3, 5, 9, 8, 12},
      {2, 0, 3, 6, 7, 9, 13}
    },
    { {0, 1, 4, 3, 6, 13, 9, 12, 15},
      {1, 2, 5, 4, 7, 14, 10, 13, 16},
      {2, 0, 3, 5, 8, 12, 11, 14, 17},
      {0, 2, 1, 8, 7, 6, 18},
      {3, 4, 5, 9, 10, 11, 19}
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
  
  /**
   * Number of elements obtained with one refinement
  **/
  const unsigned NRE[N_GEOM_ELS] = {8, 8, 8, 4, 4, 2};

  
  const unsigned referenceElementDirection[N_GEOM_ELS][3][2] = { //Endpoint1, Endpoint2 =rEED[elemem type][direction][0,1]
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

// === Geometric Element, Single - END =================
    
    
    
// === Geometric Element, FE, Single - BEGIN =================
  public:

      const unsigned GetNVE(const unsigned& elementType, const unsigned& doftype) const;
      
      const unsigned GetNFACENODES(const unsigned& elementType, const unsigned& jface, const unsigned& dof) const;
      
  private:
    
    
  /**
   * Number of degrees of freedom per geometric element per FE family:
   * linear, quadratic, biquadratic, piecewise costant, piecewise linear discontinuous
  **/
  const unsigned NVE[N_GEOM_ELS][NFE_FAMS] = {
    {8, 20, 27, 1, 4}, //hex
    {4, 10, 15, 1, 4}, //tet
    {6, 15, 21, 1, 4}, //wedge
    {4, 8, 9, 1, 3}, //quad
    {3, 6, 7, 1, 3}, //tri
    {2, 3, 3, 1, 2}  //line
  };
  
    
  const unsigned NFACENODES[ N_GEOM_ELS ][ MAXIMUM_NUMBER_OF_FACES_PER_GEOM_EL ][ NFE_FAMS_C_ZERO_LAGRANGE ] = {
    { {4, 8, 9}, // Hex
      {4, 8, 9},
      {4, 8, 9},
      {4, 8, 9},
      {4, 8, 9},
      {4, 8, 9}
    },
    { {3, 6, 7}, // Tet
      {3, 6, 7},
      {3, 6, 7},
      {3, 6, 7}
    },
    { {4, 8, 9}, // Wedge
      {4, 8, 9},
      {4, 8, 9},
      {3, 6, 7},
      {3, 6, 7}
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
    
// === Geometric Element, FE, Single - END =================





// === Mesh, Subdomains - BEGIN =================
  public:
    
      void SetElementOffsets(const std::vector < unsigned > & elementOffset, const unsigned &iproc, const unsigned &nprocs) {
        _elementOffset = elementOffset;
        _elementOwned = elementOffset[iproc + 1] - elementOffset[iproc];
        _iproc = iproc;
        _nprocs = nprocs;
      }
      
    private:

      unsigned _iproc;
      unsigned _nprocs;
      
      /** @todo Same as in Mesh, see if we can avoid duplication */
      std::vector < unsigned > _elementOffset;
      unsigned _elementOwned;

// === Mesh, Subdomains - END =================


// === Basic, Dimension - BEGIN =================
  public:
      /** To be Added */
      unsigned GetDimension() const { return _dim; }
      
  private:

      /** Dimension of the underlying Mesh */
      unsigned _dim;
// === Basic, Dimension - END =================

// === Basic, Level - BEGIN =================
  public:
     
      
  private:

      /** Pointer to the list of coarser elements */
      elem* _coarseElem;
      
      /** level of refinement of this list of elements */
      unsigned _level;
           
      
// === Basic, Level - END =================


// === Elements, Numbers - BEGIN =================
    public:
  
      /**    Return the number of elements of a certain shape, specified as an input string. Otherwise return the number of all elements */
      unsigned GetElementNumber(const char* name = "All") const;

      /** To be Added */
      unsigned GetRefinedElementNumber() const {
        return _nelr;
      };

      /** To be Added */
      void SetRefinedElementNumber(const unsigned& value) {
        _nelr = value;
      };

      /** To be Added */
      void AddToElementNumber(const unsigned& value, const char name[]);

      /** To be Added */
      void AddToElementNumber(const unsigned& value, short unsigned ielt);

    private:
      
      /** Number of elements of the Mesh */
      unsigned _nel;
      /** Number of elements of the Mesh for each Geometric type */
      unsigned _nelt[N_GEOM_ELS];
      /** Number of refined elements */
      unsigned _nelr;
// === Elements, Numbers - END =================


// === Elements, Type - BEGIN =================
  public:
    
      /** To be Added */
      unsigned GetIndex(const char name[]) const;

      /** To be Added */
      short unsigned GetElementType(const unsigned& iel);

      /** To be Added */
      MyVector< short unsigned > & GetElementTypeArray() { return _elementType; }
      
      /** To be Added */
      void SetElementType(const unsigned& iel, const short unsigned& value);
      
    
  private:
    
      MyVector< short unsigned> _elementType;
    
// === Elements, Type - END =================

      
// === Elements, Level - BEGIN =================
  public:
      void SetElementLevel(const unsigned& iel, const short unsigned& level) {
        _elementLevel[iel] = level;
      }
      
      short unsigned GetElementLevel(const unsigned &jel) {
        return _elementLevel[jel];
      }
      
    
  private:
    
      MyVector< short unsigned> _elementLevel;
      
// === Elements, Level - END =================

      
// === Elements, Groups - BEGIN =================

  public:

      /** To be Added */
      short unsigned GetElementGroup(const unsigned& iel);

      /** To be Added */
      void SetElementGroup(const unsigned& iel, const short unsigned& value);

      /** Number of groups  */
      unsigned GetElementGroupNumber() const;

      /** Number of groups  */
      void SetElementGroupNumber(const unsigned& value);
      
  private:

      /** @todo Number of groups of elements - in Gambit it is explicit at the top of the file */
      unsigned _ngroup;
      
      MyVector< short unsigned> _elementGroup;
      
// === Elements, Groups - END =================


// === Elements, Materials - BEGIN =================

  public:
      
      /** To be Added */
      void SetElementMaterial(const unsigned& iel, const short unsigned& value);

      /** To be Added */
      short unsigned GetElementMaterial(const unsigned& iel);

      void SetMaterialElementCounter( std::vector<unsigned> materialElementCounter){
        _materialElementCounter = materialElementCounter;
      }
      
      std::vector<unsigned> GetMaterialElementCounter(){
        return _materialElementCounter;
      }

  private:
    
      MyVector< short unsigned> _elementMaterial;
      /** Volume elements, 3 types of material */
      std::vector<unsigned> _materialElementCounter;
    
// === Elements, Materials - END =================



// === Elements, Elements with a common Face to the current element - BEGIN =================
   public:
     
      /** To be Added */
      void SetFaceElementIndex(const unsigned& iel, const unsigned& iface, const int& value);

      /** To be Added */
      int GetFaceElementIndex(const unsigned& iel, const unsigned& iface);

      int GetBoundaryIndex(const unsigned& iel, const unsigned& iface);
     
      /** To be added */
      void BuildElementNearFace();
      
      void ShrinkToFitElementNearFace();
      void ScatterElementNearFace();
      void LocalizeElementNearFace(const unsigned& jproc);
      void FreeLocalizedElementNearFace();

      MyMatrix <int> &  GetElementNearFaceArray() { return _elementNearFace; } 


   private:
     
      /** For every element, it is initialized to -1, and if there is a near element attached to a face, it stores that value. It is used for BCs as well */
      MyMatrix <int> _elementNearFace;

// === Elements, Elements with a common Face to the current element - END =================
 

// === Elements, all elements near the current element, including those with a common vertex - BEGIN =================
  public:
      
      void BuildElementNearElement();

      const unsigned GetElementNearElementSize(const unsigned& iel, const unsigned &layers)  {
        return (layers == 0) ? 1 : _elementNearElement.end(iel);
      };

      const unsigned GetElementNearElement(const unsigned& iel, const unsigned &j)  {
        return _elementNearElement[iel][j];
      };

      
   private:
     
      /** For each element, it gives the elements that are near the given element, including those that are only touching a common vertex
        @todo I think this should be scattered, or maybe not, if it is only used temporarily */
      MyMatrix <unsigned> _elementNearElement;
// === Elements, all elements near the current element, including those with a common vertex - END =================



// === Elements, Refinement - BEGIN =================
  public:
      
   private:
// === Elements, Refinement - END =================


      
// === Nodes - BEGIN =================
  public:

      /** To be Added */
      unsigned GetNodeNumber()const;

      /** To be Added */
      void SetNodeNumber(const unsigned& value);

  
  private:
      
      /** Number of nodes of the Mesh */
      unsigned _nvt;
// === Nodes - END =================


// === Nodes, Elements having that Node as a vertex - BEGIN =================
  public:
      
      /** To be Added */
      void BuildElementNearVertex();

      void DeleteElementNearVertex();

      /** To be Added */
      unsigned GetElementNearVertexNumber(const unsigned& inode);

      /** To be Added */
      unsigned GetElementNearVertex(const unsigned& inode, const unsigned& jnode);

   private:
      
      /** For each Node, it gives the list of elements having that Node as a vertex 
        @todo I think this should be scattered, or maybe not, if it is only used temporarily */
      MyMatrix <unsigned> _elementNearVertex;
      
// === Nodes, Elements having that Node as a vertex - END =================      


// === DOF, Local->Global (element-based) Dofmap for 1 scalar variable - BEGIN =================
  public:
      
      /** To be Added */
      unsigned GetFaceVertexIndex(const unsigned& iel, const unsigned& iface, const unsigned& inode);

      void ShrinkToFitElementDof();
      void ScatterElementDof();
      void LocalizeElementDof(const unsigned &jproc);
      void FreeLocalizedElementDof();

      /** To be Added */
      unsigned GetElementDofNumber(const unsigned& iel, const unsigned& type);

      /** Return the local->global node number */
      unsigned GetElementDofIndex(const unsigned& iel, const unsigned& inode);

      /** To be Added */
      void SetElementDofIndex(const unsigned& iel, const unsigned& inode, const unsigned& value);
      
  private:
     /** For each element, gives the conversion from local node index to global node index */
      MyMatrix <unsigned> _elementDof;
// === DOF, Local->Global (element-based) Dofmap for 1 scalar variable - END =================

   private:
     
      /** This is only going to all levels except the finest one  @todo I think it contains for each element the list of its child elements */
      MyMatrix <unsigned> _childElem;
      /** This is only going to all levels except the finest one */
      MyMatrix <unsigned> _childElemDof;


      
// === Refinement, AMR - BEGIN =================
    public:
      
      void GetAMRRestriction(Mesh *msh);
// === Refinement, AMR - END =================
      

  };
  






} //end namespace femus




#endif


// ******************** GEOMETRIC ELEMENTS, NODE NUMBERING - BEGIN **************************

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
//      |   6   \
//      |         \
//      0-----3----1


//
//	0-----2-----1
//
*/
// ******************** GEOMETRIC ELEMENTS, NODE NUMBERING - END **************************