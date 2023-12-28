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
#include "GeomElemBase.hpp"
#include "FElemTypeEnum_list.hpp"
#include "MyVector.hpp"
#include "MyMatrix.hpp"

#include "MeshGeneration.hpp"

#include <vector>
#include <map>
#include <cmath>



namespace femus {

  //Forward declarations  
  class Mesh;

  /**
   * The elem class: it contains the list of all Mesh Geometric Elements, along with several Element-based and also Node-based properties
   * @todo I believe it would even be more linear if this class did not have any function at all involving the Mesh pointer, there are very few in any case
   * @todo Some GeomEl information has to be moved to the basic classes
   * The idea is that the elem class does not modify the Mesh class, but only the other way around
  */
  class elem {

// === Friend functions and classes - BEGIN ===============

friend class Mesh;

// generation BEGIN
friend class MeshTools::Generation;
friend class MED_IO;
friend class GambitIO;
// generation END

// partitioning BEGIN
// friend class MeshPartitioning;
friend class MeshMetisPartitioning;
// partitioning END


// refinement BEGIN
friend class MeshRefinement;
// refinement END

// === Friend functions and classes - END =================


// === Constructors / Destructor - BEGIN =================
    private:

      /** constructors */
      elem(const unsigned& other_nel, const unsigned dim_in);

      //elem(elem* elc, const unsigned refindex, const std::vector < double >& coarseAmrLocal, const std::vector < double >& localizedElementType);
      elem(elem* elc, const unsigned dim_in, const unsigned refindex, const std::vector < double >& coarseAmrLocal);
// === Constructors / Destructor - END =================


// === Geometric Element, Single - BEGIN =================
  public:
    
      /** To be Added */
      const unsigned GetElementFaceNumber(const unsigned& iel, const unsigned& type = GeomElemBase::_index_for_all_faces) const   {    return GetNFC(_elementType[iel] , type);  }
  
      const unsigned GetNFC(const unsigned& elementType, const unsigned& type) const  {    return NFC[elementType][type];  }
  
      const unsigned GetIG(const unsigned& elementType, const unsigned& iface, const unsigned& jnode) const   {    return ig[elementType][iface][jnode];  }
  
      const unsigned GetReferenceElementDirection(const unsigned& elementType, const unsigned dir, const unsigned node) const {
        return directions_of_reference_element[elementType][dir][node];
      }
      
  private:

  /**
   * Number of FACES(3D), edges(2D) or point-extrema(1D) for each considered element
   * The 1st number is the quadrilaterals
   * The 2nd number is the total number of faces, "quadrilaterals + non-quadrilaterals", such that the difference "2nd - 1st" is the number of non-quadrilateral faces
   **/
  const unsigned NFC[N_GEOM_ELS][ GeomElemBase::_n_face_types_max ] = {
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
  
  
  const unsigned directions_of_reference_element[N_GEOM_ELS][3][2] = { //Endpoint1, Endpoint2 =rEED[elemem type][direction][0,1]
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

// const unsigned origin_of_reference_element[6]={26,0,12/*?*/,8,0,2};

// === Geometric Element, Single - END =================


// === Geometric Element, Single, REFINEMENT - BEGIN =================
  public:

    /** It would private if we put Solution as another friend of Mesh, but maybe it is too much for now */
    ///8 elements from refining 1 HEX, TET, WEDGE; 4 elements from refining 1 QUAD, TRI; 2 elements from refining 1 LINE
    /** MESH, REF: 8 elements from refining 1 HEX, TET, WEDGE; 4 elements from refining 1 QUAD TRI; 2 elements from refining 1 LINE // 8*DIM[2]+4*DIM[1]+2*DIM[0]; */
    const unsigned GetRefIndex(const unsigned dim) const {
      return pow(2, dim);
    }


  private:

    /** MESH, REF: 4 faces from refining 1 QUAD TRI; 2 faces from refining 1 LINE; 1 face from refining 1 point // 4*DIM[2]+2*DIM[1]+1*DIM[0]; */
    const unsigned GetRefFaceIndex(const unsigned dim) const {
      return pow(2, dim -1u);
    }

  /** Unused now */
  const unsigned GetNRE(const unsigned& elementType) const  { return NRE[elementType]; }

  /**
   * Number of elements obtained with one refinement
  **/
  const unsigned NRE[N_GEOM_ELS] = {8, 8, 8, 4, 4, 2};

// === Geometric Element, Single, REFINEMENT - END =================


// === Mesh, Basic, Dimension - BEGIN =================
  private:

      /** To be Added */
      unsigned GetDimension() const { return _dim; }

      /** Dimension of the underlying Mesh */
      unsigned _dim;
// === Mesh, Basic, Dimension - END =================

      
      

// === Mesh, Elements - BEGIN ===========================================================
  public:

      /// @todo These are called by some Writers, but I don't want to make them friends here
      void LocalizeElement_Level_Type_Group_Material(const unsigned &lproc) {
        _elementLevel.broadcast(lproc);
        _elementType.broadcast(lproc);
        _elementGroup.broadcast(lproc);
        _elementMaterial.broadcast(lproc);
      }

      /// @todo These are called by some Writers, but I don't want to make them friends here
      void FreeLocalizedElement_Level_Type_Group_Material() {
        _elementLevel.clearBroadcast();
        _elementType.clearBroadcast();
        _elementGroup.clearBroadcast();
        _elementMaterial.clearBroadcast();
      }


    private:
      
      // reorder the element according to the new element mapping
      void ReorderMeshElement_Type_Level_Group_Material___NearFace_rows_ChildElem_columns(const std::vector < unsigned >& elementMapping);

      void ScatterElement_Level_Type_Group_Material___NearFace();
      
      void BuildElem_NearFace_NearElem_using_NearVertex();
      

      void ResizeElement_Level_Type_Group_Material(const unsigned nel_in, const unsigned level_in) {
       _elementLevel.resize(nel_in, level_in);
       _elementType.resize(nel_in);
       _elementGroup.resize(nel_in);
       _elementMaterial.resize(nel_in);
      }
      
      void ScatterElement_Level_Type_Group_Material() {
        _elementLevel.scatter(_elementOffset);
        _elementType.scatter(_elementOffset);
        _elementGroup.scatter(_elementOffset);
        _elementMaterial.scatter(_elementOffset);
      }


// === Elements, Numbers - BEGIN =================
    private:
      
      /**    Return the number of elements of a certain shape, specified as an input string. Otherwise return the number of all elements */
      unsigned GetElementNumber(const std::string name = "All") const;

      /** To be Added */
      unsigned GetRefinedElementNumber() const {
        return _nelr;
      };

      static  unsigned int  InitializeNumberOfElementsFromCoarseList(elem* elc, const unsigned refindex);

      /** To be Added */
      void AddToElementNumber(const unsigned& value, const std::string name);

      /** To be Added */
      void AddToElementNumber(const unsigned& value, short unsigned ielt);

      /** To be Added */
      void SetRefinedElementNumber(const unsigned& value) {
        _nelr = value;
      };

      void InitializeNumberOfElementsPerGeomType() {
            for (unsigned g = 0; g < N_GEOM_ELS ; g++) { _nelt[g] = 0; }
      }
      
      /** Number of elements of the Mesh */
      unsigned _nel;
      /** Number of elements of the Mesh for each Geometric type */
      unsigned _nelt[N_GEOM_ELS];
      /** Number of refined elements */
      unsigned _nelr;
// === Elements, Numbers - END =================


// === Elements, Type - BEGIN =================
  public:
    
      /** @todo this is not const because it does broadcast at one point */
      MyVector< short unsigned > & GetElementTypeArray() { return _elementType; }
      
  private:
    
      /** To be Added */
      unsigned GetIndex(const std::string name) const;

      /** To be Added */
      const short unsigned GetElementType(const unsigned& iel) const;

      /** To be Added */
      void SetElementType(const unsigned& iel, const short unsigned& value);

      MyVector< short unsigned> _elementType;
    
// === Elements, Type - END =================


// === Elements, Subdomains - BEGIN =================
    private:

      unsigned _iproc;
      unsigned _nprocs;
      
      /** @todo Same as in Mesh, see if we can avoid duplication */
      std::vector < unsigned > _elementOffset;
      unsigned _elementOwned;
      
      void SetElementOffsets(const std::vector < unsigned > & elementOffset, const unsigned &iproc, const unsigned &nprocs) {

        _elementOffset = elementOffset;
        _elementOwned = elementOffset[iproc + 1] - elementOffset[iproc];

        SetProcs(iproc, nprocs);

      }

      void SetProcs(const unsigned &iproc, const unsigned &nprocs) {
        _iproc = iproc;
        _nprocs = nprocs;
      }

     unsigned get_n_elements_owned() const { return _elementOwned; }

// === Elements, Subdomains - END =================


      
// === Elements, Level - BEGIN =================
  public:
    
      const short unsigned GetElementLevel(const unsigned &jel) const {
        return _elementLevel[jel];
      }
      
      const bool GetIfElementCanBeRefined(const unsigned& iel) const {
        return (_elementLevel[iel] == _level) ? true : false;
      }
      
  private:
    
      const unsigned GetLevelOfRefinementForList() const {
        return _level;
      }

      const bool GetIfFatherHasBeenRefined(const unsigned& iel) const {
        return GetIfElementCanBeRefined(iel);
      }

      const MyVector<short unsigned> &  GetElementLevelArray() const { return _elementLevel; }

      MyVector<short unsigned> _elementLevel;
      
      /** Pointer to the list of coarser elements */
      elem* _coarseElem;
      
      /** level of refinement of this list of elements */
      unsigned _level;
      
      void SetElementLevel(const unsigned& iel, const short unsigned& level) {
        _elementLevel[iel] = level;
      }

// === Elements, Level - END =================

      
// === Elements, Groups - BEGIN =================

  public:

      /** To be Added */
      short unsigned GetElementGroup(const unsigned& iel) const;

      /** @todo If it wasn't for a nonlocal assembly, it would be private */
      void SetElementGroup(const unsigned& iel, const short unsigned& value);

  private:

      /** Number of groups  */
      unsigned GetElementGroupNumber() const;

      /** Number of groups  */
      void SetElementGroupNumber(const unsigned& value);

      /** @todo Number of groups of elements - in Gambit it is explicit at the top of the file */
      unsigned _ngroup;
      
      MyVector< short unsigned> _elementGroup;
      
// === Elements, Groups - END =================


// === Elements, Materials - BEGIN =================

  public:
      
      /** To be Added */
      short unsigned GetElementMaterial(const unsigned& iel) const;


  private:
    
      /** To be Added */
      void SetElementMaterial(const unsigned& iel, const short unsigned& value);

      MyVector< short unsigned> _elementMaterial;
    
// === Elements, Materials - END =================



// === Elements, for Each Element give the Elements with a common Face to the current element - BEGIN =================
   public:
 
      /** To be Added */
      int GetFaceElementIndex(const unsigned& iel, const unsigned& iface) const;

      int GetBoundaryIndex(const unsigned& iel, const unsigned& iface) const;
     
   private:
     
      void ShrinkToFitElementNearFace();
      void LocalizeElementNearFace(const unsigned& jproc);
      void FreeLocalizedElementNearFace();

      const MyMatrix <int> &  GetElementNearFaceArray() const { return _elementNearFace; }

      MyMatrix <int> &  GetElementNearFaceArray() { return _elementNearFace; }

      /** To be Added */
      void SetFaceElementIndex(const unsigned& iel, const unsigned& iface, const int& value);

      void ReorderElementNearFace_rows(const std::vector < unsigned >& elementMapping);

      void ScatterElementNearFace();
      
      /** To be added */
      void BuildElementNearFace();
      
      /** For every element, it is initialized to -1, and if there is a near element attached to a face, it stores that value. It is used for BCs as well */
      MyMatrix <int> _elementNearFace;

// === Elements, for Each Element give the Elements with a common Face to the current element - END =================
 

// === Elements, for Each Element give all elements near the current element, including those with a common vertex - BEGIN =================
  public:
      
      const unsigned GetElementNearElementSize(const unsigned& iel, const unsigned &layers) const {
        return (layers == 0) ? 1 : _elementNearElement.end(iel);
      };

      const unsigned GetElementNearElement(const unsigned& iel, const unsigned &j) const {
        return _elementNearElement[iel][j];
      };

      
   private:
     
      void BuildElementNearElement();

      /** For each element, it gives the elements that are near the given element, including those that are only touching a common vertex
       * It is used for Domain Decomposition and for AMR
        @todo why is this not scattered??? */
      MyMatrix <unsigned> _elementNearElement;
// === Elements, for Each Element give all elements near the current element, including those with a common vertex - END =================


// === Elements, for Each Element gives the children elements - BEGIN =================
   private:
     
      /** To be Added */
      unsigned GetChildElement(const unsigned& iel, const unsigned& json) const;

      /** To be Added */
      void SetChildElement(const unsigned& iel, const unsigned& json, const unsigned& value);

      void AllocateChildrenElement(const unsigned int& refindex, const Mesh* msh);

      void ReorderChildElement_OnCoarseElem_columns(const std::vector < unsigned >& elementMapping);

      /** This is only going to all levels except the finest one
       *  It contains for each element the list of its child elements */
      MyMatrix <unsigned> _childElem;
      
// === Elements, for Each Element gives the children elements - END =================

      
      
// === Mesh, Elements - END ===========================================================

      
      
// === Mesh, Nodes - BEGIN ===========================================================
      
      
// === Nodes, Number - BEGIN =================
  private:
      
      /** To be Added */
      unsigned GetNodeNumber() const;

      /** To be Added */
      void SetNodeNumber(const unsigned& value);

      /** Number of nodes of the Mesh */
      unsigned _nvt;
// === Nodes, Number - END =================


// === Nodes, for Each Node give the Elements having that Node as a vertex (temporary then deleted) - BEGIN =================
   private:
      
      /** To be Added */
      unsigned GetElementNearVertexNumber(const unsigned& inode) const;

      /** To be Added */
      unsigned GetElementNearVertex(const unsigned& inode, const unsigned& jnode) const;

      /** To be Added */
      void BuildElementNearVertex();

      void DeleteElementNearVertex();

      /** For each Node, it gives the list of elements having that Node as a vertex 
        @todo It is used in mesh construction and then it is deleted, so beware of not using it afterwards */
      MyMatrix <unsigned> _elementNearVertex;
      
// === Nodes, for Each Node give the Elements having that Node as a vertex (temporary then deleted) - END =================      


// === Mesh, Nodes - END ===========================================================


     
// ========= Previously, it was all info of geometric elements. From now on, there is also FE information - BEGIN ==========



// === Geometric Element, FE, Single (Local) - BEGIN =================
  public:

  /**
   * Return the number of vertices(type=0) + midpoints(type=1) + facepoints(type=2) + interiorpoits(type=2)
   **/
      const unsigned GetElementDofNumber(const unsigned& iel, const unsigned& type) const      {    return  GetNVE(_elementType[iel], type);   }
      
      const unsigned GetNVE(const unsigned& elementType, const unsigned& doftype) const        {    return NVE[elementType][doftype];     }
      
      const unsigned GetNFACENODES(const unsigned& elementType, const unsigned& jface, const unsigned& dof) const  {        return NFACENODES[elementType][jface][dof];   }
      
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
    
// === Geometric Element, FE, Single (Local) - END =================

      

// === Mesh, DOF, for Each Element return the dof of 1 scalar variable  (Local->Global (element-based) Dofmap for 1 scalar variable) - BEGIN =================
  public:

      /// @todo These are called by some Writers, but I don't want to make them friends here
      void LocalizeElementDof(const unsigned &jproc);

      /// @todo These are called by some Writers, but I don't want to make them friends here
      void FreeLocalizedElementDof();

      /** Return the local->global node number */
      unsigned GetElementDofIndex(const unsigned& iel, const unsigned& inode) const;


  private:
    
      //Dof
      void ReorderMeshElement_Dof_stuff(const std::vector < unsigned >& elementMapping);

      //Dof
      void ReorderElementDof_rows(const std::vector < unsigned >& elementMapping);

      //Dof
      // reorder the nodes according to the new node mapping
      void ReorderElementDof_columns_Using_node_mapping(const std::vector < unsigned >& nodeMapping);

      void ScatterElementDof();

      void ShrinkToFitElementDof();

      /** To be Added */
      void SetElementDofIndex(const unsigned& iel, const unsigned& inode, const unsigned& value);

      /** To be Added */
      unsigned GetFaceVertexIndex(const unsigned& iel, const unsigned& iface, const unsigned& inode) const;

     /** For each element, gives the conversion from local node index to global node index */
      MyMatrix <unsigned> _elementDof;
// === Mesh, DOF, for Each Element return the dof of 1 scalar variable  (Local->Global (element-based) Dofmap for 1 scalar variable) - END =================

      
// === Mesh, DOF, for Each Element return the dofs of all its children (or only of itself if it is not a refined element), for 1 scalar variable - BEGIN =================
   private:
     
      unsigned GetChildElementDof(const unsigned& iel, const unsigned& i0, const unsigned i1) const;

      void AllocateChildrenElementDof(const unsigned int& refindex, const Mesh* msh);

      /** To be Added */
      void SetChildElementDof(elem* elf);

      /** This is only going to all levels except the finest one */
      MyMatrix <unsigned> _childElemDof;
// === Mesh, DOF, for Each Element return the dofs of all its children (or only of itself if it is not a refined element), for 1 scalar variable - END =================


// ========= Previously, it was all info of geometric elements. From now on, there is also FE information - END ==========      

  };
  


} //end namespace femus




#endif


