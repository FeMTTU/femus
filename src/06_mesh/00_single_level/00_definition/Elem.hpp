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
   * @todo Some GeomEl information has to be moved to the basic classes
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


// === Geometric Element, Single - BEGIN =================
  public:
    
      /** To be Added */
      const unsigned GetElementFaceNumber(const unsigned& iel, const unsigned& type = 1) const   {    return GetNFC(_elementType[iel] , type);  }
  
      const unsigned GetNFC(const unsigned& elementType, const unsigned& type) const  {    return NFC[elementType][type];  }
  
      const unsigned GetIG(const unsigned& elementType, const unsigned& iface, const unsigned& jnode) const   {    return ig[elementType][iface][jnode];  }
  
      const unsigned GetNRE(const unsigned& elementType) const  { return NRE[elementType]; }
      
      const unsigned GetReferenceElementDirection(const unsigned& elementType, const unsigned dir, const unsigned node) const { 
        return referenceElementDirection[elementType][dir][node]; }
      
  private:

  /**
   * Number of FACES(3D), edges(2D) or point-extrema(1D) for each considered element
   * The 1st number is the quadrilaterals
   * The 2nd number is the total number of faces, such that the difference "2nd - 1st" is the number of triangular faces
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
    


// === Basic, Dimension - BEGIN =================
  public:
      /** To be Added */
      unsigned GetDimension() const { return _dim; }
      
  private:

      /** Dimension of the underlying Mesh */
      unsigned _dim;
// === Basic, Dimension - END =================

      
      

// === Elements - BEGIN ===========================================================
    public:
      
      //Elem
      // reorder the element according to the new element mapping
      void ReorderMeshElement_Type_Level_Group_Material___NearFace_rows_ChildElem_columns(const std::vector < unsigned >& elementMapping);

      void ScatterElement_Level_Type_Group_Material___NearFace();
      
      void BuildElem_NearFace_NearElem_using_NearVertex();
      
      
  private:
      
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
      
  public:
      
      void LocalizeElement_Level_Type_Group_Material(const unsigned &lproc) {
        _elementLevel.broadcast(lproc);
        _elementType.broadcast(lproc);
        _elementGroup.broadcast(lproc);
        _elementMaterial.broadcast(lproc);
      }
      
      void FreeLocalizedElement_Level_Type_Group_Material() {
        _elementLevel.clearBroadcast();
        _elementType.clearBroadcast();
        _elementGroup.clearBroadcast();
        _elementMaterial.clearBroadcast();
      }


      
      
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

     static  unsigned int  InitializeNumberOfElementsFromCoarseList(elem* elc, const unsigned refindex);
      
    private:
      
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
    
      /** To be Added */
      unsigned GetIndex(const char name[]) const;

      /** To be Added */
      const short unsigned GetElementType(const unsigned& iel) const;

      /** To be Added */
      MyVector< short unsigned > & GetElementTypeArray() { return _elementType; }
      
      /** To be Added */
      void SetElementType(const unsigned& iel, const short unsigned& value);
      
    
  private:
    
      MyVector< short unsigned> _elementType;
    
// === Elements, Type - END =================


// === Elements, Subdomains - BEGIN =================
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

// === Elements, Subdomains - END =================


      
// === Elements, Level - BEGIN =================
  public:
      void SetElementLevel(const unsigned& iel, const short unsigned& level) {
        _elementLevel[iel] = level;
      }
      
      const short unsigned GetElementLevel(const unsigned &jel) const {
        return _elementLevel[jel];
      }
      
      const bool GetIfElementCanBeRefined(const unsigned& iel) const {
        return (_elementLevel[iel] == _level) ? true : false;
      }
      
      const bool GetIfFatherHasBeenRefined(const unsigned& iel) const {
        return GetIfElementCanBeRefined(iel);
      }
    
  private:
    
      MyVector< short unsigned> _elementLevel;
      
      /** Pointer to the list of coarser elements */
      elem* _coarseElem;
      
      /** level of refinement of this list of elements */
      unsigned _level;
      
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



// === Elements, for Each Element give the Elements with a common Face to the current element - BEGIN =================
   public:
     
     //Elem
      void ReorderElementNearFace_rows(const std::vector < unsigned >& elementMapping);
      
 
      /** To be Added */
      void SetFaceElementIndex(const unsigned& iel, const unsigned& iface, const int& value);

      /** To be Added */
      int GetFaceElementIndex(const unsigned& iel, const unsigned& iface);

      int GetBoundaryIndex(const unsigned& iel, const unsigned& iface);
     
      void ShrinkToFitElementNearFace();
      void LocalizeElementNearFace(const unsigned& jproc);
      void FreeLocalizedElementNearFace();

      MyMatrix <int> &  GetElementNearFaceArray() { return _elementNearFace; } 


   private:
     
      void ScatterElementNearFace();
      
      /** To be added */
      void BuildElementNearFace();
      
      /** For every element, it is initialized to -1, and if there is a near element attached to a face, it stores that value. It is used for BCs as well */
      MyMatrix <int> _elementNearFace;

// === Elements, for Each Element give the Elements with a common Face to the current element - END =================
 

// === Elements, for Each Element give all elements near the current element, including those with a common vertex - BEGIN =================
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
       * It is used for Domain Decomposition and for AMR
        @todo why is this not scattered??? */
      MyMatrix <unsigned> _elementNearElement;
// === Elements, for Each Element give all elements near the current element, including those with a common vertex - END =================


// === Elements, for Each Element gives the children elements - BEGIN =================
  public:

      //Elem
      void ReorderChildElement_OnCoarseElem_columns(const std::vector < unsigned >& elementMapping);
     
      void AllocateChildrenElement(const unsigned int& refindex, const Mesh* msh);

      /** To be Added */
      void SetChildElement(const unsigned& iel, const unsigned& json, const unsigned& value);

      /** To be Added */
      unsigned GetChildElement(const unsigned& iel, const unsigned& json);
      
   private:
     
      /** This is only going to all levels except the finest one  
       *  It contains for each element the list of its child elements */
      MyMatrix <unsigned> _childElem;
      
// === Elements, for Each Element gives the children elements - END =================

      
      
// === Elements - END ===========================================================

      
      
// === Nodes - BEGIN ===========================================================
      
      
// === Nodes, Number - BEGIN =================
  public:

      /** To be Added */
      unsigned GetNodeNumber()const;

      /** To be Added */
      void SetNodeNumber(const unsigned& value);

  
  private:
      
      /** Number of nodes of the Mesh */
      unsigned _nvt;
// === Nodes, Number - END =================


// === Nodes, for Each Node give the Elements having that Node as a vertex (temporary then deleted) - BEGIN =================
  public:
      
      /** To be Added */
      void BuildElementNearVertex();

      /** To be Added */
      unsigned GetElementNearVertexNumber(const unsigned& inode);

      /** To be Added */
      unsigned GetElementNearVertex(const unsigned& inode, const unsigned& jnode);

   private:
      
      void DeleteElementNearVertex();

      /** For each Node, it gives the list of elements having that Node as a vertex 
        @todo It is used in mesh construction and then it is deleted, so beware of not using it afterwards */
      MyMatrix <unsigned> _elementNearVertex;
      
// === Nodes, for Each Node give the Elements having that Node as a vertex (temporary then deleted) - END =================      


// === Nodes - END ===========================================================


// =========       
// ========= Previously, it was all info of geometric elements. From now on, there is also FE information ==========      
// =========      



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

      

// === DOF, for Each Element return the dof of 1 scalar variable  (Local->Global (element-based) Dofmap for 1 scalar variable) - BEGIN =================
  public:
      
      /** To be Added */
      unsigned GetFaceVertexIndex(const unsigned& iel, const unsigned& iface, const unsigned& inode);

      void ShrinkToFitElementDof();
      void LocalizeElementDof(const unsigned &jproc);
      void FreeLocalizedElementDof();

      /** Return the local->global node number */
      unsigned GetElementDofIndex(const unsigned& iel, const unsigned& inode);

      /** To be Added */
      void SetElementDofIndex(const unsigned& iel, const unsigned& inode, const unsigned& value);
      
      void ScatterElementDof();
      
      //Dof
      void ReorderMeshElement_Dof_stuff(const std::vector < unsigned >& elementMapping);
      
      //Dof
      void ReorderElementDof_rows(const std::vector < unsigned >& elementMapping);
      
      //Dof
      // reorder the nodes according to the new node mapping
      void ReorderElementDof_columns_Using_node_mapping(const std::vector < unsigned >& nodeMapping);

      

     
  private:
    
     /** For each element, gives the conversion from local node index to global node index */
      MyMatrix <unsigned> _elementDof;
// === DOF, for Each Element return the dof of 1 scalar variable  (Local->Global (element-based) Dofmap for 1 scalar variable) - END =================

      
// === DOF, for Each Element return the dofs of all its children (or only of itself if it is not a refined element), for 1 scalar variable - BEGIN =================
  public:

      void AllocateChildrenElementDof(const unsigned int& refindex, const Mesh* msh);
      
      /** To be Added */
      void SetChildElementDof(elem* elf);

      unsigned GetChildElementDof(const unsigned& iel, const unsigned& i0, const unsigned i1);
      
   private:
     
      /** This is only going to all levels except the finest one */
      MyMatrix <unsigned> _childElemDof;
// === DOF, for Each Element return the dofs of all its children (or only of itself if it is not a refined element), for 1 scalar variable - END =================

      


      
// === Refinement, AMR - BEGIN =================
    public:
      
      void GetAMRRestriction(Mesh *msh) const;
// === Refinement, AMR - END =================
      

  };
  






} //end namespace femus




#endif


