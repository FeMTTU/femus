/*=========================================================================

 Program: FEMuS
 Module: Mesh
 Authors: Eugenio Aulisa

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_mesh_Mesh_hpp__
#define __femus_mesh_Mesh_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Elem.hpp"
#include "Solution.hpp"
#include "ElemType.hpp"
#include "ElemTypeEnum.hpp"
#include "ParallelObject.hpp"

#include <cassert>
#include <vector>
#include <map>



namespace femus {



class Solution;

class elem;

/**
 * The mesh class
*/

class Mesh : public ParallelObject {

// =========================
// === CONSTR-DESTR =================
// =========================
public:

    /** Constructor */
    explicit
    Mesh();

    /** destructor */
    ~Mesh();

    
    
// =========================
// === BASIC =================
// =========================
public:
    
    /** Print the mesh info for this level */
    void PrintInfo() const;

    /** MESH: Get the dimension of the problem (1D, 2D, 3D) */
    const unsigned GetDimension() const {
      return Mesh::_dimension;
    }

    /** MESH: Set the dimension of the problem (1D, 2D, 3D) */
    void SetDimension(const unsigned &dim) {
      Mesh::_dimension = dim;
    }

    /** Set the number of nodes */
    void SetNumberOfNodes(const unsigned &nnodes) {
      _nnodes = nnodes;
    };

    /** Get the number of nodes */
    unsigned GetNumberOfNodes() const {
      return _nnodes;
    }

    /** Set the number of element */
    void SetNumberOfElements(const unsigned &nelem) {
      _nelem = nelem;
    }

    /** Get the number of element */
    unsigned GetNumberOfElements() const {
      return _nelem;
    }


private:
    
    /** MESH: dimension of the problem */
    static unsigned _dimension;
    /** MESH: number of elements */
    int _nelem;
    /** MESH: number of nodes */
    unsigned _nnodes;

    /** MESH: node coordinates for each space dimension  @todo beware: this is only filled at coarse reading, then use _topology for the coordinates! */
    std::vector < std::vector < double > > _coords;

    void PrintInfoElements() const;

// =========================
// === BASIC, ELEM =================
// =========================
public:
    
    /** Also DOFMAP: Only for parallel */
    unsigned GetElementDofNumber(const unsigned &iel, const unsigned &type) const;

    /** Also DOFMAP: Only for parallel */
    unsigned GetElementFaceDofNumber(const unsigned &iel, const unsigned jface, const unsigned &type) const;

    /** Get element group*/
    short unsigned GetElementGroup(const unsigned &iel) const;
    /** Get element material*/
    short unsigned GetElementMaterial(const unsigned &iel) const;
    
    /** Get element type*/
    short unsigned GetElementType(const unsigned &iel) const;
    /** Only for parallel */
    const unsigned GetElementFaceType(const unsigned &kel, const unsigned &jface) const;
    
    /** Only for parallel */
    unsigned GetLocalFaceVertexIndex(const unsigned &iel, const unsigned &iface, const unsigned &jnode) const;

    /**  */
    unsigned GetLocalFaceVertexIndex_PassElemType(const short unsigned & el_type, const unsigned& iface, const unsigned& jnode) const;

    /** Only for parallel */
    unsigned GetElementFaceNumber(const unsigned &iel, const unsigned &type = 1) const;  ///@todo why is the default like this

    unsigned GetElementFaceNumber_PassElemType(const short unsigned & el_type, const unsigned& type = 1) const;
    
    void BuildMeshElemStructures();
    
    /** To be added */
    void BuildElementNearFace();

    /** MESH */
    elem * GetElementArray() const {
      return el;
    }
    
    /** MESH: list of all elements */
    elem *el;
    
    /** MESH: Number of elements per processor (incremental count) */
    std::vector < unsigned > _elementOffset;
 
  
    
// =========================
// === BASIC, CharacteristicLength =================
// =========================
public:
    
    /** MESH */
    void SetCharacteristicLength(const double & cLength){
      _cLength = cLength;
    }
    
    /** MESH */
    double GetCharacteristicLength(){
      return _cLength;
    };

    void ComputeCharacteristicLength();
    
    
private:

    /** Order of the domain size */
    double _cLength;

    

/// =========================
/// === FE for single elem =================
/// =========================
public:

    const elem_type * GetFiniteElement(const unsigned geom_elem, const unsigned fe_soltype) const {
        return _finiteElement[geom_elem][fe_soltype];
    }
    
    /** To be Added */
    void SetFiniteElementPtr(/*const*/ elem_type* otheFiniteElement[N_GEOM_ELS][5]);
    
    
    /** FE: Finite Element families, for each Geometric Element @todo this one day should be private */
    const elem_type *_finiteElement[N_GEOM_ELS][5];
    
    
    basis *GetBasis(const short unsigned &ielType, const short unsigned &solType);


    
    
    

    
// =========================
// === COARSE MESH GENERATION =================
// =========================
public:

    /** Only file reading */
    void ReadCoarseMeshFile (const std::string& name, const double Lref, std::vector<bool>& type_elem_flag, const bool read_groups, const bool read_boundary_groups);

      /** This function generates the coarse mesh level, $l_0$, from an input mesh file */
    void ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool> &_finiteElement_flag);

    /** This function generates the coarse mesh level, $l_0$, from an input mesh file, with option to not read groups */
    void ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool> &_finiteElement_flag, const bool read_groups, const bool read_boundary_groups);

    void ReadCoarseMeshBeforePartitioning(const std::string& name, const double Lref, std::vector<bool>& type_elem_flag, const bool read_groups, const bool read_boundary_groups);
  
    /** This function generates a coarse box mesh */
    void GenerateCoarseBoxMesh(const unsigned int nx,
                               const unsigned int ny,
                               const unsigned int nz,
                               const double xmin, const double xmax,
                               const double ymin, const double ymax,
                               const double zmin, const double zmax,
                               const ElemType type, 
                               std::vector<bool> &type_elem_flag);


    void AddBiquadraticNodesNotInMeshFile();
    
    /** Boundary names for faces, I think only used for Box mesh so far */
    std::map<unsigned int, std::string> _boundaryinfo;

private:
    
    
    /** Weights used to build the baricentric coordinate to compute the missing biquadratic nodes **/
    static const double _baricentricWeight[N_GEOM_ELS][5][18];
    
    static const unsigned _numberOfMissedBiquadraticNodes[N_GEOM_ELS];
    
    
    
// =========================
// === PARTITIONING =================
// =========================
public:

    void Partition();
    
    void PartitionForElements(std::vector < unsigned > & partition);
    

// =========================
// === MESH REFINEMENT =================
// =========================
public:

    /** MESH: Set the grid number */
    void SetLevel(const unsigned &i) {
        _level = i;
    };

    /** MESH: Get the grid number */
    unsigned GetLevel() const {
      return _level;
    }

    /** MESH: Set the coarser mesh from which this mesh is generated */
    void SetCoarseMesh( Mesh* otherCoarseMsh ){
      _coarseMsh = otherCoarseMsh;
    }

    void SetRefinementCellAndFaceIndices(const unsigned &dim) {

      Mesh::_ref_index  = pow(2, dim);     //8 elements from refining 1 HEX, TET, WEDGE; 4 elements from refining 1 QUAD TRI; 2 elements from refining 1 LINE
      Mesh::_face_index = pow(2, dim -1u);
    }

    /** MESH */
    const unsigned GetRefIndex() const {
      return Mesh::_ref_index;
    }

    /** MESH */
    const unsigned GetFaceIndex() const {
      return Mesh::_face_index;
    }
    
    /** Get if element is refined*/
    short unsigned GetRefinedElementIndex(const unsigned &iel) const;
    
    
private:
    
    void PrintInfoLevel() const;
    
    /** MESH: level of mesh in the multilevel hierarchy */
    unsigned _level;
    
    /** Pointer to the coarser mesh from which this mesh is generated, it equals NULL if _level = 0 */
    Mesh* _coarseMsh;
    
    /** MESH, REF: 8 elements from refining 1 HEX, TET, WEDGE; 4 elements from refining 1 QUAD TRI; 2 elements from refining 1 LINE */
    static unsigned _ref_index;
    /** MESH, REF: 4 faces from refining 1 QUAD TRI; 2 faces from refining 1 LINE; 1 face from refining 1 point */
    static unsigned _face_index;




// =========================
// === FE DOFMAP =================
// =========================
public:
    
    /** FE: DofMap carriers */
    void initialize_elem_offsets();
    
    void build_elem_offsets_and_reorder_mesh_elem_quantities(const std::vector <unsigned> & partition);
    
    void set_elem_counts();
    
    void mesh_reorder_node_quantities(const std::vector <unsigned> & mapping);
    
    void set_node_counts();

    void deallocate_node_mapping(std::vector < unsigned > & node_mapping) const;
    
    /** Get the dof number for the element -type- */
    unsigned GetTotalNumberOfDofs(const unsigned &type) const {
      return _dofOffset[type][_nprocs];
    }

    unsigned GetSolutionDof(const unsigned &i, const unsigned &iel, const short unsigned &solType) const;

    unsigned GetSolutionDof(const unsigned &i0,const unsigned &i1, const unsigned &ielc, const short unsigned &solType, const Mesh* mshc) const;

    /** Performs a bisection search to find the processor of the given dof */
    unsigned IsdomBisectionSearch(const unsigned &dof, const short unsigned &solType) const;

    /** FE: DofMap: Here is where the element and node global orderings are changed based on the partitioning */
    void FillISvector(vector < unsigned > &partition);

    void FillISvectorElemOffsets(std::vector < unsigned >& partition);
  
    void FillISvectorNodeOffsets();
    
    void dofmap_all_fe_families_initialize_dof_offsets();
    
    void dofmap_all_fe_families_clear_ghost_dof_list_other_procs();
    
    void dofmap_Element_based_dof_offsets_build();

    std::vector <unsigned> dofmap_Node_based_dof_offsets_Compute_Node_mapping_and_Node_ownSize();
    
    void dofmap_Node_based_dof_offsets_build_biquadratic();
    
    void dofmap_Node_based_dof_offsets_ghost_nodes_search();
    
    void dofmap_Node_based_dof_offsets_build_linear_quadratic();
    
    /** FE: DofMap: Number of owned nodes per FE family and per processor (count, non-incremental) */
    std::vector < unsigned > _ownSize[5];
    /** FE: DofMap: Number of nodes per FE family and per processor (incremental count) */
    std::vector < unsigned > _dofOffset[5];
    /** FE: DofMap: Number of ghost nodes per FE family and per processor (count, non-incremental) */
    std::vector< std::vector < int > > _ghostDofs[5];
    
private:

    
    /** FE: DofMap  k = 0, 1 */
    std::map < unsigned, unsigned > _ownedGhostMap[2];
    /** FE: DofMap  k = 0, 1 */ 
    std::vector < unsigned > _originalOwnSize[2];
    /** print node-based dofOffset counts */
    void PrintInfoNodes() const;
   
// =========================
// === FE DOFMAP & PROJECTION at SAME LEVEL (needed for node-based printing) =================
// =========================
public:
    
    /**  FE: Get the projection matrix between Lagrange FEM at the same level mesh*/
    SparseMatrix* GetQitoQjProjection(const unsigned& itype, const unsigned& jtype);

    
private:
    
    /** FE: Build the projection matrix between Lagrange FEM at the same level mesh*/
    void BuildQitoQjProjection(const unsigned& itype, const unsigned& jtype);

    /** FE: The projection matrix between Lagrange FEM at the same level mesh */
    SparseMatrix* _ProjQitoQj[3][3];

   
// =========================
// === FE DOFMAP & REFINEMENT =================
// =========================
public:
    
    /**  FE: Get the coarse to the fine projection matrix and use it to restrict only on coarse nodes i.e. projection*/
    SparseMatrix* GetCoarseToFineProjectionRestrictionOnCoarse(const unsigned& solType);

    /**  FE: Get the coarse to the fine projection matrix*/
    SparseMatrix* GetCoarseToFineProjection(const unsigned& solType);

private:
    /** FE: Build the coarse to the fine projection matrix */
    void BuildCoarseToFineProjection(const unsigned& solType, const char el_dofs[]);
    
    /** FE: The coarse to the fine projection matrix */
    SparseMatrix* _ProjCoarseToFine[5];

    
    
    
// =========================
// === TOPOLOGY: Coordinates, AMR, SolidMark (a bit of everything) =================
// =========================
public:
    /** MESH: Coordinates and other stuff */
    Solution* _topology;
    
    /** MESH: Topology */
    const unsigned GetXIndex()          const { return _xIndex; }
    const unsigned GetYIndex()          const { return _yIndex; }
    const unsigned GetZIndex()          const { return _zIndex; }
    const unsigned GetAmrIndex()        const { return _amrIndex; }
    const unsigned GetSolidMarkIndex()  const { return _solidMarkIndex; }
    
    void BuildTopologyStructures();
    
    void Topology_InitializeAndFillCoordinates();
    
    void Topology_InitializeAMR();
    
    void Topology_InitializeAndFillSolidNodeFlag();
    
    /** FSI: Allocate memory for adding fluid or solid mark */
    void AllocateAndMarkStructureNode();
    
    /** Only for parallel */
    bool GetSolidMark(const unsigned &inode) const;
    
    void GetElementNodeCoordinates(std::vector < std::vector <double > > &xv, const unsigned &iel, const unsigned &solType = 2);


  
private:
    
    // indices of the topology parallel vectors
    static const unsigned _xIndex = 0;
    static const unsigned _yIndex = 1;
    static const unsigned _zIndex = 2;
    static const unsigned _amrIndex = 3;
    static const unsigned _solidMarkIndex = 4;



// =========================
// === AMR =================
// =========================
public:
    
    /** AMR */
    static bool (* _SetRefinementFlag)(const std::vector < double >& x, const int &ElemGroupNumber, const int &level);
    static bool _IsUserRefinementFunctionDefined;
    
    /** AMR */
    bool GetIfHomogeneous(){
      return _meshIsHomogeneous;
    }

    /** AMR */
    void SetIfHomogeneous(const bool &value){
      _meshIsHomogeneous = value ;
    }

    /** AMR */
    std::vector < std::map < unsigned,  std::map < unsigned, double  > > >& GetAmrRestrictionMap() {
      return _amrRestriction;
    }
    
    /** AMR */
    std::vector < std::map < unsigned, bool > > & GetAmrSolidMark(){
      return _amrSolidMark;
    }

private:
    
    /** AMR */
    bool _meshIsHomogeneous;
    
    /** AMR: restriction map (vector of 3 FE families: linear, quadratic, biquadratic) */
    std::vector < std::map < unsigned,  std::map < unsigned, double  > > > _amrRestriction;
    /** AMR: solid mark map (vector of 3 FE families: linear, quadratic, biquadratic) */
    std::vector < std::map < unsigned, bool > > _amrSolidMark;


    
};

} //end namespace femus



#endif
