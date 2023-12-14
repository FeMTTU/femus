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
#include "ParallelObject.hpp"
#include "GeomElemBase.hpp"
#include "Elem.hpp"

#include "ElemTypeEnum.hpp"

#include "ElemType.hpp"
#include "FElemTypeEnum_list.hpp"

#include "fe_prolongation_matrices.hpp"

#include "Solution.hpp"

#include "MeshGeneration.hpp"



#include <cassert>
#include <vector>
#include <map>



namespace femus {


class elem;
class elem_type;

class Solution;

class NumericVector;
class SparseMatrix;


/**
 * The mesh class
*/

class Mesh : public ParallelObject {

// === Friend functions and classes - BEGIN =================
/// These classes can access stuff that otherwise is protected/private
friend   void MeshTools::Generation::BuildBox ( Mesh& mesh,
		      std::vector < std::vector < double> > &coords,
                      const unsigned int nx,
                      const unsigned int ny,
                      const unsigned int nz,
                      const double xmin, const double xmax,
                      const double ymin, const double ymax,
                      const double zmin, const double zmax,
                      const ElemType type,
                      std::vector<bool> &type_elem_flag );

friend class MED_IO;
friend class GambitIO;

friend class MeshPartitioning;
friend class MeshMetisPartitioning;

friend class MeshRefinement;

friend class MultiLevelMesh;
// === Friend functions and classes - END =================


// === Constructors / Destructor - BEGIN =================

public:

    /** Constructor */
    explicit
    Mesh();

    /** destructor */
    ~Mesh();

    
// === Constructors / Destructor - END =================

    

// === Geometric Element, Single - BEGIN =================
 public:
   
    /** */
    unsigned GetElementFaceNumber(const unsigned &iel, const unsigned &type = GeomElemBase::_index_for_all_faces) const;

    unsigned GetElementFaceNumber_PassElemType(const short unsigned & el_type, const unsigned& type = GeomElemBase::_index_for_all_faces) const;
    
    /**  */
    unsigned GetLocalFaceVertexIndex(const unsigned &iel, const unsigned &iface, const unsigned &jnode) const;

    /**  */
    unsigned GetLocalFaceVertexIndex_PassElemType(const short unsigned & el_type, const unsigned& iface, const unsigned& jnode) const;
    
    
// === Geometric Element, Single - END =================

    
// === Geometric Element, Single, REFINEMENT - BEGIN =================
 public:
    
    /** It would private if we put Solution as another friend of Mesh, but maybe it is too much for now */
    const unsigned GetRefIndex() const {
      return _ref_index;
    }

    
 private:
   
    /** MESH */
    const unsigned GetRefFaceIndex() const {
      return _ref_face_index;
    }
    
     void SetRefinementCellAndFaceIndices(const unsigned &dim) {
      _ref_index  = pow(2, dim);     //8 elements from refining 1 HEX, TET, WEDGE; 4 elements from refining 1 QUAD, TRI; 2 elements from refining 1 LINE
      _ref_face_index = pow(2, dim -1u);
    }

    /** MESH, REF: 8 elements from refining 1 HEX, TET, WEDGE; 4 elements from refining 1 QUAD TRI; 2 elements from refining 1 LINE // 8*DIM[2]+4*DIM[1]+2*DIM[0]; */
    unsigned _ref_index;
    
    /** MESH, REF: 4 faces from refining 1 QUAD TRI; 2 faces from refining 1 LINE; 1 face from refining 1 point // 4*DIM[2]+2*DIM[1]+1*DIM[0]; */
    unsigned _ref_face_index;

    
// === Geometric Element, Single, REFINEMENT - END =================
    


// === Mesh, BASIC, Debug - BEGIN =================

public:
    
    /** Print the mesh info */
    void PrintInfo() const;

private:
    
    void PrintInfoElements() const;
// === Mesh, BASIC, Debug - END =================


// === Mesh, BASIC, Dimension - BEGIN =================
 public:
   
    /** MESH: Get the dimension of the problem (1D, 2D, 3D) */
    const unsigned GetDimension() const {
      return Mesh::_dimension;
    }

 private:
   
    /** MESH: Set the dimension of the problem (1D, 2D, 3D) */
    void SetDimension(const unsigned &dim) {
      Mesh::_dimension = dim;
    }

    /** MESH: dimension of the problem @todo I would make it not static */
    static unsigned _dimension;
// === Mesh, BASIC, Dimension - END =================


// === Mesh, BASIC, CharacteristicLength - BEGIN =================

public:
    
    /// @todo This soon will be private
    void ComputeCharacteristicLength();
    
private:

    double GetCharacteristicLength() const {
      return _cLength;
    };

    void SetCharacteristicLength(const double & cLength){
      _cLength = cLength;
    }
    
    /** Order of the domain size */
    double _cLength;

    
// === Mesh, BASIC, CharacteristicLength - END =================



// === Mesh, Elements - BEGIN ====================================================

    
// === Elements, List - BEGIN =================
 public:
   
    /** MESH */
    inline elem * GetMeshElements() const {
      return el;
    }
    
    /** MESH: list of all elements @todo should be private */
    elem *el;
   
// === Elements, List - END =================

    
// === Elements, Number - BEGIN =================
 public:
   
    /** Get the number of element */
    unsigned GetNumberOfElements() const {
      return _nelem;
    }
    
 private:
    
    /** Set the number of element */
    void SetNumberOfElements(const unsigned &nelem) {
      _nelem = nelem;
    }

    /** MESH: number of elements */
    int _nelem;
// === Elements, Number - END =================


// === Elements, Type - BEGIN =================
 public:
    
    /** Get element type*/
    short unsigned GetElementType(const unsigned &iel) const;
    /**  */
    const unsigned GetElementFaceType(const unsigned &kel, const unsigned &jface) const;
    
// === Elements, Type - END =================

// === Elements, Subdomains - BEGIN =================
 public:
    
    /** MESH: Number of elements per processor (incremental count)  @todo should be private */
    std::vector < unsigned > _elementOffset;
 
// === Elements, Subdomains - END =================
   
   
// === Elements, Group - BEGIN =================
 public:
    
    /** Get element group*/
    short unsigned GetElementGroup(const unsigned &iel) const;
    
// === Elements, Group - END =================

   
// === Elements, Material - BEGIN =================
 public:
   
    /** Get element material @todo this should be in a separate FSI environment */
    short unsigned GetElementMaterial(const unsigned &iel) const;
    
// === Elements, Material - END =================
    
// === Mesh, Elements - END ====================================================



// === Mesh, Nodes - BEGIN ====================================================

    
// === Nodes, number - BEGIN =================
 public:
   
    /** Get the number of nodes */
    unsigned GetNumberOfNodes() const {
      return _nnodes;
    }

 private:
   
    /** Set the number of nodes */
    void SetNumberOfNodes(const unsigned &nnodes) {
      _nnodes = nnodes;
    };

    /** MESH: number of nodes */
    unsigned _nnodes;
    
// === Nodes, number - END =================
    
    
    
    
// === Nodes, coordinates read from coarse file (temporary, then use GetTopology()) - BEGIN =================

 private:
   
    /** MESH: node coordinates for each space dimension  @todo beware: this is only filled at coarse reading, then use GetTopology() for the coordinates! */
    std::vector < std::vector < double > > _coords;
// === Nodes, coordinates read from coarse file (temporary, then use GetTopology()) - END =================
    
// === Mesh, Nodes - END ====================================================
    
    
    
    
// === Mesh, COARSE GENERATION - BEGIN =================
    
// === COARSE MESH GENERATION, from FILEs - BEGIN =================

 public:

      /** This function generates the coarse mesh level, $l_0$, from an input mesh file */
    void ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool> &_finiteElement_flag);

    /** This function generates the coarse mesh level, $l_0$, from an input mesh file, with option to not read groups */
    void ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool> &_finiteElement_flag, const bool read_groups, const bool read_boundary_groups);

private:
  
    /** Only file reading */
    void ReadCoarseMeshFile (const std::string& name, const double Lref, std::vector<bool>& type_elem_flag, const bool read_groups, const bool read_boundary_groups);

    void ReadCoarseMeshBeforePartitioning(const std::string& name, const double Lref, std::vector<bool>& type_elem_flag, const bool read_groups, const bool read_boundary_groups);
  


// === COARSE MESH GENERATION, from FILEs - END =================


// === COARSE MESH GENERATION, from function - BEGIN =================
 public:
   
    /** This function generates a coarse box mesh */
    void GenerateCoarseBoxMesh(const unsigned int nx,
                               const unsigned int ny,
                               const unsigned int nz,
                               const double xmin, const double xmax,
                               const double ymin, const double ymax,
                               const double zmin, const double zmax,
                               const ElemType type, 
                               std::vector<bool> &type_elem_flag);

 public:
   
  std::map<unsigned int, std::string> GetBoundaryInfo() const {   
     return _boundaryinfo;
   }

 private:
   
  void SetBoundaryInfo(const unsigned int flag_int, const std::string flag_string) {   
    _boundaryinfo.insert(std::pair<unsigned int, std::string>(flag_int, flag_string) );
   }
   
    /** Boundary names for faces, this is only filled in the above function so far */
    std::map<unsigned int, std::string> _boundaryinfo;
    
// === COARSE MESH GENERATION, from function - END =================


// === COARSE MESH GENERATION, for all - BEGIN =================
    
 private:
    
    void AddBiquadraticNodesNotInMeshFile();
    
    /** Weights used to build the baricentric coordinate to compute the missing biquadratic nodes **/
    static const double _baricentricWeight[N_GEOM_ELS][NFE_FAMS][MAXIMUM_NUMBER_OF_NON_BIQUADRATIC_NODES_TO_USE_IN_WEIGHTING];
    
    static const unsigned _numberOfMissedBiquadraticNodes[N_GEOM_ELS];
    
// === COARSE MESH GENERATION, for all - END =================

// === Mesh, COARSE GENERATION - END =================

    

// === Mesh, Level, Current (when the Mesh is part of a multilevel hierarchy) - BEGIN  =================

public:

    /** MESH: Get the grid number */
    unsigned GetLevel() const {
      return _level;
    }

    
 private:
    
    void PrintInfoLevel() const;
    
    /** MESH: Set the grid number */
    void SetLevel(const unsigned &i) {
        _level = i;
    };

    /** MESH: level of mesh in the multi-level hierarchy */
    unsigned _level;
    
// === Mesh, Level, Current (when the Mesh is part of a multilevel hierarchy) - END  =================

    

// === FE stuff - BEGIN ====================================================

    
 
// === Geometric Element, FE, Single (FE for single geometric element) - BEGIN =================
 public:
    
    /** Also DOFMAP: Only for parallel */
    unsigned GetElementDofNumber(const unsigned &iel, const unsigned &type) const;
    
    /** Also DOFMAP: Only for parallel */
    unsigned GetElementFaceDofNumber(const unsigned &iel, const unsigned jface, const unsigned &type) const;


    
    const elem_type * GetFiniteElement(const unsigned geom_elem, const unsigned fe_soltype) const {
        return _finiteElement[geom_elem][fe_soltype];
    }
    
    /** To be Added */
    void SetFiniteElementPtr(/*const*/ elem_type* otheFiniteElement[N_GEOM_ELS][NFE_FAMS]);
    
    
    /** FE: Finite Element families, for each Geometric Element @todo this one day should be private */
    const elem_type *_finiteElement[N_GEOM_ELS][NFE_FAMS];

// === Geometric Element, FE, Single (FE for single geometric element) - END =================
    

    

// === PARTITIONING, and FE DOFMAP (from Mesh to Unknowns) - BEGIN =================

private:
  
  void PartitionElements_and_FillDofMapAllFEFamilies();
    
  std::vector < unsigned > PartitionForElements() const;
    
  std::vector < unsigned > PartitionForElements_refinement(const bool AMR, const Mesh* mshc) const;

// === PARTITIONING, and FE DOFMAP (from Mesh to Unknowns) - END =================



// === FE DOFMAP - BEGIN =================
// =========================
public:
    
    /** FE: DofMap carriers */
    
    /** Elements - BEGIN */
    void initialize_elem_offsets();
    
    void build_elem_offsets(const std::vector <unsigned> & partition);
    
    void mesh_reorder_elem_quantities();
    
    void set_elem_counts_per_subdomain();
    
    std::vector < unsigned >  elem_offsets();
    /** Elements -  END */
  
    
    /**  Nodes - BEGIN */
    std::vector < unsigned >  node_offsets();
  
    void mesh_reorder_node_quantities(const std::vector <unsigned> & mapping);
    
    void set_node_counts();

    void deallocate_node_mapping(std::vector < unsigned > & node_mapping) const;
    /**  Nodes - END */
    
    
    /** Get the dof number for the element -type- */
    unsigned GetTotalNumberOfDofs(const unsigned &type) const {
      return _dofOffset[type][_nprocs];
    }

    unsigned GetSolutionDof(const unsigned &i, const unsigned &iel, const short unsigned &solType) const;

    unsigned GetSolutionDof(const unsigned &i0,const unsigned &i1, const unsigned &ielc, const short unsigned &solType, const Mesh* mshc) const;

    /** Performs a bisection search to find the processor of the given dof */
    unsigned BisectionSearch_find_processor_of_dof(const unsigned &dof, const short unsigned &solType) const;

    /** FE: DofMap: Here is where the element and node global orderings are changed based on the partitioning */
    void FillISvectorDofMapAllFEFamilies(std::vector < unsigned > &partition);

    void FillISvectorElemOffsets(std::vector < unsigned >& partition);
  
    void FillISvectorNodeOffsets();
    
    
    void dofmap_all_fe_families_initialize();
    
    void dofmap_all_fe_families_clear_ghost_dof_list_for_other_procs();
    
    void dofmap_Element_based_dof_offsets_build();

    std::vector <unsigned> dofmap_Node_based_dof_offsets_Compute_Node_mapping_and_Node_ownSize();
    
    void dofmap_Node_based_dof_offsets_Continue_biquadratic();
    
    void dofmap_Node_based_dof_offsets_Ghost_nodes_search_Complete_biquadratic();
    
    void dofmap_Node_based_dof_offsets_Complete_linear_quadratic();
    
    inline const std::vector < unsigned > * dofmap_get_dof_offset_array() const { return  _dofOffset; }
    
    unsigned  dofmap_get_dof_offset(const unsigned soltype, const unsigned proc_id) const { return  _dofOffset[soltype][proc_id]; }
    
    unsigned  dofmap_get_own_size(const unsigned soltype, const unsigned proc_id) const { return  _ownSize[soltype][proc_id]; }
    
    std::vector < int > dofmap_get_ghost_dofs(const unsigned soltype, const unsigned proc_id) const { return  _ghostDofs[soltype][proc_id]; }

    /** FE: DofMap: Number of dofs per FE family and per processor (incremental count) @todo this should be private - long task */
    std::vector < unsigned > _dofOffset[NFE_FAMS];
    
private:

    
    /** FE: DofMap: Number of owned dofs per FE family and per processor (count, non-incremental) */
    std::vector < unsigned > _ownSize[NFE_FAMS];
    
    /** FE: DofMap: Number of ghost dofs per FE family and per processor (count, non-incremental) */
    std::vector< std::vector < int > > _ghostDofs[NFE_FAMS];
    
    /** FE: DofMap  k = 0, 1 */
    std::map < unsigned, unsigned > _ownedGhostMap[2];
    
    /** FE: DofMap  k = 0, 1 */ 
    std::vector < unsigned > _originalOwnSize[2];
    
    /** print node-based dofOffset counts */
    void PrintInfoNodes() const;

// === FE DOFMAP - END =================
    
    
// === FE DOFMAP, TOPOLOGY - BEGIN =================
public:
  
  Solution * GetTopologyToModify() { return _topology; }
  
  const Solution * GetTopology() const { return _topology; }
  
    void dofmap_build_all_fe_families_and_elem_and_node_structures();
    
    void BuildElementAndNodeStructures();

private:
  
    /** MESH: Coordinates and other stuff */
    Solution* _topology;

// === FE DOFMAP, TOPOLOGY - END =================



// === FE DOFMAP, TOPOLOGY: Coordinates - BEGIN =================

private:
    /** MESH: Topology */
    const unsigned GetXIndex()          const { return _xIndex; }
    const unsigned GetYIndex()          const { return _yIndex; }
    const unsigned GetZIndex()          const { return _zIndex; }
    
    const std::string Get_X_Name() const { return _x_name; }
    const std::string Get_Y_Name() const { return _y_name; }
    const std::string Get_Z_Name() const { return _z_name; }
    
    void GetElementNodeCoordinates(std::vector < std::vector <double > > &xv, const unsigned &iel, const unsigned &solType) const;

    void Topology_InitializeCoordinates();
    
    void Topology_FillCoordinates();
    
    // indices of the topology parallel vectors
    static const unsigned _xIndex = 0;
    static const unsigned _yIndex = 1;
    static const unsigned _zIndex = 2;
    static const std::string _x_name;
    static const std::string _y_name;  
    static const std::string _z_name;  
// === FE DOFMAP, TOPOLOGY: Coordinates - END =================


// === FE DOFMAP, TOPOLOGY: Refinement, AMR - BEGIN =================
public:
    
    const unsigned GetAmrIndex()        const { return _amrIndex; }
    const std::string GetAmrIndexName() const { return _amrIndex_name; }

private:
    void Topology_InitializeAMR();

    static const unsigned _amrIndex = 3;
    static const std::string _amrIndex_name;
    
// === FE DOFMAP, TOPOLOGY: Refinement, AMR - END =================

// === FE DOFMAP, TOPOLOGY: SolidMark  - BEGIN =================
public:

    /** Only for parallel @todo this should be in a separate FSI environment */
    bool GetSolidMark(const unsigned &inode) const;
    
    
private:

    const unsigned GetSolidMarkIndex()  const { return _solidMarkIndex; }

    /** FSI:  */
    void Topology_InitializeSolidNodeFlag();
    
    /** FSI: Allocate memory for adding fluid or solid mark @todo this should be in a separate FSI environment */
    void Topology_FillSolidNodeFlag();
    
    static const unsigned _solidMarkIndex = 4;
    static const std::string _solidMark_name;
    
// === FE DOFMAP, TOPOLOGY: SolidMark  - END =================


    
// === FE DOFMAP, REFINEMENT - BEGIN =================
 public:

 FE_Prolongation_Matrices & get_prol_matrices() { return _fe_prol_matrices; }

 private:
  
 FE_Prolongation_Matrices   _fe_prol_matrices;
 
// === FE DOFMAP, REFINEMENT - END =================
    


// === FE DOFMAP, REFINEMENT, AMR - BEGIN =================
// =========================
public:
    
    /** AMR @todo I would make it not static */
    static bool (* _SetRefinementFlag)(const std::vector < double >& x, const int &ElemGroupNumber, const int &level);
    
    /** AMR */
    static void set_is_refinement_function_defined(const bool truth_in) {
      _IsUserRefinementFunctionDefined = truth_in;
    }
    static bool get_is_refinement_function_defined() {
      return _IsUserRefinementFunctionDefined;
    }
    
    /** AMR */
    bool GetIfHomogeneous() const {
      return _meshIsHomogeneous;
    }

    /** AMR */
    const std::vector < std::map < unsigned,  std::map < unsigned, double  > > >& GetAmrRestrictionMap() const {
      return _amrRestriction;
    }
    
    /** Get if element is refined*/
    short unsigned GetRefinedElementIndex(const unsigned &iel) const;
    
private:
    
    void InitializeAndPossiblyFillAmrRestriction(const bool amr);
  
    /** AMR */
    std::vector < std::map < unsigned,  std::map < unsigned, double  > > >& GetAmrRestrictionMap() {
      return _amrRestriction;
    }
    
    /** AMR */
    void SetIfHomogeneous(const bool &value) {
      _meshIsHomogeneous = value ;
    }

    static bool _IsUserRefinementFunctionDefined;
    
    /** AMR */
    bool _meshIsHomogeneous;
    
    /** AMR: restriction map (vector of 3 FE families: linear, quadratic, biquadratic) */
    std::vector < std::map < unsigned,  std::map < unsigned, double  > > > _amrRestriction;
    
    void GetAMRRestrictionAndAMRSolidMark(Mesh *msh, const elem *el_in);
// === FE DOFMAP, REFINEMENT, AMR - END =================
    

// === FE DOFMAP, REFINEMENT, AMR, FSI - BEGIN =================
    
public:

    /** AMR */
    const std::vector < std::map < unsigned, bool > > & GetAmrSolidMark() const {
      return _amrSolidMark;
    }

private:
  
    /** AMR */
    std::vector < std::map < unsigned, bool > > & GetAmrSolidMark() {
      return _amrSolidMark;
    }
    
    /** AMR: solid mark map (vector of 3 FE families: linear, quadratic, biquadratic) */
    std::vector < std::map < unsigned, bool > > _amrSolidMark;

// === FE DOFMAP, REFINEMENT, AMR, FSI - END =================


// === FE stuff - END ====================================================


    
};






} //end namespace femus



#endif
