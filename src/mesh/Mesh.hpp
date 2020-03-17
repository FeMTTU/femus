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
#include <assert.h>

#include "vector"
#include "map"



namespace femus {



using std::vector;
class Solution;

class elem;

/**
 * The mesh class
*/

class Mesh : public ParallelObject {

public:

    /** Constructor */
    explicit
    Mesh();

    /** destructor */
    ~Mesh();

    /** Print the mesh info for this level */
    void PrintInfo();

    /** Get the dof number for the element -type- */
    unsigned GetTotalNumberOfDofs(const unsigned &type) const {
      return _dofOffset[type][_nprocs];
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
    };

    /** Get the number of element */
    unsigned GetNumberOfElements() const {
      return _nelem;
    }

    elem * GetElementArray() const {
      return el;
    }
    
    /** Get if element is refined*/
    short unsigned GetRefinedElementIndex(const unsigned &iel) const;

    /** Get element group*/
    short unsigned GetElementGroup(const unsigned &iel) const;
    /** Get element material*/
    short unsigned GetElementMaterial(const unsigned &iel) const;
    /** Get element type*/
    short unsigned GetElementType(const unsigned &iel) const;

    /** Only for parallel */
    bool GetSolidMark(const unsigned &inode) const;

    /** Only for parallel */
    unsigned GetElementDofNumber(const unsigned &iel, const unsigned &type) const;

    /** Only for parallel */
    const unsigned GetElementFaceType(const unsigned &kel, const unsigned &jface) const;

    /** Only for parallel */
    unsigned GetLocalFaceVertexIndex(const unsigned &iel, const unsigned &iface, const unsigned &jnode) const;


    /** Only for parallel */
    unsigned GetElementFaceDofNumber(const unsigned &iel, const unsigned jface, const unsigned &type) const;

    /** Only for parallel */
    unsigned GetElementFaceNumber(const unsigned &iel, const unsigned &type=1) const;

    void GetElementNodeCoordinates(std::vector < std::vector <double > > &xv, const unsigned &iel, const unsigned &solType = 2);

    /** Set the grid number */
    void SetLevel(const unsigned &i) {
        _level=i;
    };

    /** Get the grid number */
    unsigned GetLevel() const {
      return _level;
    }

    /** Set the dimension of the problem (1D, 2D, 3D) */
    void SetDimension(const unsigned &dim) {
      Mesh::_dimension = dim;
      Mesh::_ref_index = pow(2,Mesh::_dimension);  // 8*DIM[2]+4*DIM[1]+2*DIM[0];
      Mesh::_face_index = pow(2,Mesh::_dimension-1u);
    }


    /** Get the dimension of the problem (1D, 2D, 3D) */
    const unsigned GetDimension() const {
      return Mesh::_dimension;
    }

    /** To be added*/
    const unsigned GetRefIndex() const {
      return Mesh::_ref_index;
    }

    unsigned GetSolutionDof(const unsigned &i, const unsigned &iel, const short unsigned &solType) const;

    unsigned GetSolutionDof(const unsigned &i0,const unsigned &i1, const unsigned &ielc, const short unsigned &solType, const Mesh* mshc) const ;

    /** Performs a bisection search to find the processor of the given dof */
    unsigned IsdomBisectionSearch(const unsigned &dof, const short unsigned &solType) const;

    /** To be added */
    const unsigned GetFaceIndex() const {
      return Mesh::_face_index;
    }

    /** Allocate memory for adding fluid or solid mark */
    void AllocateAndMarkStructureNode();

    /** To be Added */
    void SetFiniteElementPtr(const elem_type* otheFiniteElement[6][5]);

    void Partition();
    
    void InitializeTopologyStructures();
  
    /** Only file reading */
    void ReadCoarseMeshFile (const std::string& name, const double Lref, std::vector<bool>& type_elem_flag, const bool read_groups, const bool read_boundary_groups);

      /** This function generates the coarse mesh level, $l_0$, from an input mesh file */
    void ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool> &_finiteElement_flag);

    /** This function generates the coarse mesh level, $l_0$, from an input mesh file, with option to not read groups */
    void ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool> &_finiteElement_flag, const bool read_groups, const bool read_boundary_groups);

    /** This function generates a coarse box mesh */
    void GenerateCoarseBoxMesh(const unsigned int nx,
                               const unsigned int ny,
                               const unsigned int nz,
                               const double xmin, const double xmax,
                               const double ymin, const double ymax,
                               const double zmin, const double zmax,
                               const ElemType type, std::vector<bool> &type_elem_flag);

    /** To be added */
    void FillISvector(vector < unsigned > &partition);

    /** To be added */
    void Buildkel();
    
    void BiquadraticNodesNotInGambit();
    
    std::vector < std::map < unsigned,  std::map < unsigned, double  > > >& GetAmrRestrictionMap(){
      return _amrRestriction;
    }
    
    std::vector < std::map < unsigned, bool > > & GetAmrSolidMark(){
      return _amrSolidMark;
    }

    basis *GetBasis(const short unsigned &ielType, const short unsigned &solType);
    
    // member data
    Solution* _topology;
    const elem_type *_finiteElement[6][5];

    vector < unsigned > _elementOffset;
    vector < unsigned > _ownSize[5];
    vector < unsigned > _dofOffset[5];
    vector< vector < int > > _ghostDofs[5];

    elem *el;  // topology object - list of all elements
    static bool (* _SetRefinementFlag)(const std::vector < double >& x,
                                       const int &ElemGroupNumber,const int &level);
    static bool _IsUserRefinementFunctionDefined;
    std::map<unsigned int, std::string> _boundaryinfo;

    /** Get the projection matrix between Lagrange FEM at the same level mesh*/
    SparseMatrix* GetQitoQjProjection(const unsigned& itype, const unsigned& jtype);

    /** Get the coarse to the fine projection matrix and use it to restrict only on coarse nodes i.e. projection*/
    SparseMatrix* GetCoarseToFineProjectionRestrictionOnCoarse(const unsigned& solType);

    /** Get the coarse to the fine projection matrix*/
    SparseMatrix* GetCoarseToFineProjection(const unsigned& solType);

    /** Set the coarser mesh from which this mesh is generated */
    void SetCoarseMesh( Mesh* otherCoarseMsh ){
      _coarseMsh = otherCoarseMsh;
    };

    bool GetIfHomogeneous(){
      return _meshIsHomogeneous;
    }

    void SetIfHomogeneous(const bool &value){
      _meshIsHomogeneous = value ;
    }

    void SetCharacteristicLength(const double & cLength){
      _cLenght = cLength;
    }
    
    double GetCharacteristicLength(){
      return _cLenght;
    };
    
    const unsigned GetXIndex()          const { return _xIndex; };
    const unsigned GetYIndex()          const { return _yIndex; };
    const unsigned GetZIndex()          const { return _zIndex; };
    const unsigned GetAmrIndex()        const { return _amrIndex; };
    const unsigned GetSolidMarkIndex()  const { return _solidMarkIndex; };

private:
    /** Coarser mesh from which this mesh is generated, it equals NULL if _level = 0 */
    Mesh* _coarseMsh;

    /** The projection matrix between Lagrange FEM at the same level mesh */
    SparseMatrix* _ProjQitoQj[3][3];

    /** The coarse to the fine projection matrix */
    SparseMatrix* _ProjCoarseToFine[5];

    /** Build the projection matrix between Lagrange FEM at the same level mesh*/
    void BuildQitoQjProjection(const unsigned& itype, const unsigned& jtype);

    /** Build the coarse to the fine projection matrix */
    void BuildCoarseToFineProjection(const unsigned& solType, const char el_dofs[]);

    /** Weights used to build the baricentric coordinate **/
    static const double _baricentricWeight[6][5][18];
    static const unsigned _numberOfMissedBiquadraticNodes[6];
    
    //member-data
    int _nelem;                                //< number of elements
    unsigned _nnodes;                          //< number of nodes
    unsigned _level;                           //< level of mesh in the multilevel hierarchy
    static unsigned _dimension;                //< dimension of the problem
    static unsigned _ref_index;
    static unsigned _face_index;
    
    std::map < unsigned, unsigned > _ownedGhostMap[2];
    vector < unsigned > _originalOwnSize[2];

    static const unsigned _END_IND[5];
    vector < vector < double > > _coords;

    bool _meshIsHomogeneous;
    // indices of the topology parallel vectors
    static const unsigned _xIndex = 0;
    static const unsigned _yIndex = 1;
    static const unsigned _zIndex = 2;
    static const unsigned _amrIndex = 3;
    static const unsigned _solidMarkIndex = 4;
    
    std::vector < std::map < unsigned,  std::map < unsigned, double  > > > _amrRestriction;
    std::vector < std::map < unsigned, bool > > _amrSolidMark;

    double _cLenght;
};

} //end namespace femus



#endif
