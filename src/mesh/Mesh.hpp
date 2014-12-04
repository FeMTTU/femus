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

#ifndef __mesh_hpp__
#define __mesh_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Elem.hpp"
#include "vector"
#include "map"
#include "metis.h"
#include "Solution.hpp"
#include "ElemType.hpp"
#include "ElemTypeEnum.hpp"
#include "ParallelObject.hpp"


namespace femus {



using std::vector;
class Solution;

/**
 * The mesh class
*/

class Mesh : public ParallelObject {

public:

    /** Constructor */
    explicit
    Mesh() {};

    /** destructor */
    ~Mesh();

    /** Generate mesh functions */
    
    /** This function generates the coarse mesh level, $l_0$, from an input mesh file */
    void ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool> &_finiteElement_flag);
    
    /** This function generates a coarse box mesh */
    void GenerateCoarseBoxMesh(const unsigned int nx,
                                  const unsigned int ny,
                                  const unsigned int nz,
                                  const double xmin, const double xmax,
                                  const double ymin, const double ymax,
                                  const double zmin, const double zmax,
                                  const ElemType type, std::vector<bool> &type_elem_flag);

    /** Print the mesh info for this level */
    void PrintInfo();

    /** Get the dof number for the element -type- */
    unsigned GetDofNumber(const unsigned type) const;
    
    /** Set the number of nodes */
    void SetNumberOfNodes(const unsigned nnodes) {
      _nnodes = nnodes; 
    };
    
    /** Get the number of nodes */
    unsigned GetNumberOfNodes() const {
      return _nnodes;
    }
    
    /** Set the number of element */
    void SetElementNumber(const unsigned nelem) {
      _nelem = nelem; 
    };

    /** Get the number of element */
    unsigned GetElementNumber() const {
      return _nelem;
    }

    /** Set the grid number */
    void SetGridNumber(const unsigned i) {
        _grid=i;
    };

    /** Get the grid number */
    unsigned GetGridNumber() const {
      return _grid;
    }

    /** Set the dimension of the problem (1D, 2D, 3D) */
    void SetDimension(const unsigned dim) {
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

    /** Get the metis dof from the gambit dof */
    const unsigned GetMetisDof(unsigned inode, short unsigned SolType) const {
      return IS_Gmt2Mts_dof[SolType][inode];
    }

    /** To be added */
    const unsigned GetEndIndex(const unsigned i) const {
      return _END_IND[i];
    }
    
    /** Allocate memory for adding fluid or solid mark */
    void AllocateAndMarkStructureNode();
    
    /** Get the material of the kel element */
    const unsigned GetElementMaterial(unsigned &kel) const;
    
    /** To be Added */
    void SetFiniteElementPtr(const elem_type* otheFiniteElement[6][5]);
    
    
    /** Refinement functions */
    
    /** This function generates a finer mesh level, $l_i$, from a coarser mesh level $l_{i-1}$, $i>0$ */
    void RefineMesh(const unsigned &igrid, Mesh *mshc, const elem_type* otheFiniteElement[6][5]);
    
    /** Flag all the elements to be refined */
    void FlagAllElementsToBeRefined();

    /** Flag all the even elements to be refined */
    void FlagOnlyEvenElementsToBeRefined();

    /** Flag the elements to be refined in according to a user-defined function */
    void FlagElementsToBeRefinedByUserDefinedFunction();

    /** Flag the elements to be refined in according to AMR criteria */
    void FlagElementsToBeRefinedByAMR();
    
    /** Partition Functions */
    
    /** To be added */
    void GenerateVankaPartitions_FAST( const unsigned &block_size, vector < vector< unsigned > > &blk_elements,
				       vector <unsigned> &block_element_type);
    
    /** To be added */
    void GenerateVankaPartitions_FSI( const unsigned &block_size, vector < vector< unsigned > > &block_elements,
				      vector <unsigned> &block_element_type);
    
    /** To be added */
    void GenerateVankaPartitions_FSI1( const unsigned *block_size, vector < vector< unsigned > > &block_elements,
					vector <unsigned> &block_type_range);
    
    /** To be added */    
    void GenerateVankaPartitions_METIS( const unsigned &block_size, vector < vector< unsigned > > &blk_elements);
 
    
    // member data
    Solution* _coordinate;
    const elem_type *_finiteElement[6][5];
    vector <unsigned> IS_Mts2Gmt_elem_offset;
    vector <unsigned> IS_Mts2Gmt_elem;
    vector <unsigned> own_size[5];
    vector <vector <unsigned> > MetisOffset;
    vector< vector<unsigned> > ghost_nd[5];
    vector< vector<int> > ghost_nd_mts[5];
    vector <unsigned> ghost_size[5];
    elem *el;  //< elements
    idx_t *epart;
    idx_t *npart;
    idx_t nsubdom;
    static bool (* _SetRefinementFlag)(const double &x, const double &y, const double &z,
                                       const int &ElemGroupNumber,const int &level);
    static bool _TestSetRefinementFlag;
    std::map<unsigned int, std::string> _boundaryinfo;
    

private:
  
    /** To be added */
    void copy_elr(vector <unsigned> &other_vec) const;
  
    /** To be added */
    void BuildAdjVtx();

    /** To be added */
    void Buildkmid();

    /** To be added */
    void Buildkel();
    
    /** Renumber nodes in the following order: vertices, face, center */
    void RenumberNodes(vector < vector < double> > &coords);
    
    /** Partition functions */
  
    /** Partition the mesh using the METIS partitioner */
    void GenerateMetisMeshPartition();

    
    //member-data
    int _nelem;                                   //< number of elements
    unsigned _nnodes;                              //< number of nodes
    unsigned _grid;                            //< level of mesh in the multilevel hierarchy
    static unsigned _dimension;                //< dimension of the problem
    static unsigned _ref_index;
    static unsigned _face_index;
    vector <unsigned> IS_Gmt2Mts_dof[5];        //< dof map
    vector <unsigned> IS_Gmt2Mts_dof_offset[5]; //< map offset
    static const unsigned _END_IND[5];

};

} //end namespace femus



#endif
