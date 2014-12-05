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
    void SetNumberOfElements(const unsigned nelem) {
      _nelem = nelem; 
    };

    /** Get the number of element */
    unsigned GetNumberOfElements() const {
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



} //end namespace femus



#endif
