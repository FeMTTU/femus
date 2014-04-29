/*=========================================================================

 Program: FEMUS
 Module: Mesh
 Authors: Eugenio Aulisa
 
 Copyright (c) FEMTTU
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
#include "metis.h"
#include "Solution.hpp"
#include "ElemType.hpp"
#include "ElemTypeEnum.hpp"


namespace femus {



using std::vector;
class Solution;

/** 
 * The mesh class 
*/

class mesh {
  
public:  
  Solution* _coordinate;
  
  
private:
  unsigned nvt,grid;
  static unsigned _dimension;
  static unsigned _face_index;
  vector <unsigned> IS_Gmt2Mts_dof[5];  // dof map 
  vector <unsigned> IS_Gmt2Mts_dof_offset[5]; // map offset 
  
  /** Local to global map for built-in mesh generator */
  unsigned int idx(const ElemType type, const unsigned int nx, const unsigned int i, const unsigned int j);
  
public:
  static const unsigned _END_IND[5];  
  vector< vector<unsigned> > ghost_nd[5];
  vector< vector<int> > ghost_nd_mts[5];
  vector <unsigned> ghost_size[5];
  int nel;
  static unsigned _ref_index;
  int _nprocs;
  int _iproc;
  idx_t *epart;
  idx_t *npart;
  idx_t nsubdom;
  vector <unsigned> IS_Mts2Gmt_elem_offset;
  vector <unsigned> IS_Mts2Gmt_elem;  
  vector <unsigned> own_size[5];
  vector <vector <unsigned> > MetisOffset;
  
  static bool (* _SetRefinementFlag)(const double &x, const double &y, const double &z, 
				   const int &ElemGroupNumber,const int &level);
      
  //unsigned grid;
  elem *el;  // elements
  
  /** Constructor */
  mesh();
  
  /** destructor */
  ~mesh();

  // Reading 
  
  /**
  *  This function generates the coarse mesh level, $l_0$, from an input mesh file 
  **/
  void ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool> &_type_elem_flag);
  
  void Read1D(const char infile [], vector < vector < double> > &vt);
  
  /** Read the coarse-mesh from a neutral Gambit File */
  void ReadGambit(const std::string& name, vector < vector < double> > &vt,const double Lref,std::vector<bool> &type_elem_flag);
  
  /** Built-in cube-structured mesh generator */ 
  void BuildBrick ( const unsigned int nx,
	            const unsigned int ny,
	            const unsigned int nz,
		    const double xmin, const double xmax,
		    const double ymin, const double ymax,
		    const double zmin, const double zmax,
		    const ElemType type,
		    std::vector<bool> &type_elem_flag );
  
  /**
  *  This function generates a finer mesh level, $l_i$, from a coarser mesh level $l_{i-1}$, $i>0$
  **/
  void RefineMesh(const unsigned &igrid, mesh *mshc, const elem_type* type_elem[6][5]);
  
  
  // partitioning 
  
  void generate_metis_mesh_partition();
    
  
  unsigned GetProcID() const;
  unsigned GetNumProcs() const;
  
  void BuildAdjVtx();
  void Buildkmid();
  void Buildkel();

  void copy_elr(vector <unsigned> &other_vec) const;
  unsigned GetDofNumber(const unsigned type) const;
  unsigned GetElementNumber() const;
  unsigned GetGridNumber() const;
  unsigned SetGridNumber(unsigned i){ grid=i; }; 
  void AllocateAndMarkStructureNode();
  unsigned GetDimension();
  unsigned GetRefIndex();
  
  unsigned GetMetisDof(unsigned inode, short unsigned SolType) const;
  unsigned GetEndIndex(const unsigned i) const;
  
  
};

/*----------------------- functions ----------------------------------*/
//------------------------------------------------------------------------------------------------------
inline unsigned mesh::GetMetisDof(unsigned inode, short unsigned SolType) const{
//   if(SolType<=5)
    return IS_Gmt2Mts_dof[SolType][inode];
//   else 
//     return inode;
}

inline unsigned mesh::GetProcID() const {
  return _iproc;
}

inline unsigned mesh::GetNumProcs() const {
  return _nprocs; 
}


} //end namespace femus



#endif