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

class mesh : public ParallelObject {

public:

    /** Constructor */
    mesh();

    /** destructor */
    ~mesh();

    /** This function generates the coarse mesh level, $l_0$, from an input mesh file */
    void ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool> &_type_elem_flag);

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

    /** This function generates a finer mesh level, $l_i$, from a coarser mesh level $l_{i-1}$, $i>0$ */
    void RefineMesh(const unsigned &igrid, mesh *mshc, const elem_type* type_elem[6][5]);

    /** Partition the mesh using the METIS partitioner */
    void generate_metis_mesh_partition();

    /** Print the mesh info for this level */
    void print_info();

    /** To be added */
    void BuildAdjVtx();

    /** To be added */
    void Buildkmid();

    /** To be added */
    void Buildkel();

    /** To be added */
    void copy_elr(vector <unsigned> &other_vec) const;

    /** Get the dof number for the element -type- */
    unsigned GetDofNumber(const unsigned type) const;

    /** Get the number of element */
    unsigned GetElementNumber() const;

    /** Get the grid number */
    unsigned GetGridNumber() const;

    /** Set the grid number */
    unsigned SetGridNumber(unsigned i) {
        grid=i;
    };

    /** Allocate memory for adding fluid or solid mark */
    void AllocateAndMarkStructureNode();

    /** Get the dimension of the problem (1D, 2D, 3D) */
    unsigned GetDimension();

    /** To be added*/
    unsigned GetRefIndex();

    /** Get the metis dof from the gambit dof */
    unsigned GetMetisDof(unsigned inode, short unsigned SolType) const;

    /** To be added */
    unsigned GetEndIndex(const unsigned i) const;

    /** Flag all the elements to be refined */
    void FlagAllElementsToBeRefined();

    /** Flag all the even elements to be refined */
    void FlagOnlyEvenElementsToBeRefined();

    /** Flag the elements to be refined in according to a user-defined function */
    void FlagElementsToBeRefinedByUserDefinedFunction();

    // member data
    Solution* _coordinate;
    vector <unsigned> IS_Mts2Gmt_elem_offset;
    vector <unsigned> IS_Mts2Gmt_elem;
    vector <unsigned> own_size[5];
    vector <vector <unsigned> > MetisOffset;
    vector< vector<unsigned> > ghost_nd[5];
    vector< vector<int> > ghost_nd_mts[5];
    vector <unsigned> ghost_size[5];
    elem *el;  //< elements
    int nel;
    static unsigned _ref_index;
    idx_t *epart;
    idx_t *npart;
    idx_t nsubdom;
    static const unsigned _END_IND[5];
    static bool (* _SetRefinementFlag)(const double &x, const double &y, const double &z,
                                       const int &ElemGroupNumber,const int &level);
    static bool _TestSetRefinementFlag;
    std::map<unsigned int, std::string> _boundaryinfo;

private:

    //member-data
    unsigned nvt,grid;
    static unsigned _dimension;
    static unsigned _face_index;
    vector <unsigned> IS_Gmt2Mts_dof[5];        //< dof map
    vector <unsigned> IS_Gmt2Mts_dof_offset[5]; //< map offset

    /** 2D local-to-global map for built-in mesh generator */
    unsigned int idx(const ElemType type, const unsigned int nx, const unsigned int i, const unsigned int j);

    /** 3D local-to-global map for built-in mesh generator */
    unsigned int idx(const ElemType type, const unsigned int nx, const unsigned int ny, const unsigned int i,
                     const unsigned int j, const unsigned int k);

protected:
    int _nprocs;
};

/*----------------------- functions ----------------------------------*/
//------------------------------------------------------------------------------------------------------
inline unsigned mesh::GetMetisDof(unsigned inode, short unsigned SolType) const {
    return IS_Gmt2Mts_dof[SolType][inode];
}

} //end namespace femus



#endif
