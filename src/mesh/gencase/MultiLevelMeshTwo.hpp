/*=========================================================================

 Program: FEMUS
 Module: MultiLevelMeshTwo
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_meshGencase_MultiLevelMeshTwo_hpp__
#define __femus_meshGencase_MultiLevelMeshTwo_hpp__

//C++ 
#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "hdf5.h"

#include "Typedefs.hpp"
#include "FemusConfig.hpp"
#include "FemusInputParser.hpp"
#include "ElemSto.hpp"
#include "VBTypeEnum.hpp"
#include "FETypeEnum.hpp"
#include "GeomElTypeEnum.hpp"

namespace femus {


class Files;
class Domain;


  
class MultiLevelMeshTwo  {

public:

//===== Constructors/ Destructor ===========
     MultiLevelMeshTwo () : _dim(0) { };
     MultiLevelMeshTwo (const unsigned nolevels, const unsigned dim, const GeomElType geomel_type, const std::string mesh_file_in);
    void clear ();

    //======= mesh generation functions ====
    void ElemChildToFather();
    void ComputeElemOffsetsBySubdLevel();

    void FillNodeSto();
    void ComputeNodeOffsetsBySubdLevel();
    void ComputeMaxElXNode();
    void ComputeNodeMapExtLevels();

//attributes    ************************************

    const uint _dim;               ///< spatial dimension
// ==== PARALLEL ===
    uint _iproc;       /// current subdomain
    uint _NoSubdom;    /// Number of subdomains (subdomain P)

// ==== MULTIGRID ====
    uint _NoLevels;          /// Number of Levels (L)
   
// ===== ABSTRACT GEOMEL =====
    uint      _type_FEM[VB];          ///  [VB] //just for check
    uint _elnodes[VB][QL];            ///  [VB]
    short unsigned _eltype_flag[VB];  ///  [VB]
    std::vector < std::string >  _geomelem_id;       ///   [DIM]
    std::vector < short unsigned >  _geomelem_flag;  ///   [DIM]
    
// ==== ELEMENTS ====
    uint** _n_elements_vb_lev; ///  [VB][L]      Number of Elements                                   
    int**  _off_el;            ///  [VB][L+P*_NoLevels]  offset of element numbers       
    uint** _el_map;            ///  [VB][L+P*_NoLevels]  Connectivity map     
    int**  _el_bdry_to_vol;
// ELEMENTS =============
    int ** _n_elements_sl_vb;        ///  [VB][]
    int    _n_elements_sum_levs[VB]; ///  [VB]   //of the WHOLE REFINEMENT! //LMFILLS 
    int ** _el_child_to_fath;        /// [LEV][IELEM]      //for every level, it gives you the father
    ElemStoVol**  _el_sto;                 //FILLED ACCORDING TO "how the Libmesh mesh iterator runs" which may not be id in general i think...
    ElemStoBdry** _el_sto_b;               //FILLED with OUR ORDERING, "as we find them during the volume elem loop"

// ===== NODES ======
    double* _xyz;                  //ONLY VV  ///< node coordinates
    uint*   _NoNodesXLev;          //ONLY VV  /// [LEV] ///VOLUME, at the EXTENDED LEVELS
    uint**  _Qnode_lev_Qnode_fine; //ONLY VV   ///ONLY USED FOR THE EQUATION /// FROM DOF TO NODE: Q or L dofs by using an add aux level [_NoLevels+1][NoNodes[L]]  //NO PROC info, these are DOFS in the "SERIALIZED NUMBERING"
    int**  _Qnode_fine_Qnode_lev;  //ONLY VV   ///Un nodo fine puo' o non appartenere ad un certo livello... quindi devo continuare a mantenere i -1, altrimenti dovrei fare un search che fa perdere tempo
    int**   _off_nd;               //ONLY VV  [QL][L+P*_NoLevels]  ///ONLY USED FOR THE EQUATION /// FROM DOF TO NODE: offsets for dofs      --->   // on the fine node numbering, the nodes corresponding to linear dofs are numbered FIRST... or not?                  
                                   
// NODES ===============
    int    _n_nodes;              //ONLY VV //of the WHOLE REFINEMENT! i.e. the FINE ones! //LMFILLS 
    int ** _n_nodes_sl_ql;        //ONLY VV //[QL]
    int *  _elxnode;              //ONLY VV //number of elements per node of the fine mesh
    double *    _nd_coords_libm;  //ONLY VV //node coordinates  //FILLED ACCORDING TO LIBMESH NODE ID ORDERING; then I'll print them according to my FEMUS ordering
    NodeSto**   _nd_sto;          //ONLY VV //FILLED ACCORDING TO LIBMESH NODE ID ORDERING
    int    _maxelxnode;
    
   public:    
// ====== DOMAIN SHAPE ( optional => pointer)
    Domain* _domain;
    Domain* GetDomain() const;
    void    SetDomain(Domain* );
    
    //get functions
    inline const double get_Lref() const {return _Lref;}
    inline const uint   get_dim()  const {return _dim;}

    //set functions
    inline void SetLref(const double lref_in) { _Lref = lref_in; }

// FROM LIBMESH TO FEMUS ===============
    std::vector< std::pair<int,int> > _nd_fm_libm; //from FINE FEMUS NODE ORDERING to FINE LIBMESH NODE ORDERING
    std::vector< std::pair<int,int> > _el_fm_libm; //because the EQUATION needs it for the SPARSITY PATTERN
    std::vector< std::pair<int,int> > _el_fm_libm_b;
    std::vector<int>                  _nd_libm_fm; //from FINE LIBMESH NODE ORDERING to FINE FEMUS NODE ORDERING  //TODO this is the one that is not correctly filled in debug mode
    std::vector<int>                  _el_libm_fm;  //TODO in dbg mode they do not survive in GENCASE... extremely weird thing... let us try to convert from int* to std::vector<int>
    
  protected:
    
    std::string _mesh_file;    //mesh file name from the mesh generator
 
   private:   
     
//attributes    ************************************
    double _Lref;          ///Reference length for non-dimensionalization
     
 };




} //end namespace femus



#endif
