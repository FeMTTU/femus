#ifndef __mesh123D_h
#define __mesh123D_h

//C++ 
#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "hdf5.h"

#include "Typedefs.hpp"
#include "FEMTTUConfig.h"
#include "RunTimeMap.hpp"
#include "GeomEl.hpp"
#include "ElemSto.hpp"


namespace femus {


class Files;
class Domain;


  
class Mesh  {

public:
  
    Files& _files; 
    RunTimeMap<double>  _mesh_rtmap;
    uint _dim;               ///< spatial dimension
    uint _meshVB;          /// Number of FEM "manifold" Families (VB)  F

// ===== ABSTRACT GEOMEL(S) =====
    GeomEl    _GeomEl;
    const uint _n_GeomEl;         //number of geom elem types
    uint*      _type_FEM;         //just for check
    uint _elnodes[VB][QL];

// ====== DOMAIN SHAPE (TODO optional => pointer) ----- //if I put it as reference I'd have to initialize it
    Domain* _domain;      //TODO You must remember to ALLOCATE this POINTER BEFORE USING IT!
    Domain* GetDomain() const;
    void    SetDomain(Domain* );
    void TransformElemNodesToRef(const uint vb,const double* xx_qnds,double* refbox_xyz) const;

// ==== PARALLEL ===
    uint _iproc;       /// current subdomain
    uint _NoSubdom;    /// Number of subdomains (subdomain P)

// ==== MULTIGRID ====
    uint _NoLevels;          /// Number of Levels (L)

// ==== ELEMENTS ====
    uint** _n_elements_vb_lev;        /// Number of Elements        [VB][L]                                     
    int**  _off_el;            /// offset of element numbers [VB][L+P*_NoLevels]         
    uint** _el_map;            /// Connectivity map         [VB][L+P*_NoLevels]
    int**  _el_bdry_to_vol;

// ===== NODES ======
    uint*   _NoNodesXLev;            ///VOLUME, at the EXTENDED LEVELS
    double* _xyz;                    ///< node coordinates
    uint**  _Qnode_lev_Qnode_fine;   ///ONLY USED FOR THE EQUATION /// FROM DOF TO NODE: Q or L dofs by using an add aux level [_NoLevels+1][NoNodes[L]] //ONLY VOLUME  //NO PROC info, these are DOFS in the "SERIALIZED NUMBERING"
    int**  _Qnode_fine_Qnode_lev;   ///Un nodo fine puo' o non appartenere ad un certo livello... quindi devo continuare a mantenere i -1, altrimenti dovrei fare un search che fa perdere tempo
    int**   _off_nd;                 ///ONLY USED FOR THE EQUATION /// FROM DOF TO NODE: offsets for dofs  [QL][L+P*_NoLevels]                     //ONLY VOLUME
                                     // on the fine node numbering, the nodes corresponding to linear dofs are numbered FIRST... or not?

//===== Constructors/ Destructor ===========
     Mesh (Files& files_in, RunTimeMap<double>& map_in, const double Lref);
    ~Mesh ();
    void clear ();

    //======= Print/read functions =======
    void ReadMeshFile();
    void PrintForVisualizationAllLEVAllVB() const;
    void PrintSubdomFlagOnLinCells(std::string filename) const;
    void PrintMultimeshXdmf() const;
    void PrintXDMFAllLEVAllVB() const;
    void PrintConnLinAllLEVAllVB() const;
    void PrintXDMFGridVB(std::ofstream& out, std::ostringstream& top_file, std::ostringstream& geom_file,const uint Level, const uint vb) const;
    void PrintConnLinVB(hid_t file, const uint Level, const uint vb) const; 
    void PrintMeshFile(const std::string & namefile) const;
    void PrintMeshHDF5() const;
    void PrintElemVB( hid_t file, const uint vb , int* v_inv_nd , ElemStoBase** elem_sto, std::vector<std::pair<int,int> > v_el   ) const;
    void PrintSubdomFlagOnQuadrCells(const int vb, const int Level,std::string filename) const;

    //======= mesh generation functions ====
    void ElemChildToFather();
    void ReorderElementBySubdLev_VV();
    void ReorderElementBySubdLev_BB();
    void ComputeElemOffsetsBySubdLevel();

    void FillNodeSto();
    void ReorderNodesBySubdLev();
    void ComputeNodeOffsetsBySubdLevel();
    void ComputeMaxElXNode();
    void ComputeNodeMapExtLevels();
  
//must be public right now
    std::vector< std::pair<int,int> > _el_fm_libm; //because the EQUATION needs it for the SPARSITY PATTERN
    int *                             _el_libm_fm;
    ElemStoVol**  _el_sto;                 //FILLED ACCORDING TO "how the Libmesh mesh iterator runs" which may not be id in general i think...
    int *                             _nd_libm_fm; //from FINE LIBMESH NODE ORDERING to FINE FEMUS NODE ORDERING
    int    _maxelxnode;
    
    
protected:

    const double _Lref;          ///Reference length for non-dimensionalization
    
    std::string _mesh_file;    //mesh file name from the mesh generator

    // HDF5 FIELDS ===============
    std::string _nodes_name; //name for the HDF5 dataset
    std::string _elems_name;   //name for the HDF5 dataset
//     std::string _nd_coord_folder;  //TODO why seg fault if I use them?!?
//     std::string _el_pid_name;
//     std::string _nd_map_FineToLev;

//====================================
//filled by child only  (gencase)
//====================================

    // NODES ===============
    int    _n_nodes;       //of the WHOLE REFINEMENT! i.e. the FINE ones! //LMFILLS 
    int ** _n_nodes_sl_ql; 
    int *  _elxnode;              //number of elements per node of the fine mesh
    double *    _nd_coords_libm;  //node coordinates  //FILLED ACCORDING TO LIBMESH NODE ID ORDERING; then I'll print them according to my FEMUS ordering
    NodeSto**   _nd_sto;                       //FILLED ACCORDING TO LIBMESH NODE ID ORDERING

    std::vector< std::pair<int,int> > _nd_fm_libm; //from FINE FEMUS NODE ORDERING to FINE LIBMESH NODE ORDERING

    
    // ELEMENTS =============
    int ** _n_elements_sl_vb;
    int    _n_elements_sum_levs[VB];    //of the WHOLE REFINEMENT! //LMFILLS 
    int ** _el_child_to_fath;              //for every level, it gives you the father
    ElemStoBdry** _el_sto_b;               //FILLED with OUR ORDERING, "as we find them during the volume elem loop"
    
    std::vector< std::pair<int,int> > _el_fm_libm_b;


 };




} //end namespace femus



#endif