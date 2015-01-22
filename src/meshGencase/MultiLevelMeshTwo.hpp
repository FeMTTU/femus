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
#include "FemusInputParser.hpp"
#include "GeomEl.hpp"
#include "ElemSto.hpp"
#include "VBTypeEnum.hpp"


namespace femus {


class Files;
class Domain;


  
class MultiLevelMeshTwo  {

public:

//===== Constructors/ Destructor ===========
     MultiLevelMeshTwo (const Files& files_in, const FemusInputParser<double>& map_in, const std::string mesh_file_in);
//     ~MultiLevelMeshTwo ();
    void clear ();

    //======= Print/read functions =======
    void ReadMeshFileAndNondimensionalize();
    void PrintForVisualizationAllLEVAllVB() const;
    void PrintSubdomFlagOnLinCells(std::string filename) const;
    void PrintMultimeshXdmf() const;
    void PrintXDMFAllLEVAllVB() const;
    void PrintConnLinAllLEVAllVB() const;
    void PrintXDMFGridVB(std::ofstream& out, std::ostringstream& top_file, std::ostringstream& geom_file,const uint Level, const uint vb) const;
    void PrintConnLinVB(hid_t file, const uint Level, const uint vb) const; 
    void PrintSubdomFlagOnQuadrCells(const int vb, const int Level,std::string filename) const;

    //======= mesh generation functions ====
    void ElemChildToFather();
    void ComputeElemOffsetsBySubdLevel();

    void FillNodeSto();
    void ComputeNodeOffsetsBySubdLevel();
    void ComputeMaxElXNode();
    void ComputeNodeMapExtLevels();

//attributes    ************************************
    
    ElemStoVol**  _el_sto;                 //FILLED ACCORDING TO "how the Libmesh mesh iterator runs" which may not be id in general i think...

    const uint _dim;               ///< spatial dimension
    const uint _mesh_order;

// ===== ABSTRACT GEOMEL(S) =====
    uint*      _type_FEM;         //just for check
    uint _elnodes[VB][QL];
    
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
// ELEMENTS =============
    int ** _n_elements_sl_vb;
    int    _n_elements_sum_levs[VB];    //of the WHOLE REFINEMENT! //LMFILLS 
    int ** _el_child_to_fath;              //for every level, it gives you the father
    ElemStoBdry** _el_sto_b;               //FILLED with OUR ORDERING, "as we find them during the volume elem loop"

// ===== NODES ======
    uint*   _NoNodesXLev;            ///VOLUME, at the EXTENDED LEVELS
    double* _xyz;                    ///< node coordinates
    uint**  _Qnode_lev_Qnode_fine;   ///ONLY USED FOR THE EQUATION /// FROM DOF TO NODE: Q or L dofs by using an add aux level [_NoLevels+1][NoNodes[L]] //ONLY VOLUME  //NO PROC info, these are DOFS in the "SERIALIZED NUMBERING"
    int**  _Qnode_fine_Qnode_lev;   ///Un nodo fine puo' o non appartenere ad un certo livello... quindi devo continuare a mantenere i -1, altrimenti dovrei fare un search che fa perdere tempo
    int**   _off_nd;                 ///ONLY USED FOR THE EQUATION /// FROM DOF TO NODE: offsets for dofs  [QL][L+P*_NoLevels]                     //ONLY VOLUME
                                   // on the fine node numbering, the nodes corresponding to linear dofs are numbered FIRST... or not?
// NODES ===============
    int    _n_nodes;       //of the WHOLE REFINEMENT! i.e. the FINE ones! //LMFILLS 
    int ** _n_nodes_sl_ql; 
    int *  _elxnode;              //number of elements per node of the fine mesh
    double *    _nd_coords_libm;  //node coordinates  //FILLED ACCORDING TO LIBMESH NODE ID ORDERING; then I'll print them according to my FEMUS ordering
    NodeSto**   _nd_sto;                       //FILLED ACCORDING TO LIBMESH NODE ID ORDERING
    int    _maxelxnode;
    
    // HDF5 FIELDS ===============
    std::string _nodes_name; //name for the HDF5 dataset
    std::string _elems_name;   //name for the HDF5 dataset
//     std::string _nd_coord_folder;  //TODO why seg fault if I use them?!?
//     std::string _el_pid_name;
//     std::string _nd_map_FineToLev;
  
    
   public:    
// ====== DOMAIN SHAPE (TODO optional => pointer) ----- //if I put it as reference I'd have to initialize it
    Domain* _domain;      //TODO You must remember to ALLOCATE this POINTER BEFORE USING IT!
    Domain* GetDomain() const;
    void    SetDomain(Domain* );
    void TransformElemNodesToRef(const uint elem_dim,const double* xx_qnds,double* refbox_xyz) const;
    
    //get functions
    inline const double get_Lref() const {return _Lref;}
    inline const uint   get_dim()  const {return _dim;}
    inline const GeomEl GetGeomEl(const uint dim, const uint order) const {   return _GeomEl[dim][order]; }
    inline const FemusInputParser<double>  GetRuntimeMap()  const { return _mesh_rtmap; }

    //set functions
    inline void SetLref(const double lref_in) { _Lref = lref_in; }
    
  protected:
    
    const Files& _files; 
    const FemusInputParser<double> & _mesh_rtmap;
    std::string _mesh_file;    //mesh file name from the mesh generator
 
   private:   
     
//attributes    ************************************
    
    double _Lref;          ///Reference length for non-dimensionalization
     
    std::vector< std::vector<GeomEl> >  _GeomEl;   //[DIM][QL_NODES] 
    

 };




} //end namespace femus



#endif