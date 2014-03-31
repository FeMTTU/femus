#ifndef __gencase__
#define __gencase__

//C++
#include <string>
// HDF5
#include "hdf5.h"
// libmesh
#include "libmesh/boundary_mesh.h"
#include "libmesh/mesh.h"
// FEMuS
#include "Typedefs.hpp"
#include "VBTypeEnum.hpp"
#include "FETypeEnum.hpp"
#include "RunTimeMap.hpp"
#include "ElemSto.hpp"

// Forwarding classes
class Utils;
class GeomEl;
class Domain;
class FEElemBase;



class GenCase {

private:

    // Basic data ==========
    Utils& _utils;
    int _NoSubdom;       ///< number of subdomain (one for each processor)
    int _NoLevels;       ///< number of levels (multigrid)

    // Element ===========
    const uint _dimension;
    GeomEl&  _GeomEl; ///<  GeomEl class pointer
    std::vector<FEElemBase*>& _feelems;
    uint _elnodes[VB][QL];

    // NODES ===============
    int    _n_nodes;       //of the WHOLE REFINEMENT! i.e. the FINE ones!
    int *  _NoNodesXLev;
    int ** _n_nodes_sl_ql; 
    int ** _off_nd;               //node offsets, QL (Quadratic and Linear)
    int ** _Qnode_fine_Qnode_lev; // was _gindexL  //node map  // from QUADRATIC NODE FINE, to QUADRATIC NODE of that LEVEL
    int ** _Qnode_lev_Qnode_fine;                                 
    int *  _elxnode;              //number of elements per node of the fine mesh
    int    _maxelxnode;
    NodeSto**   _nd_sto;                       //FILLED ACCORDING TO LIBMESH NODE ID ORDERING
    double *    _nd_coords_libm;  //node coordinates  //FILLED ACCORDING TO LIBMESH NODE ID ORDERING; then I'll print them according to my FEMUS ordering

    std::vector< std::pair<int,int> > _nd_fm_libm; //from FINE FEMUS NODE ORDERING to FINE LIBMESH NODE ORDERING
    int *                             _nd_libm_fm; //from FINE LIBMESH NODE ORDERING to FINE FEMUS NODE ORDERING
    
    // ELEMENTS =============
    int    _n_elements_sum_levs[VB];    //of the WHOLE REFINEMENT!
    int ** _n_elements_lev_vb;
    int ** _n_elements_sl_vb;
    int ** _off_el;
    int ** _el_child_to_fath;              //for every level, it gives you the father
    ElemStoVol**  _el_sto;                 //FILLED ACCORDING TO "how the Libmesh mesh iterator runs" which may not be id in general i think...
    ElemStoBdry** _el_sto_b;               //FILLED with OUR ORDERING, "as we find them during the volume elem loop"
    
    std::vector< std::pair<int,int> > _el_fm_libm;
    std::vector< std::pair<int,int> > _el_fm_libm_b;
    int *                             _el_libm_fm;

    // HDF5 FIELDS ===============
    std::string _nodes_name; //name for the HDF5 dataset
    std::string _elems_name;   //name for the HDF5 dataset
//     std::string _nd_coord_folder;  //TODO why seg fault if I use them?!?
//     std::string _el_pid_name;
//     std::string _nd_map_FineToLev;

public:

     GenCase(Utils& mgutils_in,GeomEl& mggeomel_in,std::vector<FEElemBase*>& feelems_in); ///<  Constructor
    ~GenCase();                                       ///<  Destructor

    void ComputeMatrix();
    void ComputeProl();
    void ComputeRest();
    void PrintOneVarMatrixHDF5(std::string name, std::string groupname, int** n_nodes_all, int count,int* Mat,int* len,int* len_off,int type1, int type2, int* FELevel );
    void PrintOneVarMGOperatorHDF5(std::string filename, std::string groupname, int* n_dofs_lev, int count,int* Rest,double* values,int* len,int* len_off, int FELevel, int FELevel2, int fe);

    void GenerateCase();
    void GenerateCoarseMesh(libMesh::Mesh* msh_coarse);
    void RefineMesh(libMesh::Mesh* msh_all_levs);
    void GenerateBoundaryMesh(libMesh::BoundaryMesh* bd_msht, libMesh::Mesh* msh_all_levs);
    void GrabMeshinfoFromLibmesh(libMesh::BoundaryMesh* bd_mesht, libMesh::Mesh* msht, libMesh::Mesh* msh);
    void CreateStructuresLevSubd();

    void PrintElemVB( hid_t file, const uint vb , int* v_inv_nd , ElemStoBase** elem_sto, std::vector<std::pair<int,int> > v_el   ) ;
    void PrintSubdomFlagOnQuadrCells(const int vb, const int Level,std::string filename) const;
    void PrintMultimeshXdmf();
    void PrintMeshHDF5();

    void ElemChildToFather();
    void ReorderElementBySubdLev_VV();
    void ReorderElementBySubdLev_BB();
    void ComputeElemOffsetsBySubdLevel();

    void FillNodeSto();
    void ReorderNodesBySubdLev();
    void ComputeNodeOffsetsBySubdLevel();
    void ComputeMaxElXNode();
    void ComputeNodeMapExtLevels();
    
};


#endif
