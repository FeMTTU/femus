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
#include "MeshTwo.hpp"

namespace femus {

// Forwarding classes
class Files;
class GeomEl;
class Domain;
class FEElemBase;



class GenCase : public Mesh {

public:

     GenCase(Files& files_in,RunTimeMap<double> map_in, const double lref, std::string mesh_file);
    ~GenCase();                                      

    void ComputeMatrix();
    void ComputeProl();
    void ComputeRest();
    void PrintOneVarMatrixHDF5(std::string name, std::string groupname, uint** n_nodes_all, int count,int* Mat,int* len,int* len_off,int type1, int type2, int* FELevel ) const;
    void PrintOneVarMGOperatorHDF5(std::string filename, std::string groupname, uint* n_dofs_lev, int count,int* Rest,double* values,int* len,int* len_off, int FELevel, int FELevel2, int fe) const;

    void GenerateCase();
    void GenerateCoarseMesh(libMesh::Mesh* msh_coarse);
    void RefineMesh(libMesh::Mesh* msh_all_levs);
    void GenerateBoundaryMesh(libMesh::BoundaryMesh* bd_msht, libMesh::Mesh* msh_all_levs);
    void GrabMeshinfoFromLibmesh(libMesh::BoundaryMesh* bd_mesht, libMesh::Mesh* msht, libMesh::Mesh* msh);
    void CreateStructuresLevSubd();

    void ElemChildToFather();
    void ReorderElementBySubdLev_VV();
    void ReorderElementBySubdLev_BB();
    void ComputeElemOffsetsBySubdLevel();

    void FillNodeSto();
    void ReorderNodesBySubdLev();
    void ComputeNodeOffsetsBySubdLevel();
    void ComputeMaxElXNode();
    void ComputeNodeMapExtLevels();
  
private:

    // Basic data ==========
    std::string _mesh_file;

    // Element ===========
    std::vector<FEElemBase*> _feelems; //these are basically used only for the embedding matrix

    // NODES ===============
    int ** _n_nodes_sl_ql; 
    int *  _elxnode;              //number of elements per node of the fine mesh
    int    _maxelxnode;
    
    
    // ELEMENTS =============
    int ** _n_elements_sl_vb;
    
};



} //end namespace femus



#endif
