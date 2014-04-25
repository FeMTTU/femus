#ifndef __gencase__
#define __gencase__

//C++
#include <string>

// HDF5
#include "hdf5.h"

// FEMuS
#include "FEMTTUConfig.h"
#include "Typedefs.hpp"
#include "VBTypeEnum.hpp"
#include "FETypeEnum.hpp"
#include "RunTimeMap.hpp"
#include "ElemSto.hpp"
#include "MeshTwo.hpp"

// libmesh
#ifdef HAVE_LIBMESH
#include "libmesh/boundary_mesh.h"
#include "libmesh/mesh.h"
#endif

namespace femus {

// Forwarding classes
class FEElemBase;


class GenCase : public Mesh {

public:

     GenCase(const Files& files_in, const RunTimeMap<double> & map_in, const double lref,const std::string mesh_file);
    ~GenCase();                 
    
    
    void PrintMeshHDF5() const;
    void PrintElemVB( hid_t file, uint vb,const std::vector<int> & v_inv_nd , ElemStoBase** elem_sto, const std::vector<std::pair<int,int> >  v_el ) const;
    void ElemChildToFather();
    void ReorderElementBySubdLev_VV();
    void ReorderElementBySubdLev_BB();
    void ReorderNodesBySubdLev();
    void ComputeMaxElXNode();

    void ComputeMatrix();
    void ComputeProl();
    void ComputeRest();
    void PrintOneVarMatrixHDF5(const std::string & name, const std::string & groupname, uint** n_nodes_all, int count,int* Mat,int* len,int* len_off,int type1, int type2, int* FELevel ) const;
    void PrintOneVarMGOperatorHDF5(const std::string & filename, const std::string & groupname, uint* n_dofs_lev, int count,int* Rest,double* values,int* len,int* len_off, int FELevel, int FELevel2, int fe) const;

    void GenerateCase();
    void CreateStructuresLevSubd();
    
#ifdef HAVE_LIBMESH
    void GenerateCoarseMesh(libMesh::Mesh* msh_coarse) const;
    void RefineMesh(libMesh::Mesh* msh_all_levs) const;
    void GenerateBoundaryMesh(libMesh::BoundaryMesh* bd_msht, libMesh::Mesh* msh_all_levs) const;
    void GrabMeshinfoFromLibmesh(libMesh::Mesh* msht, libMesh::Mesh* msh);
#endif
    
private:

    // Element ===========
    std::vector<FEElemBase*> _feelems; //these are basically used only for the embedding matrix
    std::string _mesh_file;    //mesh file name from the mesh generator
    


    std::vector< std::pair<int,int> > _nd_fm_libm; //from FINE FEMUS NODE ORDERING to FINE LIBMESH NODE ORDERING
    std::vector< std::pair<int,int> > _el_fm_libm; //because the EQUATION needs it for the SPARSITY PATTERN
    std::vector< std::pair<int,int> > _el_fm_libm_b;
    std::vector<int>                  _nd_libm_fm; //from FINE LIBMESH NODE ORDERING to FINE FEMUS NODE ORDERING  //TODO this is the one that is not correctly filled in debug mode
    std::vector<int>                  _el_libm_fm;  //TODO in dbg mode they do not survive in GENCASE... extremely weird thing... let us try to convert from int* to std::vector<int>



};



} //end namespace femus



#endif
