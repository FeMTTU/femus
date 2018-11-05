#ifndef __femus_meshGencase_GenCase_hpp__
#define __femus_meshGencase_GenCase_hpp__

//C++
#include <string>

// HDF5
#include "hdf5.h"

// FEMuS
#include "FemusConfig.hpp"
#include "Typedefs.hpp"
#include "VBTypeEnum.hpp"
#include "FETypeEnum.hpp"
#include "FemusInputParser.hpp"
#include "ElemSto.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "SystemTwo.hpp"

// libmesh
#ifdef HAVE_LIBMESH
#include "libmesh/boundary_mesh.h"
#include "libmesh/mesh.h"
#endif

namespace femus {

// Forwarding classes
class GeomElemBase;


class GenCase : public MultiLevelMeshTwo {

public:

     GenCase(const unsigned nolevels, const unsigned dim, const GeomElType geomel_type, const std::string mesh_file);
    ~GenCase();                 
    
    void ElemChildToFather();
    void ReorderElementBySubdLev_VV();
    void ReorderElementBySubdLev_BB();
    void ReorderNodesBySubdLev();
    void ComputeMaxElXNode();

    void ComputeAndPrintMGOperators(const std::string output_path);
    void ComputeAndPrintMatrix(const std::string output_path);
    void ComputeAndPrintProl(const std::string output_path);
    void ComputeAndPrintRest(const std::string output_path);
    
    static void ReadMGOps(const std::string output_path, SystemTwo * mysys);
    static void ReadMatrix(const std::string& name,  SystemTwo * mysys); 
    static void ReadProl(const std::string& name, SystemTwo * mysys);   
    static void ReadRest(const std::string& name, SystemTwo * mysys);
    
    void CreateMeshStructuresLevSubd(const std::string output_path);
    void Delete();

//functions using libmesh
    void GenerateCase(const std::string output_path);
    void GenerateCoarseMesh() const;
    void RefineMesh() const;
    void GenerateBoundaryMesh() const;
    void GrabMeshinfoFromLibmesh();

    // Element ===========
    std::vector<GeomElemBase*> _feelems; //these are basically used only for the embedding matrix

private:

#ifdef HAVE_LIBMESH
  libMesh::Mesh* _msh_coarse;
  libMesh::Mesh* _msh_all_levs;
  libMesh::BoundaryMesh* _bd_msht;
 #endif
  


};



} //end namespace femus



#endif
