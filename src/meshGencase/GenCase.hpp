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
#include "FemusInputParser.hpp"
#include "ElemSto.hpp"
#include "MultiLevelMeshTwo.hpp"

// libmesh
#ifdef HAVE_LIBMESH
#include "libmesh/boundary_mesh.h"
#include "libmesh/mesh.h"
#endif

namespace femus {

// Forwarding classes
class FEElemBase;


class GenCase : public MultiLevelMeshTwo {

public:

     GenCase(const FemusInputParser<double> & map_in, const std::string mesh_file);
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
    void CreateMeshStructuresLevSubd(const std::string output_path);
    void Delete();

//functions using libmesh
    void GenerateCase(const std::string output_path);
    void GenerateCoarseMesh() const;
    void RefineMesh() const;
    void GenerateBoundaryMesh() const;
    void GrabMeshinfoFromLibmesh();

    // Element ===========
    std::vector<FEElemBase*> _feelems; //these are basically used only for the embedding matrix

private:

#ifdef HAVE_LIBMESH
  libMesh::Mesh* _msh_coarse;
  libMesh::Mesh* _msh_all_levs;
  libMesh::BoundaryMesh* _bd_msht;
 #endif
  


};



} //end namespace femus



#endif
