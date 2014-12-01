// C++
#include <iostream>

// FEMuS
#include "FEMTTUConfig.h"
#include "FemusDefault.hpp"
#include "Files.hpp"
#include "GeomEl.hpp"
#include "Box.hpp"

#include "FEElemBase.hpp"
#include "FEQuad4.hpp"
#include "FEQuad9.hpp"
#include "FEHex8.hpp"
#include "FEHex27.hpp"
#include "FETri3.hpp"
#include "FETri6.hpp"
#include "FETet4.hpp"
#include "FETet10.hpp"
#include "FETypeEnum.hpp"
#include "RunTimeMap.hpp"
#include "TimeLoop.hpp"
#include "CmdLine.hpp"

// Libmesh
#ifdef HAVE_LIBMESH
#include "libmesh/libmesh.h" //for libmeshinit

using namespace libMesh;
#endif

//application
#include "GenCase.hpp"

using namespace femus;


// ==========================================
//         GENCASE
// ==========================================
//The goal of the gencase program is to put some files 
// into the input/ directory of the corresponding Application

//If you want to generate the mesh with libmesh, set libmesh_gen=1

//If you want to read it from an external file: 
// - set libmesh_gen=0,
// - put the external file (e.g. mesh.gam) into the CONFIG directory
// - and set the file name in F_MESH_READ in the runtime file
// Only, beware that "dimension" is the same as the mesh dimension in the file
// (anyway we have a check for that)

//This is a special application in the sense all that its configuration
//is dependent on the App that we are considering,
//so it is a sort of "service application"


// Considering instead the configuration of the classes, the classes that are involved 
// are the FE because you must know how the dofs are arranged, the GeomEl which gives you 
// the type of Geometric element of the class. Strangely, Mesh is not involved. This is because
//Mesh only "dialogues" with the output files of gencase. But clearly, Mesh depends on GeomEl,
//so there cannot be TWO DIFFERENT GeomEl for Gencase and for the App, as it happens to be the case of course.

//TODO can we make the gencase output in parallel?!
// The fact of the parallel output is what refrains me from thinking of some alternative 
// to redirecting std output INSIDE C++ instead of FROM SHELL... maybe there is a shell variable
// that holds the MPI Rank... the problem is that it would exist after a command
// of the type "mpiexec  > file$MYPROC.txt" would be launched...
//TODO see what happens with libmesh in debug mode
//TODO so far Gencase is only reading the tags <Mesh> <Box> and <Files>. Probably soon we'll remove <Files>
// and perhaps merge <Mesh> and <Box>

int main(int argc, char** argv) {

#ifdef HAVE_LIBMESH  
  LibMeshInit init(argc, argv);   // Initialize libMesh
#else
  std:: cout << "Cannot run gencase without libmesh right now" << std::endl; abort();  
#endif

  CmdLine::parse(argc,argv);
  std::string chosen_app = CmdLine::get("--app"); 
  std::cout << "**************************************************************" << std::endl;
  std::cout << "***** The application I will generate the case for is  ****** " << chosen_app << std::endl;
  std::cout << "**************************************************************" << std::endl;

  // ======= Files =====
  Files::CheckDirOrAbort("../",chosen_app);
  std::string basepath = "../" + chosen_app + "/";
  Files files(basepath);
   std::cout << "******The basepath starting from the gencase directory is ***** " << files._app_path << std::endl;

  //======== first of all, check for restart ======  //======== don't rerun gencase in that case, to avoid spending time on rebuilding the operators and so on ======      
  files.ConfigureRestart();
  if (files._restart_flag == true) { std::cout << "Do not rerun gencase in case of restart" << std::endl; abort(); }
  
  files.CheckDir(files._app_path,DEFAULT_CASEDIR); //here, we must check if the input directory where gencase writes is there  //if not, we make it

  RunTimeMap<std::string> files_map("Files",files._app_path);

// ========= GenCase =====
  RunTimeMap<double> mesh_map("Mesh",files._app_path);
  GenCase gencase(files,mesh_map,1.,files_map.get("F_MESH_READ"));
  
// if domain is a Box, instantiate it
// After that, the generation of the coarse mesh DOES depend on the DOMAIN SHAPE. in fact, the parameters of the box contain info for generating the COARSE mesh
//but now, imagine you read from mesh files...
// so, if you have a domain shape, then it should be involved in the coarse mesh generation
// if you DONT have a domain shape, then the generation of the coarse mesh should only be file-dependent  
//Ok so even though we have a shape, its mesh may come from a mesh file
// also, if you want to rotate/translate ecc, you should not do it  AFTER refining the mesh, but right after coarse generation
// i should put things together with FindCmake
// anyway, the idea is mesh.generateCoarseMesh(OptionalShape);  // if the optionalshape is there, use its function to generate the coarse mesh
                                                                // otherwise
  
          gencase.GenerateCase();


  
  std::cout << "=======End of GenCase========" << std::endl;

  return 0;  
}