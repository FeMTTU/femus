// C++
#include <iostream>

// Libmesh
#include "libmesh/libmesh.h" //for libmeshinit

// FEMuS
#include "FemusDefault.hpp"
#include "Files.hpp"
#include "Utils.hpp"
#include "GeomEl.hpp"
#include "Box.hpp"

#include "GenCase.hpp"

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


// ==========================================
//         GENCASE
// ==========================================
//The goal of the gencase program is to put some files 
// into the input/ directory of the corresponding Application

//If you want to generate the mesh with libmesh, set libmesh_gen=1

//If you want to read it from an external file: 
// - set libmesh_gen=0,
// - put the external file (e.g. mesh.gam) into the CONFIG directory
// - and set the file name in F_MESH_READ in "param_files.in"
// Only, beware that "dimension" is the same as the mesh dimension in the file
// (anyway we have a check for that)

//This is a special application in the sense all that its configuration
//is dependent on the App that we are considering,
//so it is a sort of "service application"
//in fact, it links with the libfemus as compiled for the corresponding application
//May be possible for us to ISOLATE the parameters that are 
//necessary ONLY to GENCASE?
//Well there are some parameters that are required for both, like "dimension", NDOF_FEM
//But i think we could make things much more ordered

//For instance, of the three parameter files (which could actually be put inside the class
// subdirectory in config/, as each of them is closely related to the classes Utils Files and Physics)
// we have that:
//param_phys is not involved for the gencase
//param_utils is involved
//param_files is involved

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


int main(int argc, char** argv) {

  LibMeshInit init(argc, argv);   // Initialize libMesh

  CmdLine::parse(argc,argv);
  std::string chosen_app = CmdLine::get("--app"); 
  std::cout << "**************************************************************" << std::endl;
  std::cout << "***** The application I will generate the case for is  ****** " << chosen_app << std::endl;
  std::cout << "**************************************************************" << std::endl;
  std::string basepath = "../" + chosen_app + "/";

  // ======= Files =====
  Files files(basepath);
   std::cout << "******The basepath starting from the gencase directory is ***** " << files.get_basepath() << std::endl;
  files.CheckDirOrAbort(files.get_basepath(),DEFAULT_CONFIGDIR); 
  files.get_frtmap().read();
  files.get_frtmap().print();
  files.CheckDir(files.get_basepath(),files.get_frtmap().get("INPUT_DIR")); //here, we must check if the input directory where gencase writes is there  //if not, we make it
  
  // ======= Utils =====
  Utils utils(files);
  utils._urtmap.read();
  utils._urtmap.print();

  // =======GeomEl =====
  uint geomel_type = (uint) utils._urtmap.get("geomel_type");
  uint dimension   = (uint) utils._urtmap.get("dimension");
  GeomEl geomel(dimension,geomel_type);

  // =======FEElems =====
std::vector<FEElemBase*> FEElements(QL); //these are basically used only for the embedding matrix
  for (int fe=0; fe<QL; fe++) FEElements[fe] = FEElemBase::build(&geomel,fe);

  // ======= Domain =====
  //here we miss the "service" nature of the gencase program
  //we must instantiate the domain explicitly
  //we should do the ifdefs here for all types of Domains...
  //the Domain is not related to the MESH ORDER.
  //it may be related to the MESH MAPPING if we have curved elements or not
  //but the fact of having curved elements might be more dependent on 
  //the Mesh than the Doamin shape.
  //anyway, now we do like this
      
//======== check for restart ======      
     RunTimeMap<double> timemap("TimeLoop",files.get_basepath());  //here you don't need to instantiate a TimeLoop object, but only to read its RUNTIME MAP
     timemap.read();
     timemap.print();
     TimeLoop::check_time_par(timemap);
     const uint restart      = (uint) timemap.get("restart");    // restart param
     if (restart)  {std::cout << "No GenCase because of restart flag !=0 " << std::endl; abort();}

    
  // ========= GenCase =====
      GenCase gencase(utils,geomel,FEElements);
      gencase._rtmap.read();
      gencase._rtmap.print();
      gencase.GenerateCase();

  std::cout << "=======End of GenCase========" << std::endl;

//clean
  //remember that the CHILD destructor is not called here. This is because 
  //these things are of FEElemBase type
  //TODO Do I have to put a delete [] or only delete? Do standard vectors have an overloading of the "delete" operator?
  for (int fe=0; fe<QL; fe++) delete FEElements[fe];
  
  return 0;  
}

//=======you dont need the solver library for gencase
//======= neither FEMUS_HAVE_PETSC nor HAVE_MPI
// gencase needs mpi somewhere