#include "FemusDefault.hpp"
#include "FemusInit.hpp"
#include "MultiLevelMesh.hpp"

using namespace femus;

// Test for SalomeIO reading


int main(int argc,char **args) {

  FemusInit init(argc,args,MPI_COMM_WORLD);
  
//   std::string med_file = "Mesh_1.med";
//   std::string infile = "./" + DEFAULT_CONFIGDIR + med_file;
//  
//   //Adimensional
//   double Lref = 1.;
//   
//   //MultiLevelMesh ml_msh(nm,nr,infile,"seventh",Lref,SetRefinementFlag); 
//   MultiLevelMesh ml_msh;
//   ml_msh.ReadCoarseMesh(infile,"seventh",Lref);
  
  return 0;
}
