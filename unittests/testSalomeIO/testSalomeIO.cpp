#include <sstream>
#include "FemusDefault.hpp"
#include "FemusInit.hpp"
#include "MultiLevelMesh.hpp"

using namespace femus;

// Test for SalomeIO reading


int main(int argc,char **args) {

  FemusInit init(argc,args,MPI_COMM_WORLD);
  
  std::string med_file = "two_meshes.med";
  std::ostringstream mystream; mystream << "./" << DEFAULT_CONFIGDIR << "/" << med_file;
  const std::string infile = mystream.str();
 
  //Adimensional
  double Lref = 1.;
  
  MultiLevelMesh ml_msh;
  ml_msh.ReadCoarseMesh(infile.c_str(),"seventh",Lref);
  
  return 0;
}
