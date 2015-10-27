#include <sstream>
#include "FemusDefault.hpp"
#include "FemusInit.hpp"
#include "MultiLevelMesh.hpp"
#include "WriterEnum.hpp"

using namespace femus;

// Test for SalomeIO reading


int main(int argc,char **args) {

  FemusInit init(argc,args,MPI_COMM_WORLD);
  
  std::string med_file = "GroupsANDMeshes.med";
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << med_file;
  const std::string infile = mystream.str();
 
  //Adimensional
  double Lref = 1.;
  
  MultiLevelMesh ml_msh;
  ml_msh.ReadCoarseMesh(infile.c_str(),"fifth",Lref);
  
  ml_msh.SetWriter(XDMF);
  ml_msh.GetWriter()->write(DEFAULT_OUTPUTDIR,"biquadratic");
  ml_msh.SetWriter(GMV);
  ml_msh.GetWriter()->write(DEFAULT_OUTPUTDIR,"biquadratic");
  ml_msh.SetWriter(VTK);
  ml_msh.GetWriter()->write(DEFAULT_OUTPUTDIR,"biquadratic");

  return 0;
}
