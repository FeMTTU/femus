#include <sstream>
#include "FemusDefault.hpp"
#include "FemusInit.hpp"
#include "MultiLevelMesh.hpp"
#include "WriterEnum.hpp"
#include "MultiLevelSolution.hpp"

using namespace femus;

// Test for SalomeIO reading


int main(int argc,char **args) {

  FemusInit init(argc,args,MPI_COMM_WORLD);

  std::string med_file = "FourQuad9_boundaries.med";
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << med_file;
  const std::string infile = mystream.str();

  //Adimensional
  double Lref = 1.;

  MultiLevelMesh ml_msh;
  ml_msh.ReadCoarseMesh(infile.c_str(),"fifth",Lref);
  ml_msh.PrintInfo();
  
  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution ml_sol(&ml_msh);

  // add variables to ml_sol
  ml_sol.AddSolution("u", LAGRANGE, FIRST);
  ml_sol.Initialize("All"); 
  
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");

  
  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);  //false: only Sol; true: adds EpsSol, ResSol, BdcSol
  ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR,"biquadratic",variablesToBePrinted);
  ml_sol.SetWriter(GMV);
  ml_sol.GetWriter()->SetDebugOutput(true);  //false: only Sol; true: adds EpsSol, ResSol, BdcSol
  ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR,"biquadratic",variablesToBePrinted);
  ml_sol.SetWriter(XDMF); 
  ml_sol.GetWriter()->SetDebugOutput(true);  //false: only Sol; true: adds EpsSol, ResSol, BdcSol
  ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR,"biquadratic",variablesToBePrinted);

  
  return 0;
}

