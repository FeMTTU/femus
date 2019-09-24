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

  //std::string input_file = "turek_FSI1.neu";
  std::string input_file = "testovaso.med";
//   std::string input_file = "Quad9_Four_boundaries_groups.med";
//   std::string input_file = "Quad9_Nine_without_groups.med";
//   std::string input_file = "Tri6_Two_boundaries.med";
//   std::string input_file = "Hex27_One_boundaries_groups.med";
//   std::string input_file = "Tet10_Twelve_boundaries.med"; ///@todo there seems to be an error in the output computation of biquadratic nodes
//   std::string input_file = "OneTet10.med";
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_file;
  const std::string infile = mystream.str();

  //Nondimensional
  double Lref = 1.;

  std::string fe_quad_rule("fifth");

  MultiLevelMesh ml_msh;
//   const int n_sub = 4;
//   ml_msh.GenerateCoarseBoxMesh(n_sub,n_sub,0,0.,1.,0.,1.,0.,0.,QUAD9,fe_quad_rule.c_str());
  ml_msh.ReadCoarseMesh(infile.c_str(),fe_quad_rule.c_str(),Lref);
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
  ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR,"linear",variablesToBePrinted);
  ml_sol.SetWriter(XDMF);
  ml_sol.GetWriter()->SetDebugOutput(true);  //false: only Sol; true: adds EpsSol, ResSol, BdcSol
  ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR,"linear",variablesToBePrinted);

// recent versions of Paraview do not read the GMV format
//   ml_sol.SetWriter(GMV);
//   ml_sol.GetWriter()->SetDebugOutput(true);  //false: only Sol; true: adds EpsSol, ResSol, BdcSol
//   ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR,"biquadratic",variablesToBePrinted);


  return 0;
}
