#include <sstream>
#include "FemusDefault.hpp"
#include "FemusInit.hpp"
#include "MultiLevelMesh.hpp"
#include "WriterEnum.hpp"
#include "MultiLevelSolution.hpp"

using namespace femus;

// Test for mesh file reading


int main(int argc,char **args) {

  // ======= Init ========================
  FemusInit init(argc,args,MPI_COMM_WORLD);

  // ======= Files ========================
  Files files; 
        files.CheckIODirectories();
        files.RedirectCout();
        
// fsi 3d - one layer
// volumes: 66
// faces:  66*2 + 42 = 132 + 42 = 174
// edges: 42*2 + 4 = 88
// 
// Total mesh:
// volumes: 264
// faces:  66*2  + 42*4 =  132 + 168 = 300     

 std::vector< std::string >  input_files;
 input_files.push_back("turek_FSI1.neu");
 input_files.push_back("turek_FSI1.med");
 input_files.push_back("turek_FSI1_3d.med");
 input_files.push_back("turek_FSI1_coarsest_not_yet_expanded_at_inflow.med");
//  input_files.push_back("turek_FSI1_no_bc.neu");
//    std::string input_file = "cyl.med";
//    std::string input_file = "horse2.med";
//    std::string input_file = "knot.neu";
//    std::string input_file = "dome_tri.med";
//    std::string input_file = "dome_quad.med";
//   std::string input_file = "Quad9_Four_boundaries_groups.med";
//   std::string input_file = "Quad9_Nine_without_groups.med";
//   std::string input_file = "Tri6_Two_boundaries.med"; 
//   std::string input_file = "Hex27_One_boundaries_groups.med";
//   std::string input_file = "Tet10_Twelve_boundaries.med"; ///@todo there seems to be an error in the output computation of biquadratic nodes
//   std::string input_file = "OneTet10.med";
 
        for(unsigned m = 0; m < input_files.size(); m++) {
            
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_files[m];
  const std::string infile = mystream.str();

  //Nondimensional
  double Lref = 1.;
  
  std::string fe_quad_rule("fifth");

  MultiLevelMesh ml_msh;
  
//   const int n_sub = 1;
//   ml_msh.GenerateCoarseBoxMesh(n_sub,n_sub,0,0.,1.,0.,1.,0.,0.,TRI6,fe_quad_rule.c_str());
  const bool read_groups = true;
  const bool read_boundary_groups = true;
  ml_msh.ReadCoarseMesh(infile.c_str(),fe_quad_rule.c_str(),Lref,read_groups, read_boundary_groups);
  
  ml_msh.PrintInfo();
  
  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution ml_sol(&ml_msh); //I think I cannot read the mesh alone, I need to attach at least one solution

  // add variables to ml_sol
  ml_sol.AddSolution("u", LAGRANGE, FIRST);
  ml_sol.Initialize("All"); 
  
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
  
//       std::vector < std::string > surfaceVariables;
//       surfaceVariables.push_back("X");
//       surfaceVariables.push_back("Y");
//       surfaceVariables.push_back("Z");
// 
//     ml_sol.GetWriter()->SetSurfaceVariables(surfaceVariables);
  
  const std::string output_dir = files.GetOutputPath();
//   const std::string output_dir = DEFAULT_OUTPUTDIR;
  
  
  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);  //false: only Sol; true: adds EpsSol, ResSol, BdcSol


  ml_sol.GetWriter()->Write(input_files[m], output_dir, "linear", variablesToBePrinted);
  ml_sol.GetWriter()->Write(input_files[m], output_dir, "quadratic", variablesToBePrinted);
  ml_sol.GetWriter()->Write(input_files[m], output_dir, "biquadratic", variablesToBePrinted);
  
//   ml_sol.SetWriter(XDMF); 
//   ml_sol.GetWriter()->SetDebugOutput(true);  //false: only Sol; true: adds EpsSol, ResSol, BdcSol
//   ml_sol.GetWriter()->Write(output_dir, "linear", variablesToBePrinted);
//   ml_sol.GetWriter()->Write(output_dir, "quadratic", variablesToBePrinted);
//   ml_sol.GetWriter()->Write(output_dir, "biquadratic", variablesToBePrinted);

// recent versions of Paraview do not read the GMV format
//   ml_sol.SetWriter(GMV);  
//   ml_sol.GetWriter()->SetDebugOutput(true);  //false: only Sol; true: adds EpsSol, ResSol, BdcSol
//   ml_sol.GetWriter()->Write(output_dir,"biquadratic",variablesToBePrinted);  
            
        }
        


  
  return 0;
}

