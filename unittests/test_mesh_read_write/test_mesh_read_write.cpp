#include "FemusDefault.hpp"
#include "FemusInit.hpp"
#include "MultiLevelMesh.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "WriterEnum.hpp"
#include "VTKWriter.hpp"

#include <sstream>

using namespace femus;


//====================================================
// Here various layers of an application can be activated


#define FEMUS_TEST_INIT  1

#define FEMUS_TEST_FILES  1

#define FEMUS_TEST_INPUT_PARSER  1

#define FEMUS_TEST_MESH  1

#if FEMUS_TEST_MESH != 0
  #define FEMUS_TEST_MESH_PRINT 0

  #define FEMUS_TEST_SOLUTION  1


  #if FEMUS_TEST_SOLUTION != 0
     #define FEMUS_TEST_SOLUTION_PRINT 1

     #define FEMUS_TEST_PROBLEM  1
  #endif


#endif
//====================================================




  #if FEMUS_TEST_SOLUTION != 0
double Solution_set_initial_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;

    if(!strcmp(name, "u_lag_first")) {
        value = 1.;
    }
    else if(!strcmp(name, "u_lag_serendip")) {
        value = 2.;
    }
    else if(!strcmp(name, "u_lag_second")) {
        value = 3.;
    }
    else if(!strcmp(name, "u_disc_zero")) {
        value = 4.;
    }
    else if(!strcmp(name, "u_disc_first")) {
        value = 5.;
    }

   value = x[0] ;

    return value;
}
#endif






int main(int argc,char **args) {
  
  // ======= Init - BEGIN  ========================
#if FEMUS_TEST_INIT != 0
  FemusInit init(argc, args, MPI_COMM_WORLD);
#endif
  // ======= Init - END  ========================


  // ======= Files - BEGIN  ========================
#if FEMUS_TEST_FILES != 0
  const bool use_output_time_folder = false;
  const bool redirect_cout_to_file = false;
  Files files; 
        files.CheckIODirectories(use_output_time_folder);
        files.RedirectCout(redirect_cout_to_file);
#endif
 // ======= Files - END  ========================

        
#if FEMUS_TEST_INPUT_PARSER != 0
  // it works, just pay attention that integers do not turn into unsigned
  //   files.CopyInputFiles();   // at this point everything is in the folder of the current run!!!!
  FemusInputParser< unsigned > mesh_map("Mesh", "./"/*files.GetOutputPath()*/);
  const unsigned numberOfUniformLevels =  mesh_map.get("n_levels");
#endif
  
  
  // ======= Loop over mesh files ========================
 std::vector< std::string >  input_files;

 input_files.push_back("square_0-1x0-1_divisions_2x2.med");


 
  for(unsigned m = 0; m < input_files.size(); m++) {

            
#if FEMUS_TEST_MESH != 0
    
  // ======= Mesh - BEGIN  ========================
  std::ostringstream mystream; mystream << "../../include/"  << "00_mesh_files/" << 
  "00_salome/00_parametric/with_notebook/" << input_files[m];
  const std::string infile = mystream.str();

  //Nondimensional
  double Lref = 1.;
  
  std::string fe_quad_rule("fifth");

  MultiLevelMesh ml_mesh;
  
//   const int n_sub = 1;
//   ml_mesh.GenerateCoarseBoxMesh(n_sub,n_sub,0,0.,1.,0.,1.,0.,0.,TRI6,fe_quad_rule.c_str());
  const bool read_groups = true;
  const bool read_boundary_groups = true;
  ml_mesh.ReadCoarseMesh(infile.c_str(), fe_quad_rule.c_str(), Lref, read_groups, read_boundary_groups);
//   const unsigned numberOfUniformLevels = 3;
  const unsigned erased_levels = numberOfUniformLevels - 1;
  unsigned numberOfSelectiveLevels = 0;
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
//   ml_mesh.EraseCoarseLevels(erased_levels);
  
  ml_mesh.PrintInfo();
  // ======= Mesh - END  ========================
  
//============ Solution - BEGIN ==================
  
#if FEMUS_TEST_SOLUTION != 0
  
  MultiLevelSolution ml_sol(&ml_mesh);

  #if FEMUS_TEST_PROBLEM != 0
      MultiLevelProblem   ml_prob(&ml_sol);
  #endif 
   
  const unsigned  steady_flag = 0;                //0: steady state, 2: time dependent
  const bool      is_an_unknown_of_a_pde = false; //0: not associated to any System
  ml_sol.AddSolution("u_lag_first", LAGRANGE, FIRST, steady_flag, is_an_unknown_of_a_pde);
  ml_sol.AddSolution("u_lag_serendip", LAGRANGE, SERENDIPITY, steady_flag, is_an_unknown_of_a_pde);
  ml_sol.AddSolution("u_lag_second", LAGRANGE, SECOND, steady_flag, is_an_unknown_of_a_pde);
  ml_sol.AddSolution("u_disc_zero", DISCONTINUOUS_POLYNOMIAL, ZERO, steady_flag, is_an_unknown_of_a_pde);
  ml_sol.AddSolution("u_disc_first", DISCONTINUOUS_POLYNOMIAL, FIRST, steady_flag, is_an_unknown_of_a_pde);

    for(unsigned sol = 0; sol <  ml_sol.GetSolutionSize(); sol++) {
     const std::string sol_name(ml_sol.GetSolutionName(sol));
  ml_sol.Initialize(sol_name.c_str(), Solution_set_initial_conditions, & ml_prob);
  //   ml_sol.Initialize("all"); 
    }

 
#endif

//============ Solution - END ==================

#endif
  

  
  //============ Print - BEGIN ==================
  
#if FEMUS_TEST_FILES != 0
  const std::string output_dir = files.GetOutputPath();
#endif

  std::vector < std::string > print_fe_order;
  print_fe_order.push_back("linear");
  print_fe_order.push_back("quadratic");
  print_fe_order.push_back("biquadratic");


  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
  
//       std::vector < std::string > surfaceVariables;
//       surfaceVariables.push_back("X");
//       surfaceVariables.push_back("Y");
//       surfaceVariables.push_back("Z");
// 
//     ml_sol.GetWriter()->SetSurfaceVariables(surfaceVariables);

#if FEMUS_TEST_MESH_PRINT != 0
/// The print doesn't need a solution object anymore 

  VTKWriter vtk_writer(&ml_mesh); /* cannot instantiate the father because it has pure virtual */
  vtk_writer.SetDebugOutput(true);
  VTKWriter * writer_mesh_ptr = &vtk_writer;
  
//============ Print: Loop over levels ==================
  for(unsigned l = 0; l < ml_mesh.GetNumberOfLevels(); l++) {
     for(std::vector< std::string >::iterator print_fe_it = std::begin(print_fe_order); print_fe_it != std::end(print_fe_order); ++ print_fe_it) {

        writer_mesh_ptr->Write(l+1, input_files[m], output_dir, "_only_mesh", (*print_fe_it).c_str(), variablesToBePrinted);

         }
   }
#endif


#if FEMUS_TEST_SOLUTION_PRINT != 0

  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);  //false: only Sol; true: adds EpsSol, ResSol, BdcSol
  Writer * writer_sol_ptr = ml_sol.GetWriter();
  
//============ Print: Loop over levels ==================
  for(unsigned l = 0; l < ml_mesh.GetNumberOfLevels(); l++) {
     for(std::vector< std::string >::iterator print_fe_it = std::begin(print_fe_order); print_fe_it != std::end(print_fe_order); ++ print_fe_it) {

        writer_sol_ptr->Write(l+1, input_files[m], output_dir, "", (*print_fe_it).c_str(), variablesToBePrinted);
  
//   ml_sol.SetWriter(XDMF); 
//   ml_sol.GetWriter()->SetDebugOutput(true);  //false: only Sol; true: adds EpsSol, ResSol, BdcSol
//   ml_sol.GetWriter()->Write(output_dir, *it, variablesToBePrinted);

// recent versions of Paraview do not read the GMV format
//   ml_sol.SetWriter(GMV);  
//   ml_sol.GetWriter()->SetDebugOutput(true);  //false: only Sol; true: adds EpsSol, ResSol, BdcSol
//   ml_sol.GetWriter()->Write(output_dir,"biquadratic",variablesToBePrinted);  
         }
   }

#endif

  //============ Print - END ==================


    }


  
  return 0;
}


/// @todo the print for disc_first is wrong
