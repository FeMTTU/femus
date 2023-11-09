/** tutorial/Ex3
 * This example shows how to set and solve the weak form of the nonlinear problem
 *                     -\Delta^2 u = f(x) \text{ on }\Omega,
 *            u=0 \text{ on } \Gamma,
 *      \Delta u=0 \text{ on } \Gamma,
 * on a box domain $\Omega$ with boundary $\Gamma$,
 * by using a system of second order partial differential equation.
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/



#include "FemusInit.hpp"
#include "Files.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "VTKWriter.hpp"
#include "NumericVector.hpp"

//#include "biharmonic_coupled.hpp"

#include "FE_convergence.hpp"

#include "Solution_functions_over_domains_or_mesh_files.hpp"


#include "adept.h"


#define LIBRARY_OR_USER   1 //0: library; 1: user

#if LIBRARY_OR_USER == 0
   #include "01_biharmonic_coupled.hpp"
   #define NAMESPACE_FOR_BIHARMONIC   femus
#elif LIBRARY_OR_USER == 1
   #include "biharmonic_coupled.hpp"
   #define NAMESPACE_FOR_BIHARMONIC   karthik
#endif



using namespace femus;




double GetExactSolutionValue(const std::vector < double >& x) {
  double pi = acos(-1.);
  return cos(pi * x[0]) * cos(pi * x[1]);
};


void GetExactSolutionGradient(const std::vector < double >& x, std::vector < double >& solGrad) {
  double pi = acos(-1.);
  solGrad[0]  = -pi * sin(pi * x[0]) * cos(pi * x[1]);
  solGrad[1] = -pi * cos(pi * x[0]) * sin(pi * x[1]);
};


double GetExactSolutionLaplace(const std::vector < double >& x) {
  double pi = acos(-1.);
  return -2.*pi * pi * cos(pi * x[0]) * cos(pi * x[1]);       // - pi*pi*cos(pi*x[0])*cos(pi*x[1]);
};






bool SetBoundaryCondition_bc_all_dirichlet_homogeneous(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0;
  return dirichlet;
}





int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // ======= Files - BEGIN  ========================
  const bool use_output_time_folder = false; // This allows you to run the code multiple times without overwriting. This will generate an output folder each time you run.
  const bool redirect_cout_to_file = true; // puts the output in a log file instead of the term
  Files files;
        files.CheckIODirectories(use_output_time_folder);
        files.RedirectCout(redirect_cout_to_file);

  // ======= Files - END  ========================


    // ======= System Specifics - BEGIN  ==================
  system_specifics  system_biharmonic_coupled;   //me

  system_biharmonic_coupled._mesh_files.push_back("square_-0p5-0p5x-0p5-0p5_divisions_2x2.med");
  const std::string relative_path_to_build_directory =  "../../../../";
  const std::string mesh_file = relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/02_2d/square/minus0p5-plus0p5_minus0p5-plus0p5/";  system_biharmonic_coupled._mesh_files_path_relative_to_executable.push_back(mesh_file);

  system_biharmonic_coupled._system_name = "Biharmonic";
  system_biharmonic_coupled._assemble_function = NAMESPACE_FOR_BIHARMONIC :: biharmonic_equation :: AssembleBilaplaceProblem_AD;

  system_biharmonic_coupled._boundary_conditions_types_and_values             = SetBoundaryCondition_bc_all_dirichlet_homogeneous;

  Domains::square_01by01::Function_Zero_on_boundary_5<>   system_biharmonic_coupled_function_zero_on_boundary_1;
  Domains::square_01by01::Function_Zero_on_boundary_5_Laplacian<>   system_biharmonic_coupled_function_zero_on_boundary_1_laplacian;
  system_biharmonic_coupled._assemble_function_for_rhs   = & system_biharmonic_coupled_function_zero_on_boundary_1_laplacian; //this is the RHS for the auxiliary variable v = -Delta u
  system_biharmonic_coupled._true_solution_function      = & system_biharmonic_coupled_function_zero_on_boundary_1;
  ///@todo if this is not set, nothing happens here. It is used to compute absolute errors
    // ======= System Specifics - END ==================



  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  const std::string mesh_file_total = system_biharmonic_coupled._mesh_files_path_relative_to_executable[0] + "/" + system_biharmonic_coupled._mesh_files[0];
  mlMsh.ReadCoarseMesh(mesh_file_total.c_str(), "seventh", scalingFactor);

  unsigned maxNumberOfMeshes = 5;

  std::vector < std::vector < double > > l2Norm;
  l2Norm.resize(maxNumberOfMeshes);

  std::vector < std::vector < double > > semiNorm;
  semiNorm.resize(maxNumberOfMeshes);

    std::vector<FEOrder> feOrder;
    feOrder.push_back(FIRST);
    feOrder.push_back(SERENDIPITY);
    feOrder.push_back(SECOND);





  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

    unsigned numberOfUniformLevels = i + 1;
    unsigned numberOfSelectiveLevels = 0;
    mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

    // erase all the coarse mesh levels
    mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

    // print mesh info
    mlMsh.PrintInfo();

    l2Norm[i].resize( feOrder.size() );
    semiNorm[i].resize( feOrder.size() );

    for (unsigned j = 0; j < feOrder.size(); j++) {   // loop on the FE Order

      // define the multilevel solution and attach the mlMsh object to it
      MultiLevelSolution mlSol(&mlMsh);

      // add variables to mlSol
      mlSol.AddSolution("u", LAGRANGE, feOrder[j]);
      mlSol.AddSolution("v", LAGRANGE, feOrder[j]);
      mlSol.Initialize("All");


      // define the multilevel problem attach the mlSol object to it
      MultiLevelProblem ml_prob(&mlSol);

      ml_prob.set_app_specs_pointer(& system_biharmonic_coupled);
      // ======= Problem, Files ========================
      ml_prob.SetFilesHandler(&files);

      // attach the boundary condition function and generate boundary data
      mlSol.AttachSetBoundaryConditionFunction( system_biharmonic_coupled._boundary_conditions_types_and_values );
      mlSol.GenerateBdc("u", "Steady", & ml_prob);
      mlSol.GenerateBdc("v", "Steady", & ml_prob);

      // add system Biharmonic in ml_prob as a Linear Implicit System
      NonLinearImplicitSystem& system = ml_prob.add_system < NonLinearImplicitSystem > (system_biharmonic_coupled._system_name);

      // add solution "u" to system
      system.AddSolutionToSystemPDE("u");
      system.AddSolutionToSystemPDE("v");

      // attach the assembling function to system
      system.SetAssembleFunction( system_biharmonic_coupled._assemble_function );

      // initialize and solve the system
      system.init();

      system.MGsolve();

      // convergence for u
      std::pair< double , double > norm = GetErrorNorm_L2_H1_with_analytical_sol(&mlSol, "u", GetExactSolutionValue, GetExactSolutionGradient );
      // // convergence for v
      // std::pair< double , double > norm = GetErrorNorm_L2_H1_with_analytical_sol(&mlSol, "v", LaplaceGetExactSolutionValue, LaplaceGetExactSolutionGradient );


      l2Norm[i][j]  = norm.first;
      semiNorm[i][j] = norm.second;


      // print solutions
      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back("All");

      VTKWriter vtkIO(&mlSol);
      vtkIO.Write(Files::_application_output_directory, "biquadratic", variablesToBePrinted, i);

    }
  }


  // FE_convergence::output_convergence_order();


  // ======= L2 - BEGIN  ========================
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < feOrder.size(); j++) {
      std::cout << l2Norm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < feOrder.size(); j++) {
        std::cout << log(l2Norm[i][j] / l2Norm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }
  // ======= L2 - END  ========================



  // ======= H1 - BEGIN  ========================

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "SEMINORM ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < feOrder.size(); j++) {
      std::cout << semiNorm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < feOrder.size(); j++) {
        std::cout << log(semiNorm[i][j] / semiNorm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }

  // ======= H1 - END  ========================


  return 0;
}






