/** tutorial/Ex2
 * This example shows how to set and solve the weak form of the Poisson problem
 *                    $$ \Delta u = \Delta u_exact \text{ on }\Omega, $$
 *                    $$ u=0 \text{ on } \Gamma, $$
 * on a square domain $\Omega$ with boundary $\Gamma$;
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "Files.hpp"
#include "MultiLevelSolution.hpp"
#include "LinearImplicitSystem.hpp"

#include "NumericVector.hpp"



#include "00_poisson_eqn_with_all_dirichlet_bc_AD_or_nonAD_separate.hpp"

#include "FE_convergence.hpp"


#include "../tutorial_common.hpp"

#include "../all_mesh_generation_methods.hpp"



using namespace femus;








double GetExactSolutionValue(const std::vector < double >& x) {
  double pi = acos(-1.);
  return cos(pi * x[0]) * cos(pi * x[1]);
}


void GetExactSolutionGradient(const std::vector < double >& x, vector < double >& solGrad) {
  double pi = acos(-1.);
  solGrad[0]  = -pi * sin(pi * x[0]) * cos(pi * x[1]);
  solGrad[1] = -pi * cos(pi * x[0]) * sin(pi * x[1]);
}


double GetExactSolutionLaplace(const std::vector < double >& x) {
  double pi = acos(-1.);
  return -pi * pi * cos(pi * x[0]) * cos(pi * x[1]) - pi * pi * cos(pi * x[0]) * cos(pi * x[1]);
}




void AssemblePoissonProblem_old_fe_quadrature_nonAD_interface(MultiLevelProblem& ml_prob) {
          AssemblePoissonProblem_old_fe_quadrature_nonAD(ml_prob, GetExactSolutionLaplace, NULL);
}

void AssemblePoissonProblem_old_fe_quadrature_AD_interface(MultiLevelProblem& ml_prob) {
          AssemblePoissonProblem_old_fe_quadrature_AD(ml_prob, GetExactSolutionLaplace, NULL);
}




bool square_bc_all_dirichlet(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {
  
  //this sets Dirichlet over all faces
  
  bool dirichlet = true; //dirichlet
  value = 0.;

  return dirichlet;
}




int main(int argc, char** args) {

  // ======= Init ========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // ======= Problem  ==================
  MultiLevelProblem ml_prob;

  
  // ======= Files - BEGIN  ========================
  const bool use_output_time_folder = false;
  const bool redirect_cout_to_file = true;
  Files files; 
        files.CheckIODirectories(use_output_time_folder);
        files.RedirectCout(redirect_cout_to_file);

  // ======= Problem, Files ========================
  ml_prob.SetFilesHandler(&files);
  // ======= Files - END  ========================

  
  
  // ======= Quad Rule - BEGIN ========================
   std::string fe_quad_rule("seventh");
    // ======= Quad Rule - END ========================
 
   
  // ======= Mesh file types (function, salome, gambit) - BEGIN ========================
    for (unsigned  mesh_file_type = 0; mesh_file_type < 3; mesh_file_type++) {  
  
  // define multilevel mesh
  MultiLevelMesh ml_mesh;

  // read coarse level mesh
   std::string mesh_name = Domain_square_m05p05::quad9_all_mesh_generation_methods_Structured(mesh_file_type, ml_mesh, fe_quad_rule);


   
  unsigned dim = ml_mesh.GetDimension();
  unsigned maxNumberOfMeshes;

  if (dim == 1 || dim == 2) {
    maxNumberOfMeshes = 6;
  } else {
    maxNumberOfMeshes = 4;
  }

  vector < vector < double > > l2Norm;
  l2Norm.resize(maxNumberOfMeshes);

  vector < vector < double > > semiNorm;
  semiNorm.resize(maxNumberOfMeshes);

      
  // ======= Assemble functions (AD or NON-AD) - BEGIN ========================
    std::vector < std::pair <femus::System::AssembleFunctionType, std::string> >  assemble_pointer_vec(2);
    
    assemble_pointer_vec[0].first = AssemblePoissonProblem_old_fe_quadrature_nonAD_interface;
    assemble_pointer_vec[0].second = "non-automatic_diff";
    
    assemble_pointer_vec[1].first = AssemblePoissonProblem_old_fe_quadrature_AD_interface;
    assemble_pointer_vec[1].second = "automatic_diff";

    
    
   for (unsigned func = 0; func < assemble_pointer_vec.size(); func++) {

  std::cout << std::endl;
  std::cout << std::endl;
     std::cout << "===================================" << std::endl;
     std::cout << "===================================" << std::endl;
     std::cout << assemble_pointer_vec[func].second << std::endl;
     std::cout << "===================================" << std::endl;
     std::cout << "===================================" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;


  
  // ======= Mesh refinements - BEGIN ========================
  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

    unsigned numberOfUniformLevels = i + 1;
    unsigned numberOfSelectiveLevels = 0;
    ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

    // erase all the coarse mesh levels
    ml_mesh.EraseCoarseLevels(numberOfUniformLevels - 1);

    // print mesh info
    ml_mesh.PrintInfo();

    std::vector< Unknown > unknowns = systems__generate_list_of_scalar_unknowns_for_each_FE_family_lagrangian();
    
    l2Norm[i].resize(unknowns.size());
    semiNorm[i].resize(unknowns.size());

      
      
  // ======= FE SPACES - BEGIN ========================
    for (unsigned int u = 0; u < unknowns.size(); u++) {

        // define the multilevel solution and attach the ml_mesh object to it
      MultiLevelSolution ml_sol(&ml_mesh);

      ml_sol.SetWriter(VTK);
      ml_sol.GetWriter()->SetDebugOutput(true);


      // add variables to ml_sol
      ml_sol.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order, unknowns[u]._time_order, unknowns[u]._is_pde_unknown);

      ml_sol.Initialize("All");

      // attach the boundary condition function and generate boundary data
      ml_sol.AttachSetBoundaryConditionFunction(square_bc_all_dirichlet);
      ml_sol.GenerateBdc( unknowns[u]._name.c_str() );

      // ======= Problem, Mesh and Solution  ==================
      // attach the ml_sol object
       ml_prob.SetMultiLevelMeshAndSolution(& ml_sol);
       

       ml_prob.get_systems_map().clear();
       ml_prob.set_current_system_number(0/*u*/);               //way to communicate to the assemble function, which doesn't belong to any class
      
    // ======= System - BEGIN ========================
       std::string sys_name = "Poisson" + unknowns[u]._name;
      // add system Poisson in ml_prob as a Linear Implicit System
      LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > (sys_name);


        // add solution "u" to system
      system.AddSolutionToSystemPDE( unknowns[u]._name.c_str() );

        // set unknown list
        std::vector< Unknown > unknowns_vec(1);
        unknowns_vec[0] = unknowns[u]; //need to turn this into a vector

        system.set_unknown_list_for_assembly(unknowns_vec); //way to communicate to the assemble function, which doesn't belong to any class


        // attach the assembling function to system
      system.SetAssembleFunction( assemble_pointer_vec[func].first );

      // initialize and solve the system
      system.init();
      
      system.SetOuterSolver(PREONLY);
      system.MGsolve();
    // ======= System - END ========================

  // ======= Errors - BEGIN  ========================
      std::pair< double , double > norm = GetErrorNorm_L2_H1_with_analytical_sol(& ml_sol, unknowns_vec, GetExactSolutionValue, GetExactSolutionGradient);
      l2Norm[i][u]  = norm.first;
      semiNorm[i][u] = norm.second;
  // ======= Errors - END ========================
      

  // ======= Print - BEGIN  ========================
      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back("All");
            
      std::string  output_name = mesh_name + "_" + assemble_pointer_vec[func].second + "_" + unknowns[u]._name;
//       ml_sol.GetWriter()->SetGraphVariable ("u");
      ml_sol.GetWriter()->Write(output_name.c_str(), DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, i);
  // ======= Print - END  ========================

    } //end FE  space
  // ======= FE SPACES - END ========================

      
      
} //end mesh level
  // ======= Mesh refinements - END ========================


  // ======= H1 - BEGIN  ========================
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "SEMINORM ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < semiNorm[i].size(); j++) {
      std::cout << semiNorm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < semiNorm[i].size(); j++) {
        std::cout << log(semiNorm[i][j] / semiNorm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }
  // ======= H1 - END  ========================
  
  
  // ======= L2 - BEGIN  ========================
  std::cout << std::endl;
  std::cout << std::endl;

  // print the seminorm of the error and the order of convergence between different levels
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < l2Norm[i].size(); j++) {
      std::cout << l2Norm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < l2Norm[i].size(); j++) {
        std::cout << log(l2Norm[i][j] / l2Norm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }
  // ======= L2 - END  ========================

    } //end assemble func loop
  // ======= Assemble functions (AD or NON-AD) - END ========================

    
    }
  // ======= Mesh file types (function, salome, gambit) - END ========================
    
    
  return 0;
}


