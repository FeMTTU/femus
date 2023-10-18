/** tutorial/Ex2
 * This example shows how to set and solve the weak form of the Poisson problem
 *                    $$ \Delta u = \Delta u_exact \text{ on }\Omega, $$
 *                    $$ u=0 \text{ on } \Gamma, $$
 * on a square domain $\Omega$ with boundary $\Gamma$;
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/

#include "adept.h"


#include "FemusInit.hpp"
#include "Files.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NumericVector.hpp"
#include "LinearImplicitSystem.hpp"

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
          AssemblePoissonProblem_old_fe_quadrature_nonAD(ml_prob, GetExactSolutionLaplace, GetExactSolutionValue);
}

void AssemblePoissonProblem_old_fe_quadrature_AD_interface(MultiLevelProblem& ml_prob) {
          AssemblePoissonProblem_old_fe_quadrature_AD(ml_prob, GetExactSolutionLaplace, GetExactSolutionValue);
}



bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {
  
  bool dirichlet = true; //dirichlet
  value = 0.;

  return dirichlet;
}





int main(int argc, char** args) {

  // ======= Init ========================
  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // ======= Files - BEGIN  ========================
  const bool use_output_time_folder = false;
  const bool redirect_cout_to_file = false;
  Files files; 
        files.CheckIODirectories(use_output_time_folder);
        files.RedirectCout(redirect_cout_to_file);
  // ======= Files - END  ========================

 
  // ======= Quad Rule ========================
  std::string fe_quad_rule("seventh");
 /* "seventh" is the order of accuracy that is used in the gauss integration scheme
    In the future it is not going to be an argument of the mesh function   */

  // ======= Mesh - BEGIN ========================
 // define multilevel mesh
  MultiLevelMesh ml_mesh;

  // read coarse level mesh
  const unsigned mesh_file_type = 1; //salome   
  const std::string mesh_name = Domain_square_m05p05::quad9_all_mesh_generation_methods_Structured(mesh_file_type, ml_mesh, fe_quad_rule);
  
  MultiLevelMesh ml_mesh_finest;
  const std::string mesh_name_finest = Domain_square_m05p05::quad9_all_mesh_generation_methods_Structured(mesh_file_type, ml_mesh_finest, fe_quad_rule);

 
  unsigned dim = ml_mesh.GetDimension();
  unsigned maxNumberOfMeshes;

  if (dim == 1 || dim == 2) {
    maxNumberOfMeshes = 6;
  } else {
    maxNumberOfMeshes = 4;
  }
  // ======= Mesh - END ========================

  const unsigned gap = 1;

  vector < vector < double > > l2Norm;
  l2Norm.resize(maxNumberOfMeshes + 1 - gap);

  vector < vector < double > > semiNorm;
  semiNorm.resize(maxNumberOfMeshes + 1 - gap);

  
  MultiLevelSolution * ml_sol_finest;
  

    
    std::vector< Unknown > unknowns = systems__generate_list_of_scalar_unknowns_for_each_FE_family_lagrangian();

    
    for (int i = maxNumberOfMeshes - 1; i >= 0; i--) {   // loop on the mesh level

    unsigned numberOfUniformLevels = i + 1;
    unsigned numberOfSelectiveLevels = 0;
    ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

    // erase all the coarse mesh levels
    ml_mesh.EraseCoarseLevels(numberOfUniformLevels - 1);

    
    if (i == maxNumberOfMeshes - 1) {
        unsigned numberOfUniformLevels_finest = numberOfUniformLevels;
        ml_mesh_finest.RefineMesh(numberOfUniformLevels_finest, numberOfUniformLevels_finest + numberOfSelectiveLevels, NULL);
//         ml_mesh_finest.EraseCoarseLevels(numberOfUniformLevels_finest - 1 - 1); //I need to keep the structures at all levels here so I can restrict every time
    }
    
    // print mesh info
    ml_mesh.PrintInfo();
    
    if (i < l2Norm.size()) {
    l2Norm[i].resize(unknowns.size());
    semiNorm[i].resize(unknowns.size());
    }
    
    for (unsigned j = 0; j < unknowns.size(); j++) {   // loop on the FE Order
      // define the multilevel solution and attach the ml_mesh object to it
      MultiLevelSolution ml_sol(&ml_mesh);
      
      ml_sol.SetWriter(VTK);
      ml_sol.GetWriter()->SetDebugOutput(true);
      
      // add variables to ml_sol
      ml_sol.AddSolution(unknowns[j]._name.c_str(), unknowns[j]._fe_family, unknowns[j]._fe_order, unknowns[j]._time_order, unknowns[j]._is_pde_unknown);
      ml_sol.Initialize("All");
      
      // attach the boundary condition function and generate boundary data
      ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
      ml_sol.GenerateBdc(unknowns[j]._name.c_str());
      
      // define the multilevel problem attach the ml_sol object to it
      MultiLevelProblem ml_prob(&ml_sol);

       ml_prob.SetFilesHandler(&files);
      
       ml_prob.get_systems_map().clear();
       ml_prob.set_current_system_number(0/*u*/);               //way to communicate to the assemble function, which doesn't belong to any class
      
      // add system Poisson in ml_prob as a Linear Implicit System
      LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("Poisson");

      // add solution "u" to system
      system.AddSolutionToSystemPDE( unknowns[j]._name.c_str() );

      // set unknown list
      std::vector< Unknown > unknowns_vec(1);
      unknowns_vec[0] = unknowns[j]; //need to turn this into a vector

      system.set_unknown_list_for_assembly(unknowns_vec); //way to communicate to the assemble function, which doesn't belong to any class

        
      // attach the assembling function to system
      system.SetAssembleFunction(AssemblePoissonProblem_old_fe_quadrature_nonAD_interface);
      // system.SetAssembleFunction(AssemblePoissonProblem_old_fe_quadrature_AD_interface);

      // initialize and solve the system
      system.init();
//       system.ClearVariablesToBeSolved();
//       system.AddVariableToBeSolved("All");
// 
//       ml_sol.SetWriter(VTK);
//       ml_sol.GetWriter()->SetDebugOutput(true);
//   
//       system.SetDebugLinear(true);
//       system.SetMaxNumberOfLinearIterations(6);
//       system.SetAbsoluteLinearConvergenceTolerance(1.e-4);

      system.SetOuterSolver(PREONLY);
      system.MGsolve();
      
     
    if (i < maxNumberOfMeshes - gap/*1*/) {
        //restrict the fine solution at the current level  ==================
//     for(unsigned short j = 0; j < ml_sol._mlMesh->GetNumberOfLevels(); j++) { //all levels
      const unsigned coarse = i;
      const unsigned fine   = coarse + gap;
    for (unsigned nf = 0; nf < gap; nf++) {
      ml_sol_finest->CoarsenSolutionByOneLevel_wrong( fine - nf);
    }
//        }

        //pass the restriction of the fine solution to the function that computes the error  ==================
      Solution* sol = ml_sol_finest->GetSolutionLevel(i);    // pointer to the solution (level) object
      std::pair< double , double > norm = GetErrorNorm_L2_H1_multiple_methods(&ml_sol, sol, unknowns_vec, GetExactSolutionValue, GetExactSolutionGradient );

      l2Norm[i][j]  = norm.first;
      semiNorm[i][j] = norm.second;
      
      
        
    }
    else if (i == maxNumberOfMeshes - 1) {
        //store the fine solution  ==================
              ml_sol_finest = new MultiLevelSolution (&ml_mesh_finest);  //with the declaration outside and a "new" inside it persists outside the loop scopes
              ml_sol_finest->SetWriter(VTK);
              ml_sol_finest->AddSolution(unknowns[j]._name.c_str(), unknowns[j]._fe_family, unknowns[j]._fe_order, unknowns[j]._time_order, unknowns[j]._is_pde_unknown);
              ml_sol_finest->Initialize("All");
              ml_sol_finest->AttachSetBoundaryConditionFunction(SetBoundaryCondition);
              ml_sol_finest->GenerateBdc( unknowns[j]._name.c_str() );
              
//       assert( ml_sol._mlMesh->GetNumberOfLevels() == ml_sol_finest->_mlMesh->GetNumberOfLevels() - 1 );

//     for(unsigned short k = 0; k < ml_sol._mlMesh->GetNumberOfLevels(); k++) { //all levels (only one)
    //       _solution[k]->CopySolutionToOldSolution();  //started from here
              
       const unsigned level_index_current = 0;
      //@todo there is a duplicate function in MLSol: GetSolutionLevel() and GetLevel()
     for(unsigned short j = 0; j <   ml_sol.GetSolutionLevel(level_index_current)->_Sol.size(); j++) {    //all variables
               *(ml_sol_finest->GetLevel(i)->_Sol[j]) = *(ml_sol.GetSolutionLevel(level_index_current)->_Sol[j]);
          }
//        }
      
      
       
       
    }
        
    
    
      // print solutions
      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back("All");

      ml_sol.GetWriter()->Write(files.GetOutputPath(), "biquadratic", variablesToBePrinted, i);
      
      ml_sol_finest->GetWriter()->Write(i+1, files.GetOutputPath(), "biquadratic"/*_finest*/, variablesToBePrinted, i+3 + maxNumberOfMeshes);

    }
  }
  
  delete ml_sol_finest;
  
  
  
  // ======= H1 - BEGIN  ========================
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "SEMINORM ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (int i = l2Norm.size() - 1; i >= 0; i--) {   // loop on the mesh level
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < unknowns.size(); j++) {
      std::cout << semiNorm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < l2Norm.size() - 2) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < unknowns.size(); j++) {
        std::cout << log(semiNorm[i][j] / semiNorm[i + 1][j]) / log(2.) << "  \t\t";
      }

      std::cout << std::endl;
    }

  }
  // ======= H1 - END  ========================



  // ======= L2 - BEGIN  ========================
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (int i = l2Norm.size() - 1; i >= 0; i--) {   // loop on the mesh level
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < unknowns.size(); j++) {
      std::cout << l2Norm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < l2Norm.size() - 2) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < unknowns.size(); j++) {
        std::cout << log(l2Norm[i][j] / l2Norm[i + 1][j]) / log(2.) << "  \t\t";
      }

      std::cout << std::endl;
    }

  }
  // ======= L2 - END  ========================

  return 0;
}



