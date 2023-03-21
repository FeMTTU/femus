#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "NumericVector.hpp"

#include "ElemType.hpp"





#include "../../../param.hpp"

using namespace femus;

#include  "../../../square_or_cube_03_cost_functional_without_regularization.hpp"

#include  "../../../opt_systems_elliptic_neumann.hpp"


double Solution_set_initial_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;

    if(!strcmp(name, "state")) {
        value = 0.;
    }
    else if(!strcmp(name, "control")) {
        value = 0.;
    }
    else if(!strcmp(name, "adjoint")) {
        value = 0.;
    }
    else if(!strcmp(name, "mu")) {
        value = 0.;
    }
    else if(!strcmp(name, "TargReg")) {
        value = femus::ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(x);
    }
    else if(!strcmp(name, "ContReg")) {
        value = femus::ctrl::  square_or_cube:: Domain_elements_containing_Gamma_control< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES > ::ControlDomainFlag_bdry(x);
    }
    else if(!strcmp(name, "act_flag")) {
        value = 0.;
    }


    return value;
}



bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet = true; //dirichlet
  value = 0;

  if(!strcmp(name,"control")) { //To enforce Neumann condition of control on Gamma_c
  if (faceName == 3)
    dirichlet = false;
  }

  if(!strcmp(name,"state")) {  //"state" corresponds to the first block row (u=q)
  if (faceName == 3)
    dirichlet = false;
  }

  if(!strcmp(name,"adjoint")) { //To enforce Neumann condition of adjoint on Gamma_c
  if (faceName == 3)
    dirichlet = false;
  }

  if(!strcmp(name,"mu")) {
//       value = 0.;
//   if (faceName == 3)
    dirichlet = false;
  }  
  
//     if(!strcmp(name,"adjoint")) {  //"adjoint" corresponds to the third block row
//   if (faceName == 3)    value = 1.;
//   }

  
  return dirichlet;
}





int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
 // ======= Problem  ==================
  MultiLevelProblem ml_prob;

  // ======= Files ========================
  Files files; 
  files.CheckIODirectories(true);
  files.RedirectCout(true);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double Lref = 1.;

  std::string input_file = "parametric_square_2x2.med";
//   std::string input_file = "parametric_square_4x5.med";
//   std::string input_file = "Mesh_3_groups_with_bdry_nodes.med";
//   std::string input_file = "Mesh_3_groups_with_bdry_nodes_coarser.med";
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_file;
  const std::string infile = mystream.str();

  const bool read_groups = true;
  const bool read_boundary_groups = true;
  

  mlMsh.ReadCoarseMesh(infile.c_str(), "seventh", Lref, read_groups, read_boundary_groups);
    
  
  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution ml_sol(&mlMsh);

  // add variables to ml_sol
  ml_sol.AddSolution("state", LAGRANGE, FIRST);
  ml_sol.AddSolution("control", LAGRANGE, FIRST);
  ml_sol.AddSolution("adjoint", LAGRANGE, FIRST);
  ml_sol.AddSolution("mu", LAGRANGE, FIRST);  

  ml_sol.AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  ml_sol.AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field

  //MU
  const bool      act_set_is_an_unknown_of_a_pde = false;
  std::vector<std::string> act_set_flag_name(1);  act_set_flag_name[0] = "act_flag";
  const unsigned int act_set_fake_time_dep_flag = 2;
  ml_sol.AddSolution(act_set_flag_name[0].c_str(), LAGRANGE, FIRST, act_set_fake_time_dep_flag, act_set_is_an_unknown_of_a_pde);
  ml_sol.Initialize(act_set_flag_name[0].c_str(), Solution_set_initial_conditions, & ml_prob);
  //MU
  
  ml_sol.Initialize("All");    // initialize all varaibles to zero

  ml_sol.Initialize("state",  Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize("control",  Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize("adjoint",  Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize("mu",  Solution_set_initial_conditions, & ml_prob);
  
  ml_sol.Initialize("TargReg",  Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize("ContReg",  Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize(act_set_flag_name[0].c_str(),  Solution_set_initial_conditions, & ml_prob);

  // attach the boundary condition function and generate boundary data
  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  ml_sol.GenerateBdc("state");
  ml_sol.GenerateBdc("control");
  ml_sol.GenerateBdc("adjoint");
  ml_sol.GenerateBdc("mu");  //we need add this to make the matrix iterations work...

 // ======= Problem, Mesh and Solution  ==================
 ml_prob.SetMultiLevelMeshAndSolution(& ml_sol);

  
  ml_prob.SetFilesHandler(&files);

 // add system  in ml_prob as a Linear Implicit System
  NonLinearImplicitSystemWithPrimalDualActiveSetMethod& system_opt = ml_prob.add_system < NonLinearImplicitSystemWithPrimalDualActiveSetMethod > ("LiftRestr");
  
  system_opt.SetActiveSetFlagName(act_set_flag_name);
//   system.SetMaxNumberOfNonLinearIterations(50);

  system_opt.AddSolutionToSystemPDE("state");  
  system_opt.AddSolutionToSystemPDE("control");  
  system_opt.AddSolutionToSystemPDE("adjoint");  
  system_opt.AddSolutionToSystemPDE("mu");  
  
  // attach the assembling function to system
  system_opt.SetAssembleFunction(elliptic::pure_boundary::assemble_elliptic_neumann_control);
  
  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);

  system_opt.SetDebugNonlinear(true);
//   system_opt.SetDebugFunction(femus::ctrl::cost_functional< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES,  femus::ctrl:: square_or_cube:: Domain_elements_containing_Gamma_control< femus::ctrl:: GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >  >::compute_cost_functional_regularization_bdry, ALPHA_CTRL_BDRY, BETA_CTRL_BDRY, QRULE_I );
   
  // initialize and solve the system
  system_opt.init();
  system_opt.MGsolve();
 
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
  ml_sol.GetWriter()->Write(files.GetOutputPath()/*DEFAULT_OUTPUTDIR*/, "biquadratic", variablesToBePrinted);

  return 0;
}


 
  

