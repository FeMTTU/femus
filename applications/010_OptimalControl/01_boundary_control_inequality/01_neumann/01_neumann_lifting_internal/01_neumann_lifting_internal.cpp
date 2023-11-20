#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "VTKWriter.hpp"
#include "LinearImplicitSystem.hpp"


#include "00_cost_functional.hpp"


#include "../../../param.hpp"


using namespace femus;

#include  "../../../square_or_cube_03_cost_functional_without_regularization.hpp"

#include  "../../../opt_systems_elliptic_neumann.hpp"





double InitialValueContReg(const std::vector < double >& x) {
  return femus::ctrl:: square_or_cube:: lifting_internal< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >::ControlDomainFlag_internal_restriction(x);
}

double InitialValueTargReg(const std::vector < double >& x) {
  return femus::ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(x);
}

double InitialValueState(const std::vector < double >& x) {
  return 0.;
}

double InitialValueAdjoint(const std::vector < double >& x) {
  return 0.;
}

double InitialValueControl(const std::vector < double >& x) {
  return 0.; 
}

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet = true; //dirichlet
  value = 0;
  
  if(!strcmp(name,"state")) {
  if (faceName == 3)
    dirichlet = false;
  }
    
   if(!strcmp(name,"adjoint")) {
  if (faceName == 3) { value = 0.;   dirichlet = false; }
  }
  
  return dirichlet;
}





int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double Lref = 1.;

   std::string input_file = "square_0-1x0-1_divisions_2x2.med";
//   std::string input_file = "Mesh_3_groups_with_bdry_nodes.med";
//   std::string input_file = "Mesh_3_groups_with_bdry_nodes_coarser.med";
  std::ostringstream mystream; mystream << "../../../" << Files::_application_input_directory << "/" << input_file;
  const std::string infile = mystream.str();

  const bool read_groups = true;
  const bool read_boundary_groups = true;
  

  mlMsh.ReadCoarseMesh(infile.c_str(), "seventh", Lref, read_groups, read_boundary_groups);
    
 
  unsigned numberOfUniformLevels = 2;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("state", LAGRANGE, FIRST);
  mlSol.AddSolution("control", LAGRANGE, FIRST);
  mlSol.AddSolution("adjoint", LAGRANGE, FIRST);

  mlSol.AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  mlSol.AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field

  
  mlSol.Initialize("All");    // initialize all variables to zero

  mlSol.Initialize("state", InitialValueState);
  mlSol.Initialize("control", InitialValueControl);
  mlSol.Initialize("adjoint", InitialValueAdjoint);
  
  mlSol.Initialize("TargReg", InitialValueTargReg);
  mlSol.Initialize("ContReg", InitialValueContReg);

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("state");
  mlSol.GenerateBdc("control");
  mlSol.GenerateBdc("adjoint");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem ml_prob(&mlSol);
  
  // ======= Problem, Quad Rule - BEGIN  ========================
    std::vector< std::string > fe_quad_rule_vec;
  fe_quad_rule_vec.push_back("seventh");
  
  ml_prob.SetQuadratureRuleAllGeomElemsMultiple(fe_quad_rule_vec);
  ml_prob.set_all_abstract_fe_AD_or_not();
  // ======= Problem, Quad Rule - END  ========================

  // add system  in ml_prob as a Linear Implicit System
  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("LiftRestr");
 
  system.AddSolutionToSystemPDE("state");  
  system.AddSolutionToSystemPDE("control");  
  system.AddSolutionToSystemPDE("adjoint");  
  
  // attach the assembling function to system
  system.SetAssembleFunction( femus::elliptic::lifting_internal::assemble_elliptic_neumann_control );

  // initilaize and solve the system
  system.init();
  system.MGsolve();
  
   const unsigned  n_components_ctrl = 1;
   const unsigned n_components_state = n_components_ctrl;

  std::vector<std::string> state_vars(n_components_state);  state_vars[0] = "state";
  std::vector<std::string> ctrl_vars(n_components_ctrl);     ctrl_vars[0] = "control";
  
   femus::ctrl::cost_functional< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES,  femus::ctrl:: square_or_cube:: lifting_internal< femus::ctrl:: GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES > ,
                                femus::ctrl::COST_FUNCTIONAL_WITHOUT_REG  >::compute_cost_functional_regularization_lifting_internal(ml_prob, 0, 0, state_vars, ctrl_vars, ALPHA_CTRL_VOL, BETA_CTRL_VOL, QRULE_I);

 
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("state");
  variablesToBePrinted.push_back("control");
  variablesToBePrinted.push_back("adjoint");
  
  variablesToBePrinted.push_back("TargReg");
  variablesToBePrinted.push_back("ContReg");

    // ******* Print solution *******
  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(Files::_application_output_directory, "biquadratic", variablesToBePrinted);

  return 0;
}


