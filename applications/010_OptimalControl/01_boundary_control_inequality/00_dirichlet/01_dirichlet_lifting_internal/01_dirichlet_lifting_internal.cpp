#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"

#include "00_cost_functional.hpp"


#include "../../../param.hpp"
#include  "../../../opt_systems_elliptic_dirichlet.hpp"


using namespace femus;





 // Unknowns - BEGIN  ==================
double Solution_set_initial_conditions_Unknowns(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;
    
        return value;
}



bool Solution_set_boundary_conditions_Unknowns(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet = false; //dirichlet
  value = 0.;
  
 
   if(!strcmp(name,"state")) {
       
        dirichlet = true;
        
         value = 0.;
       
    }
  
   else if(!strcmp(name, "control")) {



     ctrl:: NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS < ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >::ctrl_or_state_set_dirichlet_flags(ml_prob, faceName, x, dirichlet);
     ctrl:: NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS < ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >::ctrl_or_state_set_dirichlet_fixed_values(ml_prob, faceName, x, value);
                  
                    
                    
   }
  
    else if(!strcmp(name,"adjoint")) {
        dirichlet = true;
        value = 0.;
    }

  
    //MU
    else if(!strcmp(name,"mu")) {
    dirichlet = false;
//       value = 0.;
  }
  
  return dirichlet;
}
 // Unknowns - END  ==================


    
 // Not Unknowns - BEGIN  ==================
double Solution_set_initial_conditions_Not_Unknowns(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;
    
    if(!strcmp(name,"TargReg")) {
        value = ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(x);
    }
    else if(!strcmp(name,"ContReg")) {
        value = ctrl:: square_or_cube:: lifting_internal< ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >::ControlDomainFlag_internal_restriction(x);
    }
    else if(!strcmp(name,"act_flag")) {
        value = 0.;
    }


    return value;
}
 // Not Unknowns - END  ==================



int main(int argc, char** args) {

  // ======= Init ========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
  // ======= Problem - BEGIN ========================
  MultiLevelProblem ml_prob;
  
  // ======= Problem, Files - BEGIN  ========================
  Files files; 
        files.CheckIODirectories(ctrl::use_output_time_folder);
        files.RedirectCout(ctrl::redirect_cout_to_file);

  // ======= Problem, Files ========================
  ml_prob.SetFilesHandler(&files);
  // ======= Problem, Files - END  ========================

  
  // ======= Problem, Mesh - BEGIN ==================

  // ======= Mesh, Coarse reading - BEGIN ==================
  MultiLevelMesh ml_mesh;

  const std::string infile = Files::get_input_file_with_prefix(femus::ctrl::mesh_input, "../../../");

  
  ml_mesh.set_mesh_filename( ctrl::mesh_input );
  
  const double Lref = 1.;

  const bool read_groups = true;
  const bool read_boundary_groups = true;
    
  ml_mesh.ReadCoarseMesh(infile, Lref, read_groups, read_boundary_groups);

  // ======= Mesh, Coarse reading - END ==================

  
  // ======= Mesh: Refinement - BEGIN ==================
  unsigned numberOfUniformLevels = N_UNIFORM_LEVELS;
  unsigned numberOfSelectiveLevels = 0;
  const unsigned erased_levels = N_ERASED_LEVELS;
  
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  // ======= Mesh: Refinement - END ==================

  // ======= Mesh: Coarse erasing - BEGIN  ========================
  ml_mesh.EraseCoarseLevels(erased_levels/*numberOfUniformLevels - 1*/);
  ml_mesh.PrintInfo();
  // ======= Mesh: Coarse erasing - END  ========================

  // ======= Problem, Mesh - END ==================

  
  // ======= Problem, Solution - BEGIN ==================

  // ======= Solution, Construction - BEGIN ==================
  MultiLevelSolution ml_sol(&ml_mesh);

  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);
  
 // ======= Problem, Mesh and Solution  ==================
 ml_prob.SetMultiLevelMeshAndSolution(& ml_sol);
  // ======= Solution, Construction - END ==================


  // ======= Solutions that are Unknowns - BEGIN ==================
  std::vector< Unknown > unknowns = elliptic :: lifting_internal :: provide_list_of_unknowns( ml_mesh.GetDimension() );
 
  for (unsigned int u = 0; u < unknowns.size(); u++) { ml_sol.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order, unknowns[u]._time_order, unknowns[u]._is_pde_unknown); }
   
 for (unsigned int u = 0; u < unknowns.size(); u++)  { ml_sol.Initialize(unknowns[u]._name.c_str(), Solution_set_initial_conditions_Unknowns, & ml_prob); }
  
  ml_sol.AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions_Unknowns);
   for (unsigned int u = 0; u < unknowns.size(); u++)  {  ml_sol.GenerateBdc(unknowns[u]._name.c_str(), (unknowns[u]._time_order == 0) ? "Steady" : "Time_dependent", & ml_prob);  }
  // ======= Solutions that are Unknowns - END ==================

 
  // ======= Solutions that are not Unknowns - BEGIN  ==================
  const unsigned  steady_flag = 0;
  const bool      is_an_unknown_of_a_pde = false;
  
  // ******** targ reg - BEGIN 
  ml_sol.AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO, steady_flag, is_an_unknown_of_a_pde); //this variable is not solution of any eqn, it's just a given field
  ml_sol.Initialize("TargReg",     Solution_set_initial_conditions_Not_Unknowns, & ml_prob);
  // ******** targ reg - END

  // ******** cont reg - BEGIN 
  ml_sol.AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO, steady_flag, is_an_unknown_of_a_pde); //this variable is not solution of any eqn, it's just a given field
  ml_sol.Initialize("ContReg",     Solution_set_initial_conditions_Not_Unknowns, & ml_prob);
  // ******** cont reg - END

  // ******** active flag - BEGIN 
  //MU
  const unsigned int  n_components_ctrl = 1;
  
  const unsigned int act_set_fake_time_dep_flag = 2;  //this is needed to be able to use _SolOld
  const bool         act_flag_is_an_unknown_of_a_pde = false;

  unsigned int index_control = 0;
    for (unsigned int u = 0; u < unknowns.size(); u++) {
        if ( !(unknowns[u]._name.compare("control")) ) index_control = u;
    }
   std::vector<std::string> act_set_flag_name(n_components_ctrl);
   act_set_flag_name[0] = "act_flag";

   ml_sol.AddSolution(act_set_flag_name[0].c_str(), unknowns[index_control]._fe_family, unknowns[index_control]._fe_order, act_set_fake_time_dep_flag, act_flag_is_an_unknown_of_a_pde);               
   ml_sol.Initialize(act_set_flag_name[0].c_str(), Solution_set_initial_conditions_Not_Unknowns, & ml_prob);
  // ******** active flag - END 

  // ******** state plus cont - BEGIN 
   std::vector<std::string> state_plus_ctrl(n_components_ctrl);
   state_plus_ctrl[0] = "state_plus_ctrl";
  
  const unsigned int state_plus_ctrl_fake_time_dep_flag = 0;
  const bool         state_plus_ctrl_is_an_unknown_of_a_pde = false;
   ml_sol.AddSolution(state_plus_ctrl[0].c_str(), unknowns[index_control]._fe_family, unknowns[index_control]._fe_order, state_plus_ctrl_fake_time_dep_flag, state_plus_ctrl_is_an_unknown_of_a_pde);               
   ml_sol.Initialize(state_plus_ctrl[0].c_str(), Solution_set_initial_conditions_Unknowns, & ml_prob);

  // ******** state plus cont - END 
   
   
  // ======= Solutions that are not Unknowns - END  ==================

  
  //== Solution: CHECK SOLUTION FE TYPES between Unknowns and Not Unknowns - BEGIN  ==
  if ( ml_sol.GetSolutionType("control") != ml_sol.GetSolutionType(act_set_flag_name[0].c_str())) abort();
  //== Solution: CHECK SOLUTION FE TYPES between Unknowns and Not Unknowns - END ==
 
  // ======= Problem, Solution - END ==================
  
  
  // ======= Problem, Quad Rule - BEGIN  ========================
    std::vector< std::string > fe_quad_rule_vec;
  fe_quad_rule_vec.push_back("seventh");
  
  ml_prob.SetQuadratureRuleAllGeomElemsMultiple(fe_quad_rule_vec);
  ml_prob.set_all_abstract_fe_AD_or_not();
  // ======= Problem, Quad Rule - END  ========================
  


  // ======= Problem, System - BEGIN ========================
  NonLinearImplicitSystemWithPrimalDualActiveSetMethod & system_opt = ml_prob.add_system < NonLinearImplicitSystemWithPrimalDualActiveSetMethod > ("LiftRestr"); //MU
  
  system_opt.SetActiveSetFlagName(act_set_flag_name); //MU

  for (unsigned int u = 0; u < unknowns.size(); u++)  { 
  system_opt.AddSolutionToSystemPDE(unknowns[u]._name.c_str());
  }
  
  
  system_opt.SetAssembleFunction(elliptic::lifting_internal::assemble_elliptic_dirichlet_control);
  
// *****************
  const unsigned n_components_state = n_components_ctrl;
  std::vector<std::string> state_vars(n_components_state);  state_vars[0] = "state";
  std::vector<std::string> ctrl_vars(n_components_ctrl);     ctrl_vars[0] = "control";
  
  system_opt.set_state_vars(state_vars);
  system_opt.set_ctrl_vars(ctrl_vars);
  
  system_opt.SetDebugNonlinear(true);
//   system_opt.SetDebugFunction(ctrl::cost_functional::compute_cost_functional_regularization_lifting_internal);
// *****************
//   system_opt.SetMaxNumberOfNonLinearIterations(2);

  // initialize and solve the system
  system_opt.init();
  
  system_opt.MGsolve();
//   system.assemble_call_before_boundary_conditions(1);
  
  // ======= Problem, System  - END ========================

  // ======= Post-processing - BEGIN ========================

  // ======= Post-processing, Computations - BEGIN ========================
  
  ml_sol.add_solution( ml_sol.GetIndex(state_vars[0].c_str()),   ml_sol.GetIndex(state_plus_ctrl[0].c_str()) );
  ml_sol.add_solution( ml_sol.GetIndex(ctrl_vars[0].c_str()),   ml_sol.GetIndex(state_plus_ctrl[0].c_str()) );
  


  
  femus::ctrl::cost_functional< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES,
                                femus::ctrl:: square_or_cube:: lifting_internal< femus::ctrl:: GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES > ,
                                femus::ctrl::COST_FUNCTIONAL_WITHOUT_REG >
                                ::compute_cost_functional_regularization_lifting_internal(ml_prob, 0, 0, state_vars, ctrl_vars, ALPHA_CTRL_VOL, BETA_CTRL_VOL, QRULE_I);
  
  femus::ctrl::cost_functional< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES,
                                femus::ctrl:: square_or_cube:: Domain_elements_containing_Gamma_control< femus::ctrl:: GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >,
                                femus::ctrl::COST_FUNCTIONAL_WITHOUT_REG >
  ::compute_cost_functional_regularization_bdry(ml_prob, 0, 0, state_plus_ctrl, ctrl_vars, ALPHA_CTRL_VOL, BETA_CTRL_VOL, QRULE_I);
  // ======= Post-processing, Computations - END ========================
  
  
  // ======= Post-processing, Print - BEGIN  ========================
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
  ml_sol.GetWriter()->Write(files.GetOutputPath(), "biquadratic", variablesToBePrinted);
  // ======= Post-processing, Print - END  ========================

  // ======= Post-processing - END ========================

  // ======= Problem - END  ==================

  return 0;
  
}



