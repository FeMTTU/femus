// Solving Navier-Stokes problem using automatic differentiation and/or Picards method
// boundary conditions were set in 2D as, no slip in left,right of the box and top to bottom gravity is enforced
// therefore, U=V=0 on left and right, U=0 on top and bottom, V is free 

#include "adept.h"

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "Files.hpp"
#include "Parameter.hpp"
#include "Fluid.hpp"
#include "NumericVector.hpp"
#include "SparseMatrix.hpp"
#include "Parallel.hpp"//to get iproc HAVE_MPI is inside here

#include "Assemble_jacobian.hpp"

#include "00_cost_functional.hpp"
#include "03_opt_system_inequalities.hpp"



#include   "../manufactured_solutions.hpp"


#include   "../nsopt_params.hpp"

#include  "../../opt_systems_ns_dirichlet.hpp"


#include  "../../square_or_cube_03_cost_functional_without_regularization.hpp"




  
  

using namespace femus;



 
 
bool Solution_set_boundary_conditions_Unknowns(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char SolName[], double& value, const int faceName, const double time) {
  //1: bottom  //2: right  //3: top  //4: left  (2D square)
  //1: bottom  //2: top    //3: side            (3D cylinder)
    
  
  bool dirichlet = true;
   value = 0.;

#if exact_sol_flag == 0
// b.c. for lid-driven cavity problem, wall u_top = 1 = shear_force, v_top = 0 and u=v=0 on other 3 walls ; rhs_f = body_force = {0,0}

   if (faceName == FACE_FOR_CONTROL)  {
       
        if (x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] > SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN - 1.e-5 && 
            x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] < SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX + 1.e-5)  {
            
       if (!strcmp(SolName, "ctrl_0"))    { dirichlet = false; }
  else if (!strcmp(SolName, "ctrl_1"))    { dirichlet = false; } 
  else if (!strcmp(SolName, "ctrl_2"))    { dirichlet = false; }
  
              }
              else {
       if (!strcmp(SolName, "ctrl_0"))    { dirichlet = true; }
  else if (!strcmp(SolName, "ctrl_1"))    { dirichlet = true; } 
  else if (!strcmp(SolName, "ctrl_2"))    { dirichlet = true; } 
              }
      }
      
   
       if (!strcmp(SolName, "mu_0"))    { dirichlet = false; }
  else if (!strcmp(SolName, "mu_1"))    { dirichlet = false; } 
  else if (!strcmp(SolName, "mu_2"))    { dirichlet = false; } 
     
      
#endif

#if exact_sol_flag == 1
  //b.c. for manufactured lid driven cavity
// TOP ==========================  
   double pi = acos(-1.);
     if (faceName == FACE_FOR_CONTROL) {
       if (!strcmp(SolName, "ctrl_0"))    { value =   sin(pi* x[0]) * sin(pi* x[0]) * cos(pi* x[1]) - sin(pi* x[0]) * sin(pi* x[0]);} //lid - driven
  else if (!strcmp(SolName, "ctrl_1"))    { value = - sin(2. * pi * x[0]) * sin(pi* x[1]) + pi * x[1] * sin(2. * pi * x[0]);} 
  	
      }
#endif

      
  return dirichlet;
}


double Solution_set_initial_conditions_Unknowns(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;

    return value;
}


double Solution_set_initial_conditions_Not_Unknowns(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;

     if(!strcmp(name,"TargReg")) {
        value = femus::ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(x);
    }
    else if(!strcmp(name,"ContReg")) {
        value = femus::ctrl:: square_or_cube:: lifting_internal< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >::ControlDomainFlag_internal_restriction(x);
    }

    return value;
}












int main(int argc, char** args) {



  // ======= Init ========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

   // ======= Problem  ==================
  MultiLevelProblem ml_prob;

  
  // ======= Files - BEGIN  ========================
  Files files; 
        files.CheckIODirectories(femus::ctrl::use_output_time_folder);
	    files.RedirectCout(femus::ctrl::redirect_cout_to_file);
 
  // ======= Problem, Files ========================
  ml_prob.SetFilesHandler(&files);
  // ======= Files - END  ========================

  
  // ======= Parameters - BEGIN  ==================
  double Lref = 1.;
  double Uref = 1.;
  // add fluid material
  Parameter parameter(Lref, Uref);
  
       // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  Fluid fluid(parameter, 1, FLUID_DENSITY, "Newtonian");
  std::cout << "Fluid properties: " << std::endl;
  std::cout << fluid << std::endl;
  
  // ======= Problem, Parameters ========================
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;
  // ======= Parameters - END  ==================

  
  // ======= Problem, Quad Rule - BEGIN  ========================
    std::vector< std::string > fe_quad_rule_vec;
  fe_quad_rule_vec.push_back("seventh");
  
  ml_prob.SetQuadratureRuleAllGeomElemsMultiple(fe_quad_rule_vec);
  ml_prob.set_all_abstract_fe_AD_or_not();
  // ======= Problem, Quad Rule - END  ========================
  
  
  // ======= Mesh, Coarse reading - BEGIN ==================
  MultiLevelMesh ml_mesh;
 
    std::string mesh_folder_file = "input/";
  const std::string input_file = femus::ctrl::mesh_input;

//   std::string input_file = "square_0-1x0-1_divisions_1x2.med";
//   std::string input_file = "cyl.med"; // "fifth"
//   std::string input_file = "Mesh_3_groups_with_bdry_nodes_coarser.med";
  std::ostringstream mystream; mystream << "./" << mesh_folder_file << input_file;
  const std::string infile = mystream.str();

  
  const bool read_groups = true;
  const bool read_boundary_groups = true;
    
  ml_mesh.ReadCoarseMesh(infile, Lref, read_groups, read_boundary_groups);
  // ======= Mesh, Coarse reading - END ==================

  
  // ======= Convergence Rate, Preparation - BEGIN  ==================

  unsigned dim = ml_mesh.GetDimension();
  unsigned maxNumberOfMeshes;

  if (dim == 2) {
    maxNumberOfMeshes = N_UNIFORM_LEVELS;
  } else {
    maxNumberOfMeshes = N_UNIFORM_LEVELS;
  }

 
 
   // ======= Solutions that are Unknowns - BEGIN  ==================
  std::vector< Unknown > unknowns = navier_stokes::lifting_internal:: provide_list_of_unknowns( dim );
   // ======= Solutions that are Unknowns - END  ==================

  

#if compute_conv_flag == 1
     double comp_conv[maxNumberOfMeshes][ lifting_internal_norms::no_of_l2_norms + lifting_internal_norms::no_of_h1_norms ];
 
  // ======= Mesh ==================
  MultiLevelMesh ml_mesh_all_levels;
   ml_mesh_all_levels.ReadCoarseMesh(infile.c_str(),fe_quad_rule_vec[0].c_str(),Lref);
  
        unsigned numberOfUniformLevels_finest = maxNumberOfMeshes;
        ml_mesh_all_levels.RefineMesh(numberOfUniformLevels_finest, numberOfUniformLevels_finest, NULL);
//      ml_mesh_all_levels.EraseCoarseLevels(numberOfUniformLevels - 2);  // need to keep at least two levels to send u_(i-1) projected(prolongated) into next refinement
        
  // ======= Solution  ==================
            MultiLevelSolution * ml_sol_all_levels;
            ml_sol_all_levels = new MultiLevelSolution (& ml_mesh_all_levels);  //with the declaration outside and a "new" inside it persists outside the loop scopes

   // ======= Solutions that are Unknowns - BEGIN  ==================
  for (unsigned int u = 0; u < unknowns.size(); u++)  { 
      ml_sol_all_levels->AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order, unknowns[u]._time_order, unknowns[u]._is_pde_unknown);
   }
   
 for (unsigned int u = 0; u < unknowns.size(); u++)  { 
      ml_sol_all_levels->Initialize(unknowns[u]._name.c_str(), Solution_set_initial_conditions_Unknowns, & ml_prob);
  }
  
      ml_sol_all_levels->AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions_Unknowns);
   for (unsigned int u = 0; u < unknowns.size(); u++)  { 
     ml_sol_all_levels->GenerateBdc(unknowns[u]._name.c_str(), (unknowns[u]._time_order == 0) ? "Steady" : "Time_dependent", & ml_prob);
  }
   // ======= Solutions that are Unknowns - END  ==================


  // ======= Solutions that are not Unknowns - BEGIN  ==================
     ml_sol_all_levels->AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO);
     ml_sol_all_levels->Initialize("TargReg",     Solution_set_initial_conditions_Not_Unknowns, & ml_prob);
     
     ml_sol_all_levels->AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO);
     ml_sol_all_levels->Initialize("ContReg",     Solution_set_initial_conditions_Not_Unknowns, & ml_prob);
  // ======= Solutions that are not Unknowns - END  ==================

#endif

  // ======= Convergence Rate, Preparation - END  ==================

     
            
         for (int i = /*0*/maxNumberOfMeshes - 1; i < /*1*/maxNumberOfMeshes; i++) {   // loop on the mesh level

  // ======= Mesh: Refinement - BEGIN ==================
  unsigned numberOfUniformLevels = i + 1; 
  unsigned numberOfSelectiveLevels = 0;
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  // ======= Mesh: Refinement - END ==================

  // ======= Mesh: COARSE ERASING - BEGIN  ========================
  ml_mesh.EraseCoarseLevels(numberOfUniformLevels - 1);
  ml_mesh.PrintInfo();
  // ======= Mesh: COARSE ERASING - END  ========================

  
  // ======= Solution - BEGIN ==================
  MultiLevelSolution ml_sol(&ml_mesh);
  
  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);
  
  // ======= Problem, Mesh and Solution  ==================
  ml_prob.SetMultiLevelMeshAndSolution(& ml_sol);
  // ======= Solution - END ==================
 

  // ======= Solutions that are Unknowns - BEGIN  ==================
  for (unsigned int u = 0; u < unknowns.size(); u++)  { 
      ml_sol.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order, unknowns[u]._time_order, unknowns[u]._is_pde_unknown);
   }
   
 for (unsigned int u = 0; u < unknowns.size(); u++)  { 
      ml_sol.Initialize(unknowns[u]._name.c_str(), Solution_set_initial_conditions_Unknowns, & ml_prob);
  }
  
  ml_sol.AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions_Unknowns);
   for (unsigned int u = 0; u < unknowns.size(); u++)  { 
     ml_sol.GenerateBdc(unknowns[u]._name.c_str(), (unknowns[u]._time_order == 0) ? "Steady" : "Time_dependent", & ml_prob);
  }
  // ======= Solutions that are Unknowns - END  ==================


  // ======= Solutions that are not Unknowns - BEGIN  ==================
  ml_sol.AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO);
  ml_sol.Initialize("TargReg",     Solution_set_initial_conditions_Not_Unknowns, & ml_prob);

  ml_sol.AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO);
  ml_sol.Initialize("ContReg",     Solution_set_initial_conditions_Not_Unknowns, & ml_prob);

  // ******** active flag - BEGIN 
  const unsigned int  n_components_ctrl = dim;
  
  const unsigned int act_set_fake_time_dep_flag = 2;  //this is needed to be able to use _SolOld  //MU
  const bool      act_flag_is_an_unknown_of_a_pde = false;

  unsigned int index_control = 0;
    for (unsigned int u = 0; u < unknowns.size(); u++) {
        if ( !(unknowns[u]._name.compare("ctrl_0")) ) index_control = u;
    }
   std::vector<std::string> act_set_flag_name(n_components_ctrl);
   
   for(unsigned int d = 0; d <  act_set_flag_name.size(); d++)  {
       act_set_flag_name[d] = "act_flag_" + std::to_string(d);
    ml_sol.AddSolution(act_set_flag_name[d].c_str(), unknowns[index_control]._fe_family, unknowns[index_control]._fe_order, act_set_fake_time_dep_flag, act_flag_is_an_unknown_of_a_pde);
    ml_sol.Initialize(act_set_flag_name[d].c_str(), Solution_set_initial_conditions_Not_Unknowns, & ml_prob);
   }
  // ******** active flag - END 
   
// ======= Solutions that are not Unknowns - END  ==================
  

  // ======= Problem, System - BEGIN ========================
  NonLinearImplicitSystemWithPrimalDualActiveSetMethod & system_opt    = ml_prob.add_system < NonLinearImplicitSystemWithPrimalDualActiveSetMethod > ("NSOpt");  ///@todo this MUST return a REFERENCE, otherwise it doesn't run!
  
  system_opt.SetActiveSetFlagName(act_set_flag_name); //MU


  for (unsigned int u = 0; u < unknowns.size(); u++)  { 
  system_opt.AddSolutionToSystemPDE(unknowns[u]._name.c_str());
  }
  
//   system_opt.SetAssembleFunction(navier_stokes::lifting_internal:: assemble_ns_dirichlet_control_lifting_internal_AD);  //AD doesn't seem to work now
  system_opt.SetAssembleFunction(navier_stokes::lifting_internal:: assemble_ns_dirichlet_control_lifting_internal_nonAD);
    

 
// *****************
  std::vector< unsigned >  vector_offsets = navier_stokes::lifting_internal::provide_state_adj_ctrl_mu_offsets(dim);


const int state_pos_begin   =  vector_offsets[navier_stokes::lifting_internal::pos_index_state];
 const int ctrl_pos_begin   =  vector_offsets[navier_stokes::lifting_internal::pos_index_ctrl];
   
   
   const unsigned n_components_state = n_components_ctrl;
  std::vector<std::string> state_vars(n_components_state);
  
   for (unsigned int u = 0; u < state_vars.size(); u++)  { 
   state_vars[u] = unknowns[state_pos_begin + u]._name;
    }
    
   std::vector<std::string> ctrl_vars(n_components_ctrl);     
   for (unsigned int u = 0; u < ctrl_vars.size(); u++)  { 
   ctrl_vars[u] = unknowns[ctrl_pos_begin + u]._name;
    }
  
  

  system_opt.set_state_vars(state_vars);
  system_opt.set_ctrl_vars(ctrl_vars);

  
  // *****************
  system_opt.SetDebugNonlinear(true);
//     system_opt.SetDebugFunction(femus::ctrl::cost_functional::compute_cost_functional_regularization_lifting_internal_vec);
// *****************
  
  
  // initialize and solve the system
  system_opt.init();
 
//   system_opt.SetMaxNumberOfNonLinearIterations(30);
//   system_opt.SetNonLinearConvergenceTolerance(1.e-15);
//   system_opt.SetDebugLinear(true);
//   system_opt.SetMaxNumberOfLinearIterations(6);
//   system_opt.SetAbsoluteLinearConvergenceTolerance(1.e-14);
//   system_opt.SetOuterSolver(PREONLY);

  system_opt.MGsolve();
  // ======= Problem, System - END ========================

  
  
#if compute_conv_flag == 1
    if ( i > 0 ) {
        
//prolongation of coarser  
      ml_sol_all_levels->RefineSolution(i);
      Solution* sol_coarser_prolongated = ml_sol_all_levels->GetSolutionLevel(i);
  
      double* norm = GetErrorNorm(ml_prob,&ml_sol,sol_coarser_prolongated);
    
      for(int j = 0; j <  lifting_internal_norms::no_of_l2_norms + lifting_internal_norms::no_of_h1_norms ; j++)       comp_conv[i-1][j] = norm[j];
  
     }

    
//store the last computed solution
// 
       const unsigned level_index_current = 0;
      //@todo there is a duplicate function in MLSol: GetSolutionLevel() and GetLevel()
       const unsigned n_vars = ml_sol.GetSolutionLevel(level_index_current)->_Sol.size();
       
        for(unsigned short j = 0; j < n_vars; j++) {  
               *(ml_sol_all_levels->GetLevel(i)->_Sol[j]) = *(ml_sol.GetSolutionLevel(level_index_current)->_Sol[j]);
        }
 #endif
       
   
  // ======= Print - BEGIN  ========================
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  ml_sol.GetWriter()->Write(files.GetOutputPath(),  "biquadratic", variablesToBePrinted, i);
  // ======= Print - END  ========================

 }


#if compute_conv_flag == 1
  std::cout << "=======================================================================" << std::endl;
  std::cout << " L2-NORM ERROR and ORDER OF CONVERGENCE:\n\n";
   std::vector< std::string > norm_names_L2 = {"u_0","u_1", "p_u", "adj_0", "adj_1", "p_adj", "ctrl_0", "ctrl_1", "p_ctrl", "Vel_X" , "Vel_Y"};

   for(int j = 0; j <  norm_names_L2.size(); j++)  {
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "LEVEL\t\t\t" << norm_names_L2[j] << "\t\t\t\torder of convergence\n"; 
   for(int i = 0; i <  maxNumberOfMeshes - 1; i++){
       FE_convergence<double>::output_convergence_rate(comp_conv[i][j], comp_conv[i + 1][j], norm_names_L2[j], maxNumberOfMeshes , i );
    }
  }
  std::cout << std::endl;
  std::cout << "=======================================================================" << std::endl;
  std::cout << " H1-NORM ERROR and ORDER OF CONVERGENCE:" << std::endl;
  std::vector< std::string > norm_names_H1 = {"u_0","u_1", "adj_0","adj_1", "ctrl_0","ctrl_1", "Vel_X" , "Vel_Y"};

   for(int j = 0; j <  norm_names_H1.size(); j++)  {
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "LEVEL\t\t\t" << norm_names_H1[j] << "\t\t\t\torder of convergence\n"; 
   for(int i = 0; i <  maxNumberOfMeshes - 1; i++){
       FE_convergence<double>::output_convergence_rate(comp_conv[i][ lifting_internal_norms::no_of_l2_norms  + j], comp_conv[i + 1][ lifting_internal_norms::no_of_l2_norms  + j], norm_names_H1[j], maxNumberOfMeshes , i );
    }
  }
  std::cout << std::endl;
  std::cout << "=======================================================================" << std::endl;
#endif
 
  return 0;
}


