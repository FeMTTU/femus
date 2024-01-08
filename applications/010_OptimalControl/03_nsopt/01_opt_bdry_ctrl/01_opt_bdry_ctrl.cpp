#include "adept.h"

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NumericVector.hpp"
#include "SparseMatrix.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "Fluid.hpp"
#include "Parameter.hpp"
#include "Files.hpp"
#include "PetscMatrix.hpp"
#include "Parallel.hpp"//to get iproc HAVE_MPI is inside here

#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"

#include "ElemType.hpp"

#include "opt_common.hpp"
#include "00_cost_functional.hpp"
#include "01_opt_system.hpp"
#include "03_opt_system_inequalities.hpp"

#include   "../manufactured_solutions.hpp"

#include  "../../opt_systems_boundary_control_eqn_sobolev_integer.hpp"


#include   "../nsopt_params.hpp"

#include  "../../opt_systems_ns_dirichlet.hpp"

#include  "../../square_or_cube_03_cost_functional_without_regularization.hpp"





using namespace femus;

 
double Solution_set_initial_conditions_Unknowns(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;

    return value;
    
}


double Solution_set_initial_conditions_Not_Unknowns(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;

     if(!strcmp(name,"TargReg")) {
        value = femus::ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(x);
    }

    return value;
    
}


bool Solution_set_boundary_conditions_Unknowns(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char SolName[], double& value, const int faceName, const double time) {
  //1: bottom  //2: right  //3: top  //4: left
  
  bool dirichlet = true;
   value = 0.;
 

                if (!strcmp(SolName, "ctrl_0"))       {
//                     if (facename == FACE_FOR_CONTROL) dirichlet = false; 
  if (faceName == FACE_FOR_CONTROL) {
     if (x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] > SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN - 1.e-5 &&
         x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] < SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX + 1.e-5)  { 
         dirichlet = false;
    }
     else { 
         dirichlet = true;  
    }
  }
  else { 
      dirichlet = true;
   }
                }
                
                
                
           else if (!strcmp(SolName, "ctrl_1"))       { 
//                     if (facename == FACE_FOR_CONTROL) dirichlet = false; 
  if (faceName == FACE_FOR_CONTROL) {
     if (x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] > SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN - 1.e-5 &&
         x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] < SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX + 1.e-5)  { 
         dirichlet = false;
    }
     else { 
         dirichlet = true;  
    }
  }
  else { 
      dirichlet = true;
   }
                }
                
                
                
           else if (!strcmp(SolName, "ctrl_2"))       { 
//                     if (facename == FACE_FOR_CONTROL) dirichlet = false; 
  if (faceName == FACE_FOR_CONTROL) {
     if (x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] > SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN - 1.e-5 &&
         x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] < SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX + 1.e-5)  { 
         dirichlet = false;
    }
     else { 
         dirichlet = true;  
    }
  }
  else { 
      dirichlet = true;
   }
                }
                
                
                
                
                else if (!strcmp(SolName, "theta"))    { dirichlet = false; }

                
                
                else if (!strcmp(SolName, "u_2"))       { 
//                     if (facename == FACE_FOR_CONTROL) dirichlet = false; 
   if (faceName == FACE_FOR_CONTROL) {
     if (x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] > SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN - 1.e-5 &&
         x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] < SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX + 1.e-5)  { 
         dirichlet = false;
    }
     else { 
         dirichlet = true;  
    }
  }                   
                }
                
                
#if INEQ_FLAG != 0
   //MU
       if (!strcmp(SolName, "mu_0"))    { dirichlet = false; }
  else if (!strcmp(SolName, "mu_1"))    { dirichlet = false; } 
  else if (!strcmp(SolName, "mu_2"))    { dirichlet = false; } 
#endif  
     

                
                
#if exact_sol_flag == 0
                else if (!strcmp(SolName, "u_0"))       { if (faceName == FACE_FOR_CONTROL) dirichlet = false; }
                else if (!strcmp(SolName, "u_1"))       { if (faceName == FACE_FOR_CONTROL) dirichlet = false; }
#endif
     
#if exact_sol_flag == 1
  //b.c. for manufactured lid driven cavity
  double pi = acos(-1.);
                else if (!strcmp(SolName, "u_0"))       { if (faceName == FACE_FOR_CONTROL) value =   sin(pi* x[0]) * sin(pi* x[0]) * cos(pi* x[1]) - sin(pi* x[0]) * sin(pi* x[0]); }
                else if (!strcmp(SolName, "u_1"))       { if (faceName == FACE_FOR_CONTROL) value = - sin(2. * pi * x[0]) * sin(pi* x[1]) + pi * x[1] * sin(2. * pi * x[0]); }
 #endif
               
  return dirichlet;

}








int main(int argc, char** args) {
  
  // ======= Init ========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
  
  
  // ======= Problem  ==================
  MultiLevelProblem ml_prob; 

  
  // ======= Files - BEGIN  ========================
  Files files; 
        files.CheckIODirectories(ctrl::use_output_time_folder);
        files.RedirectCout(ctrl::redirect_cout_to_file);
  
  // ======= Problem, Files ========================
  ml_prob.SetFilesHandler(&files);
  // ======= Files - END  ========================
  
  // ======= Messages - BEGIN  ========================
  std::cout << "************ The variable for compatibility condition only works if you put it as the LAST VARIABLE" << std::endl;
  // ======= Messages - END    ========================
  
  
  // ======= Parameters - BEGIN  ==================
   //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;
 // *** apparently needed by non-AD assemble only **********************
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
  fe_quad_rule_vec.push_back("eighth");
    
  ml_prob.SetQuadratureRuleAllGeomElemsMultiple(fe_quad_rule_vec);
  ml_prob.set_all_abstract_fe_AD_or_not();
  // ======= Problem, Quad Rule - END  ========================


  
  // ======= Mesh, Coarse reading - BEGIN ==================
  MultiLevelMesh ml_mesh;

  const std::string relative_path_to_mesh_input_folder = "../../../";
  const std::string infile = Files::get_input_file_with_prefix(ctrl::mesh_input, relative_path_to_mesh_input_folder);


  const bool read_groups = true;
  const bool read_boundary_groups = true;
    
    ml_mesh.ReadCoarseMeshFileReadingBeforePartitioning(infile, Lref, read_groups, read_boundary_groups);
    
           ml_mesh.GetLevelZero(0)->dofmap_all_fe_families_initialize();

           std::vector < unsigned > elem_partition_from_mesh_file_to_new = ml_mesh.GetLevelZero(0)->elem_offsets();
  
           std::vector < unsigned > node_mapping_from_mesh_file_to_new = ml_mesh.GetLevelZero(0)->node_offsets();
           
           ml_mesh.GetLevelZero(0)->dofmap_all_fe_families_clear_ghost_dof_list_for_other_procs();
       
           ml_mesh.GetLevelZero(0)->BuildElementAndNodeStructures();
 


  ml_mesh.BuildFETypesBasedOnExistingCoarseMeshGeomElements();
  
  ml_mesh.PrepareNewLevelsForRefinement();

  // ======= Mesh, Coarse reading - END ==================
    

  // ======= Convergence Rate, Preparation - BEGIN  ==================
    
  unsigned dim = ml_mesh.GetDimension();
  unsigned maxNumberOfMeshes;

  if (dim == 2) {
    maxNumberOfMeshes = N_UNIFORM_LEVELS;
  } else {
    maxNumberOfMeshes = 2;
  }
  
  
  
   // ======= Solutions that are Unknowns - BEGIN  ==================
  std::vector< Unknown > unknowns = navier_stokes::pure_boundary:: provide_list_of_unknowns( dim );
   // ======= Solutions that are Unknowns - END  ==================

  
  
#if compute_conv_flag == 1
  MultiLevelMesh ml_mesh_all_levels;
  
     ml_mesh_all_levels.ReadCoarseMesh(infile.c_str(),fe_quad_rule_vec[0].c_str(),Lref);
 
     double comp_conv[maxNumberOfMeshes][pure_boundary_norms::no_of_l2_norms +  pure_boundary_norms::no_of_h1_norms ];
 
  
        unsigned numberOfUniformLevels_finest = maxNumberOfMeshes;
        ml_mesh_all_levels.RefineMesh(numberOfUniformLevels_finest, numberOfUniformLevels_finest, NULL);
//      ml_mesh_all_levels.EraseCoarseLevels(numberOfUniformLevels - 2);  // need to keep at least two levels to send u_(i-1) projected(prolongated) into next refinement
        
        //store the fine solution  ==================
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
   ml_sol_all_levels->InitializeBasedOnControlFaces< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >(volume_control_region.c_str(),  Solution_set_initial_conditions_Not_Unknowns, & ml_prob);
  // ======= Solutions that are not Unknowns - END  ==================
            
            
#endif

  // ======= Convergence Rate, Preparation - END  ==================

   
   
         for (int i = /*0*/maxNumberOfMeshes - 1; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

  // ======= Mesh: Refinement - BEGIN ==================
  unsigned numberOfUniformLevels = i + 1; 
  const unsigned erased_levels = numberOfUniformLevels - 1;
  unsigned numberOfSelectiveLevels = 0;
  
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  // ======= Mesh: Refinement - END ==================

  // ======= Mesh: Group info - BEGIN ==================
   ml_mesh.set_group_info(relative_path_to_mesh_input_folder);
  // ======= Mesh: Group info - END ==================

  // ======= Solution, auxiliary; needed for Boundary of Boundary of Control region - BEFORE COARSE ERASING - BEGIN  ==================
  const std::string node_based_bdry_bdry_flag_name = femus::pure_boundary::_node_based_bdry_bdry;
  const unsigned  steady_flag = 0;
  const bool      is_an_unknown_of_a_pde = false;
  
  const FEFamily node_bdry_bdry_flag_fe_fam = LAGRANGE;
  const FEOrder node_bdry_bdry_flag_fe_ord = SECOND;
  
  MultiLevelSolution * ml_sol_bdry_bdry_flag = ctrl::Gamma_control_equation_fractional_sobolev_differentiability_index<
                femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES, 
                femus::ctrl:: square_or_cube:: pure_boundary< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >
                >
::bdry_bdry_flag(files,
                                                              ml_mesh, 
                                                              infile,
                                                              node_mapping_from_mesh_file_to_new,
                                                              node_based_bdry_bdry_flag_name,
                                                              steady_flag,
                                                              is_an_unknown_of_a_pde,
                                                              node_bdry_bdry_flag_fe_fam,
                                                              node_bdry_bdry_flag_fe_ord);
  
  // ======= Solution, auxiliary - END  ==================
  
  
  
  // ======= Mesh: COARSE ERASING - BEGIN  ========================
  ml_mesh.EraseCoarseLevels(erased_levels);
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

  const std::string volume_control_region = "ContReg";
   ml_sol.AddSolution(volume_control_region.c_str(),  DISCONTINUOUS_POLYNOMIAL, ZERO);
   ml_sol.InitializeBasedOnControlFaces< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >(volume_control_region.c_str(),  Solution_set_initial_conditions_Not_Unknowns, & ml_prob);
   
   
 ctrl::Gamma_control_equation_fractional_sobolev_differentiability_index<
                femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES, 
                femus::ctrl:: square_or_cube:: Domain_elements_containing_Gamma_control< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >
                >
::bdry_bdry_flag_copy_and_delete(ml_prob,
                                ml_sol,
                                ml_mesh, 
                                erased_levels,
                                ml_sol_bdry_bdry_flag,
                                node_based_bdry_bdry_flag_name,
                                steady_flag,
                                is_an_unknown_of_a_pde,
                                node_bdry_bdry_flag_fe_fam,
                                node_bdry_bdry_flag_fe_ord);

 
 
#if INEQ_FLAG != 0
  // ******** active flag - BEGIN   //MU
  const unsigned int  n_components_ctrl = dim;
  
  const bool      act_flag_is_an_unknown_of_a_pde = false;

  unsigned int index_control = 0;
    for (unsigned int u = 0; u < unknowns.size(); u++) {
        if ( !(unknowns[u]._name.compare("ctrl_0")) ) index_control = u;
    }
  const unsigned int act_set_fake_time_dep_flag = 2;  //this is needed to be able to use _SolOld
   std::vector<std::string> act_set_flag_name(n_components_ctrl);
   
   for(unsigned int d = 0; d <  act_set_flag_name.size(); d++)  {
       act_set_flag_name[d] = "act_flag_" + std::to_string(d);
    ml_sol.AddSolution(act_set_flag_name[d].c_str(), unknowns[index_control]._fe_family, unknowns[index_control]._fe_order, act_set_fake_time_dep_flag, act_flag_is_an_unknown_of_a_pde);
    ml_sol.Initialize(act_set_flag_name[d].c_str(), Solution_set_initial_conditions_Not_Unknowns, & ml_prob);
   }
  // ******** active flag - END 
#endif
 
  // ======= Solutions that are not Unknowns - END  ==================

 
  // ======= Problem, System - BEGIN ========================
  NonLinearImplicitSystemWithPrimalDualActiveSetMethod & system_opt    = ml_prob.add_system < NonLinearImplicitSystemWithPrimalDualActiveSetMethod > ("NSOpt");
  
#if INEQ_FLAG != 0
  system_opt.SetActiveSetFlagName(act_set_flag_name); //MU
#endif

  for (unsigned int u = 0; u < unknowns.size(); u++)  { 
  system_opt.AddSolutionToSystemPDE(unknowns[u]._name.c_str());
  }
   
  
   system_opt.SetAssembleFunction(navier_stokes::pure_boundary::assemble_ns_dirichlet_control_pure_boundary);

    
// *****************
 
  std::vector< unsigned >  vector_offsets = navier_stokes::pure_boundary:: provide_state_adj_ctrl_mu_offsets(dim);


const int state_pos_begin   =  vector_offsets[navier_stokes::pure_boundary::pos_index_state];
 const int ctrl_pos_begin   =  vector_offsets[navier_stokes::pure_boundary::pos_index_ctrl];
   
   
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
  
    system_opt.SetDebugNonlinear(true);
    
//     system_opt.SetDebugFunction(ctrl::cost_functional::compute_cost_functional_regularization_bdry_vec);
    
// *****************
    
   
  system_opt.init();
  set_dense_pattern_for_unknowns(system_opt, unknowns);

   
//   initialize and solve the system
    system_opt.init();
  
    
      //----
  std::ostringstream sp_out_base2; sp_out_base2 << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "sp_";
  sp_out_base2 << "after_second_init_";
  
  unsigned n_levels = ml_mesh.GetNumberOfLevels();
  system_opt._LinSolver[n_levels - 1]->sparsity_pattern_print_nonzeros(sp_out_base2.str(), "on");
  system_opt._LinSolver[n_levels - 1]->sparsity_pattern_print_nonzeros(sp_out_base2.str(), "off");
  //----

    
//   system_opt.SetMaxNumberOfNonLinearIterations(2);
//   system_opt.SetNonLinearConvergenceTolerance(1.e-15);
//     system_opt.SetAbsoluteLinearConvergenceTolerance(1.e-14);
//     system_opt.SetOuterSolver(PREONLY);
    
   
    system_opt.MGsolve();
//   system_opt.assemble_call_before_boundary_conditions(1);
  // ======= Problem, System  - END ========================

  
  
  // ======= Convergence Rate - BEGIN  ==================
#if compute_conv_flag == 1
    if ( i > 0 ) {
        
//prolongation of coarser  
      ml_sol_all_levels->RefineSolution(i);
      Solution* sol_coarser_prolongated = ml_sol_all_levels->GetSolutionLevel(i);
  
      double* norm = GetErrorNorm(ml_prob,&ml_sol,sol_coarser_prolongated);
    
      for(int j = 0; j <  pure_boundary_norms::no_of_l2_norms + pure_boundary_norms::no_of_h1_norms ; j++)       comp_conv[i-1][j] = norm[j];
  
     }

    
//store the last computed solution
// 
       const unsigned level_index_current = 0;
      //@todo there is a duplicate function in MLSol: GetSolutionLevel() and GetLevel()
       const unsigned n_vars_sol = ml_sol.GetSolutionLevel(level_index_current)->_Sol.size();
       
        for(unsigned short j = 0; j < n_vars_sol; j++) {  
               *(ml_sol_all_levels->GetLevel(i)->_Sol[j]) = *(ml_sol.GetSolutionLevel(level_index_current)->_Sol[j]);
        }
 #endif
  // ======= Convergence Rate - END  ==================
       
 
  // ======= Print - BEGIN  ========================
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  ml_sol.GetWriter()->Write(files.GetOutputPath(),"biquadratic", variablesToBePrinted, i);
  // ======= Print - END  ========================

  }

  // ======= Convergence Rate, Postprocess - BEGIN  ==================

#if compute_conv_flag == 1
  std::cout << "=======================================================================" << std::endl;
   std::cout << " L2-NORM ERROR and ORDER OF CONVERGENCE:\n\n";
  std::vector< std::string > norm_names_L2 = {"u_0","u_1", "u_p", "adj_0","adj_1", "adj_p", "ctrl_0","ctrl_1", "theta"};

   for(int j = 0; j <  norm_names_L2.size(); j++)  {
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "LEVEL\t\t" << norm_names_L2[j] << "\t\t\t\torder of convergence\n"; 
   for(int i = 0; i <  maxNumberOfMeshes - 1; i++){
       FE_convergence<double>::output_convergence_rate(comp_conv[i][j], comp_conv[i + 1][j], norm_names_L2[j], maxNumberOfMeshes , i );
    }
  }
  std::cout << std::endl;
  std::cout << "=======================================================================" << std::endl;
  std::cout << " H1-NORM ERROR and ORDER OF CONVERGENCE:" << std::endl;
  std::vector< std::string > norm_names_H1 = {"u_0","u_1", "adj_0","adj_1", "ctrl_0","ctrl_1"};

   for(int j = 0; j <  norm_names_H1.size(); j++)  {
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "LEVEL\t\t" << norm_names_H1[j] << "\t\t\t\torder of convergence\n"; 
   for(int i = 0; i <  maxNumberOfMeshes - 1; i++){
       FE_convergence<double>::output_convergence_rate(comp_conv[i][ pure_boundary_norms::no_of_l2_norms  + j], comp_conv[i + 1][ pure_boundary_norms::no_of_l2_norms  + j], norm_names_H1[j], maxNumberOfMeshes , i );
    }
  }
  std::cout << std::endl;
  std::cout << "=======================================================================" << std::endl;
#endif
  // ======= Convergence Rate, Postprocess - END  ==================

  
  return 0;
}
 

