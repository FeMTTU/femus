#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"


#include  "00_cost_functional.hpp"
#include "01_opt_system.hpp"


#include  "../../../param.hpp"
#include  "../../../opt_systems_elliptic_dirichlet.hpp"

using namespace femus;








///@todo do a very weak impl of Laplacian
///@todo Review the ordering for phi_ctrl_x_bdry
///@todo check computation of 2nd derivatives in elem_type_template
///@todo Implement rather fast way to add inequality constraint to a problem
///@todo merge elliptic_nonlin in here
///@todo What if I did a Point domain, could I solve ODEs in time like this? :)
///@todo Re-double check that things are fine in elem_type_template, probably remove _gauss_bdry!
///@todo See if with Petsc you can enforce Dirichlet conditions using NEGATIVE indices
///@todo Remove the prints, possible cause of slowing down (maybe do assert)
///@todo The \mu/actflag pieces are now basically separated, except for setting to zero on Omega minus Gamma_c (such as is done for control)
///@todo put assembleMatrix everywhere there is a filling of the matrix!
///@todo Give the option to provide your own name to the run folder instead of the time instant. I think I did something like this when running with the external script already


           // 1 scalar weak Galerkin variable will first have element-based nodes of a certain order.
           //I will loop over the elements and take all the node dofs either of order 1 or 2, counted with repetition
           //Then I have to take the mesh skeleton (without repetition)
           //Then for the dofs on the edges how do I do? 
           // In every subdomain I will have nelems x element nodes + n skeleton dofs in that subdomain 
           // Then, when it comes to retrieving such dofs for each element, i'll retrieve the interior element nodes + the boundary dofs
           



 // Unknowns - BEGIN  ==================
double Solution_set_initial_conditions_Unknowns(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;

    return value;
    
}

///@todo notice that even if you set Dirichlet from the mesh file, here you can override it
bool Solution_set_boundary_conditions_Unknowns(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet = false; // = true; //dirichlet
  value = 0.;

  //************************state****************************************************

  if(!strcmp(name,"state")) {  //"state" corresponds to the first block row (u = q)

     ctrl:: NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS < ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >::ctrl_or_state_set_dirichlet_flags(ml_prob, faceName, x, dirichlet);

     ctrl:: NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS < ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >::ctrl_or_state_set_dirichlet_fixed_values(ml_prob, faceName, x, value);
   }
  
  
  //************************control****************************************************

  else if(!strcmp(name,"control")) {

    ctrl:: NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS < ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >::ctrl_or_state_set_dirichlet_flags(ml_prob, faceName, x, dirichlet);

    ctrl:: NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS < ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >::ctrl_or_state_set_dirichlet_fixed_values(ml_prob, faceName, x, value);
  
  }

  //************************ adjoint ****************************************************
  
  else if(!strcmp(name,"adjoint")) {
        dirichlet = true;
        value = 0.;
  }

  //************************mu****************************************************

  else if(!strcmp(name,"mu")) {
      
    dirichlet = false;

  }

  
  return dirichlet;
}
 // Unknowns - END  ==================



 // Not Unknowns - BEGIN  ==================
double Solution_set_initial_conditions_Not_Unknowns(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {
    
    double value = 0.;
    
    if(!strcmp(name, "TargReg")) {
        value = ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(x);
    }
    else if(!strcmp(name, "act_flag")) {
        value = 0.;
    }


    return value;
}
 // Not Unknowns - END  ==================








int main(int argc, char** args) {

  // ======= Init ========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);


  // ======= Problem - BEGIN  ==================
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

  const std::string relative_path_to_mesh_input_folder = "../../../";
  const std::string infile = Files::get_input_file_with_prefix(ctrl::mesh_input, relative_path_to_mesh_input_folder);

  
  ml_mesh.set_mesh_filename( ctrl::mesh_input );
  
  
  const double Lref = 1.;
  
  const bool read_groups = true;
  const bool read_boundary_groups = true;
  

  
    ml_mesh.ReadCoarseMeshFileReadingBeforePartitioning(infile, Lref, read_groups, read_boundary_groups);
    
// // //  BEGIN FillISvectorDofMapAllFEFamilies
           ml_mesh.GetLevelZero(0)->dofmap_all_fe_families_initialize();

           std::vector < unsigned > elem_partition_from_mesh_file_to_new = ml_mesh.GetLevelZero(0)->elem_offsets();
  
           std::vector < unsigned > node_mapping_from_mesh_file_to_new = ml_mesh.GetLevelZero(0)->node_offsets();
           
           ml_mesh.GetLevelZero(0)->dofmap_all_fe_families_clear_ghost_dof_list_for_other_procs();

// // //   END FillISvectorDofMapAllFEFamilies
       
           ml_mesh.GetLevelZero(0)->BuildElementAndNodeStructures();
 
           

  ml_mesh.BuildFETypesBasedOnExistingCoarseMeshGeomElements();
//   ml_mesh.BuildFETypesBasedOnExistingCoarseMeshGeomElements(fe_quad_rule_vec[0].c_str()); 
  //doesn't need dofmap. This seems to be abstract, it can be performed right after the mesh geometric elements are read. It is needed for local MG operators, as well as for Integration of shape functions...
  //The problem is that it also performs global operations such as matrix sparsity pattern, global MG operators... And these also use _dofOffset...
  //The problem is that this class actually has certain functions which have REAL structures instead of only being ABSTRACT FE families!!!
  // So:
//   - Mesh and Multimesh are real and not abstract, and rightly so 
//   - Elem is real and rightly so, and only Geometric. However it contains some abstract Geom Element, but there seems to be no overlap with FE families
  ml_mesh.PrepareNewLevelsForRefinement();       //doesn't need dofmap

  // ======= Mesh, Coarse reading - END ==================


  // ======= Mesh: Refinement - BEGIN ==================
  const unsigned numberOfUniformLevels = N_UNIFORM_LEVELS;
  const unsigned erased_levels = N_ERASED_LEVELS;
  unsigned numberOfSelectiveLevels = 0;
  
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  //RefineMesh contains a similar procedure as ReadCoarseMesh. In particular, the dofmap at each level is filled there
  // ======= Mesh: Refinement - END ==================

  // ======= Mesh: Group info - BEGIN ==================
   ml_mesh.set_group_info(relative_path_to_mesh_input_folder);
  // ======= Mesh: Group info - END ==================

  // ======= Solutions that are not Unknowns, auxiliary; needed for Boundary of Boundary of Control region - BEFORE COARSE ERASING - BEGIN  ==================
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
  
  // ======= Solutions that are not Unknowns, auxiliary - END  ==================

  
  // ======= Mesh: Coarse erasing - BEGIN  ========================
  ml_mesh.EraseCoarseLevels(erased_levels);
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
  std::vector< Unknown > unknowns = elliptic :: pure_boundary< IS_CTRL_FRACTIONAL_SOBOLEV > :: provide_list_of_unknowns( ml_mesh.GetDimension() );

  ml_sol.AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions_Unknowns);
  
  for (unsigned int u = 0; u < unknowns.size(); u++)  { 
      ml_sol.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order, unknowns[u]._time_order, unknowns[u]._is_pde_unknown);
      ml_sol.Initialize(unknowns[u]._name.c_str(), Solution_set_initial_conditions_Unknowns, & ml_prob);
      ml_sol.GenerateBdc(unknowns[u]._name.c_str(), (unknowns[u]._time_order == 0) ? "Steady" : "Time_dependent", & ml_prob);
  }
  // ======= Solutions that are Unknowns - END ==================
  

  // ======= Solutions that are not Unknowns - BEGIN  ==================
  // ******** targ reg - BEGIN 
  ml_sol.AddSolution("TargReg", DISCONTINUOUS_POLYNOMIAL, ZERO, steady_flag, is_an_unknown_of_a_pde);
  ml_sol.Initialize("TargReg",     Solution_set_initial_conditions_Not_Unknowns, & ml_prob);
  // ******** targ reg - END
  
  // ******** cont reg - BEGIN 
  const std::string volume_control_region = "ContReg";
  ml_sol.AddSolution(volume_control_region.c_str(), DISCONTINUOUS_POLYNOMIAL, ZERO, steady_flag, is_an_unknown_of_a_pde);
  ml_sol.InitializeBasedOnControlFaces< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >(volume_control_region.c_str(),  Solution_set_initial_conditions_Not_Unknowns, & ml_prob);
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
  
  
  // ******** bdry of bdry - BEGIN 
 ctrl::Gamma_control_equation_fractional_sobolev_differentiability_index<
                femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES, 
                femus::ctrl:: square_or_cube:: pure_boundary < femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >
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
  
  // ******** bdry of bdry - END 

  // ======= Solutions that are not Unknowns - END  ==================
  
  
  //== Solution: CHECK SOLUTION FE TYPES between Unknowns and Not Unknowns - BEGIN  ==
  if ( ml_sol.GetSolutionType("control") != ml_sol.GetSolutionType(act_set_flag_name[0].c_str())) abort();
  //== Solution: CHECK SOLUTION FE TYPES between Unknowns and Not Unknowns - END ==
  
  // ======= Problem, Solution - END ==================

  

  // ======= Problem, Quad Rule - BEGIN  ========================
  std::vector< std::string > fe_quad_rule_vec;
  fe_quad_rule_vec.push_back("seventh" /*"ninth"*/);
  fe_quad_rule_vec.push_back("eighth" /*"tenth"*/);

  ml_prob.SetQuadratureRuleAllGeomElemsMultiple(fe_quad_rule_vec);
  ml_prob.set_all_abstract_fe_AD_or_not();
  // ======= Problem, Quad Rule - END  ========================

  // ======= Problem, System - BEGIN ========================
  NonLinearImplicitSystemWithPrimalDualActiveSetMethod & system_opt = ml_prob.add_system < NonLinearImplicitSystemWithPrimalDualActiveSetMethod > ("BoundaryControl");
  
       // ======= Not an Unknown, but needed in the System with PDAS ========================
  system_opt.SetActiveSetFlagName(act_set_flag_name);    //MU
       // ======= Not an Unknown, but needed in the System with PDAS ========================
  
       // ======= System Unknowns ========================
  for (unsigned int u = 0; u < unknowns.size(); u++) { system_opt.AddSolutionToSystemPDE(unknowns[u]._name.c_str());  }
       // ======= System Unknowns ========================

  system_opt.SetAssembleFunction( elliptic::pure_boundary< IS_CTRL_FRACTIONAL_SOBOLEV >::assemble_elliptic_dirichlet_control );

// *****************
  const unsigned n_components_state = n_components_ctrl;
  std::vector<std::string> state_vars(n_components_state);  state_vars[0] = "state";
  std::vector<std::string> ctrl_vars(n_components_ctrl);     ctrl_vars[0] = "control";
  
  system_opt.set_state_vars(state_vars);
  system_opt.set_ctrl_vars(ctrl_vars);
  
  system_opt.SetDebugNonlinear(true);
  
//   system_opt.SetDebugFunction(ctrl::cost_functional::compute_cost_functional_regularization_bdry);

  //   ///@todo weird error if I comment this line, I expect nothing to happen but something in the assembly gets screwed up in memory I guess
// *****************
  
 

  system_opt.init();  /// I need to put this init before, later I will remove it   /// @todo it seems like you cannot do this init INSIDE A FUNCTION... understand WHY!
 
  set_dense_pattern_for_unknowns(system_opt, unknowns);

  system_opt.init();

  //----
  std::ostringstream sp_out_base2; sp_out_base2 << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "sp_";
  sp_out_base2 << "after_second_init_";
  
  unsigned n_levels = ml_mesh.GetNumberOfLevels();
  system_opt._LinSolver[n_levels - 1]->sparsity_pattern_print_nonzeros(sp_out_base2.str(), "on");
  system_opt._LinSolver[n_levels - 1]->sparsity_pattern_print_nonzeros(sp_out_base2.str(), "off");
  //----

  system_opt.MGsolve();
//   double totalAssemblyTime = 0.;
//   system_opt.nonlinear_solve_single_level(MULTIPLICATIVE, totalAssemblyTime, 0, 0);
//   system_opt.assemble_call_before_boundary_conditions(1);
  // ======= Problem, System  - END ========================

  
  // ======= Post-processing - BEGIN ========================


  // ======= Post-processing, Computations - BEGIN ========================
  femus::ctrl::cost_functional< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES,  femus::ctrl:: square_or_cube:: Domain_elements_containing_Gamma_control< femus::ctrl:: GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >   ,femus::ctrl::COST_FUNCTIONAL_WITHOUT_REG>::compute_cost_functional_regularization_bdry(ml_prob, 0, 0, state_vars, ctrl_vars, ALPHA_CTRL_BDRY, BETA_CTRL_BDRY, QRULE_I);
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




 

