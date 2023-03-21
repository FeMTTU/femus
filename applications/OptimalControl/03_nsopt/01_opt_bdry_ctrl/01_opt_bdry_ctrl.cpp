#include "adept.h"

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NumericVector.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "Fluid.hpp"
#include "Parameter.hpp"
#include "Files.hpp"
#include "PetscMatrix.hpp"
#include "paral.hpp"//to get iproc HAVE_MPI is inside here

#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"

#include "ElemType.hpp"

#include "opt_common.hpp"
#include "00_cost_functional.hpp"
#include "03_opt_system_inequalities.hpp"

using namespace femus;

#include   "../manufactured_solutions.hpp"

#include  "../../opt_systems_boundary_control_eqn_sobolev_integer.hpp"

#include   "../nsopt_params.hpp"


#include  "../../square_or_cube_03_cost_functional_without_regularization.hpp"



//****** Mesh ********************************
#define no_of_ref N_UNIFORM_LEVELS     //mesh refinements
//**************************************



 
 
//****** Convergence ********************************
#define exact_sol_flag      0  // 1 = if we want to use manufactured solution; 0 = if we use regular convention
#define compute_conv_flag   0  // 1 = if we want to compute the convergence and error ; 0 =  no error computation

#define NO_OF_L2_NORMS 9   //U,V,P,UADJ,VADJ,PADJ,ctrl_0,ctrl_1,THETA
#define NO_OF_H1_NORMS 6    //U,V,UADJ,VADJ,ctrl_0, ctrl_1
//**************************************


 


 //Unknown definition  - BEGIN ==================
   enum pos_vector_quantities {pos_index_state = 0, pos_index_adj, pos_index_ctrl, pos_index_mu};
 
 const std::vector< unsigned >  provide_state_adj_ctrl_mu_offsets(const unsigned int dimension) {
     
     std::vector<unsigned>  opt_control_offsets(4);  
     
     opt_control_offsets[pos_index_state] = 0; 
     opt_control_offsets[pos_index_adj]   = dimension + 1; 
     
     opt_control_offsets[pos_index_mu]    = 2 * (dimension + 1); 
     
#if INEQ_FLAG != 0
  opt_control_offsets[pos_index_ctrl]  = opt_control_offsets[pos_index_mu] + dimension;
#else
  opt_control_offsets[pos_index_ctrl]  =  2 * (dimension + 1);
#endif
     
     
     return opt_control_offsets;
     
 } 
 
 
 const std::vector< Unknown >  provide_list_of_unknowns(const unsigned int dimension) {


 std::vector< Unknown >  unknowns( 3*(dimension+1) + dimension);
 
 
  std::vector< unsigned >  vector_offsets = provide_state_adj_ctrl_mu_offsets(dimension);


const int state_pos_begin   =  vector_offsets[pos_index_state];
  const int adj_pos_begin   =  vector_offsets[pos_index_adj];
 const int ctrl_pos_begin   =  vector_offsets[pos_index_ctrl];
   const int mu_pos_begin   =  vector_offsets[pos_index_mu];
   

  //--- names - BEGIN
                                        unknowns[state_pos_begin + 0]._name      = "u_0";
                                        unknowns[state_pos_begin + 1]._name      = "u_1";
  if (dimension == 3)                   unknowns[state_pos_begin + 2]._name      = "u_2";
                                unknowns[state_pos_begin + dimension]._name      = "u_p";
  
                        unknowns[adj_pos_begin + 0]._name      = "adj_0";
                        unknowns[adj_pos_begin + 1]._name      = "adj_1";
  if (dimension == 3)   unknowns[adj_pos_begin + 2]._name      = "adj_2";
                unknowns[adj_pos_begin + dimension]._name      = "adj_p";
  
                       unknowns[ctrl_pos_begin + 0]._name      = "ctrl_0";
                       unknowns[ctrl_pos_begin + 1]._name      = "ctrl_1";
  if (dimension == 3)  unknowns[ctrl_pos_begin + 2]._name      = "ctrl_2";
  
               unknowns[ctrl_pos_begin + dimension]._name      = "theta";

#if INEQ_FLAG != 0
                       unknowns[mu_pos_begin + 0]._name      = "mu_0";
                       unknowns[mu_pos_begin + 1]._name      = "mu_1";
  if (dimension == 3)  unknowns[mu_pos_begin + 2]._name      = "mu_2";
#endif  
 
   //--- names - END
 

  //--- fe family - BEGIN
                                        unknowns[state_pos_begin + 0]._fe_family      = LAGRANGE;
                                        unknowns[state_pos_begin + 1]._fe_family      = LAGRANGE;
  if (dimension == 3)                   unknowns[state_pos_begin + 2]._fe_family      = LAGRANGE;
                                unknowns[state_pos_begin + dimension]._fe_family      = LAGRANGE/*DISCONTINOUS_POLYNOMIAL*/;
  
                        unknowns[adj_pos_begin + 0]._fe_family      =  LAGRANGE;
                        unknowns[adj_pos_begin + 1]._fe_family      =  LAGRANGE;
  if (dimension == 3)   unknowns[adj_pos_begin + 2]._fe_family      =  LAGRANGE;
                unknowns[adj_pos_begin + dimension]._fe_family      =  LAGRANGE/*DISCONTINOUS_POLYNOMIAL*/;
  
                       unknowns[ctrl_pos_begin + 0]._fe_family      =   LAGRANGE;
                       unknowns[ctrl_pos_begin + 1]._fe_family      =   LAGRANGE;
  if (dimension == 3)  unknowns[ctrl_pos_begin + 2]._fe_family      =   LAGRANGE;
                                                                        
               unknowns[ctrl_pos_begin + dimension]._fe_family      =   DISCONTINUOUS_POLYNOMIAL;     //theta

#if INEQ_FLAG != 0
                       unknowns[mu_pos_begin + 0]._fe_family       =   LAGRANGE;
                       unknowns[mu_pos_begin + 1]._fe_family       =   LAGRANGE;
  if (dimension == 3)  unknowns[mu_pos_begin + 2]._fe_family       =   LAGRANGE;
#endif  
 
  //--- fe family - END


   //--- fe order - BEGIN
                                        unknowns[state_pos_begin + 0]._fe_order      = SECOND;
                                        unknowns[state_pos_begin + 1]._fe_order      = SECOND;
  if (dimension == 3)                   unknowns[state_pos_begin + 2]._fe_order      = SECOND;
                                unknowns[state_pos_begin + dimension]._fe_order      = FIRST;
  
                        unknowns[adj_pos_begin + 0]._fe_order          = SECOND;
                        unknowns[adj_pos_begin + 1]._fe_order          = SECOND;
  if (dimension == 3)   unknowns[adj_pos_begin + 2]._fe_order          = SECOND;
                unknowns[adj_pos_begin + dimension]._fe_order          = FIRST;
  
                       unknowns[ctrl_pos_begin + 0]._fe_order        = SECOND;
                       unknowns[ctrl_pos_begin + 1]._fe_order        = SECOND;
  if (dimension == 3)  unknowns[ctrl_pos_begin + 2]._fe_order        = SECOND;
                                                                        
               unknowns[ctrl_pos_begin + dimension]._fe_order      =   ZERO;     //theta

#if INEQ_FLAG != 0
                       unknowns[mu_pos_begin + 0]._fe_order          = SECOND;
                       unknowns[mu_pos_begin + 1]._fe_order          = SECOND;
  if (dimension == 3)  unknowns[mu_pos_begin + 2]._fe_order          = SECOND;
#endif  
 
  
  
  
  
   //--- fe order - END
     

  
  
     for (unsigned int u = 0; u < unknowns.size(); u++) {
         
              unknowns[u]._time_order = 0;
              unknowns[u]._is_pde_unknown = true;
              unknowns[u]._is_sparse = true;
              
     }
     
 
 for (unsigned int u = 0; u < dimension; u++) {
   unknowns[ctrl_pos_begin + u]._is_sparse = IS_CTRL_FRACTIONAL_SOBOLEV ? false: true;
     }
 
 
 unknowns[ctrl_pos_begin + dimension]._is_sparse = false;    //theta
 
 
   return unknowns;
     
}

 //Unknown definition  - END ==================

 
double Solution_set_initial_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;

     if(!strcmp(name,"TargReg")) {
        value = femus::ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(x);
    }
    else if(!strcmp(name,"ContReg")) {
        value = femus::ctrl::  square_or_cube:: Domain_elements_containing_Gamma_control< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES > :: ControlDomainFlag_bdry(x);
    }

    return value;
    
}



bool Solution_set_boundary_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char SolName[], double& value, const int faceName, const double time) {
  //1: bottom  //2: right  //3: top  //4: left
  
  bool dirichlet = true;
   value = 0.;
 

                if (!strcmp(SolName, "ctrl_0"))       {
//                     if (facename == FACE_FOR_CONTROL) dirichlet = false; 
  if (faceName == FACE_FOR_CONTROL) {
     if (x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] > GAMMA_CONTROL_LOWER - 1.e-5 &&
         x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] < GAMMA_CONTROL_UPPER + 1.e-5)  { 
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
     if (x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] > GAMMA_CONTROL_LOWER - 1.e-5 &&
         x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] < GAMMA_CONTROL_UPPER + 1.e-5)  { 
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
     if (x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] > GAMMA_CONTROL_LOWER - 1.e-5 &&
         x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] < GAMMA_CONTROL_UPPER + 1.e-5)  { 
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
     if (x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] > GAMMA_CONTROL_LOWER - 1.e-5 &&
         x[ femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES ::tangential_direction_to_Gamma_control(faceName) ] < GAMMA_CONTROL_UPPER + 1.e-5)  { 
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


void assemble_ns_dirichlet_control_pure_boundary(MultiLevelProblem& ml_prob);    


double*  GetErrorNorm(const MultiLevelProblem& ml_prob, MultiLevelSolution* ml_sol, Solution* sol_coarser_prolongated);
// ||u_h - u_(h/2)||/||u_(h/2)-u_(h/4)|| = 2^alpha, alpha is order of conv 
//i.e. ||prol_(u_(i-1)) - u_(i)|| = err(i) => err(i-1)/err(i) = 2^alpha ,implemented as log(err(i)/err(i+1))/log2



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
  ml_prob.set_all_abstract_fe_multiple();
  // ======= Problem, Quad Rule - END  ========================


  
  // ======= Mesh, Coarse reading - BEGIN ==================
  MultiLevelMesh ml_mesh;

  const std::string infile = files.get_input_file_with_prefix(ctrl::mesh_input, "../../../");


  const bool read_groups = true;
  const bool read_boundary_groups = true;
    
    ml_mesh.ReadCoarseMeshFileReadingBeforePartitioning(infile.c_str(), Lref, read_groups, read_boundary_groups);
    
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
    maxNumberOfMeshes = no_of_ref;
  } else {
    maxNumberOfMeshes = 2;
  }
  
  
  
   // ======= Solutions that are Unknowns - BEGIN  ==================
  std::vector< Unknown > unknowns = provide_list_of_unknowns( dim );
   // ======= Solutions that are Unknowns - END  ==================

  
  
#if compute_conv_flag == 1
  MultiLevelMesh ml_mesh_all_levels;
  
     ml_mesh_all_levels.ReadCoarseMesh(infile.c_str(),fe_quad_rule_vec[0].c_str(),Lref);
 
     double comp_conv[maxNumberOfMeshes][NO_OF_L2_NORMS+NO_OF_H1_NORMS];
 
  
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
      ml_sol_all_levels->Initialize(unknowns[u]._name.c_str(), Solution_set_initial_conditions, & ml_prob);
  }
  
      ml_sol_all_levels->AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions);
   for (unsigned int u = 0; u < unknowns.size(); u++)  { 
     ml_sol_all_levels->GenerateBdc(unknowns[u]._name.c_str(), (unknowns[u]._time_order == 0) ? "Steady" : "Time_dependent", & ml_prob);
  }
   // ======= Solutions that are Unknowns - END  ==================


  // ======= Solutions that are not Unknowns - BEGIN  ==================
            ml_sol_all_levels->AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO);
            ml_sol_all_levels->AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO);
   ml_sol_all_levels->Initialize("TargReg",     Solution_set_initial_conditions, & ml_prob);
   ml_sol_all_levels->Initialize("ContReg",     Solution_set_initial_conditions, & ml_prob);
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

  // ======= Solution, auxiliary; needed for Boundary of Boundary of Control region - BEFORE COARSE ERASING - BEGIN  ==================
  const std::string node_based_bdry_bdry_flag_name = NODE_BASED_BDRY_BDRY;
  const unsigned  steady_flag = 0;
  const bool      is_an_unknown_of_a_pde = false;
  
  const FEFamily node_bdry_bdry_flag_fe_fam = LAGRANGE;
  const FEOrder node_bdry_bdry_flag_fe_ord = SECOND;
  
  MultiLevelSolution * ml_sol_bdry_bdry_flag = ctrl::Gamma_control_equation_fractional_sobolev_differentiability_index<
                femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES, 
                femus::ctrl:: square_or_cube:: Domain_elements_containing_Gamma_control< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >
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
      ml_sol.Initialize(unknowns[u]._name.c_str(), Solution_set_initial_conditions, & ml_prob);
  }
  
  ml_sol.AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions);
   for (unsigned int u = 0; u < unknowns.size(); u++)  { 
     ml_sol.GenerateBdc(unknowns[u]._name.c_str(), (unknowns[u]._time_order == 0) ? "Steady" : "Time_dependent", & ml_prob);
  }
  
  // ======= Solutions that are Unknowns - END  ==================

  // ======= Solutions that are not Unknowns - BEGIN  ==================
  ml_sol.AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO);
   ml_sol.Initialize("TargReg",     Solution_set_initial_conditions, & ml_prob);

   ml_sol.AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO);
   ml_sol.Initialize("ContReg",     Solution_set_initial_conditions, & ml_prob);
   
   
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
    ml_sol.Initialize(act_set_flag_name[d].c_str(), Solution_set_initial_conditions, & ml_prob);
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
   
  
   system_opt.SetAssembleFunction(assemble_ns_dirichlet_control_pure_boundary);

    
// *****************
 
  std::vector< unsigned >  vector_offsets = provide_state_adj_ctrl_mu_offsets(dim);


const int state_pos_begin   =  vector_offsets[pos_index_state];
 const int ctrl_pos_begin   =  vector_offsets[pos_index_ctrl];
   
   
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
    
      for(int j = 0; j < NO_OF_L2_NORMS+NO_OF_H1_NORMS; j++)       comp_conv[i-1][j] = norm[j];
  
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
       FE_convergence<double>::output_convergence_rate(comp_conv[i][NO_OF_L2_NORMS + j], comp_conv[i + 1][NO_OF_L2_NORMS + j], norm_names_H1[j], maxNumberOfMeshes , i );
    }
  }
  std::cout << std::endl;
  std::cout << "=======================================================================" << std::endl;
#endif
  // ======= Convergence Rate, Postprocess - END  ==================

  
  return 0;
}
 




///@todo we have to decide if order the phi by fem families, by Unknowns or by Quantities
/// I say the best way is by Quantity, because they have all the need to have phi functions as Unknowns

        

void assemble_ns_dirichlet_control_pure_boundary(MultiLevelProblem& ml_prob) {

  // Main objects - BEGIN *******************************************
 //System
 NonLinearImplicitSystemWithPrimalDualActiveSetMethod * mlPdeSys  = & ml_prob.get_system< NonLinearImplicitSystemWithPrimalDualActiveSetMethod >("NSOpt");
  
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  //System, Equation
  const char* system_name            = mlPdeSys->name().c_str();
  LinearEquationSolver*  pdeSys	 = mlPdeSys->_LinSolver[level];   
  bool assembleMatrix = mlPdeSys->GetAssembleMatrix(); 
   
  SparseMatrix*	 JAC	= pdeSys->_KK;
  NumericVector* RES 	= pdeSys->_RES;
    

  constexpr bool print_algebra_global = true;
  constexpr bool print_algebra_local = false;

  //Mesh
  Mesh*		 msh    = ml_prob._ml_msh->GetLevel(level);
  elem*		 el	= msh->el;
  const unsigned dim 	= msh->GetDimension();
  unsigned dim2     = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned nel		= msh->GetNumberOfElements();
  unsigned igrid	= msh->GetLevel();
  unsigned iproc 	= msh->processor_id();
  unsigned   nprocs = msh->n_processors();
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));


  //Solution
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;
  Solution*	 sol  	         = ml_prob._ml_sol->GetSolutionLevel(level);
  
  // Main objects - END *******************************************


  // ======= Geometry at Dofs - BEGIN  =======
   unsigned coordXType = 2; /*BIQUADR_FE*/// get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
   unsigned solType_coords = coordXType;
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
  const unsigned int dim_offset_grad = /*dim*/  3  /*2*/    ;
 
  vector< vector < double> > coordX(dim);
  for(int i = 0; i < dim; i++) {   coordX[i].reserve(max_size);   }
  
  vector< vector < double> > coordX_bd(/*space_*/dim);
  for(int i = 0; i < dim; i++) { coordX_bd[i].reserve(max_size);  }
 
  const unsigned dim_bdry = dim - 1;
  // ======= Geometry at Dofs - END  =======


 // ======= Solutions, Unknowns - BEGIN =======

  std::vector< unsigned >  vector_offsets = provide_state_adj_ctrl_mu_offsets(dim);


const int state_pos_begin   =  vector_offsets[pos_index_state];
  const int adj_pos_begin   =  vector_offsets[pos_index_adj];
 const int ctrl_pos_begin   =  vector_offsets[pos_index_ctrl];
   const int mu_pos_begin   =  vector_offsets[pos_index_mu];
   


  
    const int press_type_pos   = dim;

  const int theta_index      = ctrl_pos_begin + dim;
  
  
    std::vector<std::string> ctrl_name;
    ctrl_name.resize(dim);
    ctrl_name[0] = "ctrl_0";
    ctrl_name[1] = "ctrl_1";
    if (dim == 3)  ctrl_name[2] = "ctrl_2";
 
  const int pos_mat_ctrl = ctrl_pos_begin;
  const int pos_sol_ctrl = ctrl_pos_begin;

    const unsigned int n_components_ctrl = dim;

  std::vector< Unknown > unknowns = provide_list_of_unknowns( dim );
  
  const int n_unknowns = unknowns.size();
    

  vector < std::string > Solname_Mat(n_unknowns);  
  vector < unsigned > SolPdeIndex(n_unknowns);
  vector < unsigned > SolIndex_Mat(n_unknowns);  
  vector < unsigned > SolFEType_Mat(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    Solname_Mat[ivar] = unknowns[ivar]._name;
    SolPdeIndex[ivar]	= mlPdeSys->GetSolPdeIndex(Solname_Mat[ivar].c_str());
    SolIndex_Mat[ivar]	= ml_sol->GetIndex        (Solname_Mat[ivar].c_str());
    SolFEType_Mat[ivar]	= ml_sol->GetSolutionType(SolIndex_Mat[ivar]);
  }

  vector < unsigned > Sol_n_el_dofs_Mat_vol(n_unknowns);
  

  const double cost_functional_coeff = COST_FUNCTIONAL_COEFF;
  const double alpha = ALPHA_CTRL_BDRY;
  const double beta  = BETA_CTRL_BDRY;
 // ======= Solutions, Unknowns - END =======
  
     
  //=== Sol (quantities, not only unknowns) - BEGIN ========================================================
  
  const unsigned int n_quantities = ml_sol->GetSolutionSize();
  
 //***************************************************
    vector < std::string > Solname_quantities(n_quantities);
    
        for(unsigned ivar=0; ivar < Solname_quantities.size(); ivar++) {
            Solname_quantities[ivar] = ml_sol->GetSolutionName(ivar);
        }
 //***************************************************
 
    vector < unsigned > SolIndex_quantities(n_quantities);      //should have Sol order
    vector < unsigned > SolFEType_quantities(n_quantities);     //should have Sol order
    vector < unsigned > Sol_n_el_dofs_quantities(n_quantities); //should have Sol order
    
    for(unsigned ivar=0; ivar < n_quantities; ivar++) {
        SolIndex_quantities[ivar]    = ml_sol->GetIndex        (Solname_quantities[ivar].c_str());
        SolFEType_quantities[ivar]   = ml_sol->GetSolutionType(SolIndex_quantities[ivar]);
    }
    
    
#if INEQ_FLAG != 0
  //MU
  //************** variables for ineq constraints: act flag ****************************   
  std::vector<unsigned int> solIndex_act_flag_sol(n_components_ctrl); 
  
  ctrl::mixed_state_or_ctrl_inequality< femus:: ctrl:: square_or_cube::mixed_state_or_ctrl_inequality >::store_act_flag_in_old(mlPdeSys, ml_sol, sol, solIndex_act_flag_sol);
  //************** variables for ineq constraints: act flag ****************************   
    

  const int ineq_flag = INEQ_FLAG;
  const double c_compl = C_COMPL;
  vector < vector < double/*int*/ > > sol_actflag(n_components_ctrl);   
  vector < vector < double > > ctrl_lower(n_components_ctrl);           
  vector < vector < double > > ctrl_upper(n_components_ctrl);           
	  for (unsigned c = 0; c < n_components_ctrl; c++) {
      sol_actflag[c].reserve(max_size);
      ctrl_lower[c].reserve(max_size);
      ctrl_upper[c].reserve(max_size);
      }

    
//MU
  std::vector<unsigned int> ctrl_index_in_mat(n_components_ctrl); 
  std::vector<unsigned int>   mu_index_in_mat(n_components_ctrl);    
      for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
           ctrl_index_in_mat[kdim] =  SolPdeIndex[ctrl_pos_begin + kdim];
             mu_index_in_mat[kdim] =  SolPdeIndex[mu_pos_begin + kdim];
      }
#endif
      
    
    
  //=== Sol (quantities, not only unknowns) - END ========================================================
  

  
  // ======= Solutions, Unknowns at dofs - BEGIN =======
  vector < vector < double > > Sol_eldofs_Mat(n_unknowns);
  vector < vector < double > > gradSol_eldofs_Mat(n_unknowns);
  
  for(int k = 0; k < n_unknowns; k++) {
    Sol_eldofs_Mat[k].reserve(max_size);
    gradSol_eldofs_Mat[k].reserve(max_size*dim);    
  }
  
  
// ****** Solutions, Unknowns at dofs, theta value from proc0 - BEGIN 
  double solTheta = femus::get_theta_value(msh->n_processors(), sol, SolIndex_Mat[theta_index]);
// ****** Solutions, Unknowns at dofs, theta value from proc0  - END

  // ======= Solutions, Unknowns at dofs - END =======


  // ======= Solutions, Unknowns at quadrature points - BEGIN =======
    vector < double > SolVAR_qp(n_unknowns);   //sol_V,P_gss_of_st,adj,ctrl_ie@quadraturepoints
    vector < vector < double > > gradSolVAR_qp(n_unknowns);
    for(int k = 0; k < n_unknowns; k++) {  gradSolVAR_qp[k].resize(dim_offset_grad /*space_dim*/);  }


  vector < vector < double > >  sol_adj_x_vol_at_bdry_gss(dim);
  for (int ldim =0; ldim < dim; ldim++) sol_adj_x_vol_at_bdry_gss[ldim].reserve(max_size);
  
  vector < double > grad_adj_dot_n_res_qp;
  vector < double > grad_adj_dot_n_jac_qp;
  grad_adj_dot_n_res_qp.reserve(max_size);
  grad_adj_dot_n_jac_qp.reserve(max_size);
  // ======= Solutions, Unknowns at quadrature points - END =======
      
 
  // ====== Geometry at Quadrature points - BEGIN ==============================================================================
     std::vector < std::vector < double > >  JacI_iqp(space_dim);
     std::vector < std::vector < double > >  Jac_iqp(dim);
    for (unsigned d = 0; d < Jac_iqp.size(); d++) {   Jac_iqp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp.size(); d++) { JacI_iqp[d].resize(dim); }
    
    double detJac_iqp;
    double AbsDetJxWeight_iqp = 0.;

     std::vector < std::vector < double > >  JacI_iqp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_iqp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_iqp_bdry.size(); d++) {   Jac_iqp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp_bdry.size(); d++) { JacI_iqp_bdry[d].resize(dim-1); }
    
    double detJac_iqp_bdry;
    double AbsDetJxWeight_iqp_bdry = 0.;
    
    std::vector<double> normal_iqp(dim_offset_grad /*space_dim*/, 0.);
  // ====== Geometry at Quadrature points - END ==============================================================================
  
     
// ======= FE at Quadrature, all - BEGIN =======
      //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all;
  ml_prob.get_all_abstract_fe_multiple(elem_all);
  
  //==========================================================================================
  vector < vector < double > > phi_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_x_gss_fe(NFE_FAMS);
  
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_gss_fe[fe].reserve(max_size);
      phi_x_gss_fe[fe].reserve(max_size * dim_offset_grad);
   }
   
   
  //boundary adjoint & ctrl shape functions  
  vector < vector < double > > phi_bd_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_x_bd_gss_fe(NFE_FAMS);

    for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_bd_gss_fe[fe].reserve(max_size);
      phi_x_bd_gss_fe[fe].reserve(max_size * dim_offset_grad);
  //bdry vol adj  evaluated at bdry points
    }
    
    
  //bdry vol adj  evaluated at bdry points
   vector < vector < double > > phi_vol_at_bdry_fe(NFE_FAMS);
   vector < vector < double > > phi_x_vol_at_bdry_fe(NFE_FAMS);

    for(int fe=0; fe < NFE_FAMS; fe++) {  
         phi_vol_at_bdry_fe[fe].reserve(max_size);
       phi_x_vol_at_bdry_fe[fe].reserve(max_size * dim_offset_grad);
    }
  //==========================================================================================

 //********************* bdry cont *******************
 //*************************************************** 
  vector <double> phi_ctrl_bdry;  
  vector <double> phi_ctrl_x_bdry; 

  phi_ctrl_bdry.reserve(max_size);
  phi_ctrl_x_bdry.reserve(max_size * space_dim);
 //*************************************************** 
// ======= FE at Quadrature, all - END =======
    
  
  // ======= Equation, local - BEGIN =======
  vector < vector < int > > L2G_dofmap_Mat(n_unknowns); 
  vector < vector < double > > Res(n_unknowns);
  vector < vector < vector < double > > > Jac(n_unknowns);
  
  for(int i = 0; i < n_unknowns; i++) {     
    L2G_dofmap_Mat[i].reserve(max_size);
       Res[i].reserve(max_size);
  }
   
  if (assembleMatrix) {
    
    for(int i = 0; i < n_unknowns; i++) {
      Jac[i].resize(n_unknowns);    
      for(int j = 0; j < n_unknowns; j++) {
	Jac[i][j].reserve(max_size * max_size);	
      }
    }

  }
  
  vector < vector < double > > Jac_outer(dim);
  vector < double > Res_outer(1);

         for(int i = 0; i < dim; i++) {  Jac_outer[i].reserve(max_size); }
  // ======= Equation, local - END =======


  // ======= Equation, global - BEGIN =======
    RES->zero();
    if(assembleMatrix) JAC->zero();
  // ======= Equation, global - END =======

 
  // ======= Parameters - BEGIN ======= 
  const double IRe = ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();

  //****** Control ********************************
 double penalty_outside_control_domain_boundary = PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY;       // penalty for zero control outside Gamma_c
 double penalty_dirichlet_bc_u_equal_q = PENALTY_DIRICHLET_BC_U_EQUAL_Q_BOUNDARY;         //penalty for u = q

 double theta_value_outside_fake_element = 0.;
 //**************************************

 // ======= Parameters - END =======
  
    
   
  // ======= Fractional - BEGIN =======
  const double s_frac = S_FRAC;

  const double check_limits = 1.;//1./(1. - s_frac); // - s_frac;
   
  //--- quadrature rules -------------------
  constexpr unsigned qrule_i = QRULE_I;
  constexpr unsigned qrule_j = QRULE_J;
  constexpr unsigned qrule_k = QRULE_K;
  //----------------------
  // ======= Fractional - END =======

  
  
   ///@todo to avoid Petsc complaint about out-of-bounds allocation *******************************
//  MatSetOption(static_cast< PetscMatrix* >(JAC)->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
   //************************************
     
  
  
    
      // ***************** ADD PART - BEGIN  *******************

    if ( IS_CTRL_FRACTIONAL_SOBOLEV ) {
  
     ctrl::Gamma_control_equation_fractional_sobolev_differentiability_index<
                femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES, 
                femus::ctrl:: square_or_cube:: Domain_elements_containing_Gamma_control< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >
                >
::control_eqn_bdry(iproc,
                   nprocs,
                    ml_prob,
                    ml_sol,
                    sol,
                    msh,
                    pdeSys,
                    //-----------
                    geom_element_iel,
                    solType_coords,
                    dim,
                    space_dim,
                    dim_bdry,
                    //-----------
                    n_unknowns,
                    Solname_Mat,
                    SolFEType_Mat,
                    SolIndex_Mat,
                    SolPdeIndex,
                    Sol_n_el_dofs_Mat_vol, 
                    Sol_eldofs_Mat,  
                    L2G_dofmap_Mat,
                    max_size,
                    //-----------
                    n_quantities,
                    SolFEType_quantities,
//                     Sol_n_el_dofs_quantities, //filled inside
                    //-----------
                    elem_all,
                     //-----------
                    Jac_iqp_bdry,
                    JacI_iqp_bdry,
                    detJac_iqp_bdry,
                    AbsDetJxWeight_iqp_bdry,
                    phi_ctrl_bdry,
                    phi_ctrl_x_bdry, 
                    //-----------
                    n_components_ctrl,
                    pos_mat_ctrl,
                    pos_sol_ctrl,
                    IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY,
                    //-----------
                    JAC,
                    RES,
                    assembleMatrix,
                    //-----------
                    alpha,
                    beta,     
                    //-----------
                    s_frac,
                    check_limits,
                    USE_Cns,
                    OP_Hhalf,
                    OP_L2,
                    RHS_ONE,
                    UNBOUNDED,
                    //-----------
                    qrule_i,
                    qrule_j,
                    qrule_k,
                    Nsplit,
                    Quadrature_split_index,
                    N_DIV_UNBOUNDED,
                    //-----------
                    print_algebra_local
                    );
                    

   }
  
  else {
  
   femus::ctrl::Gamma_control_equation_integer_sobolev_differentiability_index<
                femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES, 
                femus::ctrl:: square_or_cube:: Domain_elements_containing_Gamma_control< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >
                >
                ::control_eqn_bdry(iproc,
                    ml_prob,
                    ml_sol,
                    sol,
                    msh,
                    pdeSys,
                    //-----------
                    geom_element_iel,
                    solType_coords,
                    space_dim,
                    //-----------
                    n_unknowns,
                    Solname_Mat,
                    SolFEType_Mat,
                    SolIndex_Mat,
                    SolPdeIndex,
                    Sol_n_el_dofs_Mat_vol, 
                    Sol_eldofs_Mat,  
                    L2G_dofmap_Mat,
                    max_size,
                    //-----------
                     n_quantities,
                    SolFEType_quantities,
                    //-----------
                    elem_all,
                     //-----------
                    Jac_iqp_bdry,
                    JacI_iqp_bdry,
                    detJac_iqp_bdry,
                    AbsDetJxWeight_iqp_bdry,
                     //-----------
                    phi_ctrl_bdry,
                    phi_ctrl_x_bdry, ///@todo this should be done for all components of ctrl
                    //-----------
                    n_components_ctrl,
                    pos_mat_ctrl,
                    pos_sol_ctrl,
                    IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY,
                    //-----------
                    JAC,
                    RES,
                    assembleMatrix,
                    //-----------
                    alpha,
                    beta,
                    RHS_ONE,
                    qrule_i,
                    //-----------
                    print_algebra_local
                    ) ;
  
  }
 
 
 
 
 
 
 
 
 
 
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

  // geometry *****************************
      geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
      const short unsigned ielGeom = geom_element_iel.geom_type();
      
      geom_element_iel.set_elem_center_3d(iel, solType_coords);
  // geometry end *****************************
  
  // equation *****************************
    unsigned nDofsV = msh->GetElementDofNumber(iel, SolFEType_Mat[state_pos_begin]);
    unsigned nDofsP = msh->GetElementDofNumber(iel, SolFEType_Mat[state_pos_begin + press_type_pos]);
    
    unsigned nDofsVadj = msh->GetElementDofNumber(iel,SolFEType_Mat[adj_pos_begin]);
    unsigned nDofsPadj = msh->GetElementDofNumber(iel,SolFEType_Mat[adj_pos_begin + press_type_pos]);

    unsigned nDofsGctrl     = msh->GetElementDofNumber(iel, SolFEType_Mat[ctrl_pos_begin]);
    unsigned nDofsThetactrl = msh->GetElementDofNumber(iel, SolFEType_Mat[theta_index] );

    unsigned nDofsVP = dim * nDofsV + nDofsP; //same for state and adjoint
    unsigned nDofsVPctrl = dim * nDofsGctrl + nDofsThetactrl; //control
   
    unsigned nDofsVP_tot = 2 * nDofsVP + (nDofsVPctrl);
  // equation end *****************************
  
  //***** set target domain flag ********************************** 
   int target_flag = 0;
       target_flag = ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(geom_element_iel.get_elem_center_3d());
   //***************************************       
   
 //************ set control flag - BEGIN *********************
    int does_iel_contain_a_bdry_control_face = 0;
        does_iel_contain_a_bdry_control_face = ctrl ::  square_or_cube:: Domain_elements_containing_Gamma_control< ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES > :: ControlDomainFlag_bdry(geom_element_iel.get_elem_center_3d());
 //************ initialize control node flag: for each Volume Elem, tell me if we have a Boundary Control dof *********************
    std::vector< std::vector<int> > control_node_flag_iel_jface(n_components_ctrl);
	    for(unsigned idim=0; idim < control_node_flag_iel_jface.size(); idim++) {
	          control_node_flag_iel_jface[idim].resize(nDofsGctrl);
    std::fill(control_node_flag_iel_jface[idim].begin(), control_node_flag_iel_jface[idim].end(), 0);
	    }
 //****************************- END *********************** 
  
   //###########################  fake_iel_flag - BEGIN ########################################  
    unsigned int fake_iel_flag = 0;
    unsigned int global_row_index_bdry_constr = pdeSys->KKoffset[SolPdeIndex[theta_index]][iproc];
  for (unsigned  k = 0; k < n_unknowns; k++) {
	unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType_Mat[k]);	//nDofs_V,P_of_st,adj,ctrl
	Sol_n_el_dofs_Mat_vol[k] = ndofs_unk;
	Sol_eldofs_Mat[k].resize(ndofs_unk);	//sol_V,P_of_st,adj,ctrl
	L2G_dofmap_Mat[k].resize(ndofs_unk); 
    for (unsigned i = 0; i < ndofs_unk; i++) {
	      unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType_Mat[k]);    // global to global mapping between solution node and solution dof // via local to global solution node
	  Sol_eldofs_Mat[k][i] = (*sol->_Sol[SolIndex_Mat[k]])(solDof);      // global extraction and local storage for the solution
	  L2G_dofmap_Mat[k][i] = pdeSys->GetSystemDof(SolIndex_Mat[k], SolPdeIndex[k], i, iel);
 
    if (k == SolPdeIndex[theta_index] && L2G_dofmap_Mat[k][i] == global_row_index_bdry_constr) {       fake_iel_flag = iel;  }
    }
  }
  //############################ fake_iel_flag - END #######################################
    
  //************ fake theta flag: - BEGIN this flag tells me what degrees of freedom of the current element are fake for the theta variable  *********************
    std::vector<int>  bdry_int_constr_pos_vec(1, global_row_index_bdry_constr); /*KKoffset[SolPdeIndex[PADJ]][iproc]*/
    std::vector<int> fake_theta_flag(nDofsThetactrl, 0);
    for (unsigned i = 0; i < nDofsThetactrl; i++) {
      if ( L2G_dofmap_Mat[ SolPdeIndex[theta_index] ] [i] == bdry_int_constr_pos_vec[0]) { 	fake_theta_flag[i] = 1;       }
    }
 //************ fake theta flag - END *********************

// setting Jac and Res to zero  - BEGIN ******************************* 
    for(int ivar=0; ivar<n_unknowns; ivar++) {
              Res[SolPdeIndex[ivar]].resize(Sol_n_el_dofs_Mat_vol[ivar]);
     std::fill(Res[SolPdeIndex[ivar]].begin(), Res[SolPdeIndex[ivar]].end(), 0.);         
//       memset(&Res[SolPdeIndex[ivar]][0],0.,Sol_n_el_dofs_Mat_vol[ivar]*sizeof(double));
    }
   
    for(int ivar = 0; ivar < n_unknowns; ivar++) {
      for(int jvar = 0; jvar < n_unknowns; jvar++) {
	    if(assembleMatrix) {  //MISMATCH
                  Jac[ SolPdeIndex[ivar] ] [SolPdeIndex[jvar] ].resize( Sol_n_el_dofs_Mat_vol[ivar] * Sol_n_el_dofs_Mat_vol[jvar] );
     std::fill(Jac[ SolPdeIndex[ivar] ] [SolPdeIndex[jvar] ].begin(), Jac[ SolPdeIndex[ivar] ] [SolPdeIndex[jvar] ].end(), 0.);         
// 		  memset(&Jac[ SolPdeIndex[ivar] ] [SolPdeIndex[jvar] ][0], 0., Sol_n_el_dofs_Mat_vol[ivar]*Sol_n_el_dofs_Mat_vol[jvar]*sizeof(double));
           }
        }
     }
     
    for(int ivar = 0; ivar < dim; ivar++)     std::fill(Jac_outer[ivar].begin(), Jac_outer[ivar].end(), 0.); //did not use Jac_outer as Jac itself was placing the values as expected
    Res_outer[0] = 0.;
// setting Jac and Res to zero  - END ******************************* 

  
  
//======== BoundaryLoop - BEGIN  =====================================================================

  // Perform face loop over elements that contain some control face
  if (does_iel_contain_a_bdry_control_face == 1) {
	  
      double tau = 0.;
      vector<double> normal_iqp(dim_offset_grad, 0);
	       
	  // loop on faces of the current element

      for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
          
//-------          
       const unsigned ielGeom_bd = msh->GetElementFaceType(iel, jface);    
       
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();
//-------          

	    // look for boundary faces
            const int bdry_index = el->GetFaceElementIndex(iel, jface);
            
	    if( bdry_index < 0) {
	      unsigned int face_in_rectangle_domain = -( msh->el->GetFaceElementIndex(iel, jface) + 1);

// 	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"u_0",tau,face,0.) && tau!=0.){
          if(  face_in_rectangle_domain == FACE_FOR_CONTROL) { //control face
              
//=================================================== 
		   //we use the dirichlet flag to say: if dirichlet == true, we set 1 on the diagonal. if dirichlet == false, we put the boundary equation
		  std::vector<bool> is_bc_for_control_dirichlet_on_jface(n_components_ctrl);
		  for(unsigned idim = 0; idim < is_bc_for_control_dirichlet_on_jface.size(); idim++) {
		      is_bc_for_control_dirichlet_on_jface[idim] = ml_sol->GetBdcFunctionMLProb()(& ml_prob, geom_element_iel.get_elem_center_bdry_3d(), ctrl_name[idim].c_str(), tau, face_in_rectangle_domain, 0.);
		  }
	  
	
//========= initialize gauss quantities on the boundary ============================================
		vector < double >                      SolVAR_bd_qp(n_unknowns);
		vector < vector < double > >       gradSolVAR_bd_qp(n_unknowns);
		for(int k=0; k<n_unknowns; k++) {  gradSolVAR_bd_qp[k].resize(dim_offset_grad /*space_dim*/);  }

//========= gauss_loop boundary===============================================================
		  for(unsigned iqp_bdry = 0; iqp_bdry < ml_prob.GetQuadratureRule(ielGeom_bd).GetGaussPointsNumber(); iqp_bdry++) {
    
    elem_all[qrule_i][ielGeom_bd][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), iqp_bdry, Jac_iqp_bdry, JacI_iqp_bdry, detJac_iqp_bdry, space_dim);
	elem_all[qrule_i][ielGeom_bd][solType_coords]->compute_normal(Jac_iqp_bdry, normal_iqp);
    
    AbsDetJxWeight_iqp_bdry = detJac_iqp_bdry * ml_prob.GetQuadratureRule(ielGeom_bd).GetGaussWeightsPointer()[iqp_bdry];

    elem_all[qrule_i][ielGeom_bd][SolFEType_Mat[theta_index]] ->shape_funcs_current_elem(iqp_bdry, JacI_iqp_bdry,phi_bd_gss_fe[SolFEType_Mat[theta_index]],phi_x_bd_gss_fe[SolFEType_Mat[theta_index]], boost::none, space_dim);
    
    for (unsigned ldim = 0; ldim < dim; ldim++) {
    elem_all[qrule_i][ielGeom_bd][ SolFEType_Mat[ldim + ctrl_pos_begin]] ->shape_funcs_current_elem(iqp_bdry, JacI_iqp_bdry,
                                                                                                    phi_bd_gss_fe[   SolFEType_Mat[ldim + ctrl_pos_begin] ],
                                                                                                    phi_x_bd_gss_fe[ SolFEType_Mat[ldim + ctrl_pos_begin] ], boost::none , space_dim);
    }
                

    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv_vol_at_bdry_new(geom_element_iel.get_coords_at_dofs_3d(), iqp_bdry, jface, Jac_iqp/*not_needed_here*/, JacI_iqp, detJac_iqp/*not_needed_here*/, space_dim);
    
    for (unsigned ldim = 0; ldim < dim; ldim++) {
    elem_all[qrule_i][ielGeom][SolFEType_Mat[ldim + adj_pos_begin]]->shape_funcs_vol_at_bdry_current_elem(iqp_bdry, jface, JacI_iqp, 
                                                                                                          phi_vol_at_bdry_fe  [ SolFEType_Mat[ldim + adj_pos_begin] ],                                                                                                  phi_x_vol_at_bdry_fe[ SolFEType_Mat[ldim + adj_pos_begin] ], boost::none, space_dim);
    }
     

//========== compute gauss quantities on the boundary - BEGIN ===============================================

//=============== Boundary control - BEGIN ========================================= 
		    for (unsigned  kdim = 0; kdim < n_components_ctrl; kdim++) {
                
			    unsigned int ctrl_index = kdim + ctrl_pos_begin;
                
					  SolVAR_bd_qp[ SolPdeIndex[ctrl_index] ] = 0.;
			  for(unsigned ivar2 = 0; ivar2 < dim_offset_grad; ivar2++) {  gradSolVAR_bd_qp[ SolPdeIndex[ctrl_index] ][ivar2] = 0.; }
	  
	         const unsigned ndof_bdry = msh->GetElementFaceDofNumber(iel, jface, SolFEType_Mat[ctrl_index]);

			  for(int i_bd = 0; i_bd < ndof_bdry; i_bd++) {
		                  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bd);
                          
                       SolVAR_bd_qp[SolPdeIndex[ctrl_index]]            += phi_bd_gss_fe  [ SolFEType_Mat[ctrl_index] ][i_bd]                  * Sol_eldofs_Mat[ SolPdeIndex[ ctrl_index ] ][i_vol];

			      for(unsigned ivar2 = 0; ivar2 < dim_offset_grad; ivar2++) { 
                      gradSolVAR_bd_qp[SolPdeIndex[ctrl_index]][ivar2]  += phi_x_bd_gss_fe[ SolFEType_Mat[ctrl_index] ][i_bd * dim_offset_grad + ivar2] * Sol_eldofs_Mat[ SolPdeIndex[ ctrl_index ] ][i_vol];
                        }
                        
			  }
			  
		    }
//=============== Boundary control - END ========================================= 
 
//=============== grad dot n - BEGIN ========================================= 
//     compute gauss quantities on the boundary through VOLUME interpolation
		for(unsigned ldim = 0; ldim < dim; ldim++) {   sol_adj_x_vol_at_bdry_gss[ldim].resize(dim_offset_grad); }
		grad_adj_dot_n_res_qp.resize(dim);
		grad_adj_dot_n_jac_qp.resize(dim);
		for(unsigned ldim=0; ldim<dim; ldim++) {   std::fill(sol_adj_x_vol_at_bdry_gss[ldim].begin(), sol_adj_x_vol_at_bdry_gss[ldim].end(), 0.);  }
		
		for(unsigned ldim=0; ldim<dim; ldim++) {  
		  for (int iv = 0; iv < nDofsVadj; iv++)  {
                     for (int d = 0; d < dim_offset_grad /*space_dim*/; d++) {
			   sol_adj_x_vol_at_bdry_gss[ldim][d] += Sol_eldofs_Mat[SolPdeIndex[ldim + adj_pos_begin]][iv] * phi_x_vol_at_bdry_fe[SolFEType_Mat[ldim + adj_pos_begin]][iv * dim_offset_grad /*space_dim*/ + d];//notice that the convention of the orders x y z is different from vol to bdry
                    }
		  }  
		      
		  grad_adj_dot_n_res_qp[ldim] = 0.;
		  for(unsigned d=0; d < dim_offset_grad /*space_dim*/; d++) {
		      grad_adj_dot_n_res_qp[ldim] += sol_adj_x_vol_at_bdry_gss[ldim][d] * normal_iqp[d];  
		  }
		}
		
//=============== grad dot n - END =========================================       
		  
//========== compute gauss quantities on the boundary - END ================================================

		
//=============== construct control node flag - BEGIN =========================================    
//this is all based on jface only!!!

	      /* (control_node_flag_iel_jface)       picks nodes on \Gamma_c
	         (1 - control_node_flag_iel_jface)   picks nodes on \Omega \setminus \Gamma_c
	       */
		    for(unsigned idim = 0; idim < n_components_ctrl; idim++) {
                
                unsigned int ctrl_index = idim + ctrl_pos_begin;

			if (is_bc_for_control_dirichlet_on_jface[idim] == false) {
                
	         const unsigned ndof_bdry = msh->GetElementFaceDofNumber(iel, jface, SolFEType_Mat[ctrl_index]);
             
		for(unsigned i_bdry = 0; i_bdry < ndof_bdry; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
                control_node_flag_iel_jface[idim][i_vol] = 1;
 			    }
		    }
		
        }
//=============== construct control node flag - END =========================================    
        

        const unsigned ndofs_bdry_max = msh->GetElementFaceDofNumber(iel, jface, solType_coords);
        const unsigned nDofsGctrl_bdry = msh->GetElementFaceDofNumber(iel, jface,  SolFEType_Mat[ctrl_pos_begin]);
        

            for(unsigned i_bdry = 0; i_bdry < ndofs_bdry_max; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);

//============ Boundary-Boundary and Boundary-Volume Residuals - BEGIN ============================================================================================
		  
		      for (unsigned  kdim = 0; kdim < dim; kdim++) {
			
			    double lap_res_dctrl_ctrl_bd_kdim = 0.;
                 if (i_bdry < nDofsGctrl_bdry) {
			      for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
				  lap_res_dctrl_ctrl_bd_kdim += gradSolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim] *
				                           phi_x_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin]][i_bdry * dim_offset_grad + jdim /*i_bdry + jdim*ndofs_bdry_max*/];
			      }//jdim
                 }
			
/*delta_state row */	   if (i_vol < nDofsV)      Res[kdim]                 [i_vol]  += - control_node_flag_iel_jface[kdim][i_vol] * penalty_dirichlet_bc_u_equal_q * (Sol_eldofs_Mat[SolPdeIndex[kdim + state_pos_begin]][i_vol] - Sol_eldofs_Mat[SolPdeIndex[kdim + ctrl_pos_begin]][i_vol]);	    //u-g
/*delta_adjoint row */     if (i_vol < nDofsVadj)   Res[kdim + adj_pos_begin] [i_vol]  += 0.;	   
/*delta_control row */     if (i_vol < nDofsGctrl)  Res[kdim + ctrl_pos_begin][i_vol]  += - control_node_flag_iel_jface[kdim][i_vol] * AbsDetJxWeight_iqp_bdry * (
                                                                                          IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY * alpha * SolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * phi_bd_gss_fe[SolFEType_Mat[kdim +  ctrl_pos_begin]][i_bdry]
                                                                                        + IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY * beta * lap_res_dctrl_ctrl_bd_kdim
                                                                                        - IRe * grad_adj_dot_n_res_qp[kdim]  * phi_bd_gss_fe[SolFEType_Mat[kdim +  ctrl_pos_begin]][i_bdry]
                                                                                        - /*(*sol->_Sol[SolIndex_Mat[theta_index]])(0)*/solTheta * phi_bd_gss_fe[SolFEType_Mat[kdim +  ctrl_pos_begin]][i_bdry] * normal_iqp[kdim]      //*sol->_Sol[SolIndex_Mat[theta_index]])(0) finds the global value from KKDof pos(72, 169,etc), SolVAReldof_theta gets the value in the boundary point which will be zero. Theta is just a const
                                                                                        );

		      }//kdim  

//============ Boundary-Boundary and Boundary-Volume Residuals - END ==================================================================================================

//============ Boundary-Boundary Jacobian i-BDRY/j-BDRY - BEGIN  ==================================================================================================
		      for(unsigned j_bdry = 0; j_bdry < ndofs_bdry_max; j_bdry ++) {
			  unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, jface, j_bdry);

				std::vector < double > lap_jac_dctrl_ctrl_bd(dim, 0.);
                    if (i_bdry < nDofsGctrl_bdry && j_bdry < nDofsGctrl_bdry) {

			     for (unsigned kdim = 0; kdim < dim; kdim++) {
				for (unsigned ldim = 0; ldim < dim_offset_grad /*space_dim*/; ldim++) {
///@todo the error is in here
                    lap_jac_dctrl_ctrl_bd[kdim] += phi_x_bd_gss_fe[ SolFEType_Mat[kdim + ctrl_pos_begin] ][i_bdry * dim_offset_grad + ldim] * 
                                                   phi_x_bd_gss_fe[ SolFEType_Mat[kdim + ctrl_pos_begin] ][j_bdry * dim_offset_grad + ldim];
                                        }
                                    }
                    }
                    
                 
			  for (unsigned  kdim = 0; kdim < dim; kdim++) {
                  
//DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
			    if(i_vol < nDofsV && j_vol < nDofsV && i_vol == j_vol)       		          Jac[kdim][kdim][i_vol * nDofsV + j_vol]	  +=  (1.) * penalty_dirichlet_bc_u_equal_q * control_node_flag_iel_jface[kdim][i_vol];  //u
			 
//BLOCK delta_state - control------------------------------------------------------------------------------------
			    if(i_vol < nDofsV && j_vol < nDofsGctrl && i_vol == j_vol) 	Jac[kdim][kdim + ctrl_pos_begin][i_vol * nDofsGctrl + j_vol]  += (-1.) * penalty_dirichlet_bc_u_equal_q * control_node_flag_iel_jface[kdim][i_vol];  //-g

//DIAG BLOCK delta_control - control  --------------------------------------------------------------------------------------
			  if(i_vol < nDofsGctrl && j_vol < nDofsGctrl) {
				      Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i_vol * nDofsGctrl + j_vol] +=   control_node_flag_iel_jface[kdim][i_vol] * AbsDetJxWeight_iqp_bdry * (
                                                                                                       IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY * alpha * phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin] ][i_bdry] * phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin] ][j_bdry]
                                                                                                    + IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY * beta * lap_jac_dctrl_ctrl_bd[kdim]
                                                                                                 );
                  } 
                  
               }//kdim
			
                   
            }//end j_bdry loop
//============ Boundary-Boundary Jacobian i-BDRY/j-BDRY - END  ==================================================================================================
		    
//============ Boundary-Volume Jacobian i-BDRY/j-VOL - BEGIN  ==================================================================================================
//BLOCK delta_control - adjoint------------------------------------------------------------------------------------------------
//===================loop over j in the VOLUME (while i is in the boundary)	      
		for(unsigned j = 0; j < nDofsVadj; j ++) {
            
			//=============== grad dot n  =========================================    
		    for(unsigned ldim = 0; ldim < dim; ldim++) {
			    grad_adj_dot_n_jac_qp[ldim] = 0.;
			    for(unsigned d = 0; d < dim_offset_grad /*space_dim*/; d++) {
				grad_adj_dot_n_jac_qp[ldim] += phi_x_vol_at_bdry_fe[SolFEType_Mat[ldim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + d] * normal_iqp[d];  //notice that the convention of the orders x y z is different from vol to bdry
			    }
		    }
			  //=============== grad dot n  =========================================    

			  for (unsigned kdim = 0; kdim < dim; kdim++) {
				Jac[kdim + ctrl_pos_begin][kdim + adj_pos_begin][i_vol*nDofsVadj + j] += control_node_flag_iel_jface[kdim][i_vol] * (-1.) * (AbsDetJxWeight_iqp_bdry  * phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin]][i_bdry] * IRe * grad_adj_dot_n_jac_qp[kdim]);    		      
			  }
		} // end j loop for volume 
//============ Boundary-Volume Jacobian i-BDRY/j-VOL - END  ==================================================================================================

		    }//end i_bdry loop


		    
//============ Boundary Integral Constraint Residual - BEGIN ============================================================================================
	  for (unsigned  kdim = 0; kdim < dim; kdim++) {
// 		for(unsigned i=0; i < nDofsThetactrl; i ++) { avoid because it is an element dof
/*delta_theta row */ 	/* Res[theta_index][i]*/ Res_outer[0] +=  /*fake_theta_flag[i] **/ AbsDetJxWeight_iqp_bdry * SolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * normal_iqp[kdim] ;
// 		}  
	  }
		  
//============ Boundary Integral Constraint Residual - END ============================================================================================
		
//============ Boundary Integral Constraint Jacobian - BEGIN ============================================================================================
       		for(unsigned i_bdry = 0; i_bdry < ndofs_bdry_max; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
            
 		    for (unsigned  kdim = 0; kdim < dim; kdim++) { 
			  for(unsigned i =0; i < nDofsThetactrl; i ++) {
			    if(i_vol < nDofsGctrl) {
				double temp = AbsDetJxWeight_iqp_bdry * ( phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin]][i_bdry] * normal_iqp[kdim]);
//ROW_BLOCK delta_theta - control -- loop over i in the VOLUME (while j(/i_vol) is in the boundary) -------------------------------------------------------------------------------------------------------------
			      Jac[theta_index][ctrl_pos_begin + kdim][i*nDofsGctrl + i_vol]     += - temp; /*AbsDetJxWeight_iqp_bdry * ( phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin]][i_bdry] * normal_iqp[kdim])*/
//COLUMN_BLOCK delta_control - theta ---- loop over j in the VOLUME (while i(/i_vol) is in the boundary) ---------------------------------------------------------------------------------------------------
			      Jac[ctrl_pos_begin + kdim][theta_index][i_vol*nDofsThetactrl + i] += - control_node_flag_iel_jface[kdim][i_vol] /** phi_bd_gss_fe[SolFEType_Mat[theta_index]][i]*/*temp; /*AbsDetJxWeight_iqp_bdry * ( phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin]][i_bdry] * normal_iqp[kdim]);*/
			    }//endif
			  }// i 
		    }//kdim           
            
            
            }		    
//============ Boundary Integral Constraint Jacobian - END ============================================================================================
		    
		    
                }  //end iqp_bdry loop
	  
             }    //end if control face
	     }  //end if boundary faces
      }  // loop over element faces //jface   
  } //end if control element flag

//======== BoundaryLoop - END ==== Boundary Residuals  and Jacobians ==================	
    
    ///@todo at this point, the node flag has been filled for ALL faces
    
 
 
 
//======================= VolumeLoop with Integration (and Zero boundary control outside Gamma_c) - BEGIN =====================================================    

for(unsigned iqp = 0; iqp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); iqp++) {
	
	// *** get Jacobian and test function and test function derivatives ***
       // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), iqp, Jac_iqp, JacI_iqp, detJac_iqp, space_dim);
    AbsDetJxWeight_iqp = detJac_iqp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[iqp];
   
     for(int fe=0; fe < NFE_FAMS; fe++) {
    elem_all[qrule_i][ielGeom][fe]->shape_funcs_current_elem(iqp, JacI_iqp,phi_gss_fe[fe],phi_x_gss_fe[fe], boost::none , space_dim);
      }
 
 
//geometry eval at gauss points - BEGIN ********************************
 vector < double > coordX_gss(dim, 0.);
 	for(unsigned k = 0; k <  dim; k++) {
	  for(unsigned i = 0; i < geom_element_iel.get_coords_at_dofs_3d()[k].size(); i++) { 
         coordX_gss[k] += geom_element_iel.get_coords_at_dofs_3d()[k][i] * phi_gss_fe[ solType_coords ][i];
      }
    }
//geometry eval at gauss points - END  ********************************

//begin unknowns eval at gauss points ********************************
	for(unsigned unk = 0; unk < n_unknowns; unk++) {
	    SolVAR_qp[unk] = 0.;
	    for(unsigned ivar2 = 0; ivar2 < dim_offset_grad /*space_dim*/; ivar2++) {
           gradSolVAR_qp[unk][ivar2] = 0.; 
	    }
	  
	    for(unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[unk]; i++) {
                        SolVAR_qp[unk] += phi_gss_fe[ SolFEType_Mat[unk] ][i]             * Sol_eldofs_Mat[SolPdeIndex[unk]][i];
		for(unsigned ivar2 = 0; ivar2 < dim_offset_grad /*space_dim*/; ivar2++) {       
		    gradSolVAR_qp[unk][ivar2]  += phi_x_gss_fe[ SolFEType_Mat[unk] ][i * dim_offset_grad /*space_dim*/ + ivar2] * Sol_eldofs_Mat[SolPdeIndex[unk]][i]; 
		}
	    }//ndofsunk
	  
	} //unk 
 //end unknowns eval at gauss points ********************************
	
#if exact_sol_flag == 1
//======= computation of RHS (force and desired velocity) using MMS - BEGIN =============================================== 
//state values-------------------- //non-hom bdry
vector <double>  exact_stateVel(dim);
  mms_lid_driven::value_stateVel(coordX_gss, exact_stateVel);
vector <double>  exact_lap_stateVel(dim);
  mms_lid_driven::laplace_stateVel(coordX_gss, exact_lap_stateVel);
vector < vector < double > > exact_grad_stateVel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_stateVel[k].resize(dim);
    std::fill(exact_grad_stateVel[k].begin(), exact_grad_stateVel[k].end(), 0.);
}
  mms_lid_driven::gradient_stateVel(coordX_gss,exact_grad_stateVel);

//adjoint values--------------------//hom bdry
vector <double>  exact_adjVel(dim);
  mms_lid_driven::value_adjVel(coordX_gss, exact_adjVel);
vector <double>  exact_lap_adjVel(dim);
  mms_lid_driven::laplace_adjVel(coordX_gss, exact_lap_adjVel);
vector < vector < double > > exact_grad_adjVel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_adjVel[k].resize(dim);
    std::fill(exact_grad_adjVel[k].begin(), exact_grad_adjVel[k].end(), 0.);
}
  mms_lid_driven::gradient_adjVel(coordX_gss,exact_grad_adjVel);
vector <double> exact_grad_adjPress(dim);
  mms_lid_driven::gradient_adjPress(coordX_gss, exact_grad_adjPress);

//convection terms from delta_state-------------------------------------
vector <double>  exact_conv_u_nabla_u(dim,0.);

for (unsigned k = 0; k < dim; k++){
    for (unsigned i = 0; i < dim; i++){
    exact_conv_u_nabla_u[k] += exact_grad_stateVel[k][i] * exact_stateVel[i] ; 
    }
}

//convection terms from delta_adjoint-------------------------
vector <double>  exact_conv_u_nabla_uadj(dim,0.);
vector <double>  exact_conv_nabla_uT_uadj(dim,0.);

for (unsigned k = 0; k < dim; k++){
    for (unsigned i = 0; i < dim; i++){
    exact_conv_u_nabla_uadj[k] += exact_grad_adjVel[k][i] * exact_stateVel[i] ; 
    exact_conv_nabla_uT_uadj[k] += exact_grad_stateVel[i][k] * exact_adjVel[i];
    }
}

//force and desired velocity ---------------------------------------------
vector <double> exactForce(dim,0.);
vector <double> exactVel_d(dim,0.);
for (unsigned k = 0; k < dim; k++){
    exactForce[k] = - IRe * exact_lap_stateVel[k]  + advection_flag * exact_conv_u_nabla_u[k] + exact_grad_adjPress[k];
    exactVel_d[k] =   exact_stateVel[k] + (1./cost_functional_coeff) * (IRe * exact_lap_adjVel[k] - exact_grad_adjPress[k]) 
                    + (1./cost_functional_coeff) * advection_flag * (exact_conv_u_nabla_uadj[k] - exact_conv_nabla_uT_uadj[k]);
}
//======= computation of RHS (force and desired velocity) using MMS - END =============================================== 
#endif



//============ delta_state row - BEGIN ============================================================================================

  for (unsigned i = 0; i < nDofsV; i++) {
// FIRST ROW
	for (unsigned  kdim = 0; kdim < dim; kdim++) { // velocity block row 
                double lap_res_du_u		= 0.; 
                double adv_res_uold_nablauold 	= 0.;
	      for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		    lap_res_du_u 	       += gradSolVAR_qp[SolPdeIndex[kdim]][jdim] * phi_x_gss_fe[SolFEType_Mat[kdim]][i * dim_offset_grad /*space_dim*/ + jdim];
          }
	      for (unsigned jdim = 0; jdim < dim; jdim++) {
		   adv_res_uold_nablauold  += SolVAR_qp[SolPdeIndex[jdim]]  * gradSolVAR_qp[SolPdeIndex[kdim]][jdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i];
	      }      
	      Res[kdim][i]   +=  (           
#if exact_sol_flag == 0
                                         + force[kdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i]
 #endif                                      
 #if exact_sol_flag == 1
                                       + exactForce[kdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i]
 #endif
                                       - IRe*lap_res_du_u 
                                       - advection_flag * adv_res_uold_nablauold 
                                       + SolVAR_qp[SolPdeIndex[press_type_pos]] * phi_x_gss_fe[SolFEType_Mat[kdim]][i * dim_offset_grad /*space_dim*/ + kdim]) * AbsDetJxWeight_iqp; 
	}	    
//DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
	for (unsigned j = 0; j < nDofsV; j++) {
		      vector < double > lap_jac_du_u(dim,0.);
		      vector < double > adv_uold_nablaunew(dim,0.);
	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { 
		    lap_jac_du_u[kdim] += phi_x_gss_fe[SolFEType_Mat[kdim]][i * dim_offset_grad /*space_dim*/ + jdim]*phi_x_gss_fe[SolFEType_Mat[kdim]][j * dim_offset_grad/*space_dim*/ + jdim];
            }
          }
        for (unsigned  kdim = 0; kdim < dim; kdim++) { 
         for (unsigned  jdim = 0; jdim < dim; jdim++) { //diagonal blocks only
		    adv_uold_nablaunew[kdim] 	 += SolVAR_qp[SolPdeIndex[jdim]] * phi_x_gss_fe[ SolFEType_Mat[kdim] ][j * dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i];
                }  //jdim
	      }
	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
		Jac[kdim][kdim][i*nDofsV + j] += (   IRe * lap_jac_du_u[kdim] 
                                            + advection_flag * adv_uold_nablaunew[kdim] 		// c(u_old, u_new, delta_lambda)
                                            + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i]	 // c(u_new,u_old,delta_lambda) diagonal blocks  ..... unew_nablauold
                                         ) * AbsDetJxWeight_iqp; 
              unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
		Jac[kdim][off_kdim][i*nDofsV + j] += (	advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[off_kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][off_kdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i]	// c(u_new,u_old,delta_lambda) off-diagonal blocks  ..... unew_nablauold
                                             ) * AbsDetJxWeight_iqp;
	      }
	} //j_du_u loop
     
//BLOCK Pressure
      for (unsigned j = 0; j < nDofsP; j++) {
	    for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim][press_type_pos][i*nDofsP + j] += -( phi_gss_fe[SolFEType_Mat[press_type_pos]][j] * phi_x_gss_fe[SolFEType_Mat[kdim]][i * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
	    }
      }//j_press loop
   }//i_state loop

//DIV_state
  for (unsigned i = 0; i < nDofsP; i++) {
		    double div_u_du_qp =0.;
      for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      div_u_du_qp += gradSolVAR_qp[SolPdeIndex[kdim]][kdim] ;
      }
      Res[press_type_pos][i]  +=  ( (div_u_du_qp) * phi_gss_fe[SolFEType_Mat[press_type_pos]][i] ) * AbsDetJxWeight_iqp;
      for (unsigned j = 0; j < nDofsV; j++) {
	  for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      Jac[press_type_pos][kdim][i*nDofsV + j] += -( phi_gss_fe[SolFEType_Mat[press_type_pos]][i] * phi_x_gss_fe[SolFEType_Mat[kdim]][j * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
	  }
      } //j loop
   }//i_div_state
//============ delta_state row - END ============================================================================================


    
//============ delta_adjoint row - BEGIN =============================================================================================
  
  for (unsigned i = 0; i < nDofsVadj; i++) {
// SECOND ROW
     for (unsigned kdim = 0; kdim < dim; kdim++) { 
		    double lap_res_dadj_adj 			= 0.;
		    double adv_res_phiadj_nablauold_uadjold 	= 0.;
		    double adv_res_uold_nablaphiadj_uadjold 	= 0.;
	   for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		lap_res_dadj_adj 		             += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]][jdim]*phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim];
       }
	   for (unsigned jdim = 0; jdim < dim; jdim++) {
		adv_res_phiadj_nablauold_uadjold     += phi_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i] * gradSolVAR_qp[SolPdeIndex[jdim]][kdim] 			* SolVAR_qp[SolPdeIndex[jdim + adj_pos_begin]];
		adv_res_uold_nablaphiadj_uadjold     += SolVAR_qp[SolPdeIndex[jdim]]		       * phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]];
	   }
	  Res[kdim + adj_pos_begin][i] += ( 
#if exact_sol_flag == 0
                            - cost_functional_coeff * target_flag * ctrl::square_or_cube :: cost_functional_without_regularization::DesiredTargetVec()[kdim] 			      * phi_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i]
 #endif                                      
 #if exact_sol_flag == 1
                            - cost_functional_coeff * target_flag * exactVel_d[kdim] 			      * phi_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i]
 #endif
                                        + cost_functional_coeff * target_flag * SolVAR_qp[SolPdeIndex[kdim]] * phi_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i]
                                        - IRe*lap_res_dadj_adj
                                        - advection_flag * adv_res_phiadj_nablauold_uadjold
                                        - advection_flag * adv_res_uold_nablaphiadj_uadjold
                                        + SolVAR_qp[SolPdeIndex[press_type_pos + adj_pos_begin]] * phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim]
                                        ) * AbsDetJxWeight_iqp;
      }
      
//BLOCK delta_adjoint - state------------------------------------------------------------------------------------------
     for (unsigned j = 0; j < nDofsV; j++) {
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + adj_pos_begin][kdim][i*nDofsV + j] += ( - cost_functional_coeff * target_flag * phi_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i] * phi_gss_fe[SolFEType_Mat[kdim]][j] 
                                                             + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i]    * phi_x_gss_fe[ SolFEType_Mat[kdim] ][j*dim_offset_grad /*space_dim*/ + kdim] 		* SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]   //c(delta_u, u_new, lambda_old)  diagonal blocks  ......phiadj_nablaunew_uadjold 
                                                             + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[kdim] ][j] 			* phi_x_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	 //c(u_new, delta_u, lambda_old) diagonal blocks  ......unew_nablaphiadj_uadjold
                                                            ) * AbsDetJxWeight_iqp;
              unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
		Jac[kdim + adj_pos_begin][off_kdim][i*nDofsV + j] += (  advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i] * phi_x_gss_fe[ SolFEType_Mat[off_kdim] ][j*dim_offset_grad /*space_dim*/ + kdim]		      * SolVAR_qp[SolPdeIndex[off_kdim + adj_pos_begin]]   //c(delta_u, u_new, lambda_old)  off-diagonal blocks  ......phiadj_nablaunew_uadjold 
                                                              + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[off_kdim] ][j] 		  * phi_x_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + off_kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	   //c(u_new, delta_u, lambda_old) off-diagonal blocks  ......unew_nablaphiadj_uadjold
                                                             ) * AbsDetJxWeight_iqp;
	  }            
     }//j_dadj_u loop


//DIAG BLOCK delta_adjoint - adjoint---------------------------------------------------------------------------------
     for (unsigned j = 0; j < nDofsVadj; j++) {
		    vector < double > lap_jac_dadj_adj(dim,0.);
		    vector < double > adv_uold_nablaphiadj_uadjnew(dim, 0.);
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { 
		  lap_jac_dadj_adj[kdim] += phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim] * phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + jdim];
            }
          }
       for (unsigned  kdim = 0; kdim < dim; kdim++) { 
	   for (unsigned jdim = 0; jdim < dim; jdim++) { //diagonal blocks only
	     adv_uold_nablaphiadj_uadjnew[kdim]     += SolVAR_qp[SolPdeIndex[jdim]]  * phi_x_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][j] ;
	   }
	  }
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + adj_pos_begin][kdim + adj_pos_begin][i*nDofsVadj + j] += (   IRe*lap_jac_dadj_adj[kdim] 
                                                                                + advection_flag * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim] * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][j]   //c(delta_u, u_old, lambda_new)  diagonal blocks  ......phiadj_nablauold_uadjnew  
                                                                                + advection_flag * adv_uold_nablaphiadj_uadjnew[kdim] 	//c(u_old, delta_u, lambda_new)
                                                                               ) * AbsDetJxWeight_iqp;
               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
		  Jac[kdim + adj_pos_begin][off_kdim + adj_pos_begin][i*nDofsVadj + j] += ( advection_flag * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[off_kdim]][kdim] * phi_gss_fe[ SolFEType_Mat[off_kdim + adj_pos_begin] ][j]   //c(delta_u, u_old, lambda_new)  off-diagonal blocks  ......phiadj_nablauold_uadjnew   
                                                                                  ) * AbsDetJxWeight_iqp;
	  }
    } //j_dadj_adj loop
      
//BLOCK Pressure_adj
    for (unsigned j = 0; j < nDofsPadj; j++) {
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + adj_pos_begin][press_type_pos + adj_pos_begin][i*nDofsPadj + j] += -( phi_gss_fe[SolFEType_Mat[press_type_pos + adj_pos_begin]][j] * phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
	  }
    }//j_press_adj loop
  }//i_adj loop

//DIV_adj
  for (unsigned i = 0; i < nDofsPadj; i++) {
		double div_adj_dadj_qp = 0.;
      for (unsigned kdim = 0; kdim < dim; kdim++) {
	    div_adj_dadj_qp += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin ]][kdim] ;
      }
      Res[press_type_pos + adj_pos_begin][i] += ( (div_adj_dadj_qp) * phi_gss_fe[SolFEType_Mat[press_type_pos + adj_pos_begin]][i] ) * AbsDetJxWeight_iqp;
      for (unsigned j = 0; j < nDofsVadj; j++) {
        for (unsigned kdim = 0; kdim < dim; kdim++) {
            Jac[press_type_pos + adj_pos_begin][kdim + adj_pos_begin][i*nDofsVadj + j] += - ( phi_gss_fe[SolFEType_Mat[press_type_pos + adj_pos_begin]][i] * phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
        }
      }//j loop
  }//i_div_adj

//============ delta_adjoint row - END =============================================================================================

//============ delta_control row - BEGIN  ==================================================================================================
// delta_control
    for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) {
        
         for (unsigned i = 0; i < nDofsGctrl; i++) {
       Res[kdim + ctrl_pos_begin][i] += - penalty_outside_control_domain_boundary * ( (1 - control_node_flag_iel_jface[kdim][i]) *
                             (  Sol_eldofs_Mat[SolPdeIndex[kdim + ctrl_pos_begin]][i] -  PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY)  );              //enforce control zero outside the control boundary


// //DIAG BLOCK delta_control - control--------------------------------------------------------------------------------------
     for (unsigned j = 0; j < nDofsGctrl; j++) {
	    if (i == j) {
		Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsGctrl + j] += penalty_outside_control_domain_boundary * (1 - control_node_flag_iel_jface[kdim][i]);              //enforce control zero outside the control boundary
                  } //end i==j
        }//j_dctrl_ctrl loop
     }//i_ctrl loop
  
      }  //kdim

//============ delta_control row - END  ==================================================================================================
 

 
//============ delta_mu row - BEGIN  ============================================================================================

#if INEQ_FLAG != 0
  //MU
//************ Residual, BEGIN *********************
  for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
          
  for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim]; i++) {
      
       Res[mu_pos_begin + kdim][i]  +=  (- penalty_outside_control_domain_boundary) *  (1 - control_node_flag_iel_jface[kdim][i]) * (Sol_eldofs_Mat[mu_pos_begin + kdim][i] - 0.);
      
     }
  }
//************ Residual, END *********************

  //MU
//************ Jacobian, BEGIN *********************
  for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
    for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim]; i++) {
      for (unsigned j = 0; j < Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim]; j++) {
            if (i == j) {
               Jac[mu_pos_begin + kdim][mu_pos_begin + kdim][i * Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim] + j]  +=  penalty_outside_control_domain_boundary * (1 - control_node_flag_iel_jface[kdim][i]);
            }
         }
      }
  }
//************ Jacobian, END *********************
#endif

//============ delta_mu row - END  ============================================================================================
 
 
 
      }  // end quadrature point loop
      
//======================= VolumeLoop with Integration (and Zero boundary control outside Gamma_c) - END =====================================================    

    
    
//======================= Loop without Integration For Fake Dofs related to Compatibility Condition - BEGIN =====================================================    
    
        //============ delta_theta - theta row ==================================================================================================
  for (unsigned i = 0; i < nDofsThetactrl; i++) {
             /* if ( fake_theta_flag[i] != 1 ) */             Res[ theta_index ][i]    = - (1 - fake_theta_flag[i]) * ( theta_value_outside_fake_element - Sol_eldofs_Mat[SolPdeIndex[theta_index]][i]);  // Res_outer for the exact row (i.e. when fakeflag=1 , res =0(use Res_outer) and if not 1 this loop) and this is to take care of fake placement for the rest of dofs of theta values as 8
     for (unsigned j = 0; j < nDofsThetactrl; j++) {
			         if(i == j)  Jac[ theta_index ][ theta_index ][i*nDofsThetactrl + j] = (1 - fake_theta_flag[i]) * 1.; //likewise Jac_outer (actually Jac itself works in the correct placement) for bdry integral and this is for rest of dofs
             }//j_theta loop
        }//i_theta loop
   
 //============ delta_theta row ==================================================================================================
//======================= Loop without Integration For Fake Dofs related to Compatibility Condition - END =====================================================    



 //======================= From local to global - BEGIN =====================================================
    // FIRST ALL THE BLOCKS WITHOUT THETA ROW OR COLUMN - BEGIN 
    for(unsigned i_unk = 0; i_unk < theta_index; i_unk++) {
      RES->add_vector_blocked(Res[SolPdeIndex[i_unk]], L2G_dofmap_Mat[i_unk]);
        for(unsigned j_unk = 0; j_unk < theta_index; j_unk++) {
	  if(assembleMatrix) JAC->add_matrix_blocked( Jac[ SolPdeIndex[i_unk] ][ SolPdeIndex[j_unk] ], L2G_dofmap_Mat[i_unk], L2G_dofmap_Mat[j_unk]);
        }
    }
    // FIRST ALL THE BLOCKS WITHOUT THETA ROW OR COLUMN - END 
    
    // THEN THE BLOCKS WITH THETA ROW OR COLUMN - BEGIN
	/*delta_theta-theta*/    JAC->add_matrix_blocked( Jac[ SolPdeIndex[ theta_index ] ][ SolPdeIndex[ theta_index ] ], L2G_dofmap_Mat[ theta_index ], L2G_dofmap_Mat[ theta_index ]);
	    
     if (does_iel_contain_a_bdry_control_face == 1) {
	      for (unsigned kdim = 0; kdim < dim; kdim++) {
// // //                           /*delta_control*/       RES->add_vector_blocked(Res[SolPdeIndex[n_unknowns-2-kdim]], L2G_dofmap_Mat[n_unknowns-2-kdim]); ///@todo why was this here?!?
		if(assembleMatrix) {
                          /*delta_theta-control*/ JAC->add_matrix_blocked( Jac[ SolPdeIndex[ theta_index ] ][ SolPdeIndex[theta_index -1 -kdim] ], bdry_int_constr_pos_vec, L2G_dofmap_Mat[theta_index -1 -kdim]);
                          /*delta_control-theta*/ JAC->add_matrix_blocked( Jac[ SolPdeIndex[theta_index -1 -kdim] ][ SolPdeIndex[ theta_index ] ], L2G_dofmap_Mat[theta_index -1 -kdim], bdry_int_constr_pos_vec); 
		}
      }  //kdim
     }  //add control boundary element contributions
     
     
      if (does_iel_contain_a_bdry_control_face == 1) {
          /*delta_theta(bdry constr)*/         RES->add_vector_blocked(Res_outer, bdry_int_constr_pos_vec);
	  }
	  
     /* if (L2G_dofmap_Mat[n_unknowns-1][0] != bdry_int_constr_pos_vec[0]) */ /*delta_theta(fake)*/          RES->add_vector_blocked( Res[ SolPdeIndex[ theta_index ]],       L2G_dofmap_Mat[ theta_index ]);
    // THEN THE BLOCKS WITH THETA ROW OR COLUMN - END
	  

  
     if (print_algebra_local) {
         
         //extract only the ctrl range
         std::vector <unsigned>   Sol_n_el_dofs_Mat_vol_only_ctrl(1, 0);
	      for (unsigned kdim = 0; kdim < dim; kdim++) {
	       Sol_n_el_dofs_Mat_vol_only_ctrl[0] = Sol_n_el_dofs_Mat_vol[kdim + ctrl_pos_begin];
         assemble_jacobian<double,double>::print_element_residual(iel, Res[kdim + ctrl_pos_begin], Sol_n_el_dofs_Mat_vol_only_ctrl, 10, 5);
         assemble_jacobian<double,double>::print_element_jacobian(iel, Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin], Sol_n_el_dofs_Mat_vol_only_ctrl, 10, 5);
          }
    }
     
 //======================= From local to global - END =====================================================    
   
  } //end list of elements loop for each subdomain
  

#if INEQ_FLAG != 0
  //MU in res ctrl - BEGIN  ***********************************
ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::add_one_times_mu_res_ctrl(iproc,
                               ineq_flag,
                               ctrl_index_in_mat,
                               mu_index_in_mat,
                               SolIndex_Mat,
                               sol,
                               mlPdeSys,
                               pdeSys,
                               RES);
  //MU in res ctrl - END ***********************************
#endif
  
  
  
RES->close();
if (assembleMatrix) JAC->close();  /// This is needed for the parallel, when splitting the add part from the insert part!!!
      // ***************** ADD PART - END  *******************



#if INEQ_FLAG != 0
//   ***************** INSERT PART - BEGIN (must go AFTER the sum, clearly) *******************
    
     //MU

   for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
       
// -------
   geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
      
   geom_element_iel.set_elem_center_3d(iel, solType_coords);
// -------
   
// -------
    el_dofs_unknowns_vol(sol, msh, pdeSys, iel,
                         SolFEType_Mat,
                         SolIndex_Mat,
                         SolPdeIndex,
                         Sol_n_el_dofs_Mat_vol, 
                         Sol_eldofs_Mat, 
                         L2G_dofmap_Mat);
// -------

	if ( ctrl:: square_or_cube:: Domain_elements_containing_Gamma_control< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >:: volume_elem_contains_a_Gamma_control_face( geom_element_iel.get_elem_center_3d() ) ) {


    	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {

       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);

                
       if(  femus::face_is_a_Gamma_control_face< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES > ( el, iel, iface) ) {

       ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::update_active_set_flag_for_current_nonlinear_iteration_bdry
   (msh, sol,
    iel, iface,
    geom_element_iel.get_coords_at_dofs_bdry_3d(), 
    Sol_eldofs_Mat, 
    Sol_n_el_dofs_Mat_vol, 
    c_compl,
    mu_index_in_mat,
    ctrl_index_in_mat,
    solIndex_act_flag_sol,
    ctrl_lower,
    ctrl_upper,
    sol_actflag);
 

  ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::node_insertion_bdry(iel, iface, 
                      msh,
                      L2G_dofmap_Mat,
                      mu_index_in_mat, 
                      ctrl_index_in_mat,
                      Sol_eldofs_Mat,
                      Sol_n_el_dofs_Mat_vol,
                      sol_actflag, 
                      ctrl_lower, ctrl_upper,
                      ineq_flag,
                      c_compl,
                      JAC, 
                      RES,
                      assembleMatrix
                      );
  
             }
             
       }
     }


     //============= delta_ctrl-delta_mu row ===============================
 if (assembleMatrix) {
       for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
         JAC->matrix_set_off_diagonal_values_blocked(L2G_dofmap_Mat[ctrl_index_in_mat[kdim]],  L2G_dofmap_Mat[mu_index_in_mat[kdim]], ineq_flag * 1.);
       }
}

   }
   
   
  // ***************** INSERT PART - END *******************
#endif




  RES->close();
  if (assembleMatrix) JAC->close();
   



  
  
  
  print_global_residual_jacobian(print_algebra_global,
                                 ml_prob,
                                 mlPdeSys,
                                 pdeSys,
                                 RES,
                                 JAC,
                                 iproc,
                                 assembleMatrix);
  

  
}





double*  GetErrorNorm(const MultiLevelProblem & ml_prob, MultiLevelSolution* ml_sol, Solution* sol_coarser_prolongated) {
  
    static double ErrorNormArray[NO_OF_L2_NORMS+NO_OF_H1_NORMS];
    
  unsigned level = ml_sol->_mlMesh->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*     msh = ml_sol->_mlMesh->GetLevel(level);    // pointer to the mesh (level) object
  elem*     el  = msh->el;  // pointer to the elem object in msh (level)
  Solution* sol = ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)
  
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

 // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //geometry *******************************
  vector < vector < double > > coordX(dim);    // local coordinates

  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) { 
    coordX[k].reserve(max_size);
  }
   
  //geometry *******************************

 // solution variables *******************************************
  const int n_vars_state = dim+1;
  const int n_unknowns = 3*n_vars_state; //(2.*dim)+1; //state , adjoint of velocity terms and one pressure term
  const int press_type_pos = dim;
  const int state_pos_begin = 0;
  const int adj_pos_begin   = dim+1;
  const int ctrl_pos_begin   = 2*(dim+1);
  const int theta_index = press_type_pos + ctrl_pos_begin;
  
  vector < std::string > Solname_Mat(n_unknowns);  // const char Solname_Mat[4][8] = {"u_0","u_1","u_2","u_p"};
  Solname_Mat              [state_pos_begin+0] =                "u_0";
  Solname_Mat              [state_pos_begin+1] =                "u_1";
  if (dim == 3) Solname_Mat[state_pos_begin+2] =                "u_2";
  Solname_Mat              [state_pos_begin + press_type_pos] = "u_p";
  
  Solname_Mat              [adj_pos_begin + 0] =              "adj_0";
  Solname_Mat              [adj_pos_begin + 1] =              "adj_1";
  if (dim == 3) Solname_Mat[adj_pos_begin + 2] =              "adj_2";
  Solname_Mat              [adj_pos_begin + press_type_pos] = "adj_p";

  Solname_Mat              [ctrl_pos_begin + 0] =              "ctrl_0";
  Solname_Mat              [ctrl_pos_begin + 1] =              "ctrl_1";
  if (dim == 3) Solname_Mat[ctrl_pos_begin + 2] =              "ctrl_2";
  Solname_Mat              [ctrl_pos_begin + press_type_pos] = "theta";
  
  vector < unsigned > SolIndex_Mat(n_unknowns);  
  vector < unsigned > SolFEType_Mat(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolIndex_Mat[ivar]	= ml_sol->GetIndex        (Solname_Mat[ivar].c_str());
    SolFEType_Mat[ivar]	= ml_sol->GetSolutionType(SolIndex_Mat[ivar]);
  }

  vector < double > Sol_n_el_dofs(n_unknowns);
  
  //==========================================================================================
  // velocity ************************************
  vector < vector < double > > phi_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_x_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_xx_gss_fe(NFE_FAMS);
 
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_gss_fe[fe].reserve(max_size);
      phi_x_gss_fe[fe].reserve(max_size*dim);
     phi_xx_gss_fe[fe].reserve(max_size*(3*(dim-1)));
   }
  
  //=================================================================================================
  
  // quadratures ********************************
  double AbsDetJxWeight_iqp;
  
  
  //----------- dofs ------------------------------
  vector < vector < double > > Sol_eldofs_Mat(n_unknowns);
  vector < vector < double > > gradSol_eldofs_Mat(n_unknowns);
  
  vector < vector < double > > SolVAR_coarser_prol_eldofs(n_unknowns);
  vector < vector < double > > gradSolVAR_coarser_prol_eldofs(n_unknowns);


  for(int k = 0; k < n_unknowns; k++) {
    Sol_eldofs_Mat[k].reserve(max_size);
    gradSol_eldofs_Mat[k].reserve(max_size*dim); 
    
    SolVAR_coarser_prol_eldofs[k].reserve(max_size);
    gradSolVAR_coarser_prol_eldofs[k].reserve(max_size*dim);    
  }

  //------------ at quadrature points ---------------------
  vector < double > SolVAR_qp(n_unknowns);
  vector < double > SolVAR_coarser_prol_qp(n_unknowns);
  vector < vector < double > > gradSolVAR_qp(n_unknowns);
  vector < vector < double > > gradSolVAR_coarser_prol_qp(n_unknowns);
  for(int k = 0; k < n_unknowns; k++) {
      gradSolVAR_qp[k].reserve(max_size);  
      gradSolVAR_coarser_prol_qp[k].reserve(max_size);  
  }
      
   vector  < double > l2norm (NO_OF_L2_NORMS,0.);
  vector  < double > seminorm (NO_OF_H1_NORMS,0.);

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    
  // geometry *****************************
    short unsigned ielGeom = msh->GetElementType(iel);
    
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs

    for (unsigned  k = 0; k < dim; k++) {       coordX[k].resize(nDofsX);    }
  
    for (unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }
    
      // elem average point 
    vector < double > elem_center(dim);   
    for (unsigned j = 0; j < dim; j++) {  elem_center[j] = 0.;  }
  for (unsigned j = 0; j < dim; j++) {  
      for (unsigned i = 0; i < nDofsX; i++) {
         elem_center[j] += coordX[j][i];
       }
    }
    
   for (unsigned j = 0; j < dim; j++) { elem_center[j] = elem_center[j]/nDofsX; }
  //*************************************** 
  
  // geometry end *****************************
  
  
 // equation *****************************
    unsigned nDofsV = msh->GetElementDofNumber(iel, SolFEType_Mat[state_pos_begin]);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, SolFEType_Mat[state_pos_begin + press_type_pos]);    // number of solution element dofs
    
    unsigned nDofsVadj = msh->GetElementDofNumber(iel,SolFEType_Mat[adj_pos_begin]);    // number of solution element dofs
    unsigned nDofsPadj = msh->GetElementDofNumber(iel,SolFEType_Mat[adj_pos_begin + press_type_pos]);    // number of solution element dofs

   unsigned nDofsGctrl = msh->GetElementDofNumber(iel,SolFEType_Mat[ctrl_pos_begin]);
    unsigned nDofsThetactrl = msh->GetElementDofNumber(iel,SolFEType_Mat[theta_index] );

    unsigned nDofsVP = dim * nDofsV + nDofsP; //same for state and adjoint
    unsigned nDofsVPctrl = dim * nDofsGctrl + nDofsThetactrl; //control
   
    unsigned nDofsVP_tot = 2*nDofsVP + (nDofsVPctrl);
  // equation end *****************************


   //STATE###################################################################  
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType_Mat[k]);
	Sol_n_el_dofs[k]=ndofs_unk;
       Sol_eldofs_Mat[k].resize(ndofs_unk);
       SolVAR_coarser_prol_eldofs[k].resize(ndofs_unk);
    for (unsigned i = 0; i < ndofs_unk; i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType_Mat[k]);    // global to global mapping between solution node and solution dof // via local to global solution node
       Sol_eldofs_Mat[k][i] = (*sol->_Sol[SolIndex_Mat[k]])(solDof);      // global extraction and local storage for the solution
       SolVAR_coarser_prol_eldofs[k][i] = (*sol_coarser_prolongated->_Sol[SolIndex_Mat[k]])(solDof);      // global extraction and local storage for the solution
      }
    }
  //CTRL###################################################################

 
      // ********************** Gauss point loop *******************************
      for(unsigned i_qp = 0; i_qp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); i_qp++) {
	
 
      for(int fe=0; fe < NFE_FAMS; fe++) {
	msh->_finiteElement[ielGeom][fe]->Jacobian(coordX, i_qp, AbsDetJxWeight_iqp,phi_gss_fe[fe],phi_x_gss_fe[fe],phi_xx_gss_fe[fe]);
      }
         //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
  	msh->_finiteElement[ielGeom][BIQUADR_FE]->Jacobian(coordX,i_qp,AbsDetJxWeight_iqp,phi_gss_fe[BIQUADR_FE],phi_x_gss_fe[BIQUADR_FE],phi_xx_gss_fe[BIQUADR_FE]);

 //begin unknowns eval at gauss points ********************************
	for(unsigned unk = 0; unk < n_unknowns; unk++) {
	  SolVAR_qp[unk] = 0.;
	  SolVAR_coarser_prol_qp[unk] = 0.;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSolVAR_qp[unk][ivar2] = 0.; 
	    gradSolVAR_coarser_prol_qp[unk][ivar2] = 0.; 
	  }
    }
	  
	for(unsigned unk = 0; unk <  n_unknowns; unk++) {
	  for(unsigned i = 0; i < Sol_n_el_dofs[unk]; i++) {
	    SolVAR_qp[unk] += phi_gss_fe[ SolFEType_Mat[unk] ][i] * Sol_eldofs_Mat[unk][i];
	    SolVAR_coarser_prol_qp[unk] += phi_gss_fe[ SolFEType_Mat[unk] ][i] * SolVAR_coarser_prol_eldofs[unk][i];
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      gradSolVAR_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType_Mat[unk] ][i*dim+ivar2] * Sol_eldofs_Mat[unk][i]; 
	      gradSolVAR_coarser_prol_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType_Mat[unk] ][i*dim+ivar2] * SolVAR_coarser_prol_eldofs[unk][i]; 
	    }
	  }
	  
	}  
 //end unknowns eval at gauss points ********************************


	for(unsigned unk = 0; unk < n_unknowns; unk++) {
        l2norm[unk] += ( SolVAR_qp[unk] - SolVAR_coarser_prol_qp[unk] ) * ( SolVAR_qp[unk] - SolVAR_coarser_prol_qp[unk] ) * AbsDetJxWeight_iqp ; 
        
     }
    
    
	for(unsigned unk = 0; unk < dim; unk++) {
        for(int j = 0; j < dim; j++){
        seminorm[unk] += (gradSolVAR_qp[unk][j] - gradSolVAR_coarser_prol_qp[unk][j] ) * ( gradSolVAR_qp[unk][j] - gradSolVAR_coarser_prol_qp[unk][j] ) * AbsDetJxWeight_iqp ;
        seminorm[unk + dim] += (gradSolVAR_qp[unk + adj_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + adj_pos_begin][j] ) * ( gradSolVAR_qp[unk + adj_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + adj_pos_begin][j] ) * AbsDetJxWeight_iqp ;
        seminorm[unk + 2*dim] += (gradSolVAR_qp[unk + ctrl_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j] ) * ( gradSolVAR_qp[unk + ctrl_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j] ) * AbsDetJxWeight_iqp ;
        
        }
     }
    
    } // end gauss point loop
  } //end element loop for each process


    // add the norms of all processes
  NumericVector* norm_vec_inexact;
  norm_vec_inexact = NumericVector::build().release();
  norm_vec_inexact->init(msh->n_processors(), 1 , false, AUTOMATIC);

	for(unsigned unk = 0; unk < NO_OF_L2_NORMS; unk++) {
        norm_vec_inexact->set(iproc, l2norm[unk]);
        norm_vec_inexact->close();
        l2norm[unk] = norm_vec_inexact->l1_norm();
    }

	for(unsigned unk = 0; unk < NO_OF_H1_NORMS; unk++) {
        norm_vec_inexact->set(iproc, seminorm[unk]);
        norm_vec_inexact->close();
        seminorm[unk] = norm_vec_inexact->l1_norm();
    }


  delete norm_vec_inexact;
  
 
	for(unsigned unk = 0; unk < NO_OF_L2_NORMS; unk++) {
        ErrorNormArray[unk] = sqrt(l2norm[unk]);
    }
 	for(unsigned unk = 0; unk < NO_OF_H1_NORMS; unk++) {
        ErrorNormArray[unk + NO_OF_L2_NORMS] = sqrt(seminorm[unk]);
    }
  
   return ErrorNormArray;
  
  
}
