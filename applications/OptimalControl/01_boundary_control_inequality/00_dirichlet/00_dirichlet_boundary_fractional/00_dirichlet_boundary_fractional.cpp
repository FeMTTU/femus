#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "NumericVector.hpp"
#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"
#include "Assemble_unknown.hpp"

#include "ElemType.hpp"

//for reading additional fields from MED file (based on MED ordering)
#include "MED_IO.hpp"
//for reading additional fields from MED file (based on MED ordering)




#include  "../../../param.hpp"


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
           




 //Unknown definition  ==================
 const std::vector< Unknown >  provide_list_of_unknowns(const unsigned int dimension) {
     
     
  std::vector< FEFamily > feFamily;
  std::vector< FEOrder >   feOrder;

                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
 
                        feOrder.push_back(/*FIRST*/SECOND);  //same
                        feOrder.push_back(/*FIRST*/SECOND);  //same
                        feOrder.push_back(/*FIRST*/SECOND);
                        feOrder.push_back(/*FIRST*/SECOND);  //same
 

  assert( feFamily.size() == feOrder.size() );
 
 std::vector< Unknown >  unknowns(feFamily.size());

   unknowns[0]._name      = "state";
   unknowns[1]._name      = "control";
   unknowns[2]._name      = "adjoint";
   unknowns[3]._name      = "mu";

   unknowns[0]._is_sparse = true;
   unknowns[1]._is_sparse = IS_CTRL_FRACTIONAL_SOBOLEV ? false: true;
   unknowns[2]._is_sparse = true;
   unknowns[3]._is_sparse = true;
   
     for (unsigned int u = 0; u < unknowns.size(); u++) {
         
              unknowns[u]._fe_family  = feFamily[u];
              unknowns[u]._fe_order   = feOrder[u];
              unknowns[u]._time_order = 0;
              unknowns[u]._is_pde_unknown = true;
              
     }
 
   return unknowns;
     
}



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
        value = cost_functional::ElementTargetFlag(x);
    }
    else if(!strcmp(name, "ContReg")) {
        value = ctrl::Gamma_control::ControlDomainFlag_bdry(x);
    }
    else if(!strcmp(name, "act_flag")) {
        value = 0.;
    }


    return value;
}



///@todo notice that even if you set Dirichlet from the mesh file, here you can override it
bool Solution_set_boundary_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet = false; // = true; //dirichlet
  value = 0.;

  //************************state****************************************************

  if(!strcmp(name,"state")) {  //"state" corresponds to the first block row (u = q)

     boundary_conditions::Gamma_c_double_adjacent::ctrl_or_state_set_dirichlet_flags(ml_prob, faceName, x, dirichlet);

     boundary_conditions::Gamma_c_double_adjacent::ctrl_or_state_set_dirichlet_fixed_values(ml_prob, faceName, x, value);  //function that return value

   }
  
  
  //************************control****************************************************

  else if(!strcmp(name,"control")) {

    boundary_conditions::Gamma_c_double_adjacent::ctrl_or_state_set_dirichlet_flags(ml_prob, faceName, x, dirichlet);


    value = PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

  
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



void assemble_elliptic_dirichlet_control_pure_boundary(MultiLevelProblem& ml_prob);




int main(int argc, char** args) {

  // ======= Init ========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // ======= Problem - BEGIN  ==================
  MultiLevelProblem ml_prob;

  // ======= Problem, Files - BEGIN  ========================
  Files files; 
        files.CheckIODirectories(use_output_time_folder);
        files.RedirectCout(redirect_cout_to_file);

  // ======= Problem, Files ========================
  ml_prob.SetFilesHandler(&files);
  // ======= Problem, Files - END  ========================

  
  // ======= Problem, Mesh - BEGIN ==================

  // ======= Mesh, Coarse reading - BEGIN ==================
  MultiLevelMesh ml_mesh;

  
  const std::string input_file = mesh::input;
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_file;
  const std::string infile = mystream.str();
  const double Lref = 1.;
  
  const bool read_groups = true;
  const bool read_boundary_groups = true;
  

  
    ml_mesh.ReadCoarseMeshFileReadingBeforePartitioning(infile.c_str(), Lref, read_groups, read_boundary_groups);
    
// // //  BEGIN FillISvector
           ml_mesh.GetLevelZero(0)->dofmap_all_fe_families_initialize();

           std::vector < unsigned > elem_partition_from_mesh_file_to_new = ml_mesh.GetLevelZero(0)->elem_offsets();
  
           std::vector < unsigned > node_mapping_from_mesh_file_to_new = ml_mesh.GetLevelZero(0)->node_offsets();
           
           ml_mesh.GetLevelZero(0)->dofmap_all_fe_families_clear_ghost_dof_list_for_other_procs();

// // //   END FillISvector
       
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

  // ======= Solution, auxiliary; needed for Boundary of Boundary of Control region - BEFORE COARSE ERASING - BEGIN  ==================
  const std::string node_based_bdry_bdry_flag_name = NODE_BASED_BDRY_BDRY;
  const unsigned  steady_flag = 0;
  const bool      is_an_unknown_of_a_pde = false;
  
  const FEFamily node_bdry_bdry_flag_fe_fam = LAGRANGE;
  const FEOrder node_bdry_bdry_flag_fe_ord = SECOND;
  
  MultiLevelSolution * ml_sol_bdry_bdry_flag = ctrl::Gamma_control_equation_fractional::bdry_bdry_flag(files,
                                                              ml_mesh, 
                                                              infile,
                                                              node_mapping_from_mesh_file_to_new,
                                                              node_based_bdry_bdry_flag_name,
                                                              steady_flag,
                                                              is_an_unknown_of_a_pde,
                                                              node_bdry_bdry_flag_fe_fam,
                                                              node_bdry_bdry_flag_fe_ord);
  
  // ======= Solution, auxiliary - END  ==================

  
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
  std::vector< Unknown > unknowns = provide_list_of_unknowns( ml_mesh.GetDimension() );

  ml_sol.AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions);
  
  for (unsigned int u = 0; u < unknowns.size(); u++)  { 
      ml_sol.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order, unknowns[u]._time_order, unknowns[u]._is_pde_unknown);
      ml_sol.Initialize(unknowns[u]._name.c_str(), Solution_set_initial_conditions, & ml_prob);
      ml_sol.GenerateBdc(unknowns[u]._name.c_str(), (unknowns[u]._time_order == 0) ? "Steady" : "Time_dependent", & ml_prob);
  }
  // ======= Solutions that are Unknowns - END ==================
  

  // ======= Solutions that are not Unknowns - BEGIN  ==================
  ml_sol.AddSolution("TargReg", DISCONTINUOUS_POLYNOMIAL, ZERO, steady_flag, is_an_unknown_of_a_pde);
  ml_sol.Initialize("TargReg",     Solution_set_initial_conditions, & ml_prob);
  
  
  ml_sol.AddSolution("ContReg", DISCONTINUOUS_POLYNOMIAL, ZERO, steady_flag, is_an_unknown_of_a_pde);
  ml_sol.Initialize("ContReg",     Solution_set_initial_conditions, & ml_prob);
  
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
   ml_sol.Initialize(act_set_flag_name[0].c_str(), Solution_set_initial_conditions, & ml_prob);
  // ******** active flag - END 
  
  
 ctrl::Gamma_control_equation_fractional::bdry_bdry_flag_copy_and_delete(ml_prob,
                                ml_sol,
                                ml_mesh, 
                                erased_levels,
                                ml_sol_bdry_bdry_flag,
                                node_based_bdry_bdry_flag_name,
                                steady_flag,
                                is_an_unknown_of_a_pde,
                                node_bdry_bdry_flag_fe_fam,
                                node_bdry_bdry_flag_fe_ord);
  

  // ======= Solutions that are not Unknowns - END  ==================
  
  
  //== Solution: CHECK SOLUTION FE TYPES between Unknowns and Not Unknowns - BEGIN  ==
  if ( ml_sol.GetSolutionType("control") != ml_sol.GetSolutionType(act_set_flag_name[0].c_str())) abort();
  //== Solution: CHECK SOLUTION FE TYPES between Unknowns and Not Unknowns - END ==
  
  // ======= Problem, Solution - END ==================

  

  // ======= Problem, Quad Rule - BEGIN  ========================
  std::vector< std::string > fe_quad_rule_vec;
  fe_quad_rule_vec.push_back("seventh");
  fe_quad_rule_vec.push_back("eighth");

  ml_prob.SetQuadratureRuleAllGeomElemsMultiple(fe_quad_rule_vec);
  ml_prob.set_all_abstract_fe_multiple();
  // ======= Problem, Quad Rule - END  ========================

  // ======= Problem, System - BEGIN ========================
  NonLinearImplicitSystemWithPrimalDualActiveSetMethod & system_opt = ml_prob.add_system < NonLinearImplicitSystemWithPrimalDualActiveSetMethod > ("BoundaryControl");
  
       // ======= Not an Unknown, but needed in the System with PDAS ========================
  system_opt.SetActiveSetFlagName(act_set_flag_name);    //MU
       // ======= Not an Unknown, but needed in the System with PDAS ========================
  
       // ======= System Unknowns ========================
  for (unsigned int u = 0; u < unknowns.size(); u++) { system_opt.AddSolutionToSystemPDE(unknowns[u]._name.c_str());  }
       // ======= System Unknowns ========================

  system_opt.SetAssembleFunction(assemble_elliptic_dirichlet_control_pure_boundary);

// *****************
  const unsigned n_components_state = n_components_ctrl;
  std::vector<std::string> state_vars(n_components_state);  state_vars[0] = "state";
  std::vector<std::string> ctrl_vars(n_components_ctrl);     ctrl_vars[0] = "control";
  
  system_opt.set_state_vars(state_vars);
  system_opt.set_ctrl_vars(ctrl_vars);
  
  system_opt.SetDebugNonlinear(true);
  system_opt.SetDebugFunction(ctrl::cost_functional::compute_cost_functional_regularization_bdry);
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
  ctrl::cost_functional::compute_cost_functional_regularization_bdry(ml_prob, 0, 0, state_vars, ctrl_vars);
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




  
  

//This Opt system is characterized by the following ways of setting matrix values:
// Add_values (Mat or Vec) in the volume loop
// Add_values (Mat or Vec) in the boundary loop
// Insert_values (Mat or Vec) in the boundary loop
// Insert_values (Mat or Vec) in the volume loop
// Insert_values (Mat or Vec) outside all loops
// We're going to split the two parts and add a close() at the end of each


void assemble_elliptic_dirichlet_control_pure_boundary(MultiLevelProblem & ml_prob) {
    
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = & ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("BoundaryControl");
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->el;
  
  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*            JAC = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors();

  constexpr bool print_algebra_global = false;
  constexpr bool print_algebra_local = false;
  
  

  //=============== Geometry ========================================
   unsigned solType_coords = FE_DOMAIN;
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
  CurrentElem < double > geom_element_jel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
  
  std::vector< double > normal(space_dim, 0.);
  
 //************ geom ***************************************  
  std::vector < double > coord_at_qp_bdry(space_dim);
  
  std::vector < double > phi_coords;
  std::vector < double > phi_coords_x;

  phi_coords.reserve(max_size);
  phi_coords_x.reserve(max_size * space_dim);
  

 //********************* state *********************** 
 //*************************************************** 
  std::vector <double> phi_u;
  std::vector <double> phi_u_x;

  phi_u.reserve(max_size);
  phi_u_x.reserve(max_size * space_dim);
  
  
  //boundary state shape functions
  std::vector <double> phi_u_bdry;  
  std::vector <double> phi_u_x_bdry; 

  phi_u_bdry.reserve(max_size);
  phi_u_x_bdry.reserve(max_size * space_dim);
  
 //***************************************************  
 //***************************************************  

  
 //********************** adjoint ********************
 //*************************************************** 
  std::vector <double> phi_adj;
  std::vector <double> phi_adj_x;

  phi_adj.reserve(max_size);
  phi_adj_x.reserve(max_size * space_dim);
 

  //boundary adjoint shape functions  
  std::vector <double> phi_adj_bdry;  
  std::vector <double> phi_adj_x_bdry; 

  phi_adj_bdry.reserve(max_size);
  phi_adj_x_bdry.reserve(max_size * space_dim);

  
  //volume shape functions at boundary
  std::vector <double> phi_adj_vol_at_bdry;
  std::vector <double> phi_adj_x_vol_at_bdry;
  phi_adj_vol_at_bdry.reserve(max_size);
  phi_adj_x_vol_at_bdry.reserve(max_size * space_dim);
  
  std::vector <double> sol_adj_x_vol_at_bdry_gss(space_dim);
 //*************************************************** 
 //*************************************************** 

  
 //********************* bdry cont *******************
 //*************************************************** 
  std::vector <double> phi_ctrl_bdry;  
  std::vector <double> phi_ctrl_x_bdry; 

  phi_ctrl_bdry.reserve(max_size);
  phi_ctrl_x_bdry.reserve(max_size * space_dim);
 //*************************************************** 

    const unsigned int n_components_ctrl = 1;
    const unsigned int first_loc_comp_ctrl = 0;


    //MU
  //************** act flag ****************************   
    std::vector <unsigned int> solIndex_act_flag_sol(n_components_ctrl);

    ctrl::ctrl_inequality::store_act_flag_in_old(mlPdeSys, ml_sol, sol, solIndex_act_flag_sol);
  
  
  //********* variables for ineq constraints *****************
     std::vector <std::vector < double/*int*/ > > sol_actflag(n_components_ctrl);    //flag for active set
     std::vector <std::vector < double > > ctrl_lower(n_components_ctrl);  
     std::vector <std::vector < double > > ctrl_upper(n_components_ctrl);  
  
     for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) {
          sol_actflag[kdim].reserve(max_size);
           ctrl_lower[kdim].reserve(max_size);
           ctrl_upper[kdim].reserve(max_size);
     }
     
     const int ineq_flag = INEQ_FLAG;
     const double c_compl = C_COMPL;
  //***************************************************  
  //MU
  
  
//***************************************************
//********* WHOLE SET OF VARIABLES ******************
    const unsigned int n_unknowns = mlPdeSys->GetSolPdeIndex().size();
    const unsigned int n_quantities = ml_sol->GetSolutionSize();

//***************************************************
    enum Pos_in_matrix {pos_mat_state = 0, pos_mat_ctrl, pos_mat_adj, pos_mat_mu}; //these are known at compile-time 
                    ///@todo these are the positions in the MlSol object or in the Matrix? I'd say the matrix, but we have to check where we use it...

                    
    assert(pos_mat_state   == mlPdeSys->GetSolPdeIndex("state"));
    assert(pos_mat_ctrl    == mlPdeSys->GetSolPdeIndex("control"));
    assert(pos_mat_adj     == mlPdeSys->GetSolPdeIndex("adjoint"));
    assert(pos_mat_mu      == mlPdeSys->GetSolPdeIndex("mu"));
//***************************************************

    std::vector < std::string > Solname_Mat(n_unknowns);  //this coincides with Pos_in_matrix
    Solname_Mat[0] = "state";
    Solname_Mat[1] = "control";
    Solname_Mat[2] = "adjoint";
    Solname_Mat[3] = "mu";

 //***************************************************
   enum Pos_in_Sol {pos_sol_state = 0, pos_sol_ctrl, pos_sol_adj, pos_sol_mu, pos_sol_targreg, pos_sol_contreg, pos_sol_actflag}; //these are known at compile-time 

        assert(pos_sol_actflag == solIndex_act_flag_sol[0]);
        
        
  std::vector<unsigned int> ctrl_index_in_mat(1);  ctrl_index_in_mat[0] = pos_mat_ctrl;
  std::vector<unsigned int>   mu_index_in_mat(1);    mu_index_in_mat[0] = pos_mat_mu;

//***************************************************
    
 //***************************************************
    std::vector < std::string > Solname_quantities(n_quantities);
    
        for(unsigned ivar=0; ivar < Solname_quantities.size(); ivar++) {
            Solname_quantities[ivar] = ml_sol->GetSolutionName(ivar);
        }
 //***************************************************
        
    std::vector < unsigned > SolIndex_Mat(n_unknowns);      //should have Mat order
    std::vector < unsigned > SolFEType_Mat(n_unknowns);       //should have Mat order
    std::vector < unsigned > SolPdeIndex(n_unknowns);     //should have Mat order, of course

    std::vector < unsigned > SolIndex_quantities(n_quantities);      //should have Sol order
    std::vector < unsigned > SolFEType_quantities(n_quantities);     //should have Sol order
    std::vector < unsigned > Sol_n_el_dofs_quantities_vol(n_quantities); //should have Sol order
 
  

    for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
        SolIndex_Mat[ivar]    = ml_sol->GetIndex        (Solname_Mat[ivar].c_str());
        SolFEType_Mat[ivar]   = ml_sol->GetSolutionType(SolIndex_Mat[ivar]);
        SolPdeIndex[ivar] = mlPdeSys->GetSolPdeIndex(Solname_Mat[ivar].c_str());
    }
    
    for(unsigned ivar=0; ivar < n_quantities; ivar++) {
        SolIndex_quantities[ivar]    = ml_sol->GetIndex        (Solname_quantities[ivar].c_str());
        SolFEType_quantities[ivar]   = ml_sol->GetSolutionType(SolIndex_quantities[ivar]);
    }    

    std::vector < unsigned > Sol_n_el_dofs_Mat_vol(n_unknowns);    //should have Mat order

//***************************************************
    std::vector < std::vector < double > >  Sol_eldofs_Mat(n_unknowns);  //should have Mat order
    for(int k = 0; k < n_unknowns; k++) {        Sol_eldofs_Mat[k].reserve(max_size);    }


    //----------- quantities (at dof objects) ------------------------------
    std::vector< int >       L2G_dofmap_Mat_AllVars;
    L2G_dofmap_Mat_AllVars.reserve( n_unknowns * max_size );
    std::vector < std::vector < int > >     L2G_dofmap_Mat(n_unknowns);     //should have Mat order
    for(int i = 0; i < n_unknowns; i++) {
        L2G_dofmap_Mat[i].reserve(max_size);
    }
    
 //*************************************************** 
  std::vector < double > Res;   Res.reserve( n_unknowns*max_size);                         //should have Mat order
  std::vector < double > Jac;   Jac.reserve( n_unknowns*max_size * n_unknowns*max_size);   //should have Mat order
 //*************************************************** 

 
 //********************* DATA ************************ 
  const double u_des = cost_functional::DesiredTarget();
  const double alpha = ALPHA_CTRL_BDRY;
  const double beta  = BETA_CTRL_BDRY;
  const double penalty_outside_control_domain_boundary = PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY;       // penalty for zero control outside Gamma_c and zero mu outside Gamma_c
  const double penalty_dirichlet_bc_u_equal_q = PENALTY_DIRICHLET_BC_U_EQUAL_Q_BOUNDARY;         //penalty for u=q
 //*************************************************** 
  
  RES->zero();  
  if (assembleMatrix)  JAC->zero();

 //*************************************************** 
// ---
     std::vector < std::vector < double > >  JacI_iqp(space_dim);
     std::vector < std::vector < double > >  Jac_iqp(dim);
    for (unsigned d = 0; d < Jac_iqp.size(); d++) {   Jac_iqp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp.size(); d++) { JacI_iqp[d].resize(dim); }
    
    double detJac_iqp;
  double weight_iqp = 0.;
// ---

// ---
    std::vector < std::vector < double > >  JacI_iqp_bdry(space_dim);
    std::vector < std::vector < double > >  Jac_iqp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_iqp_bdry.size(); d++) {   Jac_iqp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp_bdry.size(); d++) { JacI_iqp_bdry[d].resize(dim-1); }
    
    double detJac_iqp_bdry;
  double weight_iqp_bdry = 0.;

// ---
    
//*************************************************** 

    //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all;
  ml_prob.get_all_abstract_fe_multiple(elem_all);
//*************************************************** 
//***************************************************
  
  const unsigned dim_bdry = dim - 1;
  
  const double s_frac = S_FRAC;

  const double check_limits = 1.;//1./(1. - s_frac); // - s_frac;

  


  //--- quadrature rules -------------------
  constexpr unsigned qrule_i = QRULE_I;
  constexpr unsigned qrule_j = QRULE_J;
  constexpr unsigned qrule_k = QRULE_K;
  //----------------------


    
  if ( IS_CTRL_FRACTIONAL_SOBOLEV ) {
  
     ctrl::Gamma_control_equation_fractional::control_eqn_bdry_fractional_sobolev_differentiability_index(iproc,
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
                    weight_iqp_bdry,
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
  
   ctrl::Gamma_control_equation::control_eqn_bdry_integer_sobolev_differentiability_index(iproc,
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
                    weight_iqp_bdry,
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
                    RHS_ONE,
                    qrule_i,
                    //-----------
                    print_algebra_local
                    ) ;
  
  }
  
  

//**************************                    
// AAA do not close this because later they will be filled with the rest of the system!!!      
//    JAC->close();
//   RES->close();

   //print JAC and RES to files
   //You print only what is being sent to the global matrix. If nothing is sent, nothing is printed                 
   // I will keep this print here for later because it highlights what positions were filled in the matrix
   // If I remove everything above here, it seems like the very last diagonal position is filled... why? from where?
//   JAC->close(); //JAC->zero();
//   const unsigned nonlin_iter = 0/*mlPdeSys->GetNonlinearIt()*/;
//     assemble_jacobian< double, double >::print_global_jacobian(assembleMatrix, ml_prob, JAC, nonlin_iter);
// //     assemble_jacobian< double, double >::print_global_residual(ml_prob, RES, nonlin_iter);
//   std::cout << "****************************" << std::endl;
//**************************                    

  
                    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

// -------
    geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    const short unsigned ielGeom = geom_element_iel.geom_type();

    geom_element_iel.set_elem_center_3d(iel, solType_coords);
// -------
    

 //***************************************************
   el_dofs_unknowns_vol(sol, msh, pdeSys, iel,
                        SolFEType_Mat,
                        SolIndex_Mat,
                        SolPdeIndex,
                        Sol_n_el_dofs_Mat_vol, 
                        Sol_eldofs_Mat,  
                        L2G_dofmap_Mat);
  
   el_dofs_quantities_vol(sol, msh, iel, SolFEType_quantities, Sol_n_el_dofs_quantities_vol); 
  //***************************************************
   
 
    unsigned int nDof_max          = ElementJacRes<double>::compute_max_n_dofs(Sol_n_el_dofs_Mat_vol);
    
    unsigned int sum_Sol_n_el_dofs = ElementJacRes<double>::compute_sum_n_dofs(Sol_n_el_dofs_Mat_vol);

    Res.resize(sum_Sol_n_el_dofs);
    std::fill(Res.begin(), Res.end(), 0.);

    Jac.resize(sum_Sol_n_el_dofs * sum_Sol_n_el_dofs);
    std::fill(Jac.begin(), Jac.end(), 0.);
    
    L2G_dofmap_Mat_AllVars.resize(0);
      for (unsigned  k = 0; k < n_unknowns; k++)     L2G_dofmap_Mat_AllVars.insert(L2G_dofmap_Mat_AllVars.end(), L2G_dofmap_Mat[k].begin(), L2G_dofmap_Mat[k].end());
 //***************************************************

      
  //************* set target domain flag **************
   int target_flag = 0;
   target_flag = cost_functional::ElementTargetFlag(geom_element_iel.get_elem_center_3d());
 //*************************************************** 
   

 //************ set control flag *********************
   std::vector< std::vector< int > > control_node_flag = 
       ctrl::Gamma_control::is_dof_associated_to_Gamma_control_equation(msh, ml_sol, & ml_prob, iel, geom_element_iel, solType_coords, Solname_Mat, SolFEType_Mat, Sol_n_el_dofs_Mat_vol, pos_mat_ctrl, n_components_ctrl);
  //*************************************************** 
 

	if ( ctrl::Gamma_control::volume_elem_contains_a_Gamma_control_face(geom_element_iel.get_elem_center_3d()) ) {
	  
	  std::vector<double> normal(space_dim, 0.);
	       
	  // loop on faces of the current element

	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
// -------
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();

       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, iface);    
// -------
       
// -------
       std::vector<unsigned int> Sol_el_n_dofs_current_face(n_unknowns); ///@todo the active flag is not an unknown! However, if I add to the quantities something that has higher order, then I will have error below when I take the max number of dofs on a face. Since the act_flag must have the same FE family as the control, then I can limit this array to n_unknowns instead of n_quantities

       for (unsigned  k = 0; k < Sol_el_n_dofs_current_face.size(); k++) {
                 if (SolFEType_quantities[k] < 3) Sol_el_n_dofs_current_face[k] = msh->GetElementFaceDofNumber(iel, iface, SolFEType_quantities[k]);  ///@todo fix this absence
       }
       
       const unsigned nDof_max_bdry = ElementJacRes<double>::compute_max_n_dofs(Sol_el_n_dofs_current_face);
// -------
       
		
	    if( ctrl::Gamma_control::face_is_a_Gamma_control_face(msh->el, iel, iface) ) {
              
 
//========= initialize gauss quantities on the boundary ============================================
                double sol_ctrl_bdry_gss = 0.;
                double sol_adj_bdry_gss = 0.;
                std::vector<double> sol_ctrl_x_bdry_gss(space_dim);   

//========= initialize gauss quantities on the boundary ============================================
		
        const unsigned n_gauss_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussPointsNumber();
        
    
		for(unsigned ig_bdry = 0; ig_bdry < n_gauss_bdry; ig_bdry++) {
    
    elem_all[qrule_i][ielGeom_bdry][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), ig_bdry, Jac_iqp_bdry, JacI_iqp_bdry, detJac_iqp_bdry, space_dim);
	elem_all[qrule_i][ielGeom_bdry][solType_coords]->compute_normal(Jac_iqp_bdry, normal);
    
    weight_iqp_bdry = detJac_iqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussWeightsPointer()[ig_bdry];

    elem_all[qrule_i][ielGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(ig_bdry, JacI_iqp_bdry, phi_ctrl_bdry, phi_ctrl_x_bdry, boost::none, space_dim);
    elem_all[qrule_i][ielGeom_bdry][SolFEType_quantities[pos_sol_state]]->shape_funcs_current_elem(ig_bdry, JacI_iqp_bdry, phi_u_bdry, phi_u_x_bdry,  boost::none, space_dim);
    elem_all[qrule_i][ielGeom_bdry][SolFEType_quantities[pos_sol_adj]]  ->shape_funcs_current_elem(ig_bdry, JacI_iqp_bdry, phi_adj_bdry, phi_adj_x_bdry,  boost::none, space_dim);


    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv_vol_at_bdry_new(geom_element_iel.get_coords_at_dofs_3d(), ig_bdry, iface, Jac_iqp/*not_needed_here*/, JacI_iqp, detJac_iqp/*not_needed_here*/, space_dim);
    elem_all[qrule_i][ielGeom][SolFEType_quantities[pos_sol_adj]]->shape_funcs_vol_at_bdry_current_elem(ig_bdry, iface, JacI_iqp, phi_adj_vol_at_bdry, phi_adj_x_vol_at_bdry, boost::none, space_dim);
     
//     msh->_finiteElement[ielGeom][SolFEType_quantities[pos_sol_adj]]->fill_volume_shape_funcs_at_boundary_quadrature_points_on_current_elem(geom_element_iel.get_coords_at_dofs(), geom_element_iel.get_coords_at_dofs_bdry_3d(), iface, ig_bdry, phi_adj_vol_at_bdry, phi_adj_x_vol_at_bdry);

		  
//========== compute gauss quantities on the boundary ===============================================
		  sol_ctrl_bdry_gss = 0.;
                  std::fill(sol_ctrl_x_bdry_gss.begin(), sol_ctrl_x_bdry_gss.end(), 0.);
                  
		      for (unsigned int i_bdry = 0; i_bdry < phi_ctrl_bdry.size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_ctrl_bdry_gss +=  Sol_eldofs_Mat[pos_mat_ctrl][i_vol] * phi_ctrl_bdry[i_bdry];
                            for (unsigned int d = 0; d < space_dim; d++) {
			      sol_ctrl_x_bdry_gss[d] += Sol_eldofs_Mat[pos_mat_ctrl][i_vol] * phi_ctrl_x_bdry[i_bdry * space_dim + d];
			    }
		      }
		      
		      
		  sol_adj_bdry_gss = 0.;
		      for (unsigned int i_bdry = 0; i_bdry <  phi_adj_bdry.size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_adj_bdry_gss  +=  Sol_eldofs_Mat[pos_mat_adj][i_vol] * phi_adj_bdry[i_bdry];
              }		      
		      
//=============== grad dot n for residual ========================================= 
//     compute gauss quantities on the boundary through VOLUME interpolation
           std::fill(sol_adj_x_vol_at_bdry_gss.begin(), sol_adj_x_vol_at_bdry_gss.end(), 0.);
		      for (unsigned int iv = 0; iv < Sol_n_el_dofs_Mat_vol[pos_mat_adj]; iv++)  {
			
         for (unsigned int d = 0; d < space_dim; d++) {
			      sol_adj_x_vol_at_bdry_gss[d] += Sol_eldofs_Mat[pos_mat_adj][iv] * phi_adj_x_vol_at_bdry[iv * space_dim + d];
			    }
		      }  
		      
    double grad_adj_dot_n_res = 0.;
        for(unsigned d=0; d < space_dim; d++) {
	  grad_adj_dot_n_res += sol_adj_x_vol_at_bdry_gss[d] * normal[d];  
	}
//=============== grad dot n  for residual =========================================       

//========== compute gauss quantities on the boundary ================================================

		  // *** phi_i loop ***
		  for(unsigned i_bdry=0; i_bdry < nDof_max_bdry; i_bdry++) {
              
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);

                 double lap_rhs_dctrl_ctrl_bdry_gss_i = 0.;
                 for (unsigned d = 0; d < space_dim; d++) {
                       if ( i_vol < Sol_n_el_dofs_Mat_vol[pos_mat_ctrl] )  lap_rhs_dctrl_ctrl_bdry_gss_i +=  phi_ctrl_x_bdry[i_bdry * space_dim + d] * sol_ctrl_x_bdry_gss[d];
                 }
                 
		 
//============ Bdry Residuals - BEGIN  ==================	
                Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_state, i_vol) ] +=
                    - control_node_flag[first_loc_comp_ctrl][i_vol] * penalty_dirichlet_bc_u_equal_q     * KEEP_ADJOINT_PUSH * (   Sol_eldofs_Mat[pos_mat_state][i_vol] - Sol_eldofs_Mat[pos_mat_ctrl][i_vol] )
                    - control_node_flag[first_loc_comp_ctrl][i_vol] *  weight_iqp_bdry * KEEP_ADJOINT_PUSH * grad_adj_dot_n_res * phi_u_bdry[i_bdry];   // u = q


                Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_ctrl, i_vol) ]  += 
                    - control_node_flag[first_loc_comp_ctrl][i_vol] *  weight_iqp_bdry *
                          (    IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY * alpha * phi_ctrl_bdry[i_bdry]  * sol_ctrl_bdry_gss
							+  IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY * beta * lap_rhs_dctrl_ctrl_bdry_gss_i 
							                            - KEEP_ADJOINT_PUSH * grad_adj_dot_n_res * phi_ctrl_bdry[i_bdry]
// 							                           -         phi_ctrl_bdry[i_bdry]*sol_adj_bdry_gss // for Neumann control
						  );  //boundary optimality condition
                Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_adj, i_vol) ]  += 0.; 
//============ Bdry Residuals - END  ==================    
		    
		    for(unsigned j_bdry=0; j_bdry < nDof_max_bdry; j_bdry ++) {
		         unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, iface, j_bdry);

//============ Bdry Jacobians - BEGIN ==================	


// FIRST BLOCK ROW
//============ u = q - BEGIN ===========================	    
                 
if ( i_vol == j_vol )  {
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_state, i_vol, j_vol) ] += 
		     penalty_dirichlet_bc_u_equal_q * KEEP_ADJOINT_PUSH * ( control_node_flag[first_loc_comp_ctrl][i_vol]) * ( 1.);
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_ctrl, i_vol, j_vol) ]  += 
		     penalty_dirichlet_bc_u_equal_q * KEEP_ADJOINT_PUSH * ( control_node_flag[first_loc_comp_ctrl][i_vol]) * (-1.);
		}
//============ u = q - END ===========================

		    

// SECOND BLOCK ROW
//=========== boundary control eqn =============	    

//========block delta_control / control ========
              double  lap_mat_dctrl_ctrl_bdry_gss = 0.;
		      for (unsigned d = 0; d < space_dim; d++) {  lap_mat_dctrl_ctrl_bdry_gss += phi_ctrl_x_bdry[i_bdry * space_dim + d] * phi_ctrl_x_bdry[j_bdry * space_dim + d];    }

          
              Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_ctrl, pos_mat_ctrl, i_vol, j_vol) ] 
			+=  control_node_flag[first_loc_comp_ctrl][i_vol] *  weight_iqp_bdry * ( IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY * alpha * phi_ctrl_bdry[i_bdry] * phi_ctrl_bdry[j_bdry] 
			                                              +  IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY * beta *  lap_mat_dctrl_ctrl_bdry_gss);   
    
		   
//============ Bdry Jacobians - END ==================	
				
	      }  //end j loop
	      
//===================loop over j in the VOLUME (while i is in the boundary)	      
	for(unsigned j=0; j < nDof_max; j ++) {
		      
  //=============== grad dot n  =========================================    
    double grad_adj_dot_n_mat = 0.;
        for(unsigned d = 0; d< space_dim; d++) {
	  grad_adj_dot_n_mat += phi_adj_x_vol_at_bdry[j * space_dim + d] * normal[d];  //notice that the convention of the orders x y z is different from vol to bdry
	}
//=============== grad dot n  =========================================    

		      
//==========block delta_control/adjoint ========
		     Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_ctrl, pos_mat_adj, i_vol, j) ]  += 
		     control_node_flag[first_loc_comp_ctrl][i_vol] * (-1.) * weight_iqp_bdry * KEEP_ADJOINT_PUSH * grad_adj_dot_n_mat * phi_ctrl_bdry[i_bdry];    		      

//==========block delta_state/adjoint ========
		     Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_adj, i_vol, j) ] += 
		     control_node_flag[first_loc_comp_ctrl][i_vol] * (1.) * weight_iqp_bdry * KEEP_ADJOINT_PUSH * grad_adj_dot_n_mat * phi_u_bdry[i_bdry];  
		      
		    }   //end loop i_bdry // j_vol
	      
	      

		  }  //end i loop
		}  //end ig_bdry loop
		
	      }    //end if control face
	      
	  }    //end loop over faces
	  
	} //end if control element flag
	

//========= gauss value quantities on the volume ==============  
	double sol_u_gss = 0.;
	double sol_adj_gss = 0.;
	std::vector<double> sol_u_x_gss(space_dim);     std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::vector<double> sol_adj_x_gss(space_dim);   std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
//=============================================== 
 
 
      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom).GetGaussPointsNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), ig, Jac_iqp, JacI_iqp, detJac_iqp, space_dim);
    weight_iqp = detJac_iqp * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom).GetGaussWeightsPointer()[ig];

    elem_all[qrule_i][ielGeom][SolFEType_quantities[pos_sol_state]]->shape_funcs_current_elem(ig, JacI_iqp, phi_u, phi_u_x, boost::none, space_dim);
    elem_all[qrule_i][ielGeom][SolFEType_quantities[pos_sol_adj]]  ->shape_funcs_current_elem(ig, JacI_iqp, phi_adj, phi_adj_x, boost::none, space_dim);

          
	sol_u_gss = 0.;
	sol_adj_gss = 0.;
	std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	
	for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[pos_mat_state]; i++) {
	                                                sol_u_gss      += Sol_eldofs_Mat[pos_mat_state][i] * phi_u[i];
                   for (unsigned d = 0; d < space_dim; d++)   sol_u_x_gss[d] += Sol_eldofs_Mat[pos_mat_state][i] * phi_u_x[i * space_dim + d];
          }
	
	for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[pos_mat_adj]; i++) {
	                                                sol_adj_gss      += Sol_eldofs_Mat[pos_mat_adj][i] * phi_adj[i];
                   for (unsigned d = 0; d < space_dim; d++)   sol_adj_x_gss[d] += Sol_eldofs_Mat[pos_mat_adj][i] * phi_adj_x[i * space_dim + d];
        }

//==========FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {

              double laplace_rhs_du_adj_i = 0.;
              double laplace_rhs_dadj_u_i = 0.;
              double laplace_rhs_du_u_i = 0.;

              for (unsigned kdim = 0; kdim < space_dim; kdim++) {
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_state] )  laplace_rhs_du_adj_i  +=  phi_u_x   [i * space_dim + kdim] * sol_adj_x_gss[kdim];
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_adj] )    laplace_rhs_dadj_u_i  +=  phi_adj_x [i * space_dim + kdim] * sol_u_x_gss[kdim];
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_state] )  laplace_rhs_du_u_i    +=  phi_u_x   [i * space_dim + kdim] * sol_u_x_gss[kdim];
	      }

//============ Volume residuals - BEGIN  ==================	    
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_state, i) ] += - weight_iqp * (
          - laplace_rhs_du_adj_i
          + target_flag *
          #if COST_FUNCTIONAL_TYPE == 0
                phi_u[i] * ( sol_u_gss - u_des)
          #elif COST_FUNCTIONAL_TYPE == 1
                laplace_rhs_du_u_i
          #endif
          );
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_ctrl, i) ]  += - penalty_outside_control_domain_boundary * ( (1 - control_node_flag[first_loc_comp_ctrl][i]) * (  Sol_eldofs_Mat[pos_mat_ctrl][i] - PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY)  );
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_adj, i) ]   += - weight_iqp * (-1.) * (laplace_rhs_dadj_u_i);
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_mu, i) ]    += - penalty_outside_control_domain_boundary * ( (1 - control_node_flag[first_loc_comp_ctrl][i]) * (  Sol_eldofs_Mat[pos_mat_mu][i] - 0.)  );  //MU
//============  Volume Residuals - END ==================	    
	      
	      
//============ Volume Jacobians - BEGIN  ==================	    
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {
                
              double laplace_mat_dadj_u = 0.;
              double laplace_mat_du_adj = 0.;
              double laplace_mat_du_u = 0.;

              for (unsigned kdim = 0; kdim < space_dim; kdim++) {
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_adj]     && j < Sol_n_el_dofs_Mat_vol[pos_mat_state] )     laplace_mat_dadj_u        +=  (phi_adj_x [i * space_dim + kdim] * phi_u_x   [j * space_dim + kdim]);
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_state]   && j < Sol_n_el_dofs_Mat_vol[pos_mat_adj] )   laplace_mat_du_adj        +=  (phi_u_x   [i * space_dim + kdim] * phi_adj_x [j * space_dim + kdim]);
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_state]   && j < Sol_n_el_dofs_Mat_vol[pos_mat_state] )     laplace_mat_du_u          += (phi_u_x  [i * space_dim + kdim] * phi_u_x  [j * space_dim + kdim]);
		
	      }

              //============ delta_state row ============================
              // BLOCK delta_state / state
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_state, i, j) ]  += weight_iqp  * target_flag *
		#if COST_FUNCTIONAL_TYPE == 0
             phi_u[i] * phi_u[j]
        #elif COST_FUNCTIONAL_TYPE == 1
             laplace_mat_du_u
        #endif
        ;

              //BLOCK delta_state / adjoint
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_adj, i, j) ]  += weight_iqp * (-1.) * laplace_mat_du_adj;
	      
	      
              //=========== delta_control row ===========================
              //enforce control zero outside the control boundary
	      if ( i==j )
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_ctrl, pos_mat_ctrl, i, j) ]  += penalty_outside_control_domain_boundary * ( (1 - control_node_flag[first_loc_comp_ctrl][i]));    /*weight * phi_adj[i]*phi_adj[j]*/
              
	      //=========== delta_adjoint row ===========================
	      // BLOCK delta_adjoint / state
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_adj, pos_mat_state, i, j) ]  += weight_iqp * (-1.) * laplace_mat_dadj_u;

	      
	      //============= delta_mu row ===============================
	        if ( i==j )   
		  Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_mu, pos_mat_mu, i, j) ]  += penalty_outside_control_domain_boundary * ( (1 - control_node_flag[first_loc_comp_ctrl][i]));    //MU
          
	         } // end phi_j loop
           } // endif assemble_matrix
//============ Volume Jacobians - END ==================	    

        } // end phi_i loop
        
      } // end gauss point loop

  

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    RES->add_vector_blocked(Res, L2G_dofmap_Mat_AllVars);

    if (assembleMatrix) {
      JAC->add_matrix_blocked(Jac, L2G_dofmap_Mat_AllVars, L2G_dofmap_Mat_AllVars);
    }
    
    
    //========== dof-based part, without summation
 
     if (print_algebra_local) {
         assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat_vol, 10, 5);
         assemble_jacobian<double,double>::print_element_jacobian(iel, Jac, Sol_n_el_dofs_Mat_vol, 10, 5);
     }
     
     
  } //end element loop for each process
  


  //MU in res ctrl - BEGIN  ***********************************
ctrl::ctrl_inequality::add_one_times_mu_res_ctrl(iproc,
                               ineq_flag,
                               ctrl_index_in_mat,
                               mu_index_in_mat,
                               SolIndex_Mat,
                               sol,
                               mlPdeSys,
                               pdeSys,
                               RES);
  //MU in res ctrl - END ***********************************

    
  // ***************** END ASSEMBLY - ADD PART *******************

RES->close();
if (assembleMatrix) JAC->close();  /// This is needed for the parallel, when splitting the add part from the insert part!!!


//   ***************** INSERT PART - BEGIN (must go AFTER the sum, clearly) *******************
 /// @todo One very important thing to consider: we have some PENALTIES that were set before during the SUMMATION part.
 // Now, if we do INSERT, we may end up OVERWRITING certain values, SUCH AS THOSE PENALTIES!!!
 // So you have to be very careful here!
    
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

	if ( ctrl::Gamma_control::volume_elem_contains_a_Gamma_control_face( geom_element_iel.get_elem_center_3d() ) ) {


    	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {

       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);

                
       if(  ctrl::Gamma_control::face_is_a_Gamma_control_face( el, iel, iface) ) {

       ctrl::ctrl_inequality::update_active_set_flag_for_current_nonlinear_iteration_bdry
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
 

  ctrl::ctrl_inequality::node_insertion_bdry(iel, iface, 
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

  RES->close();
  if (assembleMatrix) JAC->close();  ///@todo is it needed? I think so
   
  
  print_global_residual_jacobian(print_algebra_global,
                                 ml_prob,
                                 mlPdeSys,
                                 pdeSys,
                                 RES,
                                 JAC,
                                 iproc,
                                 assembleMatrix);
  

  return;
}







