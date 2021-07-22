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



#define FACE_FOR_CONTROL        2  /* 1-2 x coords, 3-4 y coords, 5-6 z coords */




#include "../../../param.hpp"


//***** Implementation-related ****************** 
#define IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY    0
//**************************************

//***** Quadrature-related ****************** 
// for integrations in the same element
#define Nsplit 0
#define Quadrature_split_index  0

//for semi-analytical integration in the unbounded domain
#define N_DIV_UNBOUNDED  10

#define QRULE_I   0
#define QRULE_J   1
#define QRULE_K   QRULE_I
//**************************************

//***** Operator-related ****************** 
  #define RHS_ONE             1.
  #define KEEP_ADJOINT_PUSH   0
#define IS_CTRL_FRACTIONAL_SOBOLEV   1
#define S_FRAC 0.5

#define NORM_GIR_RAV  0

#if NORM_GIR_RAV == 0

  #define OP_L2       0
  #define OP_H1       0
  #define OP_Hhalf    1

  #define UNBOUNDED   0

  #define USE_Cns     1

#elif NORM_GIR_RAV == 1 

  #define OP_L2       1
  #define OP_H1       0
  #define OP_Hhalf    1

  #define UNBOUNDED   0

  #define USE_Cns     0
#endif
//**************************************


//***** Domain-related ****************** 
#define EX_1        GAMMA_CONTROL_LOWER
#define EX_2        GAMMA_CONTROL_UPPER
#define EY_1        0.
#define EY_2        1.   ///@todo  see here

#define DOMAIN_EX_1 0
#define DOMAIN_EX_2 1
//**************************************


#define FE_DOMAIN  2 //with 0 it only works in serial, you must put 2 to make it work in parallel...: that's because when you fetch the dofs from _topology you get the wrong indices

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

using namespace femus;



 //Unknown definition  ==================
 const std::vector< Unknown >  provide_list_of_unknowns(const unsigned int dimension) {
     
     
  std::vector< FEFamily > feFamily;
  std::vector< FEOrder >   feOrder;

                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
 
                        feOrder.push_back(/*FIRST*/SECOND);
                        feOrder.push_back(/*FIRST*/SECOND);
                        feOrder.push_back(/*FIRST*/SECOND);
                        feOrder.push_back(/*FIRST*/SECOND);
 

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
        value = ElementTargetFlag(x);
    }
    else if(!strcmp(name, "ContReg")) {
        value = ControlDomainFlag_bdry(x);
    }
    else if(!strcmp(name, "act_flag")) {
        value = 0.;
    }


    return value;
}



///@todo notice that even if you set Dirichlet from the mesh file, here you can override it
bool Solution_set_boundary_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet; // = true; //dirichlet
  value = 0.;

  if(!strcmp(name,"control")) {
      
  if (faceName == FACE_FOR_CONTROL) {
     if (x[ axis_direction_Gamma_control(faceName) ] > GAMMA_CONTROL_LOWER - 1.e-5 && x[ axis_direction_Gamma_control(faceName) ] < GAMMA_CONTROL_UPPER + 1.e-5)  { 
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

  else if(!strcmp(name,"state")) {  //"state" corresponds to the first block row (u = q)
      
  if (faceName == FACE_FOR_CONTROL) {
      
     if (x[ axis_direction_Gamma_control(faceName) ] > GAMMA_CONTROL_LOWER - 1.e-5 && x[ axis_direction_Gamma_control(faceName) ] < GAMMA_CONTROL_UPPER + 1.e-5) { 
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

  else if(!strcmp(name,"mu")) {
      
    dirichlet = false;

  }
  
  else { dirichlet = true; }
  
//     if(!strcmp(name,"adjoint")) { 
//     dirichlet = false;
//   }

  
  return dirichlet;
}


void ComputeIntegral(const MultiLevelProblem& ml_prob);

void AssembleOptSys(MultiLevelProblem& ml_prob);


int main(int argc, char** args) {

  // ======= Init ========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
  // ======= Files ========================
  const bool use_output_time_folder = false;
  const bool redirect_cout_to_file = true;
  Files files; 
        files.CheckIODirectories(use_output_time_folder);
        files.RedirectCout(redirect_cout_to_file);

  // ======= Quad Rule ========================
  //right now only one quadrature rule is used, so there is no possibility of quadrature point offset to try to avoid numerical cancellation
  //quadr rule order
  /*const*/ std::vector< std::string > fe_quad_rule_vec;
  fe_quad_rule_vec.push_back("seventh");
  fe_quad_rule_vec.push_back("eighth");

  // ======= Mesh  ==================
  MultiLevelMesh ml_mesh;

  
//   std::string input_file = "parametric_square_1x1.med";
//   std::string input_file = "parametric_square_1x2.med";
//   std::string input_file = "parametric_square_2x2.med";
//   std::string input_file = "parametric_square_4x5.med";
//   std::string input_file = "Mesh_3_groups_with_bdry_nodes.med";
  std::string input_file = "Mesh_3_groups_with_bdry_nodes_coarser.med";
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_file;
  const std::string infile = mystream.str();
  const double Lref = 1.;
  
  const bool read_groups = true;
  const bool read_boundary_groups = true;
  
// // // =================================================================  
// // // ================= Mesh: UNPACKING ReadCoarseMesh - BEGIN ================================================  
// // // =================================================================  
//   ml_mesh.ReadCoarseMesh(infile.c_str(), fe_quad_rule_vec[0].c_str(), Lref, read_groups, read_boundary_groups);

//   ml_mesh.ReadCoarseMeshOnlyFileReading(infile.c_str(), Lref, read_groups, read_boundary_groups);
  
    ml_mesh.ReadCoarseMeshOnlyFileReadingBeforePartitioning(infile.c_str(), Lref, read_groups, read_boundary_groups);
//     ml_mesh.GetLevelZero(0)->Partition();
       std::vector < unsigned > partition;
       ml_mesh.GetLevelZero(0)->PartitionForElements(partition);
//        ml_mesh.GetLevelZero(0)->FillISvector(partition);
       
           ml_mesh.GetLevelZero(0)->initialize_elem_dof_offsets();
           std::vector < unsigned > mapping;
           ml_mesh.GetLevelZero(0)->build_elem_offsets_and_dofs_element_based(partition, mapping);
           ml_mesh.GetLevelZero(0)->from_mesh_file_to_femus_node_partition_mapping_ownSize(partition, mapping);
           ml_mesh.GetLevelZero(0)->end_building_dof_offset_biquadratic_and_coord_reordering(mapping);
           ml_mesh.GetLevelZero(0)->ghost_nodes_search();
           ml_mesh.GetLevelZero(0)->complete_dof_offsets();
       
       partition.resize(0);
    ml_mesh.GetLevelZero(0)->ReadCoarseMeshAfterPartitioning();

  ml_mesh.BuildElemType(fe_quad_rule_vec[0].c_str());
  ml_mesh.AllocateAllLevels();
// // // =================================================================  
// // // ================= Mesh: UNPACKING ReadCoarseMesh - END ===============================================  
// // // =================================================================
  


  // ======= Mesh: REFINING ========================
  const unsigned numberOfUniformLevels = N_UNIFORM_LEVELS;
  const unsigned erased_levels = N_ERASED_LEVELS;
  unsigned numberOfSelectiveLevels = 0;
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // ======= Solution, auxiliary - BEFORE COARSE ERASING  ==================
  const unsigned  steady_flag = 0;
  const bool      is_an_unknown_of_a_pde = false;
  MultiLevelSolution * ml_sol_aux = new MultiLevelSolution(&ml_mesh);
  
  ml_sol_aux->SetWriter(VTK);
  ml_sol_aux->GetWriter()->SetDebugOutput(true);
  
  const std::string node_based_bdry_flag_name = "node_based_bdry_flag";
  const FEFamily node_flag_fe_fam = LAGRANGE;
  const FEOrder node_flag_fe_ord = SECOND;
  ml_sol_aux->AddSolution(node_based_bdry_flag_name.c_str(), node_flag_fe_fam, node_flag_fe_ord, steady_flag, is_an_unknown_of_a_pde);
  ml_sol_aux->Initialize(node_based_bdry_flag_name.c_str());
      // ======= COARSE READING and REFINEMENT ========================
  ml_sol_aux->GetSolutionLevel(0)->GetSolutionName(node_based_bdry_flag_name.c_str()) = MED_IO(*ml_mesh.GetLevel(0)).node_based_flag_read_from_file(infile, mapping);
  for(unsigned l = 1; l < ml_mesh.GetNumberOfLevels(); l++) {
     ml_sol_aux->RefineSolution(l);
  }
  
  std::vector < std::string > variablesToBePrinted_aux;
  variablesToBePrinted_aux.push_back("all");
  for(unsigned l = 0; l < ml_mesh.GetNumberOfLevels(); l++) {
  ml_sol_aux->GetWriter()->Write(l+1, "aux", files.GetOutputPath(), "", "biquadratic", variablesToBePrinted_aux);
   }

  
  // ======= Mesh: COARSE ERASING ========================
  ml_mesh.EraseCoarseLevels(erased_levels);
  ml_mesh.PrintInfo();
  
  
  // ======= Solution  ==================
  MultiLevelSolution ml_sol(&ml_mesh);
  
  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);

  // ======= Problem  ==================
  MultiLevelProblem ml_prob(&ml_sol);
  
  ml_prob.SetFilesHandler(&files);
  ml_prob.SetQuadratureRuleAllGeomElemsMultiple(fe_quad_rule_vec);
  ml_prob.set_all_abstract_fe_multiple();

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
  
  //MU
  const std::string act_set_flag_name = "act_flag";
  const unsigned int act_set_fake_time_dep_flag = 2;
  ml_sol.AddSolution(act_set_flag_name.c_str(), LAGRANGE, /*FIRST*/SECOND, act_set_fake_time_dep_flag, is_an_unknown_of_a_pde);
  ml_sol.Initialize(act_set_flag_name.c_str(), Solution_set_initial_conditions, & ml_prob);
  //MU
  
  //---- node_based_bdry_flag ------
  ml_sol.AddSolution(node_based_bdry_flag_name.c_str(), node_flag_fe_fam, node_flag_fe_ord, steady_flag, is_an_unknown_of_a_pde);
  ml_sol.Initialize(node_based_bdry_flag_name.c_str(), Solution_set_initial_conditions, & ml_prob);
  // copy ml_sol_aux at the non-removed levels into ml_sol
  for(unsigned l = 0; l < ml_mesh.GetNumberOfLevels(); l++) {
      *(ml_sol.GetSolutionLevel(l)->_Sol[ ml_sol.GetIndex(node_based_bdry_flag_name.c_str()) ]) =
      *(ml_sol_aux->GetSolutionLevel(l + erased_levels)->_Sol[ ml_sol_aux->GetIndex(node_based_bdry_flag_name.c_str()) ]);
  }
  delete ml_sol_aux;
  //---- node_based_bdry_flag ------

  //-- CHECK SOLUTION FE TYPES --------
  if ( ml_sol.GetSolutionType("control") != ml_sol.GetSolutionType("state")) abort();
  if ( ml_sol.GetSolutionType("control") != ml_sol.GetSolutionType("mu")) abort();
  if ( ml_sol.GetSolutionType("control") != ml_sol.GetSolutionType(act_set_flag_name.c_str())) abort();
  //-- CHECK SOLUTION FE TYPES --------
  
  

  // ======= Solutions that are not Unknowns - END  ==================

  
  // ======= System - BEGIN ========================
  NonLinearImplicitSystemWithPrimalDualActiveSetMethod & system = ml_prob.add_system < NonLinearImplicitSystemWithPrimalDualActiveSetMethod > ("BoundaryControl");
  
  system.SetAssembleFunction(AssembleOptSys);

// *****************
  system.SetDebugNonlinear(true);
  system.SetDebugFunction(ComputeIntegral);  //weird error if I comment this line, I expect nothing to happen but something in the assembly gets screwed up in memory I guess
// *****************
  
       // ======= System Unknowns ========================
  for (unsigned int u = 0; u < unknowns.size(); u++)  system.AddSolutionToSystemPDE(unknowns[u]._name.c_str());  
  
       // ======= Not an Unknown, but needed in the System with PDAS ========================
  system.SetActiveSetFlagName(act_set_flag_name);    //MU
       // ======= Not an Unknown, but needed in the System with PDAS ========================

 
   system.init();  /// I need to put this init before, later I will remove it   /// @todo it seems like you cannot do this init INSIDE A FUNCTION... understand WHY!
 
  set_dense_pattern_for_unknowns(system, unknowns);
  // ======= System  - END ========================

  //   initialize and solve the system
  system.init();

  //----
  std::ostringstream sp_out_base2; sp_out_base2 << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "sp_";
  sp_out_base2 << "after_second_init_";
  
  unsigned n_levels = ml_mesh.GetNumberOfLevels();
  system._LinSolver[n_levels - 1]->sparsity_pattern_print_nonzeros(sp_out_base2.str(), "on");
  system._LinSolver[n_levels - 1]->sparsity_pattern_print_nonzeros(sp_out_base2.str(), "off");
  //----

  system.MGsolve();
//   double totalAssemblyTime = 0.;
//   system.nonlinear_solve_single_level(MULTIPLICATIVE, totalAssemblyTime, 0, 0);
//   system.assemble_call_before_boundary_conditions(2);
  
  // ======= Print ========================
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
  ml_sol.GetWriter()->Write(files.GetOutputPath(), "biquadratic", variablesToBePrinted);

  return 0;
}




  
  

//This Opt system is characterized by the following ways of setting matrix values:
// Add_values (Mat or Vec) in the volume loop
// Add_values (Mat or Vec) in the boundary loop
// Insert_values (Mat or Vec) in the boundary loop
// Insert_values (Mat or Vec) in the volume loop
// Insert_values (Mat or Vec) outside all loops
// We're going to split the two parts and add a close() at the end of each


void AssembleOptSys(MultiLevelProblem& ml_prob) {
    
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("BoundaryControl");
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->el;
  
  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             KK = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors();

  constexpr bool print_algebra_global = true;
  constexpr bool print_algebra_local = false;
  
  

  //=============== Geometry ========================================
   unsigned solType_coords = FE_DOMAIN;
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
  CurrentElem < double > geom_element_jel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
  
  std::vector< double > normal(space_dim, 0.);
  
 //************ geom ***************************************  
  vector < double > coord_at_qp_bdry(space_dim);
  
  vector <double> phi_coords;
  vector <double> phi_coords_x;
  vector <double> phi_coords_xx; 

  phi_coords.reserve(max_size);
  phi_coords_x.reserve(max_size * space_dim);
  phi_coords_xx.reserve(max_size * dim2);
  
 //********************* bdry cont *******************
 //*************************************************** 
  vector <double> phi_coords_iqp_bdry;  
  vector <double> phi_coords_x_iqp_bdry; 

  phi_coords_iqp_bdry.reserve(max_size);
  phi_coords_x_iqp_bdry.reserve(max_size * space_dim);

 //*************************************************** 

 //********************* state *********************** 
 //*************************************************** 
  vector <double> phi_u;
  vector <double> phi_u_x;
  vector <double> phi_u_xx;

  phi_u.reserve(max_size);
  phi_u_x.reserve(max_size * space_dim);
  phi_u_xx.reserve(max_size * dim2);
  
  
  //boundary state shape functions
  vector <double> phi_u_bdry;  
  vector <double> phi_u_x_bdry; 

  phi_u_bdry.reserve(max_size);
  phi_u_x_bdry.reserve(max_size * space_dim);
  
 //***************************************************  
 //***************************************************  

  
 //********************** adjoint ********************
 //*************************************************** 
  vector <double> phi_adj;
  vector <double> phi_adj_x;
  vector <double> phi_adj_xx;

  phi_adj.reserve(max_size);
  phi_adj_x.reserve(max_size * space_dim);
  phi_adj_xx.reserve(max_size * dim2);
 

  //boundary adjoint shape functions  
  vector <double> phi_adj_bdry;  
  vector <double> phi_adj_x_bdry; 

  phi_adj_bdry.reserve(max_size);
  phi_adj_x_bdry.reserve(max_size * space_dim);

  
  //volume shape functions at boundary
  vector <double> phi_adj_vol_at_bdry;
  vector <double> phi_adj_x_vol_at_bdry;
  phi_adj_vol_at_bdry.reserve(max_size);
  phi_adj_x_vol_at_bdry.reserve(max_size * space_dim);
  vector <double> sol_adj_x_vol_at_bdry_gss(space_dim);
 //*************************************************** 
 //*************************************************** 

  
 //********************* bdry cont *******************
 //*************************************************** 
  vector <double> phi_ctrl_bdry;  
  vector <double> phi_ctrl_x_bdry; 

  phi_ctrl_bdry.reserve(max_size);
  phi_ctrl_x_bdry.reserve(max_size * space_dim);
 //*************************************************** 

  //MU
  //************** act flag ****************************   
  unsigned int solIndex_act_flag_sol; 
  unsigned int solFEType_act_flag_sol;
  store_act_flag_in_old(mlPdeSys, ml_sol, sol,
                        solIndex_act_flag_sol, //this becomes a vector
                        solFEType_act_flag_sol //remove this one, only Index
                       );
  
  
  //********* variables for ineq constraints *****************
  vector < double/*int*/ >  sol_actflag;   sol_actflag.reserve(max_size); //flag for active set
  vector < double >  ctrl_lower;   ctrl_lower.reserve(max_size);
  vector < double >  ctrl_upper;   ctrl_upper.reserve(max_size);
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

    vector < std::string > Solname_Mat(n_unknowns);  //this coincides with Pos_in_matrix
    Solname_Mat[0] = "state";
    Solname_Mat[1] = "control";
    Solname_Mat[2] = "adjoint";
    Solname_Mat[3] = "mu";

 //***************************************************
   enum Pos_in_Sol {pos_sol_state = 0, pos_sol_ctrl, pos_sol_adj, pos_sol_mu, pos_sol_targreg, pos_sol_contreg, pos_sol_actflag}; //these are known at compile-time 

        assert(pos_sol_actflag == solIndex_act_flag_sol);
//***************************************************
    
 //***************************************************
    vector < std::string > Solname_quantities(n_quantities);
    
        for(unsigned ivar=0; ivar < Solname_quantities.size(); ivar++) {
            Solname_quantities[ivar] = ml_sol->GetSolutionName(ivar);
        }
 //***************************************************
        
    vector < unsigned > SolIndex_Mat(n_unknowns);      //should have Mat order
    vector < unsigned > SolFEType_Mat(n_unknowns);       //should have Mat order
    vector < unsigned > SolPdeIndex(n_unknowns);     //should have Mat order, of course

    vector < unsigned > SolIndex_quantities(n_quantities);      //should have Sol order
    vector < unsigned > SolFEType_quantities(n_quantities);     //should have Sol order
    vector < unsigned > Sol_n_el_dofs_quantities_vol(n_quantities); //should have Sol order
 
  

    for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
        SolIndex_Mat[ivar]    = ml_sol->GetIndex        (Solname_Mat[ivar].c_str());
        SolFEType_Mat[ivar]   = ml_sol->GetSolutionType(SolIndex_Mat[ivar]);
        SolPdeIndex[ivar] = mlPdeSys->GetSolPdeIndex(Solname_Mat[ivar].c_str());
    }
    
    for(unsigned ivar=0; ivar < n_quantities; ivar++) {
        SolIndex_quantities[ivar]    = ml_sol->GetIndex        (Solname_quantities[ivar].c_str());
        SolFEType_quantities[ivar]   = ml_sol->GetSolutionType(SolIndex_quantities[ivar]);
    }    

    vector < unsigned > Sol_n_el_dofs_Mat_vol(n_unknowns);    //should have Mat order

//***************************************************
    vector < vector < double > >  sol_eldofs_Mat(n_unknowns);  //should have Mat order
    for(int k = 0; k < n_unknowns; k++) {        sol_eldofs_Mat[k].reserve(max_size);    }


    //----------- quantities (at dof objects) ------------------------------
    std::vector< int >       L2G_dofmap_Mat_AllVars;
    L2G_dofmap_Mat_AllVars.reserve( n_unknowns * max_size );
    vector < vector < int > >     L2G_dofmap_Mat(n_unknowns);     //should have Mat order
    for(int i = 0; i < n_unknowns; i++) {
        L2G_dofmap_Mat[i].reserve(max_size);
    }
    
 //*************************************************** 
  std::vector < double > Res;   Res.reserve( n_unknowns*max_size);                         //should have Mat order
  std::vector < double > Jac;   Jac.reserve( n_unknowns*max_size * n_unknowns*max_size);   //should have Mat order
 //*************************************************** 

 
 //********************* DATA ************************ 
  const double u_des = DesiredTarget();
  const double alpha = ALPHA_CTRL_BDRY;
  const double beta  = BETA_CTRL_BDRY;
  const double penalty_outside_control_boundary = 1.e50;       // penalty for zero control outside Gamma_c and zero mu outside Gamma_c
  const double penalty_strong_bdry = 1.e20;  // penalty for boundary equation on Gamma_c
  const double penalty_ctrl = 1.e10;         //penalty for u=q
 //*************************************************** 
  
  RES->zero();  
  if (assembleMatrix)  KK->zero();

 //*************************************************** 
// ---
     std::vector < std::vector < double > >  JacI_iqp(space_dim);
     std::vector < std::vector < double > >  Jac_iqp(dim);
    for (unsigned d = 0; d < Jac_iqp.size(); d++) {   Jac_iqp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp.size(); d++) { JacI_iqp[d].resize(dim); }
    
    double detJac_iqp;
// ---

// ---
    std::vector < std::vector < double > >  JacI_iqp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_iqp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_iqp_bdry.size(); d++) {   Jac_iqp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp_bdry.size(); d++) { JacI_iqp_bdry[d].resize(dim-1); }
    
    double detJac_iqp_bdry;
// ---
    
// ---
     std::vector < std::vector < double > >  JacI_jqp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_jqp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_jqp_bdry.size(); d++) {   Jac_jqp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_jqp_bdry.size(); d++) { JacI_jqp_bdry[d].resize(dim-1); }
    
    double detJac_jqp_bdry;
    
  double weight_iqp = 0.;
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

  double C_ns = 2 * (1 - USE_Cns) + USE_Cns * s_frac * pow(2, (2. * s_frac)) * tgamma((dim_bdry + 2. * s_frac) / 2.) / (pow(M_PI, dim_bdry / 2.) * tgamma(1 -  s_frac)) ;
  
//*************************************************** 
  unsigned n_max = pow(2,dim_bdry);
  std::vector < double > extremes(n_max);
  extremes[0] = EX_1;
  extremes[1] = EX_2;
  if(dim_bdry == 2){
    extremes[2] = EY_1;
    extremes[3] = EY_2;
  }
  
  std::vector < std::vector < double > > ex_control(dim); // control zone extremes. Built as: 2D) ( x1 y1 ) ( x2 y2 )    3D) ( x1 y1 z1 ) ( x2 y2 z2 ) ( x3 y3 z3 ) 
  for(unsigned d = 0; d < dim; d++) {
      ex_control[d].reserve(n_max);
  }
  int control_xyz = (FACE_FOR_CONTROL - 1) / 2;
  bool ctrl_min_max = (FACE_FOR_CONTROL - 1) % 2;
  for(unsigned d = 0; d < dim; d++) {
    for(unsigned n_e = 0; n_e < n_max; n_e++){
      if(control_xyz == d) ex_control[d][n_e] =  (ctrl_min_max)? DOMAIN_EX_2:DOMAIN_EX_1;
      else ex_control[d][n_e] = extremes[n_e];
    }
  }
//*************************************************** 


  //--- quadrature rules -------------------
  constexpr unsigned qrule_i = QRULE_I;
  constexpr unsigned qrule_j = QRULE_J;
  constexpr unsigned qrule_k = QRULE_K;
  //----------------------

    const unsigned int n_components_ctrl = 1;
    const unsigned int first_loc_comp_ctrl = 0;


    
  if ( IS_CTRL_FRACTIONAL_SOBOLEV ) {
  
     control_eqn_bdry_fractional(iproc,
                   nprocs,
                    ml_prob,
                    ml_sol,
                    sol,
                    msh,
                    pdeSys,
                    //-----------
                    geom_element_iel,
                    geom_element_jel,
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
                    sol_eldofs_Mat,  
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
                    phi_coords_iqp_bdry,
                    phi_coords_x_iqp_bdry, 
                    //-----------
                    Jac_jqp_bdry,
                    JacI_jqp_bdry,
                    detJac_jqp_bdry,
//                     weight_bdry,
//                     phi_ctrl_bdry,
//                     phi_ctrl_x_bdry, 
                    //-----------
                    n_components_ctrl,
                    pos_mat_ctrl,
                    pos_sol_ctrl,
                    IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY,
                    //-----------
                    KK,
                    RES,
                    assembleMatrix,
                    //-----------
                    alpha,
                    beta,     
                    s_frac,
                    check_limits,
                    C_ns,
                    OP_Hhalf,
                    OP_L2,
                    RHS_ONE,
                    UNBOUNDED,
                    //-----------
                    ex_control,
                    //-----------
                    qrule_i,
                    qrule_j,
                    qrule_k,
                    Nsplit,
                    Quadrature_split_index,
                    N_DIV_UNBOUNDED
                    );
                    
  }
  
  else {
  
   control_eqn_bdry(iproc,
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
                    sol_eldofs_Mat,  
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
                    KK,
                    RES,
                    assembleMatrix,
                    //-----------
                    alpha,
                    beta,
                    RHS_ONE,
                    qrule_i
                    ) ;
  
  }
  
  

//**************************                    
// AAA do not close this because later they will be filled with the rest of the system!!!      
//    KK->close();
//   RES->close();

   //print JAC and RES to files
   //You print only what is being sent to the global matrix. If nothing is sent, nothing is printed                 
   // I will keep this print here for later because it highlights what positions were filled in the matrix
   // If I remove everything above here, it seems like the very last diagonal position is filled... why? from where?
//   KK->close(); //KK->zero();
//   const unsigned nonlin_iter = 0/*mlPdeSys->GetNonlinearIt()*/;
//     assemble_jacobian< double, double >::print_global_jacobian(assembleMatrix, ml_prob, KK, nonlin_iter);
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
                        sol_eldofs_Mat,  
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
   target_flag = ElementTargetFlag(geom_element_iel.get_elem_center_3d());
 //*************************************************** 
   

 //************ set control flag *********************
   std::vector< std::vector< int > > control_node_flag = 
       is_dof_associated_to_boundary_control_equation(msh, ml_sol, & ml_prob, iel, geom_element_iel, solType_coords, Solname_Mat, SolFEType_Mat, Sol_n_el_dofs_Mat_vol, pos_mat_ctrl, n_components_ctrl);
  //*************************************************** 
 

	if ( volume_elem_contains_a_boundary_control_face(geom_element_iel.get_elem_center_3d()) ) {
	  
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
       
		
	    if( face_is_a_boundary_control_face(msh->el, iel, iface) ) {
              
 
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
			
			sol_ctrl_bdry_gss +=  sol_eldofs_Mat[pos_mat_ctrl][i_vol] * phi_ctrl_bdry[i_bdry];
                            for (unsigned int d = 0; d < space_dim; d++) {
			      sol_ctrl_x_bdry_gss[d] += sol_eldofs_Mat[pos_mat_ctrl][i_vol] * phi_ctrl_x_bdry[i_bdry * space_dim + d];
			    }
		      }
		      
		      
		  sol_adj_bdry_gss = 0.;
		      for (unsigned int i_bdry = 0; i_bdry <  phi_adj_bdry.size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_adj_bdry_gss  +=  sol_eldofs_Mat[pos_mat_adj][i_vol] * phi_adj_bdry[i_bdry];
              }		      
		      
//=============== grad dot n for residual ========================================= 
//     compute gauss quantities on the boundary through VOLUME interpolation
           std::fill(sol_adj_x_vol_at_bdry_gss.begin(), sol_adj_x_vol_at_bdry_gss.end(), 0.);
		      for (unsigned int iv = 0; iv < Sol_n_el_dofs_Mat_vol[pos_mat_adj]; iv++)  {
			
         for (unsigned int d = 0; d < space_dim; d++) {
			      sol_adj_x_vol_at_bdry_gss[d] += sol_eldofs_Mat[pos_mat_adj][iv] * phi_adj_x_vol_at_bdry[iv * space_dim + d];
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
                    - control_node_flag[first_loc_comp_ctrl][i_vol] * penalty_ctrl     * KEEP_ADJOINT_PUSH * (   sol_eldofs_Mat[pos_mat_state][i_vol] - sol_eldofs_Mat[pos_mat_ctrl][i_vol] )
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
//============ u = q ===========================	    
                 
if ( i_vol == j_vol )  {
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_state, i_vol, j_vol) ] += 
		     penalty_ctrl * KEEP_ADJOINT_PUSH * ( control_node_flag[first_loc_comp_ctrl][i_vol]);
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_ctrl, i_vol, j_vol) ]  += 
		     penalty_ctrl * KEEP_ADJOINT_PUSH * ( control_node_flag[first_loc_comp_ctrl][i_vol]) * (-1.);
		}
//============ u = q ===========================

		    

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

    elem_all[qrule_i][ielGeom][SolFEType_quantities[pos_sol_state]]->shape_funcs_current_elem(ig, JacI_iqp, phi_u, phi_u_x, phi_u_xx, space_dim);
    elem_all[qrule_i][ielGeom][SolFEType_quantities[pos_sol_adj]]  ->shape_funcs_current_elem(ig, JacI_iqp, phi_adj, phi_adj_x, phi_adj_xx, space_dim);

          
	sol_u_gss = 0.;
	sol_adj_gss = 0.;
	std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	
	for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[pos_mat_state]; i++) {
	                                                sol_u_gss      += sol_eldofs_Mat[pos_mat_state][i] * phi_u[i];
                   for (unsigned d = 0; d < space_dim; d++)   sol_u_x_gss[d] += sol_eldofs_Mat[pos_mat_state][i] * phi_u_x[i * space_dim + d];
          }
	
	for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[pos_mat_adj]; i++) {
	                                                sol_adj_gss      += sol_eldofs_Mat[pos_mat_adj][i] * phi_adj[i];
                   for (unsigned d = 0; d < space_dim; d++)   sol_adj_x_gss[d] += sol_eldofs_Mat[pos_mat_adj][i] * phi_adj_x[i * space_dim + d];
        }

//==========FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {
	  
              double laplace_rhs_du_adj_i = 0.;
              double laplace_rhs_dadj_u_i = 0.;
              for (unsigned kdim = 0; kdim < space_dim; kdim++) {
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_state] )  laplace_rhs_du_adj_i +=  phi_u_x   [i * space_dim + kdim] * sol_adj_x_gss[kdim];
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_adj] )    laplace_rhs_dadj_u_i +=  phi_adj_x [i * space_dim + kdim] * sol_u_x_gss[kdim];
	      }
	      
//============ Volume residuals - BEGIN  ==================	    
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_state, i) ] += - weight_iqp * ( target_flag * phi_u[i] * ( sol_u_gss - u_des)  - laplace_rhs_du_adj_i ); 
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_ctrl, i) ]  += - penalty_outside_control_boundary * ( (1 - control_node_flag[first_loc_comp_ctrl][i]) * (  sol_eldofs_Mat[pos_mat_ctrl][i] - 0.)  );
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_adj, i) ]   += - weight_iqp * (-1.) * (laplace_rhs_dadj_u_i);
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_mu, i) ]    += - penalty_outside_control_boundary * ( (1 - control_node_flag[first_loc_comp_ctrl][i]) * (  sol_eldofs_Mat[pos_mat_mu][i] - 0.)  );  //MU
//============  Volume Residuals - END ==================	    
	      
	      
//============ Volume Jacobians - BEGIN  ==================	    
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {
                
              double laplace_mat_dadj_u = 0.;
              double laplace_mat_du_adj = 0.;

              for (unsigned kdim = 0; kdim < space_dim; kdim++) {
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_adj]     && j < Sol_n_el_dofs_Mat_vol[pos_mat_state] )     laplace_mat_dadj_u        +=  (phi_adj_x [i * space_dim + kdim] * phi_u_x   [j * space_dim + kdim]);
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_state]   && j < Sol_n_el_dofs_Mat_vol[pos_mat_adj] )   laplace_mat_du_adj        +=  (phi_u_x   [i * space_dim + kdim] * phi_adj_x [j * space_dim + kdim]);
		
	      }

              //============ delta_state row ============================
              // BLOCK delta_state / state
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_state, i, j) ]  += weight_iqp  * target_flag *  phi_u[i] * phi_u[j];   
              //BLOCK delta_state / adjoint
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_adj, i, j) ]  += weight_iqp * (-1.) * laplace_mat_du_adj;
	      
	      
              //=========== delta_control row ===========================
              //enforce control zero outside the control boundary
	      if ( i==j )
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_ctrl, pos_mat_ctrl, i, j) ]  += penalty_outside_control_boundary * ( (1 - control_node_flag[first_loc_comp_ctrl][i]));    /*weight * phi_adj[i]*phi_adj[j]*/
              
	      //=========== delta_adjoint row ===========================
	      // BLOCK delta_adjoint / state
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_adj, pos_mat_state, i, j) ]  += weight_iqp * (-1.) * laplace_mat_dadj_u;

	      
	      //============= delta_mu row ===============================
	        if ( i==j )   
		  Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_mu, pos_mat_mu, i, j) ]  += penalty_outside_control_boundary * ( (1 - control_node_flag[first_loc_comp_ctrl][i]));    //MU
          
	         } // end phi_j loop
           } // endif assemble_matrix
//============ Volume Jacobians - END ==================	    

        } // end phi_i loop
        
      } // end gauss point loop

  

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    RES->add_vector_blocked(Res, L2G_dofmap_Mat_AllVars);

    if (assembleMatrix) {
      KK->add_matrix_blocked(Jac, L2G_dofmap_Mat_AllVars, L2G_dofmap_Mat_AllVars);
    }
    
    
    //========== dof-based part, without summation
 
     if (print_algebra_local) {
         assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat_vol, 10, 5);
         assemble_jacobian<double,double>::print_element_jacobian(iel, Jac, Sol_n_el_dofs_Mat_vol, 10, 5);
     }
     
     
  } //end element loop for each process
  

  //MU
add_one_times_mu_res_ctrl(iproc,
                               ineq_flag,
                               pos_mat_ctrl,
                               pos_mat_mu,
                               SolIndex_Mat,
                               sol,
                               mlPdeSys,
                               pdeSys,
                               RES);

    
  // ***************** END ASSEMBLY - ADD PART *******************

RES->close();
if (assembleMatrix) KK->close();  /// This is needed for the parallel, when splitting the add part from the insert part!!!


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
                        SolFEType_Mat, SolIndex_Mat, SolPdeIndex,
                        Sol_n_el_dofs_Mat_vol, sol_eldofs_Mat, L2G_dofmap_Mat);
// -------

	if ( volume_elem_contains_a_boundary_control_face( geom_element_iel.get_elem_center_3d() ) ) {


    	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {

       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);

                
       if(  face_is_a_boundary_control_face( el, iel, iface) ) {

       update_active_set_flag_for_current_nonlinear_iteration_bdry
   (msh, sol,
    iel, iface,
    geom_element_iel.get_coords_at_dofs_bdry_3d(), 
    sol_eldofs_Mat, 
    Sol_n_el_dofs_Mat_vol, 
    pos_mat_mu,               //this becomes a vector
    pos_mat_ctrl,             //this becomes a vector
    c_compl, 
    ctrl_lower, ctrl_upper,   //this becomes a vector
    sol_actflag,              //this becomes a vector
    solFEType_act_flag_sol, //remove this one, only Index
    solIndex_act_flag_sol);   //this becomes a vector
 

  node_insertion_bdry(iel, iface, 
                      msh,
                      L2G_dofmap_Mat,
                      pos_mat_mu, 
                      pos_mat_ctrl,
                      sol_eldofs_Mat,
                      Sol_n_el_dofs_Mat_vol,
                      sol_actflag, 
                      solFEType_act_flag_sol,  //remove this one, only Index
                      ineq_flag,
                      c_compl,
                      ctrl_lower, ctrl_upper,
                      KK, RES,
                      assembleMatrix
                      );
  
             }
             
       }
     }


     //============= delta_ctrl-delta_mu row ===============================
 if (assembleMatrix) { KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap_Mat[pos_mat_ctrl],  L2G_dofmap_Mat[pos_mat_mu], ineq_flag * 1.); }   //this becomes a vector

   }
   
   
  // ***************** INSERT PART - END *******************

  RES->close();
  if (assembleMatrix) KK->close();  ///@todo is it needed? I think so
    
    
  if (print_algebra_global) {
    assemble_jacobian< double, double >::print_global_jacobian(assembleMatrix, ml_prob, KK, mlPdeSys->GetNonlinearIt());
//     assemble_jacobian< double, double >::print_global_residual(ml_prob, RES,  mlPdeSys->GetNonlinearIt());

    RES->close();
    std::ostringstream res_out; res_out << ml_prob.GetFilesHandler()->GetOutputPath() << "./" << "res_" << mlPdeSys->GetNonlinearIt()  << ".txt";
    pdeSys->print_with_structure_matlab_friendly(iproc, res_out.str().c_str(), RES);

  }
     
     

  return;
}


 
  
void ComputeIntegral(const MultiLevelProblem& ml_prob)  {
  
  
  const NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("BoundaryControl");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->el;

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //=============== Geometry ========================================
   unsigned solType_coords = FE_DOMAIN;
 
  CurrentElem < double > geom_element_iel(dim, msh);
    
  constexpr unsigned int space_dim = 3;
  
  std::vector<double> normal(space_dim, 0.);
 //***************************************************

  //=============== Integration ========================================

 //***************************************************
  double alpha = ALPHA_CTRL_BDRY;
  double beta  = BETA_CTRL_BDRY;
  
 //*************** state ***************************** 
 //*************************************************** 
  vector <double> phi_u;     phi_u.reserve(max_size);
  vector <double> phi_u_x;   phi_u_x.reserve(max_size * space_dim);
  vector <double> phi_u_xx;  phi_u_xx.reserve(max_size * dim2);

 
  unsigned solIndex_u;
  solIndex_u = ml_sol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u = ml_sol->GetSolutionType(solIndex_u);    // get the finite element type for "state"

  vector < double >  sol_u; // local solution
  sol_u.reserve(max_size);
  
  double u_gss = 0.;
 //*************************************************** 
 //***************************************************

  
 //************** desired ****************************
 //***************************************************
  vector <double> phi_udes;
  vector <double> phi_udes_x;
  vector <double> phi_udes_xx;

    phi_udes.reserve(max_size);
    phi_udes_x.reserve(max_size * space_dim);
    phi_udes_xx.reserve(max_size * dim2);
 
  
//  unsigned solIndexTdes;
//   solIndexTdes = ml_sol->GetIndex("Tdes");    // get the position of "state" in the ml_sol object
//   unsigned solTypeTdes = ml_sol->GetSolutionType(solIndexTdes);    // get the finite element type for "state"

  vector < double >  sol_udes;
  sol_udes.reserve(max_size);

  double udes_gss = 0.;
 //***************************************************
 //***************************************************

 //************** cont *******************************
 //***************************************************
  vector <double> phi_ctrl_bdry;  
  vector <double> phi_ctrl_x_bdry; 

  phi_ctrl_bdry.reserve(max_size);
  phi_ctrl_x_bdry.reserve(max_size * space_dim);

  unsigned solIndex_ctrl = ml_sol->GetIndex("control");
  unsigned solType_ctrl = ml_sol->GetSolutionType(solIndex_ctrl);

   vector < double >  sol_ctrl;   sol_ctrl.reserve(max_size);
 //***************************************************
 //*************************************************** 
  
 //********** DATA *********************************** 
  double u_des = DesiredTarget();
 //*************************************************** 
  
  double integral_target = 0.;
  double integral_alpha  = 0.;
  double integral_beta   = 0.;

  

 //*************************************************** 
  constexpr unsigned qrule_i = QRULE_I;
  
     std::vector < std::vector < double > >  JacI_qp(space_dim);
     std::vector < std::vector < double > >  Jac_qp(dim);
    for (unsigned d = 0; d < Jac_qp.size(); d++) {   Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp.size(); d++) { JacI_qp[d].resize(dim); }
    
    double detJac_iqp;

     std::vector < std::vector < double > >  JacI_qp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_qp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_qp_bdry.size(); d++) {   Jac_qp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp_bdry.size(); d++) { JacI_qp_bdry[d].resize(dim-1); }
    
    double detJac_iqp_bdry;
  
  double weight_iqp; 
  double weight_iqp_bdry = 0.;
    
      //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all;
  ml_prob.get_all_abstract_fe_multiple(elem_all);
 //*************************************************** 
  
  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    const short unsigned ielGeom = geom_element_iel.geom_type();

  //************* set target domain flag **************
   geom_element_iel.set_elem_center_3d(iel, solType_coords);

   int target_flag = 0;
   target_flag = ElementTargetFlag(geom_element_iel.get_elem_center_3d());
 //***************************************************

   
 //*********** state ********************************* 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);
    sol_u    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
      unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);
    }
 //*********** state ********************************* 


 //*********** cont ********************************** 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);
    sol_ctrl    .resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);
    } 

 //*********** cont ********************************** 
 
 
 //*********** udes ********************************** 
    unsigned nDof_udes  = msh->GetElementDofNumber(iel, solType_u);
    sol_udes    .resize(nDof_udes);
    for (unsigned i = 0; i < sol_udes.size(); i++) {
            sol_udes[i] = u_des;  //dof value
    } 
 //*********** udes ********************************** 

 
 //********** ALL VARS ******************************* 
    int nDof_max    =  nDof_u;   //  TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
    if(nDof_udes > nDof_max) 
    {
      nDof_max = nDof_udes;
      }
    
    if(nDof_ctrl > nDof_max)
    {
      nDof_max = nDof_ctrl;
    }
    
 //***************************************************

  
	if ( volume_elem_contains_a_boundary_control_face( geom_element_iel.get_elem_center_3d() ) ) {
	  
	       
	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, iface);    
       const unsigned nve_bdry_ctrl = msh->GetElementFaceDofNumber(iel,iface,solType_ctrl);
       
// ----------
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();
// ----------

		
	    if( face_is_a_boundary_control_face(msh->el, iel, iface) ) {

	
		//============ initialize gauss quantities on the boundary ==========================================
                double sol_ctrl_bdry_gss = 0.;
                std::vector<double> sol_ctrl_x_bdry_gss(space_dim);
		//============ initialize gauss quantities on the boundary ==========================================
		
		for(unsigned ig_bdry = 0; ig_bdry < ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussPointsNumber(); ig_bdry++) {
		  
    elem_all[qrule_i][ielGeom_bdry][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), ig_bdry, Jac_qp_bdry, JacI_qp_bdry, detJac_iqp_bdry, space_dim);
    weight_iqp_bdry = detJac_iqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussWeightsPointer()[ig_bdry];
    elem_all[qrule_i][ielGeom_bdry][solType_ctrl] ->shape_funcs_current_elem(ig_bdry, JacI_qp_bdry, phi_ctrl_bdry, phi_ctrl_x_bdry, boost::none, space_dim);

		  
		 //========== compute gauss quantities on the boundary ===============================================
		  sol_ctrl_bdry_gss = 0.;
                  std::fill(sol_ctrl_x_bdry_gss.begin(), sol_ctrl_x_bdry_gss.end(), 0.);
		      for (int i_bdry = 0; i_bdry < nve_bdry_ctrl; i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_ctrl_bdry_gss +=  sol_ctrl[i_vol] * phi_ctrl_bdry[i_bdry];
                            for (int d = 0; d < space_dim; d++) {
			      sol_ctrl_x_bdry_gss[d] += sol_ctrl[i_vol] * phi_ctrl_x_bdry[i_bdry * space_dim + d];
			    }
		      }
		      
		      double laplace_ctrl_surface = 0.;  for (int d = 0; d < space_dim; d++) { laplace_ctrl_surface += sol_ctrl_x_bdry_gss[d] * sol_ctrl_x_bdry_gss[d]; }

                 //========= compute gauss quantities on the boundary ================================================
                  integral_alpha +=  weight_iqp_bdry * sol_ctrl_bdry_gss * sol_ctrl_bdry_gss; 
                  integral_beta  +=  weight_iqp_bdry * laplace_ctrl_surface;
                 
             }
	      } //end face
	      
	  }  // loop over element faces   
	  
	} //end if control element flag

//=====================================================================================================================  
//=====================================================================================================================  
//=====================================================================================================================  
  
  
   
      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom).GetGaussPointsNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), ig, Jac_qp, JacI_qp, detJac_iqp, space_dim);
    weight_iqp = detJac_iqp * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom).GetGaussWeightsPointer()[ig];

    elem_all[qrule_i][ielGeom][solType_u]                 ->shape_funcs_current_elem(ig, JacI_qp, phi_u, phi_u_x, phi_u_xx, space_dim);
    elem_all[qrule_i][ielGeom][solType_u/*solTypeTdes*/]  ->shape_funcs_current_elem(ig, JacI_qp, phi_udes, phi_udes_x, phi_udes_xx, space_dim);
    
	u_gss     = 0.;  for (unsigned i = 0; i < nDof_u; i++)        u_gss += sol_u[i]     * phi_u[i];
	udes_gss  = 0.;  for (unsigned i = 0; i < nDof_udes; i++)  udes_gss += sol_udes[i]  * phi_udes[i];  

               integral_target += target_flag * weight_iqp * (u_gss  - udes_gss) * (u_gss - udes_gss);
	  
      } // end gauss point loop
      
  } //end element loop

  double total_integral = 0.5 * integral_target + 0.5 * alpha * integral_alpha + 0.5 * beta * integral_beta;
  
  
  ////////////////////////////////////////
       std::cout << "integral on processor " << iproc << ": " << total_integral << std::endl;

   double J = 0.;
      MPI_Allreduce( &total_integral, &J, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );


    std::cout << "@@@@@@@@@@@@@@@@ functional value: " << J << std::endl;
  
//   std::cout << "The value of the integral_target is " << std::setw(11) << std::setprecision(10) << integral_target << std::endl;
//   std::cout << "The value of the integral_alpha  is " << std::setw(11) << std::setprecision(10) << integral_alpha << std::endl;
//   std::cout << "The value of the integral_beta   is " << std::setw(11) << std::setprecision(10) << integral_beta << std::endl;
//   std::cout << "The value of the total integral  is " << std::setw(11) << std::setprecision(10) << total_integral << std::endl;
 
return;
  
}
  
