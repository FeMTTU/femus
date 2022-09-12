// Solving Navier-Stokes problem using automatic differentiation and/or Picards method
// boundary conditions were set in 2D as, no slip in left,right of the box and top to bottom gravity is enforced
// therefore, U=V=0 on left and right, U=0 on top and bottom, V is free 

#include "adept.h"

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "Files.hpp"
#include "Parameter.hpp"
#include "Fluid.hpp"
#include "MultiLevelSolution.hpp"
#include "NumericVector.hpp"
#include "paral.hpp"//to get iproc HAVE_MPI is inside here

#include "Assemble_jacobian.hpp"


#include   "../manufactured_solutions.hpp"

#define FACE_FOR_CONTROL  2
#define FACE_FOR_TARGET    2


#include   "../nsopt_params.hpp"

  const double cost_functional_coeff = 1.;
  const double alpha = ALPHA_CTRL_VOL;
  const double beta  = BETA_CTRL_VOL;
  const double penalty_outside_control_domain = PENALTY_OUTSIDE_CONTROL_DOMAIN;         // penalty for zero control outside

  
//****** Mesh ********************************
  #define no_of_ref    N_UNIFORM_LEVELS  //mesh refinements

  
#define exact_sol_flag 0 // 1 = if we want to use manufactured solution; 0 = if we use regular convention
#define compute_conv_flag 0 // 1 = if we want to compute the convergence and error ; 0 =  no error computation

#define NO_OF_L2_NORMS 11   //U,V,P,adj_0,adj_1,PADJ,ctrl_0,ctrl_1,PCTRL,U+U0,V+V0
#define NO_OF_H1_NORMS 8    //U,V,adj_0,adj_1,ctrl_0,ctrl_1,U+U0,V+V0


using namespace femus;

 
bool Solution_set_boundary_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char SolName[], double& value, const int faceName, const double time) {
  //1: bottom  //2: right  //3: top  //4: left  (2D square)
  //1: bottom  //2: top    //3: side            (3D cylinder)
    
  
  bool dirichlet = true;
   value = 0.;

#if exact_sol_flag == 0
// b.c. for lid-driven cavity problem, wall u_top = 1 = shear_force, v_top = 0 and u=v=0 on other 3 walls ; rhs_f = body_force = {0,0}

   if (faceName == FACE_FOR_CONTROL)  {
        if (x[ axis_direction_Gamma_control(faceName) ] > GAMMA_CONTROL_LOWER - 1.e-5 && x[ axis_direction_Gamma_control(faceName) ] < GAMMA_CONTROL_UPPER + 1.e-5)  { 
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


double Solution_set_initial_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;

     if(!strcmp(name,"TargReg")) {
        value = ElementTargetFlag(x);
    }
    else if(!strcmp(name,"ContReg")) {
        value = ControlDomainFlag_internal_restriction(x);
    }

    return value;
}



 void  name_of_unknowns(std::vector< Unknown > & unknowns, const unsigned int dimension) {

  const int state_pos_begin = 0;
  const int adj_pos_begin   =      dimension + 1;
  const int ctrl_pos_begin  = 2 * (dimension + 1);
  const int mu_pos_begin    = 3 * (dimension + 1);
  
                        unknowns[state_pos_begin + 0]._name    = "u_0";
                        unknowns[state_pos_begin + 1]._name    = "u_1";
  if (dimension == 3)   unknowns[state_pos_begin + 2]._name    = "u_2";
                                unknowns[dimension]._name      = "p_u";
  
                        unknowns[adj_pos_begin + 0]._name      = "adj_0";
                        unknowns[adj_pos_begin + 1]._name      = "adj_1";
  if (dimension == 3)   unknowns[adj_pos_begin + 2]._name      = "adj_2";
                unknowns[adj_pos_begin + dimension]._name      = "p_adj";
  
                       unknowns[ctrl_pos_begin + 0]._name      = "ctrl_0";
                       unknowns[ctrl_pos_begin + 1]._name      = "ctrl_1";
  if (dimension == 3)  unknowns[ctrl_pos_begin + 2]._name      = "ctrl_2";
               unknowns[ctrl_pos_begin + dimension]._name      = "p_ctrl";

                       unknowns[mu_pos_begin + 0]._name      = "mu_0";
                       unknowns[mu_pos_begin + 1]._name      = "mu_1";
  if (dimension == 3)  unknowns[mu_pos_begin + 2]._name      = "mu_2";
 
     
}


 //Unknown definition  ==================
 const std::vector< Unknown >  provide_list_of_unknowns(const unsigned int dimension) {
     
     
  std::vector< FEFamily > feFamily;
  std::vector< FEOrder >   feOrder;

                        feFamily.push_back(LAGRANGE);   //state
                        feFamily.push_back(LAGRANGE);
  if (dimension == 3)   feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE/*DISCONTINOUS_POLYNOMIAL*/);
  
                        feFamily.push_back(LAGRANGE);   //adjoint
                        feFamily.push_back(LAGRANGE);
  if (dimension == 3)   feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE/*DISCONTINOUS_POLYNOMIAL*/);
  
                        feFamily.push_back(LAGRANGE);   //control
                        feFamily.push_back(LAGRANGE);
  if (dimension == 3)   feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
 
                        feFamily.push_back(LAGRANGE);   //mu
                        feFamily.push_back(LAGRANGE);
  if (dimension == 3)   feFamily.push_back(LAGRANGE);
 
  
  
                        feOrder.push_back(SECOND);
                        feOrder.push_back(SECOND);
  if (dimension == 3)   feOrder.push_back(SECOND);
                        feOrder.push_back(FIRST);
  
                        feOrder.push_back(SECOND);
                        feOrder.push_back(SECOND);
  if (dimension == 3)   feOrder.push_back(SECOND);
                        feOrder.push_back(FIRST);
  
                        feOrder.push_back(SECOND);
                        feOrder.push_back(SECOND);
  if (dimension == 3)   feOrder.push_back(SECOND);
                        feOrder.push_back(FIRST);
 
                        feOrder.push_back(SECOND);   //mu
                        feOrder.push_back(SECOND);
  if (dimension == 3)   feOrder.push_back(SECOND);

  assert( feFamily.size() == feOrder.size() );
 
 std::vector< Unknown >  unknowns(feFamily.size());
 
    name_of_unknowns(unknowns, dimension);
  
     for (unsigned int u = 0; u < unknowns.size(); u++) {
         
              unknowns[u]._fe_family  = feFamily[u];
              unknowns[u]._fe_order   = feOrder[u];
              unknowns[u]._time_order = 0;
              unknowns[u]._is_pde_unknown = true;
              unknowns[u]._is_sparse = true;
              
     }
     
 
   return unknowns;
     
}





void AssembleNavierStokesOpt_nonAD(MultiLevelProblem &ml_prob);
void AssembleNavierStokesOpt_AD   (MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );

void ComputeIntegral(const MultiLevelProblem& ml_prob);

double*  GetErrorNorm(const MultiLevelProblem& ml_prob, MultiLevelSolution* ml_sol, Solution* sol_coarser_prolongated);
// ||u_h - u_(h/2)||/||u_(h/2)-u_(h/4)|| = 2^alpha, alpha is order of conv 
//i.e. ||prol_(u_(i-1)) - u_(i)|| = err(i) => err(i-1)/err(i) = 2^alpha ,implemented as log(err(i)/err(i+1))/log2



int main(int argc, char** args) {



  // ======= Init ========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

   // ======= Problem  ==================
  MultiLevelProblem ml_prob;

  
  // ======= Files - BEGIN  ========================
  const bool use_output_time_folder = false;
  const bool redirect_cout_to_file = false;
  Files files; 
        files.CheckIODirectories(use_output_time_folder);
	    files.RedirectCout(redirect_cout_to_file);
 
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
  ml_prob.set_all_abstract_fe_multiple();
  // ======= Problem, Quad Rule - END  ========================
  
  
  // ======= Mesh, Coarse reading - BEGIN ==================
  MultiLevelMesh ml_mesh;
 
    std::string mesh_folder_file = "input/";
//   std::string input_file = "parametric_square_1x1.med";
//   std::string input_file = "parametric_square_1x2.med";
//   std::string input_file = "cyl.med"; // "fifth"
  std::string input_file = "Mesh_3_groups_with_bdry_nodes_coarser.med";
  std::ostringstream mystream; mystream << "./" << mesh_folder_file << input_file;
  const std::string infile = mystream.str();

  
  const bool read_groups = true;
  const bool read_boundary_groups = true;
    
  ml_mesh.ReadCoarseMeshFileReadingBeforePartitioning(infile.c_str(), Lref, read_groups, read_boundary_groups);
    
  ml_mesh.GetLevelZero(0)->build_dofmap_all_fe_families_and_elem_and_node_structures();
 

  ml_mesh.BuildFETypesBasedOnExistingCoarseMeshGeomElements();
  
  ml_mesh.PrepareNewLevelsForRefinement();


  ml_mesh.InitializeQuadratureWithFEEvalsOnExistingCoarseMeshGeomElements(fe_quad_rule_vec[0].c_str()); ///@todo keep it only for compatibility with old ElemType, because of its destructor 
  // ======= Mesh, Coarse reading - END ==================

  
  // ======= Convergence Rate, Preparation - BEGIN  ==================

  unsigned dim = ml_mesh.GetDimension();
  unsigned maxNumberOfMeshes;

  if (dim == 2) {
    maxNumberOfMeshes = no_of_ref;
  } else {
    maxNumberOfMeshes = no_of_ref /*4*/;
  }

 
 
   // ======= Solutions that are Unknowns - BEGIN  ==================
  std::vector< Unknown > unknowns = provide_list_of_unknowns( dim );
   // ======= Solutions that are Unknowns - END  ==================

  

#if compute_conv_flag == 1
     double comp_conv[maxNumberOfMeshes][NO_OF_L2_NORMS+NO_OF_H1_NORMS];
 
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
      ml_sol_all_levels->Initialize(unknowns[u]._name.c_str(), Solution_set_initial_conditions, & ml_prob);
  }
  
      ml_sol_all_levels->AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions);
   for (unsigned int u = 0; u < unknowns.size(); u++)  { 
     ml_sol_all_levels->GenerateBdc(unknowns[u]._name.c_str(), (unknowns[u]._time_order == 0) ? "Steady" : "Time_dependent", & ml_prob);
  }
   // ======= Solutions that are Unknowns - END  ==================


  // ======= Solutions that are not Unknowns - BEGIN  ==================
     ml_sol_all_levels->AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO);
     ml_sol_all_levels->Initialize("TargReg",     Solution_set_initial_conditions, & ml_prob);
     ml_sol_all_levels->AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO);
     ml_sol_all_levels->Initialize("ContReg",     Solution_set_initial_conditions, & ml_prob);
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
    ml_sol.Initialize(act_set_flag_name[d].c_str(), Solution_set_initial_conditions, & ml_prob);
   }
  // ******** active flag - END 
   
// ======= Solutions that are not Unknowns - END  ==================
  

  // ======= Problem, System - BEGIN ========================
  NonLinearImplicitSystemWithPrimalDualActiveSetMethod & system_opt    = ml_prob.add_system < NonLinearImplicitSystemWithPrimalDualActiveSetMethod > ("NSOpt");  ///@todo this MUST return a REFERENCE, otherwise it doesn't run!
  
  system_opt.SetActiveSetFlagName(act_set_flag_name); //MU


  for (unsigned int u = 0; u < unknowns.size(); u++)  { 
  system_opt.AddSolutionToSystemPDE(unknowns[u]._name.c_str());
  }
  
//   system_opt.SetAssembleFunction(AssembleNavierStokesOpt_AD);  //AD doesn't seem to work now
  system_opt.SetAssembleFunction(AssembleNavierStokesOpt_nonAD);
    
// *****************
  system_opt.SetDebugNonlinear(true);
  system_opt.SetDebugFunction(ComputeIntegral);
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
    
      for(int j = 0; j < NO_OF_L2_NORMS+NO_OF_H1_NORMS; j++)       comp_conv[i-1][j] = norm[j];
  
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
       FE_convergence<double>::output_convergence_rate(comp_conv[i][NO_OF_L2_NORMS + j], comp_conv[i + 1][NO_OF_L2_NORMS + j], norm_names_H1[j], maxNumberOfMeshes , i );
    }
  }
  std::cout << std::endl;
  std::cout << "=======================================================================" << std::endl;
#endif
 
  return 0;
}







void AssembleNavierStokesOpt_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled
std::cout << " ********************************  AD SYSTEM ******************************************** " << std::endl;
  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem * mlPdeSys   = & ml_prob.get_system<NonLinearImplicitSystem> ("NSOpt");   // pointer to the nonlinear implicit system named "NSOpt" 
   const unsigned level = mlPdeSys->GetLevelToAssemble();
 
  Mesh*          msh          	= ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         	= msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  ml_sol    = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        	= ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


   LinearEquationSolver*  pdeSys	 = mlPdeSys->_LinSolver[level];   
 SparseMatrix*    JAC         	= pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

  // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //geometry *******************************
  unsigned coordXType = 2; /*BIQUADR_FE*/// get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  unsigned solType_coords = coordXType;  //FE_DOMAIN = 0; //we do linear FE this time // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
  const unsigned int dim_offset_grad = dim  /*3*/  /*2*/    ;
 
  vector < vector < double > > coordX(dim);    // local coordinates

  for (unsigned  k = 0; k < dim; k++) { 
    coordX[k].reserve(max_size);
  }
  //geometry *******************************

//STATE######################################################################
  //velocity *******************************
  vector < unsigned > solVIndex(dim);
  solVIndex[0] = ml_sol->GetIndex("u_0");
  solVIndex[1] = ml_sol->GetIndex("u_1");

  if (dim == 3) solVIndex[2] = ml_sol->GetIndex("u_2");

  unsigned solVType = ml_sol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"
  vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("u_0");
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("u_1");

  if (dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("u_2");
  
  vector < vector < adept::adouble > >  solV(dim);    // local solution
   vector< vector < adept::adouble > > aResV(dim);    // local redidual vector
   
 for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(max_size);
    aResV[k].reserve(max_size);
  }

  
  vector <double> phiV_gss;  // local test function
  vector <double> phiV_x_gss; // local test function first order partial derivatives
  vector <double> phiV_xx_gss; // local test function second order partial derivatives

  phiV_gss.reserve(max_size);
  phiV_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);
  phiV_xx_gss.reserve(max_size * dim2);
  

  //velocity *******************************

  //pressure *******************************
  unsigned solPIndex;
  solPIndex = ml_sol->GetIndex("p_u");    // get the position of "p_u" in the ml_sol object
  unsigned solPType = ml_sol->GetSolutionType(solPIndex);    // get the finite element type for "p_u"

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("p_u");    // get the position of "p_u" in the pdeSys object

  vector < adept::adouble >  solP; // local solution
  vector< adept::adouble > aResP; // local redidual vector
  
  solP.reserve(max_size);
  aResP.reserve(max_size);
  
  double* phiP_gss;
  //pressure *******************************
//STATE######################################################################
  
//ADJOINT######################################################################
  //velocity *******************************
  vector < unsigned > solVadjIndex(dim);
  solVadjIndex[0] = ml_sol->GetIndex("adj_0");
  solVadjIndex[1] = ml_sol->GetIndex("adj_1");

  if (dim == 3) solVadjIndex[2] = ml_sol->GetIndex("adj_2");

  unsigned solVadjType = ml_sol->GetSolutionType(solVadjIndex[0]);
 vector < unsigned > solVPdeadjIndex(dim);
  solVPdeadjIndex[0] = mlPdeSys->GetSolPdeIndex("adj_0");
  solVPdeadjIndex[1] = mlPdeSys->GetSolPdeIndex("adj_1");

  if (dim == 3) solVPdeadjIndex[2] = mlPdeSys->GetSolPdeIndex("adj_2");
  
  vector < vector < adept::adouble > >  solVadj(dim);
   vector< vector < adept::adouble > > aResVadj(dim);
   
 for (unsigned  k = 0; k < dim; k++) {
    solVadj[k].reserve(max_size);
    aResVadj[k].reserve(max_size);
  }

  
  vector <double> phiVadj_gss;  // local test function
  vector <double> phiVadj_x_gss; // local test function first order partial derivatives
  vector <double> phiVadj_xx_gss; // local test function second order partial derivatives

  phiVadj_gss.reserve(max_size);
  phiVadj_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);
  phiVadj_xx_gss.reserve(max_size * dim2);
  
  //velocity *******************************

  //pressure *******************************
  unsigned solPadjIndex;
  solPadjIndex = ml_sol->GetIndex("p_adj");
  unsigned solPadjType = ml_sol->GetSolutionType(solPadjIndex);    // get the finite element type for "p_adj"

  unsigned solPPdeadjIndex;
  solPPdeadjIndex = mlPdeSys->GetSolPdeIndex("p_adj");    // get the position of "p_adj" in the pdeSys object

  vector < adept::adouble >  solPadj; // local solution
  vector< adept::adouble > aResPadj; // local redidual vector
  
  solPadj.reserve(max_size);
  aResPadj.reserve(max_size);
  
  double* phiPadj_gss;
  //pressure *******************************
//ADJOINT######################################################################

  
//CONTROL######################################################################
  //velocity *******************************
  vector < unsigned > solVctrlIndex(dim);
  solVctrlIndex[0] = ml_sol->GetIndex("ctrl_0");
  solVctrlIndex[1] = ml_sol->GetIndex("ctrl_1");

  if (dim == 3) solVctrlIndex[2] = ml_sol->GetIndex("ctrl_2");

  unsigned solVctrlType = ml_sol->GetSolutionType(solVctrlIndex[0]);
 vector < unsigned > solVPdectrlIndex(dim);
 
  solVPdectrlIndex[0] = mlPdeSys->GetSolPdeIndex("ctrl_0");
  solVPdectrlIndex[1] = mlPdeSys->GetSolPdeIndex("ctrl_1");
  if (dim == 3) solVPdectrlIndex[2] = mlPdeSys->GetSolPdeIndex("ctrl_2");
  
  vector < vector < adept::adouble > >  solVctrl(dim);    // local solution
   vector< vector < adept::adouble > > aResVctrl(dim);    // local redidual vector
   
 for (unsigned  k = 0; k < dim; k++) {
    solVctrl[k].reserve(max_size);
    aResVctrl[k].reserve(max_size);
  }

  
  vector <double> phiVctrl_gss;  // local test function
  vector <double> phiVctrl_x_gss; // local test function first order partial derivatives
  vector <double> phiVctrl_xx_gss; // local test function second order partial derivatives

  phiVctrl_gss.reserve(max_size);
  phiVctrl_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);
  phiVctrl_xx_gss.reserve(max_size * dim2);
  

  //velocity *******************************

  //pressure *******************************
  unsigned solPctrlIndex;
  solPctrlIndex = ml_sol->GetIndex("p_ctrl");    // get the position of "p_ctrl" in the ml_sol object
  unsigned solPctrlType = ml_sol->GetSolutionType(solPctrlIndex);    // get the finite element type for "p_ctrl"

  unsigned solPPdectrlIndex;
  solPPdectrlIndex = mlPdeSys->GetSolPdeIndex("p_ctrl");    // get the position of "p_ctrl" in the pdeSys object

  vector < adept::adouble >  solPctrl; // local solution
  vector< adept::adouble > aResPctrl; // local redidual vector
  
  solPctrl.reserve(max_size);
  aResPctrl.reserve(max_size);
  
  double* phiPctrl_gss;
  //pressure *******************************
//CONTROL######################################################################


  
  //Nondimensional values ******************
  double IRe 		= ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
  //Nondimensional values ******************
  
  
  
  
  vector< int > L2G_dofmap_Mat; // local to global pdeSys dofs
  L2G_dofmap_Mat.reserve(3 *(dim + 1) *max_size);

  vector< double > Res; // local redidual vector
  Res.reserve(3 *(dim + 1) *max_size);

  vector < double > Jac;
  Jac.reserve(3* (dim + 1) *max_size * 3*(dim + 1) *max_size);

  JAC->zero(); // Set to zero all the entries of the Global Matrix

//*************************************************** 
     std::vector < std::vector < double > >  JacI_iqp(space_dim);
     std::vector < std::vector < double > >  Jac_iqp(dim);
    for (unsigned d = 0; d < Jac_iqp.size(); d++) {   Jac_iqp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp.size(); d++) { JacI_iqp[d].resize(dim); }
    
   double AbsDetJxWeight_iqp;
   double detJac_iqp;

    //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
//*************************************************** 
 
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {


   // geometry *****************************
      geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
      const short unsigned ielGeom = geom_element_iel.geom_type();
  // geometry end *****************************

  // equation
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
    unsigned nDofsVP = dim * nDofsV + nDofsP;
   
    
    unsigned nDofsVadj = msh->GetElementDofNumber(iel, solVadjType);    // number of solution element dofs
    unsigned nDofsPadj = msh->GetElementDofNumber(iel, solPadjType);    // number of solution element dofs

     unsigned nDofsVctrl = msh->GetElementDofNumber(iel, solVctrlType);    // number of solution element dofs
    unsigned nDofsPctrl = msh->GetElementDofNumber(iel, solPctrlType);    // number of solution element dofs
    
    
     unsigned nDofsVP_tot = 3*nDofsVP;
     
     
    
    for (unsigned  k = 0; k < dim; k++)  {
      solV[k].resize(nDofsV);
      solVadj[k].resize(nDofsVadj);
      solVctrl[k].resize(nDofsVctrl);
    }
    solP.resize(nDofsP);
    solPadj.resize(nDofsPadj);
    solPctrl.resize(nDofsPctrl);
    

//element matrices and vectors
    // resize local arrays
    L2G_dofmap_Mat.resize(nDofsVP_tot);
    
//     Jac.resize(nDofsVP * nDofsVP);

    for (unsigned  k = 0; k < dim; k++) {
      aResV[k].resize(nDofsV);    //resize
      std::fill(aResV[k].begin(), aResV[k].end(), 0);    //set aRes to zero
    
      aResVadj[k].resize(nDofsVadj);    //resize
      std::fill(aResVadj[k].begin(), aResVadj[k].end(), 0);    //set aRes to zero
      
      aResVctrl[k].resize(nDofsVctrl);    //resize
      std::fill(aResVctrl[k].begin(), aResVctrl[k].end(), 0);    //set aRes to zero

    }

    aResP.resize(nDofsP);    //resize
    std::fill(aResP.begin(), aResP.end(), 0);    //set aRes to zero

     aResPadj.resize(nDofsPadj);    //resize
    std::fill(aResPadj.begin(), aResPadj.end(), 0);    //set aRes to zero
    
    aResPctrl.resize(nDofsPctrl);    //resize
    std::fill(aResPctrl.begin(), aResPctrl.end(), 0);    //set aRes to zero

  //*************************************** 
  
  //***** set target domain flag ********************************** 
  geom_element_iel.set_elem_center_3d(iel, solType_coords);

   int target_flag = 0;
   target_flag = ElementTargetFlag(geom_element_iel.get_elem_center_3d()/*elem_center*/);
//***************************************   
    
    
   //STATE###################################################################  
    // velocity ************
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof // via local to global solution node

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        L2G_dofmap_Mat[i + k * nDofsV] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);
      }
    }
    
    // pressure *************
    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // global to global mapping between solution node and solution dof // via local to global solution node
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      L2G_dofmap_Mat[i + dim * nDofsV] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }
//STATE###################################################################

//ADJ###################################################################
     // velocity ************
    for (unsigned i = 0; i < nDofsVadj; i++) {
      unsigned solVadjDof = msh->GetSolutionDof(i, iel, solVadjType);    // global to global mapping between solution node and solution dof // via local to global solution node

      for (unsigned  k = 0; k < dim; k++) {
        solVadj[k][i] = (*sol->_Sol[solVadjIndex[k]])(solVadjDof);      // global extraction and local storage for the solution
        L2G_dofmap_Mat[i + k * nDofsV +nDofsVP] = pdeSys->GetSystemDof(solVadjIndex[k], solVPdeadjIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }
    
    // pressure *************
    for (unsigned i = 0; i < nDofsPadj; i++) {
      unsigned solPadjDof = msh->GetSolutionDof(i, iel, solPadjType);    // global to global mapping between solution node and solution dof // via local to global solution node
      solPadj[i] = (*sol->_Sol[solPadjIndex])(solPadjDof);      // global extraction and local storage for the solution
      L2G_dofmap_Mat[i + dim * nDofsV +nDofsVP] = pdeSys->GetSystemDof(solPadjIndex, solPPdeadjIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }
//ADJ###################################################################


//CTRL###################################################################
     // velocity ************
    for (unsigned i = 0; i < nDofsVctrl; i++) {
      unsigned solVctrlDof = msh->GetSolutionDof(i, iel, solVctrlType);    // global to global mapping between solution node and solution dof // via local to global solution node

      for (unsigned  k = 0; k < dim; k++) {
        solVctrl[k][i] = (*sol->_Sol[solVctrlIndex[k]])(solVctrlDof);      // global extraction and local storage for the solution
        L2G_dofmap_Mat[i + k * nDofsV + 2*nDofsVP] = pdeSys->GetSystemDof(solVctrlIndex[k], solVPdectrlIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }
    
    // pressure *************
    for (unsigned i = 0; i < nDofsPctrl; i++) {
      unsigned solPctrlDof = msh->GetSolutionDof(i, iel, solPctrlType);    // global to global mapping between solution node and solution dof // via local to global solution node
      solPctrl[i] = (*sol->_Sol[solPctrlIndex])(solPctrlDof);      // global extraction and local storage for the solution
      L2G_dofmap_Mat[i + dim * nDofsV + 2*nDofsVP] = pdeSys->GetSystemDof(solPctrlIndex, solPPdectrlIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }
//CTRL###################################################################




      // start a new recording of all the operations involving adept::adouble variables
      s.new_recording();

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {

//STATE#############################################################################3	
        // *** get gauss point weight, test function and test function partial derivatives ***
   elem_all[ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), ig, Jac_iqp, JacI_iqp, detJac_iqp, space_dim);
   
    AbsDetJxWeight_iqp = detJac_iqp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[ig];
   
    elem_all[ielGeom][solVType]->shape_funcs_current_elem(ig, JacI_iqp, phiV_gss, phiV_x_gss, phiV_xx_gss , space_dim);

        phiP_gss = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

	
        vector < adept::adouble > solV_gss(dim, 0);
        vector < vector < adept::adouble > > gradSolV_gss(dim);
        vector < double > coordX_gss(dim, 0.);

        for (unsigned  k = 0; k < dim; k++) {
          gradSolV_gss[k].resize(dim_offset_grad);
          std::fill(gradSolV_gss[k].begin(), gradSolV_gss[k].end(), 0);
        }

          for (unsigned  k = 0; k < dim; k++) {
        for (unsigned i = 0; i < geom_element_iel.get_coords_at_dofs_3d()[k].size(); i++) {
            coordX_gss[k] += geom_element_iel.get_coords_at_dofs_3d()[k][i] * phiV_gss[i];    ///@todo change phiV into phi_coords!!!
             }
          }
          
          for (unsigned  k = 0; k < dim; k++) {
        for (unsigned i = 0; i < solV[k].size(); i++) {
            
            solV_gss[k]   += solV[k][i] * phiV_gss[i];
        

          for (unsigned j = 0; j < dim_offset_grad; j++) {
              gradSolV_gss[k][j] += solV[k][i] * phiV_x_gss[i * dim_offset_grad + j] ;
            }
            
          }
        }

        adept::adouble solP_gss = 0;

        for (unsigned i = 0; i < nDofsP; i++) {
          solP_gss += phiP_gss[i] * solP[i];
        }


//STATE###############################################################################

//ADJOINT#############################################################################3	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[ielGeom][solVadjType]->shape_funcs_current_elem(ig, JacI_iqp, phiVadj_gss, phiVadj_x_gss, phiVadj_xx_gss , space_dim);

        phiPadj_gss = msh->_finiteElement[ielGeom][solPadjType]->GetPhi(ig);

	
        vector < adept::adouble > solVadj_gss(dim, 0);
        vector < vector < adept::adouble > > gradSolVadj_gss(dim);

        for (unsigned  k = 0; k < dim; k++) {
          gradSolVadj_gss[k].resize(dim_offset_grad /*space_dim*/);
          std::fill(gradSolVadj_gss[k].begin(), gradSolVadj_gss[k].end(), 0);
        }

        for (unsigned i = 0; i < nDofsVadj; i++) {
          for (unsigned  k = 0; k < dim; k++) {
            solVadj_gss[k] += phiVadj_gss[i] * solVadj[k][i];
          }

          for (unsigned j = 0; j < dim_offset_grad /*space_dim*/; j++) {
            for (unsigned  k = 0; k < dim; k++) {
              gradSolVadj_gss[k][j] += phiVadj_x_gss[i * dim_offset_grad /*space_dim*/ + j] * solVadj[k][i];
            }
          }
        }

        adept::adouble solPadj_gss = 0;

        for (unsigned i = 0; i < nDofsPadj; i++) {
          solPadj_gss += phiPadj_gss[i] * solPadj[i];
        }

//ADJOINT###############################################################################



//CONTROL#############################################################################3	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[ielGeom][solVctrlType]->shape_funcs_current_elem(ig, JacI_iqp, phiVctrl_gss, phiVctrl_x_gss, phiVctrl_xx_gss , space_dim);

        phiPctrl_gss = msh->_finiteElement[ielGeom][solPctrlType]->GetPhi(ig);

	
        vector < adept::adouble > solVctrl_gss(dim, 0);
        vector < vector < adept::adouble > > gradSolVctrl_gss(dim);

        for (unsigned  k = 0; k < dim; k++) {
          gradSolVctrl_gss[k].resize(dim_offset_grad /*space_dim*/);
          std::fill(gradSolVctrl_gss[k].begin(), gradSolVctrl_gss[k].end(), 0);
        }

        for (unsigned i = 0; i < nDofsVctrl; i++) {
          for (unsigned  k = 0; k < dim; k++) {
            solVctrl_gss[k] += phiVctrl_gss[i] * solVctrl[k][i];
          }

          for (unsigned j = 0; j < dim_offset_grad /*space_dim*/; j++) {
            for (unsigned  k = 0; k < dim; k++) {
              gradSolVctrl_gss[k][j] += phiVctrl_x_gss[i * dim_offset_grad /*space_dim*/ + j] * solVctrl[k][i];
            }
          }
        }

        adept::adouble solPctrl_gss = 0;

        for (unsigned i = 0; i < nDofsPctrl; i++) {
          solPctrl_gss += phiPctrl_gss[i] * solPctrl[i];
        }
//CONTROL###############################################################################


#if exact_sol_flag == 1
//computation of RHS (force and desired velocity) using MMS=============================================== 
//state values--------------------
vector <double>  exact_stateVel(dim, 0.);
   mms_state_control::value_stateVel(coordX_gss, exact_stateVel);
vector < vector < double > > exact_grad_stateVel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_stateVel[k].resize(dim);
    std::fill(exact_grad_stateVel[k].begin(), exact_grad_stateVel[k].end(), 0.);
}
   mms_state_control::gradient_stateVel(coordX_gss,exact_grad_stateVel);
vector <double>  exact_lap_stateVel(dim, 0.);
   mms_state_control::laplace_stateVel(coordX_gss, exact_lap_stateVel);
vector <double> exact_grad_statePress(dim, 0.);
   mms_state_control::gradient_statePress(coordX_gss, exact_grad_statePress);

//control values-------------------------------
vector <double>  exact_ctrlVel(dim);
   mms_state_control::value_ctrlVel(coordX_gss, exact_ctrlVel);
vector < vector < double > > exact_grad_ctrlVel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_ctrlVel[k].resize(dim);
    std::fill(exact_grad_ctrlVel[k].begin(), exact_grad_ctrlVel[k].end(), 0.);
}
   mms_state_control::gradient_ctrlVel(coordX_gss,exact_grad_ctrlVel);
vector <double>  exact_lap_ctrlVel(dim);
   mms_state_control::laplace_ctrlVel(coordX_gss, exact_lap_ctrlVel);

//convection terms from delta_state-------------------------------------
vector <double>  exact_conv_u_nabla_u(dim,0.);
vector <double>  exact_conv_u_nabla_uctrl(dim,0.);
vector <double>  exact_conv_uctrl_nabla_u(dim,0.);
vector <double>  exact_conv_uctrl_nabla_uctrl(dim,0.);

for (unsigned k = 0; k < dim; k++){
    for (unsigned i = 0; i < dim; i++){
    exact_conv_u_nabla_u[k] += exact_grad_stateVel[k][i] * exact_stateVel[i] ; 
    exact_conv_u_nabla_uctrl[k] += exact_grad_ctrlVel[k][i] * exact_stateVel[i] ; 
    exact_conv_uctrl_nabla_u[k] += exact_grad_stateVel[k][i] * exact_ctrlVel[i] ; 
    exact_conv_uctrl_nabla_uctrl[k] += exact_grad_ctrlVel[k][i] * exact_ctrlVel[i] ; 
    }
}

//convection terms from delta_adjoint-------------------------
vector <double>  exact_conv_u_nabla_uadj(dim,0.);
vector <double>  exact_conv_nabla_uT_uadj(dim,0.);
vector <double>  exact_conv_nabla_uctrlT_uadj(dim,0.);
vector <double>  exact_conv_uctrl_nabla_uadj(dim,0.);

for (unsigned k = 0; k < dim; k++){
    for (unsigned i = 0; i < dim; i++){
    exact_conv_u_nabla_uadj[k] += exact_grad_stateVel[k][i] * exact_stateVel[i] ; 
    exact_conv_nabla_uT_uadj[k] += exact_grad_stateVel[i][k] * exact_stateVel[i];
    exact_conv_nabla_uctrlT_uadj[k] += exact_grad_ctrlVel[i][k] * exact_stateVel[i];  
    exact_conv_uctrl_nabla_uadj[k] += exact_grad_stateVel[k][i] * exact_ctrlVel[i] ; 
    }
}

//force and desired velocity ---------------------------------------------
vector <double> exactForce(dim,0.);
vector <double> exactVel_d(dim,0.);
for (unsigned k = 0; k < dim; k++){
    exactForce[k] = - IRe * exact_lap_stateVel[k] - IRe * exact_lap_ctrlVel[k] 
                    + advection_flag * (exact_conv_u_nabla_u[k] + exact_conv_u_nabla_uctrl[k] + exact_conv_uctrl_nabla_u[k] + exact_conv_uctrl_nabla_uctrl[k]) 
                    + exact_grad_statePress[k];
    exactVel_d[k] =   exact_stateVel[k] + exact_ctrlVel[k] 
                    + (1./cost_functional_coeff) * ( IRe * exact_lap_stateVel[k] - exact_grad_statePress[k]) 
                    + (1./cost_functional_coeff) * advection_flag * (exact_conv_u_nabla_uadj[k] - exact_conv_nabla_uT_uadj[k] - exact_conv_nabla_uctrlT_uadj[k] + exact_conv_uctrl_nabla_uadj[k]);
}

//computation of RHS (force and desired velocity) using MMS=============================================== 
#endif



        // *** phiV_i loop ***
        for (unsigned i = 0; i < nDofsV; i++) {
          vector < adept::adouble > NSV_gss(dim, 0.);
	  vector < adept::adouble > NSVadj_gss(dim, 0.);
	  vector < adept::adouble > NSVctrl_gss(dim, 0.);
	  
          for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      
          for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { //focus on single partial derivative
	    
          
	      NSV_gss[kdim]   	 	+=  IRe*phiV_x_gss[i * dim_offset_grad /*space_dim*/ + jdim]*gradSolV_gss[kdim][jdim]; 
						    // /*deformation_tensor*//*IRe * phiV_x_gss[i * dim + jdim] * 
						    // (gradSolV_gss[kdim][jdim] + gradSolV_gss[jdim][kdim])*/;  //diffusion
	      NSV_gss[kdim] 		+= IRe*phiV_x_gss[i * dim_offset_grad /*space_dim*/ + jdim]*gradSolVctrl_gss[kdim][jdim];	 //delta_state-control
	      NSVadj_gss[kdim]   	+=  IRe*phiVadj_x_gss[i * dim_offset_grad /*space_dim*/ + jdim]*gradSolVadj_gss[kdim][jdim];  
	      NSVctrl_gss[kdim]   	+=   beta * phiVctrl_x_gss[i * dim_offset_grad /*space_dim*/ + jdim] * gradSolVctrl_gss[kdim][jdim];
	      NSVctrl_gss[kdim] 	+=   - IRe*phiVctrl_x_gss[i * dim_offset_grad /*space_dim*/ + jdim]*gradSolVadj_gss[kdim][jdim];  //nabla_delta_control-nabla_adjoint
	  }  //jdim loop

          for (unsigned jdim = 0; jdim < dim; jdim++) { //focus on single partial derivative
 
          NSV_gss[kdim]  		+=  advection_flag * phiV_gss[i] * (solV_gss[jdim] * gradSolV_gss[kdim][jdim]);                                     //advection (u_hat . \nabla) u_hat
	      NSV_gss[kdim] 		+=  advection_flag * phiV_gss[i] * (solV_gss[jdim] * gradSolVctrl_gss[kdim][jdim]);                                  //advection (u_hat . \nabla) u_0
	      NSV_gss[kdim] 		+=  advection_flag * phiV_gss[i] * (solVctrl_gss[jdim] * gradSolV_gss[kdim][jdim]);                                  //advection (u_0 . \nabla) u_hat 
	      NSV_gss[kdim] 		+=  advection_flag * phiV_gss[i] * (solVctrl_gss[jdim] * gradSolVctrl_gss[kdim][jdim]);                              //advection (u_0 . \nabla) u_0
	      
	      NSVadj_gss[kdim]		+=   advection_flag * phiVadj_gss[i] * gradSolV_gss[jdim][kdim] * solVadj_gss[jdim];           // c(delta u,u,lambda)
	      NSVadj_gss[kdim]		+=   advection_flag * solV_gss[jdim] * phiVadj_x_gss[i * dim + jdim] * solVadj_gss[kdim] ;     // c(u,delta u, lambda)
	      NSVadj_gss[kdim]		+=   advection_flag * phiVadj_gss[i] * gradSolVctrl_gss[jdim][kdim] * solVadj_gss[jdim];       // c(delta u,u0,lambda)
	      NSVadj_gss[kdim]		+=   advection_flag * solVctrl_gss[jdim] * phiVadj_x_gss[i * dim + jdim] * solVadj_gss[kdim] ;       // c(u0,delta u,lambda)
	      
	      NSVctrl_gss[kdim]		+= -  advection_flag * solV_gss[jdim] * phiVctrl_x_gss[i * dim + jdim] * solVadj_gss[kdim] ;     // c(u,delta u0, lambda)
	      NSVctrl_gss[kdim]		+= -  advection_flag * phiVctrl_gss[i] * gradSolV_gss[jdim][kdim] * solVadj_gss[jdim];           // c(delta u0,u,lambda)
	      NSVctrl_gss[kdim]		+= -  advection_flag * phiVctrl_gss[i] * gradSolVctrl_gss[jdim][kdim] * solVadj_gss[jdim];       // c(delta u0,u0,lambda)
	      NSVctrl_gss[kdim]		+= -  advection_flag * solVctrl_gss[jdim] * phiVctrl_x_gss[i * dim + jdim] * solVadj_gss[kdim] ;       // c(u0,delta u0,lambda)
						  
	  }  //jdim loop
	  
 #if exact_sol_flag == 0
              NSV_gss[kdim]     += - force[kdim] * phiV_gss[i];
	      NSVadj_gss[kdim] 		+=  + cost_functional_coeff* target_flag * DesiredTargetVel()[kdim] * phiVadj_gss[i];
  	      NSVctrl_gss[kdim]   	+=  - cost_functional_coeff* target_flag * DesiredTargetVel()[kdim] * phiVctrl_gss[i];
#endif
 #if exact_sol_flag == 1
              NSV_gss[kdim]     += - exactForce[kdim] * phiV_gss[i];
	      NSVadj_gss[kdim] 		+=  + cost_functional_coeff* target_flag * exactVel_d[kdim] * phiVadj_gss[i];
 	      NSVctrl_gss[kdim]   	+=  - cost_functional_coeff* target_flag * exactVel_d[kdim] * phiVctrl_gss[i];
#endif        
              
              
          NSVadj_gss[kdim]		+=  - cost_functional_coeff * target_flag * solV_gss[kdim]*phiVadj_gss[i]; //delta_adjoint-state
	      NSVadj_gss[kdim] 		+=  - cost_functional_coeff * target_flag * solVctrl_gss[kdim]*phiVadj_gss[i]; //delta_adjoint-control
	      NSVctrl_gss[kdim] 	+=    cost_functional_coeff * target_flag             * solV_gss[kdim]*phiVctrl_gss[i]; //delta_control-state
	      NSVctrl_gss[kdim]   	+=   (cost_functional_coeff * target_flag + alpha) * solVctrl_gss[kdim] * phiVctrl_gss[i] ;
           
            //velocity-pressure block
          NSV_gss[kdim] 	    += - solP_gss * phiV_x_gss[i * dim_offset_grad /*space_dim*/ + kdim];
	      NSVadj_gss[kdim] 	    += - solPadj_gss * phiVadj_x_gss[i * dim_offset_grad /*space_dim*/ + kdim];
          NSVctrl_gss[kdim] 	+= - solPctrl_gss * phiVctrl_x_gss[i * dim_offset_grad /*space_dim*/ + kdim];
	    
	  } //kdim loop


          for (unsigned  kdim = 0; kdim < dim; kdim++) { // - (b-Ax)
            aResV[kdim][i] 		+=   NSV_gss[kdim]     * AbsDetJxWeight_iqp;
	    aResVadj[kdim][i]   	+=   NSVadj_gss[kdim]  * AbsDetJxWeight_iqp;
            aResVctrl[kdim][i]    	+=   NSVctrl_gss[kdim] * AbsDetJxWeight_iqp;
							 	    
	  }
        } // end phiV_i loop

        // *** phiP_i loop ***
        for (unsigned i = 0; i < nDofsP; i++) {
          for (int kdim = 0; kdim < dim; kdim++) {
            aResP[i] 		+= - (gradSolV_gss[kdim][kdim]) * phiP_gss[i]  * AbsDetJxWeight_iqp;
	    aResPadj[i]  	+= - (gradSolVadj_gss[kdim][kdim]) * phiPadj_gss[i]  * AbsDetJxWeight_iqp;
	    aResPctrl[i]   	+= - (gradSolVctrl_gss[kdim][kdim]) * phiPctrl_gss[i]  * AbsDetJxWeight_iqp;
	    
          }
        } // end phiP_i loop
      } // end gauss point loop
      
              



    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize(nDofsVP_tot);    //resize

    for (int i = 0; i < nDofsV; i++) {
      for (unsigned  kdim = 0; kdim < dim; kdim++) {
        Res[ i +  kdim * nDofsV ]            =  - aResV[kdim][i].value();
	Res[ i +  kdim * nDofsV + nDofsVP]   =  - aResVadj[kdim][i].value();
	Res[ i +  kdim * nDofsV + 2*nDofsVP] =  - aResVctrl[kdim][i].value();
      }
    }

    for (int i = 0; i < nDofsP; i++) {
      Res[ i + dim * nDofsV ]            = - aResP[i].value();
      Res[ i + dim * nDofsV + nDofsVP]   = - aResPadj[i].value();
      Res[ i + dim * nDofsV + 2*nDofsVP] = - aResPctrl[i].value();
    }
  

    RES->add_vector_blocked(Res, L2G_dofmap_Mat);

    //Extarct and store the Jacobian
       Jac.resize(nDofsVP_tot * nDofsVP_tot);

      // define the dependent variables
      for (unsigned  kdim = 0; kdim < dim; kdim++) { s.dependent(&aResV[kdim][0], nDofsV);}              s.dependent(&aResP[0], nDofsP);
      for (unsigned  kdim = 0; kdim < dim; kdim++) { s.dependent(&aResVadj[kdim][0], nDofsVadj); }       s.dependent(&aResPadj[0], nDofsPadj);
      for (unsigned  kdim = 0; kdim < dim; kdim++) { s.dependent(&aResVctrl[kdim][0], nDofsVctrl);  }    s.dependent(&aResPctrl[0], nDofsPctrl);
      
      // define the independent variables
      for (unsigned  kdim = 0; kdim < dim; kdim++) {  s.independent(&solV[kdim][0], nDofsV); }         s.independent(&solP[0], nDofsP);
      for (unsigned  kdim = 0; kdim < dim; kdim++) { s.independent(&solVadj[kdim][0], nDofsVadj); }    s.independent(&solPadj[0], nDofsPadj);
      for (unsigned  kdim = 0; kdim < dim; kdim++) { s.independent(&solVctrl[kdim][0], nDofsVctrl); }  s.independent(&solPctrl[0], nDofsPctrl);

      // get the and store jacobian matrix (row-major)
      s.jacobian(&Jac[0] , true);
      
      JAC->add_matrix_blocked(Jac, L2G_dofmap_Mat, L2G_dofmap_Mat);
 
      s.clear_independents();
      s.clear_dependents();
    
  } //end element loop for each process

  RES->close();
//   RES->print();

  JAC->close();

  // ***************** END ASSEMBLY *******************
}


void ComputeIntegral(const MultiLevelProblem& ml_prob) {

   const NonLinearImplicitSystem * mlPdeSys   = & ml_prob.get_system<NonLinearImplicitSystem> ("NSOpt");   // pointer to the nonlinear implicit system named "NSOpt"
 
  const unsigned level = mlPdeSys->GetLevelToAssemble();
 

  Mesh*          msh          	= ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         	= msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  ml_sol    = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        	= ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object
  
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

  // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //geometry *******************************
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  unsigned solType_coords = coordXType;  //FE_DOMAIN = 0; //we do linear FE this time // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
  const unsigned int dim_offset_grad = /*dim*/  3  /*2*/    ;
 
  vector < vector < double > > coordX(dim);    // local coordinates

  for (unsigned  k = 0; k < dim; k++) { 
    coordX[k].reserve(max_size);
  }
  
  double AbsDetJxWeight_iqp;
  
  
  //geometry *******************************

//STATE######################################################################
  vector < unsigned > solVIndex(dim);
  solVIndex[0] = ml_sol->GetIndex("u_0");
  solVIndex[1] = ml_sol->GetIndex("u_1");

  if (dim == 3) solVIndex[2] = ml_sol->GetIndex("u_2");

  unsigned solVType = ml_sol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"
  
  vector < vector < double > >  solV(dim);    // local solution
  vector <double >  V_gss(dim, 0.);    //  solution
   
 for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(max_size);
  }

  
  vector <double> phiV_gss;  // local test function
  vector <double> phiV_x_gss; // local test function first order partial derivatives

  phiV_gss.reserve(max_size);
  phiV_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);
  
//STATE######################################################################
  

//CONTROL######################################################################
  vector < unsigned > solVctrlIndex(dim);
  solVctrlIndex[0] = ml_sol->GetIndex("ctrl_0");
  solVctrlIndex[1] = ml_sol->GetIndex("ctrl_1");

  if (dim == 3) solVctrlIndex[2] = ml_sol->GetIndex("ctrl_2");

  unsigned solVctrlType = ml_sol->GetSolutionType(solVctrlIndex[0]);
  
  vector < vector < double > >  solVctrl(dim);    // local solution
  vector < double >   Vctrl_gss(dim, 0.);    //  solution
   
 for (unsigned  k = 0; k < dim; k++) {
    solVctrl[k].reserve(max_size);
  }

  
  vector <double> phiVctrl_gss;  // local test function
  vector <double> phiVctrl_x_gss; // local test function first order partial derivatives

  phiVctrl_gss.reserve(max_size);
  phiVctrl_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);
  
//CONTROL######################################################################

// Vel_desired##################################################################
  vector <double> phiVdes_gss;  // local test function
  vector <double> phiVdes_x_gss; // local test function first order partial derivatives

  phiVdes_gss.reserve(max_size);
  phiVdes_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);

  vector <double>  solVdes(dim,0.);
  vector<double> Vdes_gss(dim, 0.);  
  
// Vel_desired##################################################################




double  integral_target_alpha = 0.;
double	integral_beta   = 0.;
double	integral_gamma  = 0.;
double  integral_div_ctrl = 0.;



//*************************************************** 
     std::vector < std::vector < double > >  JacI_iqp(space_dim);
     std::vector < std::vector < double > >  Jac_iqp(dim);
    for (unsigned d = 0; d < Jac_iqp.size(); d++) {   Jac_iqp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp.size(); d++) { JacI_iqp[d].resize(dim); }
    
    double detJac_iqp;

    //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
//*************************************************** 

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

   // geometry *****************************
      geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
      const short unsigned ielGeom = geom_element_iel.geom_type();
  // geometry end *****************************
   
// equation
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);
//     unsigned nDofsVdes = msh->GetElementDofNumber(iel, solVType); 
    unsigned nDofsVctrl = msh->GetElementDofNumber(iel, solVctrlType);
    
     
    for (unsigned  k = 0; k < dim; k++)  {
      solV[k].resize(nDofsV);
      solVctrl[k].resize(nDofsVctrl);
//       solVdes[k].resize(nDofsVdes);
    }
  //*************************************** 
  
  //***** set target domain flag ********************************** 
   geom_element_iel.set_elem_center_3d(iel, solType_coords);

   int target_flag = 0;
   target_flag = ElementTargetFlag(geom_element_iel.get_elem_center_3d());
//***************************************       
    
 //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = ControlDomainFlag_internal_restriction(geom_element_iel.get_elem_center_3d());
 //*************************************************** 
   
    
 //STATE###################################################################  
    // velocity ************
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
      }
    }
//STATE###################################################################

//CONTROL###################################################################  
    // velocity ************
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVctrlDof = msh->GetSolutionDof(i, iel, solVctrlType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solVctrl[k][i] = (*sol->_Sol[solVctrlIndex[k]])(solVctrlDof);      // global extraction and local storage for the solution
      }
    }
//CONTROL###################################################################




  //DESIRED VEL###################################################################  
    // velocity ************
//     for (unsigned i = 0; i < nDofsV; i++) {
//       unsigned solVdesDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < solVdes.size() /*dim*/; k++) {
        solVdes[k]/*[i]*/ = DesiredTargetVel()[k] /*(*sol->_Sol[solVIndex[k]])(solVdesDof)*/;      // global extraction and local storage for the solution
      }
//     }
 //DESIRED VEL###################################################################

 
 
 
      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {

//STATE#############################################################################	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), ig, Jac_iqp, JacI_iqp, detJac_iqp, space_dim);
    AbsDetJxWeight_iqp = detJac_iqp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[ig];
   
    elem_all[ielGeom][solVType]->shape_funcs_current_elem(ig, JacI_iqp, phiV_gss, phiV_x_gss, boost::none , space_dim);
    elem_all[ielGeom][solVctrlType]->shape_funcs_current_elem(ig, JacI_iqp, phiVctrl_gss, phiVctrl_x_gss,  boost::none, space_dim);
    elem_all[ielGeom][solVType /*solVdes*/]->shape_funcs_current_elem(ig, JacI_iqp, phiVdes_gss, phiVdes_x_gss,  boost::none, space_dim);

	
	  vector < vector < double > > gradVctrl_gss(dim);
      for (unsigned  k = 0; k < dim; k++) {
          gradVctrl_gss[k].resize(dim_offset_grad /*space_dim*/);
          std::fill(gradVctrl_gss[k].begin(), gradVctrl_gss[k].end(), 0);
        }
	
    for (unsigned  k = 0; k < dim; k++) {
      V_gss[k]       = 0.;
      Vdes_gss[k]    = 0.;
       Vctrl_gss[k]  = 0.;
    }
    
      for (unsigned i = 0; i < nDofsV; i++) {
        for (unsigned  k = 0; k < dim; k++) {
            V_gss[k] += solV[k][i] * phiV_gss[i];
            Vdes_gss[k] += solVdes[k]/*[i]*/ * phiVdes_gss[i];
		}
      }
	
      for (unsigned i = 0; i < nDofsVctrl; i++) {
        for (unsigned  k = 0; k < dim; k++) {
            Vctrl_gss[k] += solVctrl[k][i] * phiVctrl_gss[i];
	 }
     for (unsigned j = 0; j < dim_offset_grad /*space_dim*/; j++) {
            for (unsigned  k = 0; k < dim; k++) {
              gradVctrl_gss[k][j] += phiVctrl_x_gss[i * dim_offset_grad /*space_dim*/ + j] * solVctrl[k][i];
            }
          }
      }
          
//                 for (unsigned i = 0; i < nDofsV; i++) {

      for (unsigned  k = 0; k < dim; k++) {
          integral_div_ctrl +=  AbsDetJxWeight_iqp * gradVctrl_gss[k][k] /** phiVctrl_gss[i]*/;
      }
//       }
	
      for (unsigned  k = 0; k < dim; k++) {
	 integral_target_alpha +=  target_flag * (V_gss[k] + Vctrl_gss[k] - Vdes_gss[k]) * (V_gss[k] + Vctrl_gss[k] - Vdes_gss[k]) * AbsDetJxWeight_iqp; 
	 integral_beta	+=  control_el_flag * ((Vctrl_gss[k]) * (Vctrl_gss[k]) * AbsDetJxWeight_iqp);
      }
      for (unsigned  k = 0; k < dim; k++) {
	for (unsigned  j = 0; j < dim; j++) {	
		integral_gamma	  +=  control_el_flag * gradVctrl_gss[k][j] * gradVctrl_gss[k][j] * AbsDetJxWeight_iqp;
	}
      }
   
  
      }// end gauss point loop
    } //end element loop  

       std::ostringstream filename_out; filename_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "Integral_computation"  << ".txt";

       std::ofstream intgr_fstream;
  if (paral::get_rank() == 0 ) {
      intgr_fstream.open(filename_out.str().c_str(),std::ios_base::app);
      intgr_fstream << " ***************************** Non Linear Iteration "<< mlPdeSys->GetNonlinearIt() << " *********************************** " <<  std::endl << std::endl;
      intgr_fstream << "The value of the target functional for " << "alpha " <<   std::setprecision(0) << std::scientific << cost_functional_coeff << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_target_alpha << std::endl;
      intgr_fstream << "The value of the L2 control for        " << "beta  " <<   std::setprecision(0) << std::scientific << alpha  << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_beta         << std::endl;
      intgr_fstream << "The value of the H1 control for        " << "gamma " <<   std::setprecision(0) << std::scientific << beta << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_gamma        << std::endl;
      intgr_fstream << "The value of the total integral is " << std::setw(11) << std::setprecision(10) <<  integral_target_alpha * cost_functional_coeff*0.5  + integral_beta *alpha*0.5 + integral_gamma *beta*0.5 << std::endl;
      intgr_fstream << "The value of the divergence of the control is " << std::setw(11) << std::setprecision(10) <<  integral_div_ctrl << std::endl;
      intgr_fstream <<  std::endl;
      intgr_fstream.close();  //you have to close to disassociate the file from the stream
}  

    
//     std::cout << "The value of the integral of target for alpha "<< std::setprecision(0)<< std::scientific<<  cost_functional_coeff<< " is " << std::setw(11) << std::setprecision(10) << std::fixed<< integral_target_alpha << std::endl;
//     std::cout << "The value of the integral of beta for beta "<<  std::setprecision(0)<<std::scientific<<alpha << " is " << std::setw(11) << std::setprecision(10) <<  std::fixed<< integral_beta << std::endl;
//     std::cout << "The value of the integral of gamma for gamma "<< std::setprecision(0)<<std::scientific<<beta<< " is " << std::setw(11) << std::setprecision(10) <<  std::fixed<< integral_gamma << std::endl; 
//     std::cout << "The value of the total integral is " << std::setw(11) << std::setprecision(10) <<  integral_target_alpha *(cost_functional_coeff*0.5)+ integral_beta *(alpha*0.5) + integral_gamma*(beta*0.5) << std::endl; 
   
    
    return; 
	  
  
}



void AssembleNavierStokesOpt_nonAD(MultiLevelProblem& ml_prob) {
     
 std::cout << " ********************************  NON-AD SYSTEM ******************************************** " << std::endl;

  // ======= Main objects - BEGIN =======
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
  constexpr bool print_algebra_local = true;

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
  
  // ======= Main objects - END =======
  
  
  
  // ======= Geometry at Dofs - BEGIN  =======
  unsigned coordXType = 2; /*BIQUADR_FE*/// get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  unsigned solType_coords = coordXType;  //FE_DOMAIN = 0; //we do linear FE this time // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  vector< vector < double> > coordX(dim);
  for(int i = 0; i < dim; i++) { coordX[i].reserve(max_size);  }

  constexpr unsigned int space_dim = 3;
 
  const unsigned int dim_offset_grad = /*dim*/  3  /*2*/    ;
  // ======= Geometry at Dofs - END  =======
  
  // ======= Solutions, Unknowns - BEGIN =======
  const int n_vars_state = dim + 1;
  const int vel_type_pos = 0;
  const int press_type_pos = dim;

  const int state_pos_begin = 0;               ///@todo make in agreement with Unknowns function
  const int adj_pos_begin   =     n_vars_state;
  const int ctrl_pos_begin  = 2 * n_vars_state;
  const int mu_pos_begin    = 3 * n_vars_state;

  
  std::vector< Unknown > unknowns = provide_list_of_unknowns( dim );
  
  const int n_unknowns = unknowns.size();
 
  
  vector < unsigned > SolPdeIndex(n_unknowns);
  vector < unsigned > SolIndex(n_unknowns);  
  vector < unsigned > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolPdeIndex[ivar]	= mlPdeSys->GetSolPdeIndex(unknowns[ivar]._name.c_str());
    SolIndex[ivar]	= ml_sol->GetIndex        (unknowns[ivar]._name.c_str());
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);
  }

  vector < unsigned > Sol_n_el_dofs_Mat_vol(n_unknowns);
  
  
      const unsigned int n_components_ctrl = dim;
  double penalty_outside_control_domain = 1.e20;         ///@todo  this number affects convergence or not! // penalty for zero control outside 
  // ======= Solutions, Unknowns - END =======

      
  // ======= Solutions, not Unknowns - BEGIN =======

  //MU
  //************** variables for ineq constraints: act flag ****************************   
  std::vector<unsigned int> solIndex_act_flag_sol(n_components_ctrl); 
  
  ctrl_inequality::store_act_flag_in_old(mlPdeSys, ml_sol, sol, solIndex_act_flag_sol);
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
      
  // ======= Solutions, not Unknowns - END =======
 

  
  
  // ======= Solutions, Unknowns at dofs - BEGIN =======
  vector < vector < double > > Sol_eldofs_Mat(n_unknowns);
  vector < vector < double > > gradSol_eldofs_Mat(n_unknowns);
  
  for(int k=0; k<n_unknowns; k++) {
    Sol_eldofs_Mat[k].reserve(max_size);
    gradSol_eldofs_Mat[k].reserve(max_size*dim);    
  }
  // ======= Solutions, Unknowns at dofs - END =======

  // ======= Solutions, Unknowns at quadrature points - BEGIN =======
  vector < double > SolVAR_qp(n_unknowns);
    vector < vector < double > > gradSolVAR_qp(n_unknowns);
    for(int k=0; k<n_unknowns; k++) {  gradSolVAR_qp[k].resize(dim_offset_grad /*space_dim*/);  }
  // ======= Solutions, Unknowns at quadrature points - END =======
      
  //============ Geometry at Quadrature points - BEGIN ==============================================================================
     std::vector < std::vector < double > >  JacI_iqp(space_dim);
     std::vector < std::vector < double > >  Jac_iqp(dim);
    for (unsigned d = 0; d < Jac_iqp.size(); d++) {   Jac_iqp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp.size(); d++) { JacI_iqp[d].resize(dim); }
    
    double AbsDetJxWeight_iqp = 0.;
    double detJac_iqp;
  //============ Geometry at Quadrature points - END ==============================================================================

   
// ======= FE at Quadrature, all - BEGIN =======
    //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
  
  //==========================================================================================
  vector < vector < double > > phi_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_x_gss_fe(NFE_FAMS);
 
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_gss_fe[fe].reserve(max_size);
      phi_x_gss_fe[fe].reserve(max_size * dim_offset_grad /*space_dim*/);
   }
// ======= FE at Quadrature, all - END =======


  
  // ======= Equation, local - BEGIN =======
  vector < vector < int > > L2G_dofmap_Mat(n_unknowns); 
  vector < vector < double > > Res(n_unknowns);
  vector < vector < vector < double > > > Jac(n_unknowns);
 
  for(int i = 0; i < n_unknowns; i++) {     
    L2G_dofmap_Mat[i].reserve(max_size);
      Res[i].reserve(max_size);
  }
   
  if(assembleMatrix) {
    for(int i = 0; i < n_unknowns; i++) {
      Jac[i].resize(n_unknowns);    
      for(int j = 0; j < n_unknowns; j++) {
	Jac[i][j].reserve(max_size*max_size);	
      }
    }
  }
  // ======= Equation, local - END =======
  
  // ======= Equation, global - BEGIN =======
   RES->zero();
    if(assembleMatrix) JAC->zero();
  // ======= Equation, global - END =======
  
    
    
  // ======= Parameters - BEGIN ======= 
  double IRe = ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
  // ======= Parameters - END =======

  
   
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

  // geometry *****************************
      geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
      const short unsigned ielGeom = geom_element_iel.geom_type();
  // geometry end *****************************

  
  // equation *****************************
    unsigned nDofsV = msh->GetElementDofNumber(iel, SolFEType[vel_type_pos]);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, SolFEType[state_pos_begin + press_type_pos]);    // number of solution element dofs
    
    unsigned nDofsVadj = msh->GetElementDofNumber(iel,SolFEType[adj_pos_begin]);    // number of solution element dofs
    unsigned nDofsPadj = msh->GetElementDofNumber(iel,SolFEType[adj_pos_begin + press_type_pos]);    // number of solution element dofs

    unsigned nDofsVctrl = msh->GetElementDofNumber(iel,SolFEType[ctrl_pos_begin]);    // number of solution element dofs
    unsigned nDofsPctrl = msh->GetElementDofNumber(iel,SolFEType[ctrl_pos_begin + press_type_pos] );    // number of solution element dofs

    unsigned nDofsVP = dim * nDofsV + nDofsP;
    unsigned nDofsVP_tot = 3*nDofsVP;
  // equation end *****************************

    
    
    
    
  //***** set target domain flag ********************************** 
   geom_element_iel.set_elem_center_3d(iel, solType_coords);

   int target_flag = 0;
   target_flag = ElementTargetFlag(geom_element_iel.get_elem_center_3d()/*elem_center*/);
   //***************************************       
   
   //###################################################################  
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
	Sol_n_el_dofs_Mat_vol[k]=ndofs_unk;
       Sol_eldofs_Mat[k].resize(ndofs_unk);
       L2G_dofmap_Mat[k].resize(ndofs_unk); 
    for (unsigned i = 0; i < ndofs_unk; i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);    // global to global mapping between solution node and solution dof // via local to global solution node
       Sol_eldofs_Mat[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);      // global extraction and local storage for the solution
       L2G_dofmap_Mat[k][i] = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    
    //###################################################################
    for(int ivar=0; ivar<n_unknowns; ivar++) {
      
      Res[SolPdeIndex[ivar]].resize(Sol_n_el_dofs_Mat_vol[ivar]);
      memset(&Res[SolPdeIndex[ivar]][0],0.,Sol_n_el_dofs_Mat_vol[ivar]*sizeof(double));
    }
   
    for(int ivar=0; ivar<n_unknowns; ivar++) {
      for(int jvar=0; jvar<n_unknowns; jvar++) {
      if(assembleMatrix){  //MISMATCH
	Jac[ SolPdeIndex[ivar] ][ SolPdeIndex[jvar] ].resize(Sol_n_el_dofs_Mat_vol[ivar]*Sol_n_el_dofs_Mat_vol[jvar]);
	memset(&Jac[SolPdeIndex[ivar]][SolPdeIndex[jvar]][0],0.,Sol_n_el_dofs_Mat_vol[ivar]*Sol_n_el_dofs_Mat_vol[jvar]*sizeof(double));
      }
    }
  }
  
    //=============================================================================

    


 //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = ControlDomainFlag_internal_restriction(geom_element_iel.get_elem_center_3d());

  std::vector< std::vector< int > > control_node_flag(n_components_ctrl);
       
	  for (unsigned c = 0; c < n_components_ctrl; c++) {
              control_node_flag[c].resize(Sol_n_el_dofs_Mat_vol[ctrl_pos_begin + c]);
              std::fill(control_node_flag[c].begin(), control_node_flag[c].end(), 0);   
         }

  if (control_el_flag == 1) {
	  for (unsigned c = 0; c < n_components_ctrl; c++) {
      std::fill(control_node_flag[c].begin(), control_node_flag[c].end(), 1);
         }
  }
  //*************************************************** 
   
 ///@todo if I want to restrict the control lifting, I might just set to zero the dof values here! Also, I have to remove the equations for the unneeded dofs with a PENALTY!
 /// Well, this can also be done in the Initialization function I think... 
 /// The problem there is that it is only DOF-BASED, it doesn't receive the ELEMENT INFORMATION. We should change that

//    if (control_el_flag == 0) {
// 	  for (unsigned c = 0; c < n_components_ctrl; c++) {
//       std::fill(Sol_eldofs_Mat[ctrl_pos_begin + c].begin(), Sol_eldofs_Mat[ctrl_pos_begin + c].end(), 0.);
//          }
//    }

// it seems that with this the convergence for the control variables is worse...
 
 
   
      for(unsigned iqp = 0; iqp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); iqp++) {
	
    elem_all[ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), iqp, Jac_iqp, JacI_iqp, detJac_iqp, space_dim);
    AbsDetJxWeight_iqp = detJac_iqp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[iqp];
   
     for(int fe=0; fe < NFE_FAMS; fe++) {
    elem_all[ielGeom][fe]->shape_funcs_current_elem(iqp, JacI_iqp,phi_gss_fe[fe],phi_x_gss_fe[fe], boost::none, space_dim);
      }


 //begin unknowns eval at gauss points ********************************
	for(unsigned unk = 0; unk < n_unknowns; unk++) {
	  SolVAR_qp[unk] = 0.;
	  for(unsigned ivar2=0; ivar2<dim_offset_grad /*space_dim*/; ivar2++){ 
	    gradSolVAR_qp[unk][ivar2] = 0.; 
	  }
	  
	  for(unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[unk]; i++) {
	    SolVAR_qp[unk] += phi_gss_fe[ SolFEType[unk] ][i] * Sol_eldofs_Mat[unk][i];
	    for(unsigned ivar2=0; ivar2 < dim_offset_grad /*space_dim*/; ivar2++) {
	      gradSolVAR_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType[unk] ][i*dim_offset_grad /*space_dim*/+ivar2] * Sol_eldofs_Mat[unk][i]; 
	    }
	  }
	  
	}  
 //end unknowns eval at gauss points ********************************
	
      vector < double > coordX_gss(dim, 0.);
 	for(unsigned k = 0; k <  dim; k++) {
	  for(unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[k]; i++) {
         coordX_gss[k] += coordX[k][i] * phi_gss_fe[ SolFEType[k] ][i];
      }
    }
	
//  // I x = 5 test ********************************
// 	for(unsigned i_unk = 0; i_unk<n_unknowns; i_unk++) { 
// 	    for(unsigned i_dof=0; i_dof < Sol_n_el_dofs_Mat_vol[i_unk]; i_dof++) {
// 		/*if ( i_unk!=0 && i_unk!=1 && i_unk!=2 && i_unk!=3 && i_unk!=4 && i_unk!=6 && i_unk!=7 )*/  Res[SolPdeIndex[i_unk]][i_dof] +=  (               0.* phi_gss_fe[SolFEType[i_unk]][i_dof] 
// 		                                    - SolVAR_qp[i_unk]*phi_gss_fe[SolFEType[i_unk]][i_dof] )*AbsDetJxWeight_iqp;
// 		  for(unsigned j_unk = 0; j_unk<n_unknowns; j_unk++) {
// 		  	for(unsigned j_dof=0; j_dof < Sol_n_el_dofs_Mat_vol[j_unk]; j_dof++) {
// 			  
// 		              if (i_unk==j_unk /*&& i_unk!=0 && i_unk!=1 && i_unk!=2 && i_unk!=3 && i_unk!=4 && i_unk!=6 && i_unk!=7*/)   {
// 				Jac[ SolPdeIndex[i_unk] ][ SolPdeIndex[j_unk] ][ i_dof*Sol_n_el_dofs_Mat_vol[i_unk] + j_dof ] += 
// 				        ( phi_gss_fe[SolFEType[i_unk]][i_dof]*phi_gss_fe[SolFEType[j_unk]][j_dof] )*AbsDetJxWeight_iqp;
// 			      }
// 			  
// 			} //j_dof
// 		  }  //j_unk
// 	    }  //i_dof
// 	}  //i_unk
//  // I x = 5 test ********************************
 

#if exact_sol_flag == 1
//computation of RHS (force and desired velocity) using MMS - BEGIN =============================================== 
//state values--------------------
vector <double>  exact_stateVel(dim, 0.);
   mms_state_control::value_stateVel(coordX_gss, exact_stateVel);
vector < vector < double > > exact_grad_stateVel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_stateVel[k].resize(dim);
    std::fill(exact_grad_stateVel[k].begin(), exact_grad_stateVel[k].end(), 0.);
}
   mms_state_control::gradient_stateVel(coordX_gss,exact_grad_stateVel);
vector <double>  exact_lap_stateVel(dim, 0.);
   mms_state_control::laplace_stateVel(coordX_gss, exact_lap_stateVel);
vector <double> exact_grad_statePress(dim, 0.);
   mms_state_control::gradient_statePress(coordX_gss, exact_grad_statePress);

//control values-------------------------------
vector <double>  exact_ctrlVel(dim);
   mms_state_control::value_ctrlVel(coordX_gss, exact_ctrlVel);
vector < vector < double > > exact_grad_ctrlVel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_ctrlVel[k].resize(dim);
    std::fill(exact_grad_ctrlVel[k].begin(), exact_grad_ctrlVel[k].end(), 0.);
}
   mms_state_control::gradient_ctrlVel(coordX_gss,exact_grad_ctrlVel);
vector <double>  exact_lap_ctrlVel(dim);
   mms_state_control::laplace_ctrlVel(coordX_gss, exact_lap_ctrlVel);

//convection terms from delta_state-------------------------------------
vector <double>  exact_conv_u_nabla_u(dim,0.);
vector <double>  exact_conv_u_nabla_uctrl(dim,0.);
vector <double>  exact_conv_uctrl_nabla_u(dim,0.);
vector <double>  exact_conv_uctrl_nabla_uctrl(dim,0.);

for (unsigned k = 0; k < dim; k++){
    for (unsigned i = 0; i < dim; i++){
    exact_conv_u_nabla_u[k] += exact_grad_stateVel[k][i] * exact_stateVel[i] ; 
    exact_conv_u_nabla_uctrl[k] += exact_grad_ctrlVel[k][i] * exact_stateVel[i] ; 
    exact_conv_uctrl_nabla_u[k] += exact_grad_stateVel[k][i] * exact_ctrlVel[i] ; 
    exact_conv_uctrl_nabla_uctrl[k] += exact_grad_ctrlVel[k][i] * exact_ctrlVel[i] ; 
    }
}

//convection terms from delta_adjoint-------------------------
vector <double>  exact_conv_u_nabla_uadj(dim,0.);
vector <double>  exact_conv_nabla_uT_uadj(dim,0.);
vector <double>  exact_conv_nabla_uctrlT_uadj(dim,0.);
vector <double>  exact_conv_uctrl_nabla_uadj(dim,0.);

for (unsigned k = 0; k < dim; k++){
    for (unsigned i = 0; i < dim; i++){
    exact_conv_u_nabla_uadj[k] += exact_grad_stateVel[k][i] * exact_stateVel[i] ; 
    exact_conv_nabla_uT_uadj[k] += exact_grad_stateVel[i][k] * exact_stateVel[i];
    exact_conv_nabla_uctrlT_uadj[k] += exact_grad_ctrlVel[i][k] * exact_stateVel[i];  
    exact_conv_uctrl_nabla_uadj[k] += exact_grad_stateVel[k][i] * exact_ctrlVel[i] ; 
    }
}

//force and desired velocity ---------------------------------------------
vector <double> exactForce(dim,0.);
vector <double> exactVel_d(dim,0.);
for (unsigned k = 0; k < dim; k++){
    exactForce[k] = - IRe * exact_lap_stateVel[k] - IRe * exact_lap_ctrlVel[k] 
                    + advection_flag * (exact_conv_u_nabla_u[k] + exact_conv_u_nabla_uctrl[k] + exact_conv_uctrl_nabla_u[k] + exact_conv_uctrl_nabla_uctrl[k]) 
                    + exact_grad_statePress[k];
    exactVel_d[k] =   exact_stateVel[k] + exact_ctrlVel[k] 
                    + (1./cost_functional_coeff) * ( IRe * exact_lap_stateVel[k] - exact_grad_statePress[k]) 
                    + (1./cost_functional_coeff) * advection_flag * (exact_conv_u_nabla_uadj[k] - exact_conv_nabla_uT_uadj[k] - exact_conv_nabla_uctrlT_uadj[k] + exact_conv_uctrl_nabla_uadj[k]);
}

//computation of RHS (force and desired velocity) using MMS - END =============================================== 
#endif

 
 
 
//============ delta_state row - BEGIN  ============================================================================================

//************ Residual, Velocity, BEGIN *********************

    for (unsigned  kdim = 0; kdim < dim; kdim++) { // velocity block row 
          
         for (unsigned i = 0; i < nDofsV; i++) {

	              double lap_res_du_u_kdim_i 			= 0.; 
		      double lap_res_du_ctrl_kdim_i 			= 0.;
		      double adv_res_uold_nablauold_kdim_i 		= 0.;
		      double adv_res_uold_nablauctrlold_kdim_i 	= 0.;
		      double adv_res_uctrlold_nablauold_kdim_i 	= 0.;
		      double adv_res_uctrlold_nablauctrlold_kdim_i 	= 0.;
              
	      for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		    lap_res_du_u_kdim_i  		     += gradSolVAR_qp[SolPdeIndex[kdim]][jdim]		            * phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad + jdim];
		    lap_res_du_ctrl_kdim_i 		 += gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim] * phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad + jdim];
          }
	      for (unsigned jdim = 0; jdim < dim; jdim++) {
		   adv_res_uold_nablauold_kdim_i 	     += SolVAR_qp[SolPdeIndex[jdim]] 		  * gradSolVAR_qp[SolPdeIndex[kdim]][jdim]		    * phi_gss_fe[ SolFEType[kdim] ][i];
		  adv_res_uold_nablauctrlold_kdim_i 	 += SolVAR_qp[SolPdeIndex[jdim]] 		  * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim] * phi_gss_fe[ SolFEType[kdim] ][i];
		  adv_res_uctrlold_nablauold_kdim_i	     += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]] * gradSolVAR_qp[SolPdeIndex[kdim]][jdim] 		    * phi_gss_fe[ SolFEType[kdim] ][i];
		  adv_res_uctrlold_nablauctrlold_kdim_i  += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]] * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim] * phi_gss_fe[ SolFEType[kdim] ][i];
	      }
	      
	      
	      Res[kdim][i]   +=   AbsDetJxWeight_iqp * (         
#if exact_sol_flag == 0
                                         + force[kdim] * phi_gss_fe[ SolFEType[kdim] ][i]
 #endif                                      
 #if exact_sol_flag == 1
                                       + exactForce[kdim] * phi_gss_fe[ SolFEType[kdim] ][i]
 #endif
                                           - IRe*lap_res_du_u_kdim_i 
                                           - IRe*lap_res_du_ctrl_kdim_i
                                           - advection_flag * adv_res_uold_nablauold_kdim_i 
                                           - advection_flag * adv_res_uold_nablauctrlold_kdim_i
                                           - advection_flag * adv_res_uctrlold_nablauold_kdim_i
					                       - advection_flag * adv_res_uctrlold_nablauctrlold_kdim_i
					 );
	}
      
}
	
//************ Residual, Velocity, END *********************

//************ Jacobian, Velocity, BEGIN *********************

//DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
              
  for (unsigned  kdim = 0; kdim < dim; kdim++) {
              
    for (unsigned i = 0; i < nDofsV; i++) {

	  for (unsigned j = 0; j < nDofsV; j++) {

    double lap_jac_du_u_kdim_i_j = 0.;
    double adv_uold_nablaunew_kdim_i_j = 0.;
    double adv_uctrlold_nablaunew_kdim_i_j = 0.;

        
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		    lap_jac_du_u_kdim_i_j += phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad /*space_dim*/ + jdim]	* 	phi_x_gss_fe[SolFEType[kdim]][ j * dim_offset_grad /*space_dim*/ + jdim ];
            }
          for (unsigned  jdim = 0; jdim < dim; jdim++) {
		    adv_uold_nablaunew_kdim_i_j     += SolVAR_qp[SolPdeIndex[jdim]]  * phi_x_gss_fe[ SolFEType[kdim] ][j * dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim] ][i];
		    adv_uctrlold_nablaunew_kdim_i_j += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]] * phi_x_gss_fe[ SolFEType[kdim] ][j * dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim] ][i];
            }  //jdim
                
              
        Jac[kdim][kdim][i * nDofsV + j] += (  IRe * lap_jac_du_u_kdim_i_j 
						    + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim] 		  * phi_gss_fe[ SolFEType[kdim] ][i]	 // c(u_hat_new,u_hat_old,delta_lambda) diagonal blocks  ..... unew_nablauold
						    + advection_flag * adv_uold_nablaunew_kdim_i_j 								 // c(u_hat_old, u_hat_new, delta_lambda)
						    + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][kdim] * phi_gss_fe[ SolFEType[kdim] ][i]	 // c(u_hat_new,u_0_old,delta_lambda) diagonal blocks  ..... unew_nablauctrlold
						    + advection_flag * adv_uctrlold_nablaunew_kdim_i_j 						  	 // c(u_0_old, u_hat_new, delta_lambda)
						    ) * AbsDetJxWeight_iqp; 
						    
               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
		Jac[kdim][off_kdim][i * nDofsV + j] += ( +	advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][off_kdim] 		    * phi_gss_fe[ SolFEType[kdim] ][i]	// c(u_hat_new,u_hat_old,delta_lambda) off-diagonal blocks  ..... unew_nablauold
						      + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][off_kdim] * phi_gss_fe[ SolFEType[kdim] ][i]	// c(u_hat_new,u_0_old,delta_lambda) off-diagonal blocks  ..... unew_nablauctrlold
						      ) * AbsDetJxWeight_iqp;
                              
	      
	} //j_du_u loop
	
   }//i_state

       
}
	      

//BLOCK delta_state - control------------------------------------------------------------------------------------
  for (unsigned  kdim = 0; kdim < dim; kdim++) {
              
   for (unsigned i = 0; i < nDofsV; i++) {
       
	for (unsigned j = 0; j < nDofsVctrl; j++) {
        
		       double  lap_jac_du_ctrl_kdim_i_j = 0.;
		       double  adv_uold_nablauctrlnew_kdim_i_j = 0.;
		       double  adv_uctrlold_nablauctrlnew_kdim_i_j = 0.;
              
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { 
		    lap_jac_du_ctrl_kdim_i_j += phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad /*space_dim*/ + jdim] * phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][j * dim_offset_grad /*space_dim*/ + jdim];
            }
            
            
		  for (unsigned  jdim = 0; jdim < dim; jdim++) {  //diagonal blocks only
		    adv_uold_nablauctrlnew_kdim_i_j	 += SolVAR_qp[SolPdeIndex[jdim]] 		  * phi_x_gss_fe[ SolFEType[kdim +  ctrl_pos_begin] ][j * dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim] ][i];
		    adv_uctrlold_nablauctrlnew_kdim_i_j += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]] * phi_x_gss_fe[ SolFEType[kdim +  ctrl_pos_begin] ][j * dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim] ][i];
                }  //jdim
	      


	      Jac[kdim][kdim + ctrl_pos_begin ][i*nDofsVctrl + j] += (+ IRe * lap_jac_du_ctrl_kdim_i_j 
									+ advection_flag * adv_uold_nablauctrlnew_kdim_i_j														 		 // c(u_hat_old, u_0_new, delta_lambda)
									+ advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim]		       * phi_gss_fe[ SolFEType[kdim] ][i]	 // c(u_0_new,u_hat_old,delta_lambda) diagonal blocks  ..... uctrlnew_nablauold
									+ advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][j] * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][kdim] * phi_gss_fe[ SolFEType[kdim] ][i]	 // c(u_0_new,u_0_old,delta_lambda) diagonal blocks  ..... uctrlnew_nablauctrlold
									+ advection_flag * adv_uctrlold_nablauctrlnew_kdim_i_j															 // c(u_0_old, u_0_new, delta_lambda)
									) * AbsDetJxWeight_iqp;
									
               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
		Jac[kdim][off_kdim + ctrl_pos_begin][i*nDofsVctrl + j] += (   +	advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim + ctrl_pos_begin] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][off_kdim]		     * phi_gss_fe[ SolFEType[kdim] ][i]	 // c(u_0_new,u_hat_old,delta_lambda) off-diagonal blocks  ..... uctrlnew_nablauold
									      + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim + ctrl_pos_begin] ][j] * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][off_kdim] * phi_gss_fe[ SolFEType[kdim] ][i]	 // c(u_0_new,u_0_old,delta_lambda) off-diagonal blocks  ..... uctrlnew_nablauctrlold
									  ) * AbsDetJxWeight_iqp;
	} //j_du_ctrl loop

   }//i_state loop

       
}
   
   
      
//************ Jacobian, Velocity, END *********************


//************ Residual & Jacobian, Pressure, BEGIN *********************
   
//BLOCK Pressure
   for (unsigned  kdim = 0; kdim < dim; kdim++) {
         
     for (unsigned i = 0; i < nDofsV; i++) {
       	      Res[kdim][i]   +=   AbsDetJxWeight_iqp *  SolVAR_qp[SolPdeIndex[press_type_pos]] * phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad /*space_dim*/ + kdim];
              
      for (unsigned j = 0; j < nDofsP; j++) {
         Jac[kdim][press_type_pos][i * nDofsP + j] += -( phi_gss_fe[SolFEType[press_type_pos]][j] * phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
         } //j_press loop
      
    }//i_state loop
      
  }

   
//DIV_state -----------
	    double div_u_du_qp = 0.;
      for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      div_u_du_qp += gradSolVAR_qp[SolPdeIndex[kdim]][kdim];
      }


  for (unsigned i = 0; i < nDofsP; i++) {  
      Res[press_type_pos][i]  +=  ( (div_u_du_qp) * phi_gss_fe[SolFEType[press_type_pos]][i] ) * AbsDetJxWeight_iqp;
   }
      
   for (unsigned  kdim = 0; kdim < dim; kdim++) {
     for (unsigned i = 0; i < nDofsP; i++) {
      for (unsigned j = 0; j < nDofsV; j++) {
	      Jac[press_type_pos][kdim][i * nDofsV + j] += - ( phi_gss_fe[SolFEType[press_type_pos]][i] * phi_x_gss_fe[SolFEType[kdim]][j * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
        } //j loop
	  }
	  
   }
//************ Residual & Jacobian, Pressure, END *********************


//============ delta_state row - END  ============================================================================================


    
//============ delta_adjoint row - BEGIN  =============================================================================================
  
  for (unsigned kdim = 0; kdim < dim; kdim++) { 
          
  for (unsigned i = 0; i < nDofsVadj; i++) {

		    double lap_res_dadj_adj_kdim_i 			= 0.;
		    double adv_res_phiadj_nablauold_uadjold_kdim_i 	= 0.;
		    double adv_res_uold_nablaphiadj_uadjold_kdim_i 	= 0.;
		    double adv_res_phiadj_nablauctrlold_uadjold_kdim_i = 0.;
		    double adv_res_uctrlold_nablaphiadj_uadjold_kdim_i = 0.;
            
        for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		   lap_res_dadj_adj_kdim_i    += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]][jdim]  *  phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim];
          }
          
       for (unsigned jdim = 0; jdim < dim; jdim++) {
		adv_res_phiadj_nablauold_uadjold_kdim_i     += phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i] * gradSolVAR_qp[SolPdeIndex[jdim]][kdim] 			* SolVAR_qp[SolPdeIndex[jdim + adj_pos_begin]];
		adv_res_uold_nablaphiadj_uadjold_kdim_i     += SolVAR_qp[SolPdeIndex[jdim]]		       * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]];
		adv_res_phiadj_nablauctrlold_uadjold_kdim_i += phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i] * gradSolVAR_qp[SolPdeIndex[jdim  + ctrl_pos_begin]][kdim]	* SolVAR_qp[SolPdeIndex[jdim + adj_pos_begin]];
		adv_res_uctrlold_nablaphiadj_uadjold_kdim_i += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]]  * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]];
	   }
	   
	   
	  Res[kdim + adj_pos_begin][i] += ( 
#if exact_sol_flag == 0
                            - cost_functional_coeff * target_flag * DesiredTargetVel()[kdim] 			      * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]
 #endif                                      
 #if exact_sol_flag == 1
                            - cost_functional_coeff * target_flag * exactVel_d[kdim] 			      * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]
 #endif
					    + cost_functional_coeff * target_flag * SolVAR_qp[SolPdeIndex[kdim]] 		      * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]
					    + cost_functional_coeff * target_flag * SolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]
					    - IRe * lap_res_dadj_adj_kdim_i
					    - advection_flag * adv_res_phiadj_nablauold_uadjold_kdim_i
					    - advection_flag * adv_res_uold_nablaphiadj_uadjold_kdim_i
					    - advection_flag * adv_res_phiadj_nablauctrlold_uadjold_kdim_i
					    - advection_flag * adv_res_uctrlold_nablaphiadj_uadjold_kdim_i
					    ) * AbsDetJxWeight_iqp;
      }
  }
  
//BLOCK delta_adjoint - state------------------------------------------------------------------------------------------
  for (unsigned kdim = 0; kdim < dim; kdim++) {
          
   for (unsigned i = 0; i < nDofsVadj; i++) {
 
     for (unsigned j = 0; j < nDofsV; j++) {
         
	      Jac[kdim + adj_pos_begin][kdim][i*nDofsV + j] += ( - cost_functional_coeff * target_flag * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i] * phi_gss_fe[SolFEType[kdim]][j] 
								 + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i]    * phi_x_gss_fe[ SolFEType[kdim] ][j*dim_offset_grad /*space_dim*/ + kdim] 		* SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]   //c(delta_u, u_hat_new, lambda_old)  diagonal blocks  ......phiadj_nablaunew_uadjold 
								 + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim] ][j] 			* phi_x_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	 //c(u_hat_new, delta_u, lambda_old) diagonal blocks  ......unew_nablaphiadj_uadjold
								     ) * AbsDetJxWeight_iqp;
                                     
               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
               
		Jac[kdim + adj_pos_begin][off_kdim][i * nDofsV + j] += (+ advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i] * phi_x_gss_fe[ SolFEType[off_kdim] ][j*dim_offset_grad /*space_dim*/ + kdim]		      * SolVAR_qp[SolPdeIndex[off_kdim + adj_pos_begin]]   //c(delta_u, u_hat_new, lambda_old)  off-diagonal blocks  ......phiadj_nablaunew_uadjold 
								      + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim] ][j] 		  * phi_x_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + off_kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	   //c(u_hat_new, delta_u, lambda_old) off-diagonal blocks  ......unew_nablaphiadj_uadjold
								      ) * AbsDetJxWeight_iqp;
	    
	  }//j_dadj_u loop
     }
   }
   
   
//BLOCK delta_adjoint - control-----------------------------------------------------------------------------------------
  for (unsigned kdim = 0; kdim < dim; kdim++) {
      
   for (unsigned i = 0; i < nDofsVadj; i++) {
       
     for (unsigned j = 0; j < nDofsVctrl; j++) {
         
	     Jac[kdim + adj_pos_begin][kdim + ctrl_pos_begin][i*nDofsVctrl + j] += ( - cost_functional_coeff * target_flag * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i] * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][j] 
										    + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i]  * phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][j*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]   //c(delta_u, u_0_new, lambda_old)  diagonal blocks  ......phiadj_nablauctrlnew_uadjold 
										    + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][j] * phi_x_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + kdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]   //c(u_0_new, delta_u, lambda_old) diagonal blocks  ......uctrlnew_nablaphiadj_uadjold
										    ) * AbsDetJxWeight_iqp;
                                            
               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
	     Jac[kdim + adj_pos_begin][off_kdim + ctrl_pos_begin][i*nDofsVctrl + j] += (+ advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i]      * phi_x_gss_fe[ SolFEType[off_kdim + ctrl_pos_begin] ][j*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[off_kdim + adj_pos_begin]]   //c(delta_u, u_0_new, lambda_old)  off-diagonal blocks  ......phiadj_nablauctrlnew_uadjold 
											+ advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim + ctrl_pos_begin] ][j] * phi_x_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + off_kdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	   //c(u_0_new, delta_u, lambda_old) off-diagonal blocks  ......uctrlnew_nablaphiadj_uadjold
											) * AbsDetJxWeight_iqp;
	    
	  }//j_dadj_ctrl loop
	  
     }
       
   }

//DIAG BLOCK delta_adjoint - adjoint---------------------------------------------------------------------------------
  for (unsigned  kdim = 0; kdim < dim; kdim++) {
              
   for (unsigned i = 0; i < nDofsVadj; i++) {
       
     for (unsigned j = 0; j < nDofsVadj; j++) {
         
		     double lap_jac_dadj_adj_kdim_i_j = 0.;
		     double adv_uold_nablaphiadj_uadjnew_kdim_i_j = 0.;
		     double adv_uctrlold_nablaphiadj_uadjnew_kdim_i_j = 0.;
            
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		  lap_jac_dadj_adj_kdim_i_j  +=  phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + jdim];
            }

            for (unsigned jdim = 0; jdim < dim; jdim++) { //diagonal blocks only
	     adv_uold_nablaphiadj_uadjnew_kdim_i_j    += SolVAR_qp[SolPdeIndex[jdim]] 		     * phi_x_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j] ;
	     adv_uctrlold_nablaphiadj_uadjnew_kdim_i_j  += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]] * phi_x_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j] ;
	   }
	   

	   Jac[kdim + adj_pos_begin][kdim + adj_pos_begin][i * nDofsVadj + j] += ( + IRe*lap_jac_dadj_adj_kdim_i_j 
										    + advection_flag * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim] 		  * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j]   //c(delta_u, u_hat_old, lambda_new)  diagonal blocks  ......phiadj_nablauold_uadjnew  
										    + advection_flag * adv_uold_nablaphiadj_uadjnew_kdim_i_j 	//c(u_hat_old, delta_uhat, lambda_new)
										    + advection_flag * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][kdim] * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j]   //c(delta_u, u_0_old, lambda_new)  diagonal blocks  ......phiadj_nablauctrlold_uadjnew  
										    + advection_flag * adv_uctrlold_nablaphiadj_uadjnew_kdim_i_j   //c(u_0_old, delta_uhat, lambda_new)
											) * AbsDetJxWeight_iqp;
                                            
               unsigned int off_kdim = (kdim+1) % dim; //off-diagonal blocks
		  Jac[kdim + adj_pos_begin][off_kdim + adj_pos_begin][i * nDofsVadj + j] += (+ advection_flag * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[off_kdim]][kdim]		    * phi_gss_fe[ SolFEType[off_kdim + adj_pos_begin] ][j]   //c(delta_u, u_hat_old, lambda_new)  off-diagonal blocks  ......phiadj_nablauold_uadjnew   
											   + advection_flag * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[off_kdim + ctrl_pos_begin]][kdim] * phi_gss_fe[ SolFEType[off_kdim + adj_pos_begin] ][j]   //c(delta_u, u_0_old, lambda_new)  off-diagonal blocks  ......phiadj_nablauctrlold_uadjnew  
											  ) * AbsDetJxWeight_iqp;
	    
      } //j_dadj_adj loop
      
  }//i_adj loop
  
 }

  
//************ Residual & Jacobian, Pressure Adj, BEGIN *********************
//BLOCK Pressure_adj
  
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
  for (unsigned i = 0; i < nDofsVadj; i++) {
      
      Res[kdim + adj_pos_begin][i] +=  AbsDetJxWeight_iqp * (SolVAR_qp[SolPdeIndex[press_type_pos + adj_pos_begin]] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim]);
      
    for (unsigned j = 0; j < nDofsPadj; j++) {
	      Jac[kdim + adj_pos_begin][press_type_pos + adj_pos_begin][i*nDofsPadj + j] += - (    phi_gss_fe[SolFEType[press_type_pos + adj_pos_begin]][j] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
           }//j_press_adj loop
        }//i_adj loop
	  }

//DIV_adj
		double div_adj_dadj_qp = 0.;
      for (unsigned kdim = 0; kdim < dim; kdim++) {
	    div_adj_dadj_qp += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin ]][kdim] ;
      }
      
  for (unsigned i = 0; i < nDofsPadj; i++) {
      Res[press_type_pos + adj_pos_begin][i] += ( (div_adj_dadj_qp) * phi_gss_fe[ SolFEType[press_type_pos + adj_pos_begin] ][i] ) * AbsDetJxWeight_iqp;
  }//i_div_adj
      
    for (unsigned kdim = 0; kdim < dim; kdim++) {
     for (unsigned i = 0; i < nDofsPadj; i++) {
      for (unsigned j = 0; j < nDofsVadj; j++) {
	    Jac[press_type_pos + adj_pos_begin][kdim + adj_pos_begin][i * nDofsVadj + j] += - (   phi_gss_fe[SolFEType[press_type_pos + adj_pos_begin]][i] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
          }//j loop
        }//i_div_adj
	   }
//************ Residual & Jacobian, Pressure Adj, END *********************

//============ delta_adjoint row - END  =============================================================================================


//============ delta_control row - BEGIN ==================================================================================================


//************ Residual, BEGIN *********************

      if ( control_el_flag == 1)    {
          
      for (unsigned kdim = 0; kdim < dim; kdim++) {
          
         for (unsigned i = 0; i < nDofsVctrl; i++) {
      
		    double lap_res_dctrl_ctrl_kdim_i 			 = 0.;
		    double lap_res_dctrl_adj_kdim_i 			 = 0.;
		    double adv_res_phictrl_nablauold_uadjold_kdim_i 	 = 0.;
		    double adv_res_uold_nablaphictrl_uadjold_kdim_i 	 = 0.;
		    double adv_res_phictrl_nablauctrlold_uadjold_kdim_i = 0.;
		    double adv_res_uctrlold_nablaphictrl_uadjold_kdim_i = 0.;
            
     for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		lap_res_dctrl_ctrl_kdim_i		      += gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim] 	*   phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim];
		lap_res_dctrl_adj_kdim_i 		      += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]][jdim] 	*   phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim];
         }
         
      for (unsigned jdim = 0; jdim < dim; jdim++) {
		adv_res_uold_nablaphictrl_uadjold_kdim_i     += SolVAR_qp[SolPdeIndex[jdim]]			 * phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]];
		adv_res_phictrl_nablauold_uadjold_kdim_i     += phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i] * gradSolVAR_qp[SolPdeIndex[jdim]][kdim]   			   * SolVAR_qp[SolPdeIndex[jdim + adj_pos_begin]];
		adv_res_phictrl_nablauctrlold_uadjold_kdim_i += phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i] * gradSolVAR_qp[SolPdeIndex[jdim  + ctrl_pos_begin]][kdim]	   * SolVAR_qp[SolPdeIndex[jdim + adj_pos_begin]];
		adv_res_uctrlold_nablaphictrl_uadjold_kdim_i += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]]   * phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]];
      }
      
      
      Res[kdim + ctrl_pos_begin][i] +=  AbsDetJxWeight_iqp * (
#if exact_sol_flag == 0
                     + cost_functional_coeff * target_flag * DesiredTargetVel()[kdim] * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i]
 #endif                                      
 #if exact_sol_flag == 1
                     + cost_functional_coeff * target_flag * exactVel_d[kdim] * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i]
 #endif
					 - cost_functional_coeff * target_flag * SolVAR_qp[SolPdeIndex[kdim]] * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i]
					 - cost_functional_coeff * target_flag * SolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i]
					 - alpha * SolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i]
					 - beta * lap_res_dctrl_ctrl_kdim_i
					 + IRe * lap_res_dctrl_adj_kdim_i
					+ advection_flag * adv_res_uold_nablaphictrl_uadjold_kdim_i
					+ advection_flag * adv_res_phictrl_nablauold_uadjold_kdim_i
					+ advection_flag * adv_res_phictrl_nablauctrlold_uadjold_kdim_i
					+ advection_flag * adv_res_uctrlold_nablaphictrl_uadjold_kdim_i				
					);
//       }
      
      
      }

          
    }
    
      }
      
    else if ( control_el_flag == 0) {
    
        for (unsigned kdim = 0; kdim < dim; kdim++) {
           for (unsigned i = 0; i < nDofsVctrl; i++) {

        Res[kdim + ctrl_pos_begin][i] +=       (- penalty_outside_control_domain)  *  (1 - control_node_flag[kdim][i]) * (Sol_eldofs_Mat[kdim + ctrl_pos_begin][i] - 0.);
                 }
            }
    }
  
//************ Residual, END *********************

//************ Jacobian, BEGIN *********************
      if ( control_el_flag == 1)    {
//BLOCK delta_control - state------------------------------------------------------------------------------------------------
   for (unsigned kdim = 0; kdim < dim; kdim++) {

     for (unsigned i = 0; i < nDofsVctrl; i++) {
    
      for (unsigned j = 0; j < nDofsV; j++) {
          
	      Jac[kdim + ctrl_pos_begin][kdim][i*nDofsV + j] +=  AbsDetJxWeight_iqp * (
                                 + cost_functional_coeff * target_flag * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i] * phi_gss_fe[SolFEType[kdim]][j] 
								 - advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim] ][j] 		       * phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	 //c(u_hat_new, delta_u0, lambda_old) diagonal blocks  ......unew_nablaphictrl_uadjold
								 - advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i]  * phi_x_gss_fe[ SolFEType[kdim] ][j*dim_offset_grad /*space_dim*/ + kdim] 			* SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]   //c(delta_u0, u_hat_new, lambda_old)  diagonal blocks  ......phictrl_nablaunew_uadjold 
								    );
          
               unsigned int off_kdim = (kdim+1) % dim; //off-diagonal blocks
               
	      Jac[kdim + ctrl_pos_begin][off_kdim][i * nDofsV + j] +=  AbsDetJxWeight_iqp * ( - advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim] ][j] 		   * phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i*dim_offset_grad /*space_dim*/ + off_kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	 //c(u_hat_new, delta_u0, lambda_old) off-diagonal blocks  ......unew_nablaphictrl_uadjold
								      - advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i] * phi_x_gss_fe[ SolFEType[off_kdim] ][j*dim_offset_grad /*space_dim*/ + kdim] 			* SolVAR_qp[SolPdeIndex[off_kdim + adj_pos_begin]]   //c(delta_u0, u_hat_new, lambda_old)  off-diagonal blocks  ......phictrl_nablaunew_uadjold 
								    );
	    
        }//j_dctrl_u loop
     }
  }
      
//BLOCK delta_control - adjoint------------------------------------------------------------------------------------------------
 for (unsigned  kdim = 0; kdim < dim; kdim++) {
              
   for (unsigned i = 0; i < nDofsVctrl; i++) {
      for (unsigned j = 0; j < nDofsVadj; j++) {
          
		    double  lap_jac_dctrl_adj_kdim_i_j = 0.;
		    double  adv_uold_nablaphictrl_uadjnew_kdim_i_j = 0.;
		    double  adv_uctrlold_nablaphictrl_uadjnew_kdim_i_j = 0.;
            
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { 
		lap_jac_dctrl_adj_kdim_i_j += phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]	*	phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + jdim];
            }

            for (unsigned jdim = 0; jdim < dim; jdim++) { //diagonal blocks only
	     adv_uold_nablaphictrl_uadjnew_kdim_i_j     += SolVAR_qp[SolPdeIndex[jdim]] 		     * phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i*dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j] ;
	     adv_uctrlold_nablaphictrl_uadjnew_kdim_i_j += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]]* phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i*dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j] ;
	   }
	   

	   Jac[kdim + ctrl_pos_begin][kdim + adj_pos_begin][i*nDofsVadj + j] += AbsDetJxWeight_iqp * ( - IRe * lap_jac_dctrl_adj_kdim_i_j 
										    - advection_flag * adv_uold_nablaphictrl_uadjnew_kdim_i_j 	//c(u_hat_old, delta_u0, lambda_new)
										    - advection_flag * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim] 		    * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j]   //c(delta_u0, u_hat_old, lambda_new)  diagonal blocks  ......phictrl_nablauold_uadjnew  
										    - advection_flag * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][kdim]  * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j]   //c(delta_u0, u_0_old, lambda_new)  diagonal blocks  ......phictrl_nablauctrlold_uadjnew  
										    - advection_flag * adv_uctrlold_nablaphictrl_uadjnew_kdim_i_j 	//c(u_0_old, delta_u0, lambda_new)
										      );
                                              
               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
               
	      Jac[kdim + ctrl_pos_begin][off_kdim + adj_pos_begin][i*nDofsVadj + j] += AbsDetJxWeight_iqp *  ( - advection_flag * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[off_kdim]][kdim]  		     * phi_gss_fe[ SolFEType[off_kdim + adj_pos_begin] ][j]   //c(delta_u0, u_hat_old, lambda_new)  off-diagonal blocks  ......phictrl_nablauold_uadjnew  
										         - advection_flag * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[off_kdim + ctrl_pos_begin]][kdim]  * phi_gss_fe[ SolFEType[off_kdim + adj_pos_begin] ][j]   //c(delta_u0, u_0_old, lambda_new)  off-diagonal blocks  ......phictrl_nablauctrlold_uadjnew  
										      );
	    
	     } //j_dctrl_adj loop
      }
 }
 
//DIAG BLOCK delta_control - control--------------------------------------------------------------------------------------
 for (unsigned  kdim = 0; kdim < dim; kdim++) {
    for (unsigned i = 0; i < nDofsVctrl; i++) {
      for (unsigned j = 0; j < nDofsVctrl; j++) {
          
            double  lap_jac_dctrl_ctrl_kdim_i_j = 0.;
              
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { 
		lap_jac_dctrl_ctrl_kdim_i_j += phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]	*	phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][j * dim_offset_grad /*space_dim*/ + jdim];
            }

            Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsVctrl + j] += (  + (cost_functional_coeff * target_flag + alpha) * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i] * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][j]
											+  beta * lap_jac_dctrl_ctrl_kdim_i_j 
											- advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i]    * phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][j*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]   //c(delta_u0, u_0_new, lambda_old)  diagonal blocks  ......phictrl_nablauctrlnew_uadjold 
											- advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][j] 	* phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	 //c(u_0_new, delta_u0, lambda_old) diagonal blocks  ......uctrlnew_nablaphictrl_uadjold
											) * AbsDetJxWeight_iqp;
                                            
               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
               
	      Jac[kdim + ctrl_pos_begin][off_kdim + ctrl_pos_begin][i*nDofsVctrl + j] += ( - advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i]     * phi_x_gss_fe[ SolFEType[off_kdim + ctrl_pos_begin] ][j*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[off_kdim + adj_pos_begin]]   //c(delta_u0, u_0_new, lambda_old)  off-diagonal blocks  ......phictrl_nablauctrlnew_uadjold 
											   - advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim + ctrl_pos_begin] ][j] * phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i*dim_offset_grad /*space_dim*/ + off_kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	 //c(u_0_new, delta_u0, lambda_old) off-diagonal blocks  ......uctrlnew_nablaphictrl_uadjold
											) * AbsDetJxWeight_iqp;
	    
      }//j_dctrl_ctrl loop

    }//i_ctrl loop

  }
    
   
  }
      
    else if ( control_el_flag == 0) {
        
    for (unsigned  kdim = 0; kdim < dim; kdim++) {
     for (unsigned i = 0; i < nDofsVctrl; i++) {
      for (unsigned j = 0; j < nDofsVctrl; j++) {
          if (i == j) {
          Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsVctrl + j]  +=  penalty_outside_control_domain * (1 - control_node_flag[kdim][i]);
          }
      }
    }
  }
        
    } 
 //************ Jacobian, END *********************
   
    
    
//BLOCK Pressure_ctrl
  for (unsigned kdim = 0; kdim < dim; kdim++) {
          
for (unsigned i = 0; i < nDofsVctrl; i++) {      
            Res[kdim + ctrl_pos_begin][i] +=  AbsDetJxWeight_iqp * ( SolVAR_qp[SolPdeIndex[press_type_pos + ctrl_pos_begin]] * phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim]);   //pressure

      for (unsigned j = 0; j < nDofsPctrl; j++) {
	      Jac[kdim + ctrl_pos_begin][press_type_pos + ctrl_pos_begin][i*nDofsPctrl + j] += -( phi_gss_fe[SolFEType[press_type_pos + ctrl_pos_begin]][j] * phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
      }//j_press_ctrl
      
  }//i_ctrl loop
  
 }

  
//DIV_ctrl
		  double div_ctrl_dctrl_qp = 0.;
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
		div_ctrl_dctrl_qp += gradSolVAR_qp[ SolPdeIndex[ctrl_pos_begin + kdim] ][kdim] ;
	  }
	  
  for (unsigned i = 0; i < nDofsPctrl; i++) {
	  Res[press_type_pos + ctrl_pos_begin][i] += ( (div_ctrl_dctrl_qp) * phi_gss_fe[ SolFEType[ctrl_pos_begin + press_type_pos ] ][i] ) * AbsDetJxWeight_iqp;
  }//i_div_ctrl
  
      
 for (unsigned kdim = 0; kdim < dim; kdim++) {
   for (unsigned i = 0; i < nDofsPctrl; i++) {
	  for (unsigned j = 0; j < nDofsVctrl; j++) {
		Jac[press_type_pos + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsVctrl + j] += - ( phi_gss_fe[SolFEType[press_type_pos + ctrl_pos_begin]][i] * phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][j * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
	  }//j loop
    }//i_div_ctrl
 }
 
//============ delta_control row - END ==================================================================================================
 
//============ delta_mu row - BEGIN  ============================================================================================
  //MU
//************ Residual, BEGIN *********************
  for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
          
  for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim]; i++) {
      
       Res[mu_pos_begin + kdim][i]  +=  (- penalty_outside_control_domain) *  (1 - control_node_flag[kdim][i]) * (Sol_eldofs_Mat[mu_pos_begin + kdim][i] - 0.);
      
     }
  }
//************ Residual, END *********************

  //MU
//************ Jacobian, BEGIN *********************
  for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
    for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim]; i++) {
      for (unsigned j = 0; j < Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim]; j++) {
            if (i == j) {
               Jac[mu_pos_begin + kdim][mu_pos_begin + kdim][i * Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim] + j]  +=  penalty_outside_control_domain * (1 - control_node_flag[kdim][i]);
            }
         }
      }
  }
//************ Jacobian, END *********************


//============ delta_mu row - END  ============================================================================================
 
 
 
      }  // end quadrature point loop
 

      //***************************************************************************************************************
      

    //Sum the local matrices/vectors into the Global Matrix/Vector
    for(unsigned i_unk = 0; i_unk < n_unknowns; i_unk++) {
      RES->add_vector_blocked(Res[SolPdeIndex[i_unk]],L2G_dofmap_Mat[i_unk]);
        for(unsigned j_unk=0; j_unk < n_unknowns; j_unk++) {
	  if(assembleMatrix) JAC->add_matrix_blocked( Jac[ SolPdeIndex[i_unk] ][ SolPdeIndex[j_unk] ], L2G_dofmap_Mat[i_unk], L2G_dofmap_Mat[j_unk]);
        }
    }
 
   //--------------------------------------------------------------------------------------------------------  
  } //end list of elements loop for each subdomain
  
  
  //MU in res ctrl - BEGIN  ***********************************
ctrl_inequality::add_one_times_mu_res_ctrl(iproc,
                               ineq_flag,
                               ctrl_index_in_mat,
                               mu_index_in_mat,
                               SolIndex,
                               sol,
                               mlPdeSys,
                               pdeSys,
                               RES);
  //MU in res ctrl - END ***********************************
    
    
    
     
RES->close();
if (assembleMatrix) JAC->close();  /// This is needed for the parallel, when splitting the add part from the insert part!!!
      // ***************** ADD PART - END  *******************
 
  

//   ***************** INSERT PART - BEGIN (must go AFTER the sum, clearly) *******************
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

// -------
   geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
      
   geom_element_iel.set_elem_center_3d(iel, solType_coords);
// -------
   
// -------
    el_dofs_unknowns_vol(sol, msh, pdeSys, iel,
                        SolFEType,
                        SolIndex,
                        SolPdeIndex,
                        Sol_n_el_dofs_Mat_vol,
                        Sol_eldofs_Mat,
                        L2G_dofmap_Mat);
      
    //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = ControlDomainFlag_internal_restriction(geom_element_iel.get_elem_center_3d());
 
    
    if (control_el_flag == 1) {

        
  ctrl_inequality::update_active_set_flag_for_current_nonlinear_iteration
  (msh,
   sol,
   iel,
   geom_element_iel.get_coords_at_dofs/*_3d*/(),
   Sol_eldofs_Mat,
   Sol_n_el_dofs_Mat_vol,
   c_compl,
   mu_index_in_mat,
   ctrl_index_in_mat,
   solIndex_act_flag_sol,
   ctrl_lower,
   ctrl_upper,
   sol_actflag);
  
      


    ctrl_inequality::node_insertion(iel,
                   msh,
                   L2G_dofmap_Mat,
                   mu_index_in_mat,
                   ctrl_index_in_mat,
                    Sol_eldofs_Mat,
                   Sol_n_el_dofs_Mat_vol,
                   sol_actflag,
                   ctrl_lower,
                   ctrl_upper,
                   ineq_flag,
                   c_compl,
                   JAC,
                   RES,
                   assembleMatrix
                   );
   

       }

 
    
    
 //============= delta_ctrl-delta_mu row ===============================
  if (assembleMatrix) { 
             for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
JAC->matrix_set_off_diagonal_values_blocked(L2G_dofmap_Mat[ ctrl_index_in_mat[kdim] ], L2G_dofmap_Mat[ mu_index_in_mat[kdim] ], ineq_flag * 1.);
             }
            }
  
  
  }
//   ***************** INSERT PART - END (must go AFTER the sum, clearly) *******************
  
  
  
 RES->close();
 JAC->close();
  // ***************** END ASSEMBLY *******************

    
  
  print_global_residual_jacobian(print_algebra_global,
                                 ml_prob,
                                 mlPdeSys,
                                 pdeSys,
                                 RES,
                                 JAC,
                                 iproc,
                                 assembleMatrix);
  


}
 

 
 double*  GetErrorNorm(const MultiLevelProblem& ml_prob, MultiLevelSolution* ml_sol, Solution* sol_coarser_prolongated) {
  
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
  const int n_vars = dim+1;
  const int n_unknowns = 3*n_vars; //(2.*dim)+1; //state , adjoint of velocity terms and one pressure term
  const int vel_type_pos = 0;
  const int press_type_pos = dim;
  const int state_pos_begin = 0;
  const int adj_pos_begin   = dim+1;
  const int ctrl_pos_begin   = 2*(dim+1);
  
  vector < std::string > Solname(n_unknowns);  // const char Solname[4][8] = {"u_0","u_1","u_2","p_u"};
  Solname              [state_pos_begin+0] =                "u_0";
  Solname              [state_pos_begin+1] =                "u_1";
  if (dim == 3) Solname[state_pos_begin+2] =                "u_2";
  Solname              [state_pos_begin + press_type_pos] = "p_u";
  
  Solname              [adj_pos_begin + 0] =              "adj_0";
  Solname              [adj_pos_begin + 1] =              "adj_1";
  if (dim == 3) Solname[adj_pos_begin + 2] =              "adj_2";
  Solname              [adj_pos_begin + press_type_pos] = "p_adj";

  Solname              [ctrl_pos_begin + 0] =              "ctrl_0";
  Solname              [ctrl_pos_begin + 1] =              "ctrl_1";
  if (dim == 3) Solname[ctrl_pos_begin + 2] =              "ctrl_2";
  Solname              [ctrl_pos_begin + press_type_pos] = "p_ctrl";
  
  vector < unsigned > SolIndex(n_unknowns);  
  vector < unsigned > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolIndex[ivar]	= ml_sol->GetIndex        (Solname[ivar].c_str());
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);
  }

  vector < double > Sol_n_el_dofs_Mat_vol(n_unknowns);
  
  //==========================================================================================
  // velocity ************************************
  vector < vector < double > > phi_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_x_gss_fe(NFE_FAMS);
 
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_gss_fe[fe].reserve(max_size);
      phi_x_gss_fe[fe].reserve(max_size*dim);
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
    unsigned nDofsV = msh->GetElementDofNumber(iel, SolFEType[vel_type_pos]);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, SolFEType[state_pos_begin + press_type_pos]);    // number of solution element dofs
    
    unsigned nDofsVadj = msh->GetElementDofNumber(iel,SolFEType[adj_pos_begin]);    // number of solution element dofs
    unsigned nDofsPadj = msh->GetElementDofNumber(iel,SolFEType[adj_pos_begin + press_type_pos]);    // number of solution element dofs

    unsigned nDofsVctrl = msh->GetElementDofNumber(iel,SolFEType[ctrl_pos_begin]);    // number of solution element dofs
    unsigned nDofsPctrl = msh->GetElementDofNumber(iel,SolFEType[ctrl_pos_begin + press_type_pos] );    // number of solution element dofs

    unsigned nDofsVP = dim * nDofsV + nDofsP;
    unsigned nDofsVP_tot = 3*nDofsVP;
  // equation end *****************************


   //STATE###################################################################  
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
	Sol_n_el_dofs_Mat_vol[k]=ndofs_unk;
       Sol_eldofs_Mat[k].resize(ndofs_unk);
       SolVAR_coarser_prol_eldofs[k].resize(ndofs_unk);
    for (unsigned i = 0; i < ndofs_unk; i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);    // global to global mapping between solution node and solution dof // via local to global solution node
       Sol_eldofs_Mat[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);      // global extraction and local storage for the solution
       SolVAR_coarser_prol_eldofs[k][i] = (*sol_coarser_prolongated->_Sol[SolIndex[k]])(solDof);      // global extraction and local storage for the solution
      }
    }
  //CTRL###################################################################

 
      // ********************** Gauss point loop *******************************
      for(unsigned iqp = 0;iqp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); iqp++) {
	
 
      for(int fe=0; fe < NFE_FAMS; fe++) {
	msh->_finiteElement[ielGeom][fe]->Jacobian(coordX,iqp, AbsDetJxWeight_iqp,phi_gss_fe[fe],phi_x_gss_fe[fe], boost::none);
      }
         //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
  	msh->_finiteElement[ielGeom][BIQUADR_FE]->Jacobian(coordX,iqp,AbsDetJxWeight_iqp,phi_gss_fe[BIQUADR_FE],phi_x_gss_fe[BIQUADR_FE], boost::none);

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
	  for(unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[unk]; i++) {
	    SolVAR_qp[unk] += phi_gss_fe[ SolFEType[unk] ][i] * Sol_eldofs_Mat[unk][i];
	    SolVAR_coarser_prol_qp[unk] += phi_gss_fe[ SolFEType[unk] ][i] * SolVAR_coarser_prol_eldofs[unk][i];
//         std::cout << SolVAR_qp[unk] << " \t " << SolVAR_coarser_prol_qp[unk] << std::endl;
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      gradSolVAR_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType[unk] ][i*dim+ivar2] * Sol_eldofs_Mat[unk][i]; 
	      gradSolVAR_coarser_prol_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType[unk] ][i*dim+ivar2] * SolVAR_coarser_prol_eldofs[unk][i]; 
//         std::cout << gradSolVAR_qp[unk][ivar2] << " \t " << gradSolVAR_coarser_prol_qp[unk][ivar2] << std::endl;
	    }
	  }
	  
	}  
 //end unknowns eval at gauss points ********************************


	for(unsigned unk = 0; unk < n_unknowns; unk++) {
        l2norm[unk] += ( SolVAR_qp[unk] - SolVAR_coarser_prol_qp[unk] ) * ( SolVAR_qp[unk] - SolVAR_coarser_prol_qp[unk] ) * AbsDetJxWeight_iqp ; 
        
     }
    
    for(int  i = 0; i < dim; i++) {

        l2norm[n_unknowns + i] += ( ( SolVAR_qp[vel_type_pos + i] + SolVAR_qp[ctrl_pos_begin + i] ) - ( SolVAR_coarser_prol_qp[vel_type_pos + i]  + SolVAR_coarser_prol_qp[ctrl_pos_begin + i] ) ) * ( ( SolVAR_qp[vel_type_pos + i] + SolVAR_qp[ctrl_pos_begin + i] ) - ( SolVAR_coarser_prol_qp[vel_type_pos + i] + SolVAR_coarser_prol_qp[ctrl_pos_begin + i] ) )  * AbsDetJxWeight_iqp ;
    
    }

          
	for(unsigned unk = 0; unk < dim; unk++) {
        for(int j = 0; j < dim; j++){
        seminorm[unk] += (gradSolVAR_qp[unk][j] - gradSolVAR_coarser_prol_qp[unk][j] ) * ( gradSolVAR_qp[unk][j] - gradSolVAR_coarser_prol_qp[unk][j] ) * AbsDetJxWeight_iqp ;
        seminorm[unk + dim] += (gradSolVAR_qp[unk + adj_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + adj_pos_begin][j] ) * ( gradSolVAR_qp[unk + adj_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + adj_pos_begin][j] ) * AbsDetJxWeight_iqp ;
        seminorm[unk + 2*dim] += (gradSolVAR_qp[unk + ctrl_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j] ) * ( gradSolVAR_qp[unk + ctrl_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j] ) * AbsDetJxWeight_iqp ;
        seminorm[unk + 3*dim] += ((gradSolVAR_qp[unk][j]+gradSolVAR_qp[unk + ctrl_pos_begin][j]) - (gradSolVAR_coarser_prol_qp[unk][j]+gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j])) * ((gradSolVAR_qp[unk][j]+gradSolVAR_qp[unk + ctrl_pos_begin][j]) - (gradSolVAR_coarser_prol_qp[unk][j]+gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j])) * AbsDetJxWeight_iqp ;
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
