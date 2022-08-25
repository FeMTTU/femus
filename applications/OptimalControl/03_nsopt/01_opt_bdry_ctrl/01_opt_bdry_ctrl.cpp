#include "adept.h"

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NumericVector.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "Fluid.hpp"
#include "Parameter.hpp"
#include "Files.hpp"
#include "PetscMatrix.hpp"
#include "paral.hpp"//to get iproc HAVE_MPI is inside here

#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"

#include "ElemType.hpp"

using namespace femus;

#include   "../manufactured_solutions.hpp"

#define FACE_FOR_CONTROL          1        /* 1-2 x coords, 3-4 y coords, 5-6 z coords */

#define FACE_FOR_TARGET    1


#include   "../nsopt_params.hpp"

  const double cost_functional_coeff = 1.;
  const double alpha = ALPHA_CTRL_BDRY;
  const double beta  = BETA_CTRL_BDRY;

  #define RHS_ONE             0.

  
#define QRULE_I   0

//***** Implementation-related ****************** 
#define IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY   1  // 1 internal routine; 0 external routine
// you have to be careful with the nonlinear iterations, because sometimes the algorithm restarts and the nonlinear index is set back to zero!!!
// The first matrix before the boundary conditions and the first residual seem to be exactly the same...
// Maybe what is different is the matrix and residual AFTER the Boundary Conditions...  
  
  //Maybe the problem is the normal derivative of the adjoint..., all its vector components...
//**************************************

//****** Mesh ********************************
#define no_of_ref N_UNIFORM_LEVELS     //mesh refinements



//****** Convergence ********************************
#define exact_sol_flag      0  // 1 = if we want to use manufactured solution; 0 = if we use regular convention
#define compute_conv_flag   0  // 1 = if we want to compute the convergence and error ; 0 =  no error computation

#define NO_OF_L2_NORMS 9   //U,V,P,UADJ,VADJ,PADJ,GX,GY,THETA
#define NO_OF_H1_NORMS 6    //U,V,UADJ,VADJ,GX, GY
//**************************************



//**************************************
 double penalty_outside_control_boundary = 1.e50;       // penalty for zero control outside Gamma_c
 double penalty_ctrl = 1.e10;         //penalty for u = q
 double theta_value_outside_fake_element = 0.;
 //**************************************
 
 
 
//************** how to retrieve theta from proc0 ************************************* 
const double get_theta_value(const unsigned int nprocs, const Solution * sol, const unsigned int sol_theta_index) {
    
NumericVector* local_theta_vec;

          local_theta_vec = NumericVector::build().release();
        local_theta_vec->init(*sol->_Sol[sol_theta_index]);
        sol->_Sol[sol_theta_index]->localize(*local_theta_vec);
        
        PetscScalar* values;
        VecGetArray(static_cast<PetscVector*>(local_theta_vec)->vec(), &values);
        double theta_value = values[0];
        if (nprocs == 1) {
            if ( (*sol->_Sol[sol_theta_index])(0) != theta_value) abort();
        }
        
        return theta_value;
}
//*************************************************** 


 
double Solution_set_initial_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;

     if(!strcmp(name,"TargReg")) {
        value = ElementTargetFlag(x);
    }
    else if(!strcmp(name,"ContReg")) {
        value = ControlDomainFlag_bdry(x);
    }

    return value;
    
}



bool Solution_set_boundary_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char SolName[], double& value, const int faceName, const double time) {
  //1: bottom  //2: right  //3: top  //4: left
  
  bool dirichlet = true;
   value = 0.;
 

                if (!strcmp(SolName, "GX"))       {
//                     if (facename == FACE_FOR_CONTROL) dirichlet = false; 
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
                
                
                
           else if (!strcmp(SolName, "GY"))       { 
//                     if (facename == FACE_FOR_CONTROL) dirichlet = false; 
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
                
                
                
           else if (!strcmp(SolName, "GZ"))       { 
//                     if (facename == FACE_FOR_CONTROL) dirichlet = false; 
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
                
                
                
                
                else if (!strcmp(SolName, "THETA"))    { dirichlet = false; }

                
                
                else if (!strcmp(SolName, "W"))       { 
//                     if (facename == FACE_FOR_CONTROL) dirichlet = false; 
   if (faceName == FACE_FOR_CONTROL) {
     if (x[ axis_direction_Gamma_control(faceName) ] > GAMMA_CONTROL_LOWER - 1.e-5 && x[ axis_direction_Gamma_control(faceName) ] < GAMMA_CONTROL_UPPER + 1.e-5)  { 
         dirichlet = false;
    }
     else { 
         dirichlet = true;  
    }
  }                   
                }
                
                
     
#if exact_sol_flag == 0
                else if (!strcmp(SolName, "U"))       { if (faceName == FACE_FOR_CONTROL) dirichlet = false; }
                else if (!strcmp(SolName, "V"))       { if (faceName == FACE_FOR_CONTROL) dirichlet = false; }
#endif
     
#if exact_sol_flag == 1
  //b.c. for manufactured lid driven cavity
  double pi = acos(-1.);
                else if (!strcmp(SolName, "U"))       { if (faceName == FACE_FOR_CONTROL) value =   sin(pi* x[0]) * sin(pi* x[0]) * cos(pi* x[1]) - sin(pi* x[0]) * sin(pi* x[0]); }
                else if (!strcmp(SolName, "V"))       { if (faceName == FACE_FOR_CONTROL) value = - sin(2. * pi * x[0]) * sin(pi* x[1]) + pi * x[1] * sin(2. * pi * x[0]); }
 #endif
               
  return dirichlet;

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
  
                        feFamily.push_back(DISCONTINUOUS_POLYNOMIAL);
 
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
  
                        feOrder.push_back(ZERO);
 

  assert( feFamily.size() == feOrder.size() );
 
 std::vector< Unknown >  unknowns(feFamily.size());
 
  const int adj_pos_begin   = dimension + 1;
  const int ctrl_pos_begin  = 2 * (dimension + 1);

                                        unknowns[0]._name      = "U";
                                        unknowns[1]._name      = "V";
  if (dimension == 3)                   unknowns[2]._name      = "W";
                                unknowns[dimension]._name      = "P";
                        unknowns[adj_pos_begin + 0]._name      = "UADJ";
                        unknowns[adj_pos_begin + 1]._name      = "VADJ";
  if (dimension == 3)   unknowns[adj_pos_begin + 2]._name      = "WADJ";
                unknowns[adj_pos_begin + dimension]._name      = "PADJ";
                       unknowns[ctrl_pos_begin + 0]._name      = "GX";
                       unknowns[ctrl_pos_begin + 1]._name      = "GY";
  if (dimension == 3)  unknowns[ctrl_pos_begin + 2]._name      = "GZ";
  
               unknowns[ctrl_pos_begin + dimension]._name      = "THETA";

     for (unsigned int u = 0; u < unknowns.size(); u++) {
         
              unknowns[u]._fe_family  = feFamily[u];
              unknowns[u]._fe_order   = feOrder[u];
              unknowns[u]._time_order = 0;
              unknowns[u]._is_pde_unknown = true;
              unknowns[u]._is_sparse = true;
              
     }
 
 
 unknowns[ctrl_pos_begin + dimension]._is_sparse = false;
 
   return unknowns;
     
}



void AssembleNavierStokesOpt(MultiLevelProblem& ml_prob);    

void ComputeIntegral(const MultiLevelProblem& ml_prob);

double*  GetErrorNorm(const MultiLevelProblem& ml_prob, MultiLevelSolution* ml_sol, Solution* sol_coarser_prolongated);
// ||u_h - u_(h/2)||/||u_(h/2)-u_(h/4)|| = 2^alpha, alpha is order of conv 
//i.e. ||prol_(u_(i-1)) - u_(i)|| = err(i) => err(i-1)/err(i) = 2^alpha ,implemented as log(err(i)/err(i+1))/log2



int main(int argc, char** args) {
  
  // ======= Init ========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
    // ======= Files ========================
  const bool use_output_time_folder = false;
  const bool redirect_cout_to_file = false;
  Files files; 
        files.CheckIODirectories(use_output_time_folder);
        files.RedirectCout(redirect_cout_to_file);
  
    // ======= Quad Rule ========================
  std::vector< std::string > fe_quad_rule_vec;
  fe_quad_rule_vec.push_back("seventh");
//   fe_quad_rule_vec.push_back("eighth");
    
  // ======= Parameter  ==================
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

  
  // ======= Mesh  ==================
  MultiLevelMesh ml_mesh;
	
//   std::string input_file = "square_parametric.med";
  std::string input_file = "parametric_square_1x1.med";
//   std::string input_file = "square_4x5.med";
//     std::string input_file = "Mesh_3_groups.med";
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_file;
  const std::string infile = mystream.str();
  
 ml_mesh.ReadCoarseMesh(infile.c_str(),fe_quad_rule_vec[0].c_str(),Lref);


  // ======= Problem  ==================
  MultiLevelProblem ml_prob; 

  ml_prob.SetFilesHandler(&files);
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;
  ml_prob.SetQuadratureRuleAllGeomElemsMultiple(fe_quad_rule_vec);
  ml_prob.set_all_abstract_fe_multiple();
    
  // ======= Refinement  ==================
    
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

         for (int i = /*0*/maxNumberOfMeshes - 1; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

  unsigned numberOfUniformLevels = i + 1; 
  unsigned numberOfSelectiveLevels = 0;
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  ml_mesh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  ml_mesh.PrintInfo();

  MultiLevelSolution ml_sol(&ml_mesh);
  
  // ======= Problem  ==================
  ml_prob.SetMultiLevelMeshAndSolution(& ml_sol);
 
  
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
  ml_sol.AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO);
   ml_sol.Initialize("TargReg",     Solution_set_initial_conditions, & ml_prob);
   ml_sol.Initialize("ContReg",     Solution_set_initial_conditions, & ml_prob);
  // ======= Solutions that are not Unknowns - END  ==================

 
  // ======= System - BEGIN ========================
  NonLinearImplicitSystem& system_opt    = ml_prob.add_system < NonLinearImplicitSystem > ("NSOpt");

  for (unsigned int u = 0; u < unknowns.size(); u++)  { 
  system_opt.AddSolutionToSystemPDE(unknowns[u]._name.c_str());
  }
   
  
  // attach the assembling function to system
   system_opt.SetAssembleFunction(AssembleNavierStokesOpt);

   
  system_opt.init();
  set_dense_pattern_for_unknowns(system_opt, unknowns);
  // ======= System  - END ========================

   
  // initialize and solve the system
    system_opt.init();
    system_opt.ClearVariablesToBeSolved();
    system_opt.AddVariableToBeSolved("All");
  
    ml_sol.SetWriter(VTK);
    ml_sol.GetWriter()->SetDebugOutput(true);
    
    system_opt.SetDebugNonlinear(true);
    system_opt.SetDebugFunction(ComputeIntegral);
    
//   system_opt.SetMaxNumberOfNonLinearIterations(2);
  
//   system_opt.SetNonLinearConvergenceTolerance(1.e-15);
//     system_opt.SetAbsoluteLinearConvergenceTolerance(1.e-14);
    system_opt.SetOuterSolver(PREONLY);
    
    
   
    system_opt.MGsolve();
//   system_opt.assemble_call_before_boundary_conditions(1);

  
  
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
       
 
  // ======= Print ========================
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  ml_sol.GetWriter()->Write(files.GetOutputPath()/*DEFAULT_OUTPUTDIR*/,"biquadratic", variablesToBePrinted, i);
 
  //Destroy all the new systems
//   ml_prob.clear();
  }

//  delete ml_sol_all_levels; 

#if compute_conv_flag == 1
  std::cout << "=======================================================================" << std::endl;
   std::cout << " L2-NORM ERROR and ORDER OF CONVERGENCE:\n\n";
  std::vector< std::string > norm_names_L2 = {"U","V", "P", "UADJ","VADJ", "PADJ", "GX","GY", "THETA"};

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
  std::vector< std::string > norm_names_H1 = {"U","V", "UADJ","VADJ", "GX","GY"};

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
  
  return 0;
}
 




///@todo we have to decide if order the phi by fem families, by Unknowns or by Quantities
/// I say the best way is by Quantity, because they have all the need to have phi functions as Unknowns

        

void AssembleNavierStokesOpt(MultiLevelProblem& ml_prob) {
     
  NonLinearImplicitSystem * mlPdeSys  = & ml_prob.get_system< NonLinearImplicitSystem >("NSOpt");
 const unsigned level = mlPdeSys->GetLevelToAssemble();

  bool assembleMatrix = mlPdeSys->GetAssembleMatrix(); 
   
  Solution*	 sol  	         = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearEquationSolver*  pdeSys	 = mlPdeSys->_LinSolver[level];   
  const char* pdename            = mlPdeSys->name().c_str();
  
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;
  
  Mesh*		 msh    = ml_prob._ml_msh->GetLevel(level);
  elem*		 el	= msh->el;
  SparseMatrix*	 JAC	= pdeSys->_KK;
  NumericVector* RES 	= pdeSys->_RES;
  
  
   //***** @todo to avoid Petsc complaint *******************************
//  MatSetOption(static_cast< PetscMatrix* >(JAC)->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
   //************************************

    
  //data
  const unsigned dim 	= msh->GetDimension();
  unsigned dim2     = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned igrid	= msh->GetLevel();
  unsigned iproc 	= msh->processor_id();
  unsigned   nprocs = msh->n_processors();


  constexpr bool print_algebra_global = false;
  constexpr bool print_algebra_local = false;

  
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));


   //=============== Geometry ========================================
   unsigned coordXType = 2; /*BIQUADR_FE*/// get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
   unsigned solType_coords = coordXType;  //FE_DOMAIN = 0; //we do linear FE this time // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
  const unsigned int dim_offset_grad = /*dim*/  3  /*2*/    ;
 
  std::vector<double> normal(dim_offset_grad /*space_dim*/, 0.);
 
  // geometry *******************************************
  vector< vector < double> > coordX(dim);
  vector< vector < double> > coordX_bd(/*space_*/dim);
  for(int i=0;i<dim;i++) {   
       coordX[i].reserve(max_size); 
       coordX_bd[i].reserve(max_size); 
  }
  // geometry *******************************************

// //  //***************************************************  
// //   
// //   vector < double > coord_at_qp_bdry(space_dim);
// //   
// //   vector <double> phi_coords;
// //   vector <double> phi_coords_x;
// //   vector <double> phi_coords_xx; 
// // 
// //   phi_coords.reserve(max_size);
// //   phi_coords_x.reserve(max_size * space_dim);
// //   phi_coords_xx.reserve(max_size * dim2);
// //   
// //   //boundary shape functions
// //   vector <double> phi_coords_bdry;  
// //   vector <double> phi_coords_x_bdry; 
// // 
// //   phi_coords_bdry.reserve(max_size);
// //   phi_coords_x_bdry.reserve(max_size * space_dim);
// //  //*************************************************** 
 
  // solution variables *******************************************
  
    std::vector<std::string> ctrl_name;
    ctrl_name.resize(dim);
    ctrl_name[0] = "GX";
    ctrl_name[1] = "GY";
     if (dim == 3)  ctrl_name[2] = "GZ";
 
  
  //=== Unknowns - BEGIN ========================================================
  const int n_vars_state = dim + 1;
  const int n_unknowns = 3 * n_vars_state; //(2.*dim)+1; //state , adjoint of velocity terms and one pressure term
  
  const int state_pos_begin  = 0;
  const int press_type_pos   = dim;
  const int adj_pos_begin    = press_type_pos + 1;
  const int ctrl_pos_begin   = 2 * (dim + 1);
  const int theta_index      = ctrl_pos_begin + dim;
  
  
  const int pos_mat_ctrl = ctrl_pos_begin;

  
  
  vector < std::string > Solname_Mat(n_unknowns);  // const char Solname_Mat[4][8] = {"U","V","W","P"};
  Solname_Mat              [state_pos_begin+0] =                "U";
  Solname_Mat              [state_pos_begin+1] =                "V";
  if (dim == 3) Solname_Mat[state_pos_begin+2] =                "W";
  Solname_Mat              [state_pos_begin + press_type_pos] = "P";

  Solname_Mat              [adj_pos_begin + 0] =              "UADJ";
  Solname_Mat              [adj_pos_begin + 1] =              "VADJ";
  if (dim == 3) Solname_Mat[adj_pos_begin + 2] =              "WADJ";
  Solname_Mat              [adj_pos_begin + press_type_pos] = "PADJ";

  Solname_Mat              [ctrl_pos_begin + 0] =                  ctrl_name[0];
  Solname_Mat              [ctrl_pos_begin + 1] =                  ctrl_name[1];
  if (dim == 3) Solname_Mat[ctrl_pos_begin + 2] =                  ctrl_name[2];
  Solname_Mat              [theta_index ]       =               "THETA";
  
  vector < unsigned > SolPdeIndex(n_unknowns);
  vector < unsigned > SolIndex_Mat(n_unknowns);  
  vector < unsigned > SolFEType_Mat(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolPdeIndex[ivar]	= mlPdeSys->GetSolPdeIndex(Solname_Mat[ivar].c_str());
    SolIndex_Mat[ivar]	= ml_sol->GetIndex        (Solname_Mat[ivar].c_str());
    SolFEType_Mat[ivar]	= ml_sol->GetSolutionType(SolIndex_Mat[ivar]);
  }

  vector < unsigned > Sol_n_el_dofs_Mat_vol(n_unknowns);
  //=== Unknowns - END ========================================================
  
  //=== Sol (quantities, not unknowns) - BEGIN ========================================================
  const int pos_sol_ctrl = ctrl_pos_begin;
  
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
  //=== Sol (quantities, not unknowns) - END ========================================================
  
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
   

  vector < vector < double > >  sol_adj_x_vol_at_bdry_gss(dim);
  for (int ldim =0; ldim < dim; ldim++) sol_adj_x_vol_at_bdry_gss[ldim].reserve(max_size);
  
  vector < double > grad_adj_dot_n_res;
  vector < double > grad_adj_dot_n_jac;
  grad_adj_dot_n_res.reserve(max_size);
  grad_adj_dot_n_jac.reserve(max_size);
  
 //********************* bdry cont *******************
 //*************************************************** 
  vector <double> phi_ctrl_bdry;  
  vector <double> phi_ctrl_x_bdry; 

  phi_ctrl_bdry.reserve(max_size);
  phi_ctrl_x_bdry.reserve(max_size * space_dim);
 //*************************************************** 

  //=================================================================================================
  
  
  // equation - BEGIN  ***********************************
  vector < vector < int > > L2G_dofmap_Mat(n_unknowns); 
  vector < vector < double > > Res(n_unknowns);
  vector < vector < vector < double > > > Jac(n_unknowns);
  
  vector < vector < double > > Jac_outer(dim);
  vector < double > Res_outer(1);

  
  for(int i = 0; i < n_unknowns; i++) {     
    L2G_dofmap_Mat[i].reserve(max_size);
       Res[i].reserve(max_size);
  }
   
  if (assembleMatrix) {
    
    for(int i = 0; i < n_unknowns; i++) {
      Jac[i].resize(n_unknowns);    
      for(int j = 0; j < n_unknowns; j++) {
	Jac[i][j].reserve(max_size*max_size);	
      }
    }

         for(int i = 0; i < dim; i++) {  Jac_outer[i].reserve(max_size); }
    
  }
  // equation - END  ***********************************

  
  //----------- dofs - BEGIN  ------------------------------
  vector < vector < double > > sol_eldofs_Mat(n_unknowns); //sol_V,P_of_st,adj,ctrl
  vector < vector < double > > gradsol_eldofs_Mat(n_unknowns);
  
  for(int k = 0; k < n_unknowns; k++) {
    sol_eldofs_Mat[k].reserve(max_size);
    gradsol_eldofs_Mat[k].reserve(max_size*dim);    
  }
  //----------- dofs - END  ------------------------------

  //------------ at quadrature points - BEGIN ---------------------
    vector < double > SolVAR_qp(n_unknowns);   //sol_V,P_gss_of_st,adj,ctrl_ie@quadraturepoints
    vector < vector < double > > gradSolVAR_qp(n_unknowns);
    for(int k = 0; k < n_unknowns; k++) {  gradSolVAR_qp[k].resize(dim_offset_grad /*space_dim*/);  }
      
      
 //*************************************************** 
     std::vector < std::vector < double > >  JacI_qp(space_dim);
     std::vector < std::vector < double > >  Jac_qp(dim);
    for (unsigned d = 0; d < Jac_qp.size(); d++) {   Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp.size(); d++) { JacI_qp[d].resize(dim); }
    
    double detJac_qp;

     std::vector < std::vector < double > >  JacI_iqp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_iqp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_iqp_bdry.size(); d++) {   Jac_iqp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp_bdry.size(); d++) { JacI_iqp_bdry[d].resize(dim-1); }
    
    double detJac_iqp_bdry;
    
    double weight_iqp = 0.;
    double weight_iqp_bdry = 0.;
    
      //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all;
  ml_prob.get_all_abstract_fe_multiple(elem_all);
//*************************************************** 
  //------------ at quadrature points - END ---------------------
 
  
//************** how to retrieve theta from proc0 ************************************* 
 double solTheta = get_theta_value(msh->n_processors(), sol, SolIndex_Mat[theta_index]);
//*************************************************** 

  const double IRe = ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
  
  
  // Set to zero all the global structures
    RES->zero();
    if(assembleMatrix) JAC->zero();

    
 
  constexpr unsigned qrule_i = QRULE_I;
  //----------------------

    const unsigned int n_components_ctrl = dim;
 
 
   if ( !(IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY) ) {
                      
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
                     //-----------
                    phi_ctrl_bdry,
                    phi_ctrl_x_bdry, ///@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ todo
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
                    qrule_i
                    ) ;
  
  }
 
 
 
 
 
 
 
 
 
  // ****************** element loop *******************
 
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

    unsigned nDofsGctrl = msh->GetElementDofNumber(iel,SolFEType_Mat[ctrl_pos_begin]);
    unsigned nDofsThetactrl = msh->GetElementDofNumber(iel,SolFEType_Mat[theta_index] );

    unsigned nDofsVP = dim * nDofsV + nDofsP; //same for state and adjoint
    unsigned nDofsVPctrl = dim * nDofsGctrl + nDofsThetactrl; //control
   
    unsigned nDofsVP_tot = 2 * nDofsVP + (nDofsVPctrl);
  // equation end *****************************
  
  //***** set target domain flag ********************************** 
   int target_flag = 0;
       target_flag = ElementTargetFlag(geom_element_iel.get_elem_center_3d());
   //***************************************       
   
 //************ set control flag *********************
    int does_iel_contain_a_bdry_control_face = 0;
        does_iel_contain_a_bdry_control_face = ControlDomainFlag_bdry(geom_element_iel.get_elem_center_3d());
 //************ initialize control node flag: for each Volume Elem, tell me if we have a Boundary Control dof *********************
    std::vector< std::vector<int> > control_node_flag_iel_jface(n_components_ctrl);
	    for(unsigned idim=0; idim < control_node_flag_iel_jface.size(); idim++) {
	          control_node_flag_iel_jface[idim].resize(nDofsGctrl);
    std::fill(control_node_flag_iel_jface[idim].begin(), control_node_flag_iel_jface[idim].end(), 0);
	    }
 //*************************************************** 
  
   //###################################################################  
    unsigned int fake_iel_flag = 0;
    unsigned int global_row_index_bdry_constr = pdeSys->KKoffset[SolPdeIndex[theta_index]][iproc];
  for (unsigned  k = 0; k < n_unknowns; k++) {
	unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType_Mat[k]);	//nDofs_V,P_of_st,adj,ctrl
	Sol_n_el_dofs_Mat_vol[k] = ndofs_unk;
	sol_eldofs_Mat[k].resize(ndofs_unk);	//sol_V,P_of_st,adj,ctrl
	L2G_dofmap_Mat[k].resize(ndofs_unk); 
    for (unsigned i = 0; i < ndofs_unk; i++) {
	      unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType_Mat[k]);    // global to global mapping between solution node and solution dof // via local to global solution node
	  sol_eldofs_Mat[k][i] = (*sol->_Sol[SolIndex_Mat[k]])(solDof);      // global extraction and local storage for the solution
	         L2G_dofmap_Mat[k][i] = pdeSys->GetSystemDof(SolIndex_Mat[k], SolPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
 
    if (k == SolPdeIndex[theta_index] && L2G_dofmap_Mat[k][i] == global_row_index_bdry_constr) {       fake_iel_flag = iel;  }
    }
  }
  //###################################################################
    
   //************ set fake theta flag: this flag tells me what degrees of freedom of the current element are fake for the theta variable *********************
    std::vector<int>  bdry_int_constr_pos_vec(1, global_row_index_bdry_constr); /*KKoffset[SolPdeIndex[PADJ]][iproc]*/
    std::vector<int> fake_theta_flag(nDofsThetactrl, 0);
    for (unsigned i = 0; i < nDofsThetactrl; i++) {
      if ( L2G_dofmap_Mat[ SolPdeIndex[theta_index] ] [i] == bdry_int_constr_pos_vec[0]) { 	fake_theta_flag[i] = 1;       }
    }
    
   
 //************ end set fake theta flag *********************

// setting Jac and Res to zero  - BEGIN ******************************* 
    for(int ivar=0; ivar<n_unknowns; ivar++) {
              Res[SolPdeIndex[ivar]].resize(Sol_n_el_dofs_Mat_vol[ivar]);
     std::fill(Res[SolPdeIndex[ivar]].begin(), Res[SolPdeIndex[ivar]].end(), 0.);         
//       memset(&Res[SolPdeIndex[ivar]][0],0.,Sol_n_el_dofs_Mat_vol[ivar]*sizeof(double));
    }
   
    for(int ivar=0; ivar<n_unknowns; ivar++) {
      for(int jvar=0; jvar<n_unknowns; jvar++) {
	    if(assembleMatrix){  //MISMATCH
                  Jac[ SolPdeIndex[ivar] ] [SolPdeIndex[jvar] ].resize(Sol_n_el_dofs_Mat_vol[ivar]*Sol_n_el_dofs_Mat_vol[jvar]);
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
      vector<double> normal(dim_offset_grad, 0);
	       
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

// 	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"U",tau,face,0.) && tau!=0.){
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
	elem_all[qrule_i][ielGeom_bd][solType_coords]->compute_normal(Jac_iqp_bdry, normal);
    
    weight_iqp_bdry = detJac_iqp_bdry * ml_prob.GetQuadratureRule(ielGeom_bd).GetGaussWeightsPointer()[iqp_bdry];

    elem_all[qrule_i][ielGeom_bd][SolFEType_Mat[theta_index]] ->shape_funcs_current_elem(iqp_bdry, JacI_iqp_bdry,phi_bd_gss_fe[SolFEType_Mat[theta_index]],phi_x_bd_gss_fe[SolFEType_Mat[theta_index]], boost::none, space_dim);
    
    elem_all[qrule_i][ielGeom_bd][SolFEType_Mat[ctrl_pos_begin]] ->shape_funcs_current_elem(iqp_bdry, JacI_iqp_bdry,phi_bd_gss_fe[SolFEType_Mat[ctrl_pos_begin]],phi_x_bd_gss_fe[SolFEType_Mat[ctrl_pos_begin]], boost::none , space_dim);


    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv_vol_at_bdry_new(geom_element_iel.get_coords_at_dofs_3d(), iqp_bdry, jface, Jac_qp/*not_needed_here*/, JacI_qp, detJac_qp/*not_needed_here*/, space_dim);
    elem_all[qrule_i][ielGeom][SolFEType_Mat[adj_pos_begin]]->shape_funcs_vol_at_bdry_current_elem(iqp_bdry, jface, JacI_qp, phi_vol_at_bdry_fe[SolFEType_Mat[adj_pos_begin]],phi_x_vol_at_bdry_fe[SolFEType_Mat[adj_pos_begin]], boost::none, space_dim);
     

//========== compute gauss quantities on the boundary - BEGIN ===============================================

//=============== Boundary control - BEGIN ========================================= 
		    for (unsigned  kdim = 0; kdim < n_components_ctrl; kdim++) {
                
			    unsigned int ctrl_index = kdim + ctrl_pos_begin;
                
					  SolVAR_bd_qp[ SolPdeIndex[ctrl_index] ] = 0.;
			  for(unsigned ivar2 = 0; ivar2 < dim_offset_grad; ivar2++) {  gradSolVAR_bd_qp[ SolPdeIndex[ctrl_index] ][ivar2] = 0.; }
	  
	         const unsigned ndof_bdry = msh->GetElementFaceDofNumber(iel, jface, SolFEType_Mat[ctrl_index]);

			  for(int i_bd = 0; i_bd < ndof_bdry; i_bd++) {
		                  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bd);
                          
                       SolVAR_bd_qp[SolPdeIndex[ctrl_index]]            += phi_bd_gss_fe  [ SolFEType_Mat[ctrl_index] ][i_bd]                  * sol_eldofs_Mat[ SolPdeIndex[ ctrl_index ] ][i_vol];

			      for(unsigned ivar2 = 0; ivar2 < dim_offset_grad; ivar2++) { 
                      gradSolVAR_bd_qp[SolPdeIndex[ctrl_index]][ivar2]  += phi_x_bd_gss_fe[ SolFEType_Mat[ctrl_index] ][i_bd * dim_offset_grad + ivar2] * sol_eldofs_Mat[ SolPdeIndex[ ctrl_index ] ][i_vol];
                        }
                        
			  }
			  
		    }
//=============== Boundary control - END ========================================= 
 
//=============== grad dot n - BEGIN ========================================= 
//     compute gauss quantities on the boundary through VOLUME interpolation
		for(unsigned ldim = 0; ldim < dim; ldim++) {   sol_adj_x_vol_at_bdry_gss[ldim].resize(dim_offset_grad); }
		grad_adj_dot_n_res.resize(dim);
		grad_adj_dot_n_jac.resize(dim);
		for(unsigned ldim=0; ldim<dim; ldim++) {   std::fill(sol_adj_x_vol_at_bdry_gss[ldim].begin(), sol_adj_x_vol_at_bdry_gss[ldim].end(), 0.);  }
		
		for(unsigned ldim=0; ldim<dim; ldim++) {  
		  for (int iv = 0; iv < nDofsVadj; iv++)  {
                     for (int d = 0; d < dim_offset_grad /*space_dim*/; d++) {
			   sol_adj_x_vol_at_bdry_gss[ldim][d] += sol_eldofs_Mat[SolPdeIndex[ldim + adj_pos_begin]][iv] * phi_x_vol_at_bdry_fe[SolFEType_Mat[ldim + adj_pos_begin]][iv * dim_offset_grad /*space_dim*/ + d];//notice that the convention of the orders x y z is different from vol to bdry
                    }
		  }  
		      
		  grad_adj_dot_n_res[ldim] = 0.;
		  for(unsigned d=0; d < dim_offset_grad /*space_dim*/; d++) {
		      grad_adj_dot_n_res[ldim] += sol_adj_x_vol_at_bdry_gss[ldim][d] * normal[d];  
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
        

        const unsigned nve_bdry_max = msh->GetElementFaceDofNumber(iel, jface, solType_coords);
        
  // *** phi_i loop ***
            for(unsigned i_bdry = 0; i_bdry < nve_bdry_max; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);

//============ Boundary-Boundary and Boundary-Volume Residuals - BEGIN ============================================================================================
		  
		      for (unsigned  kdim = 0; kdim < dim; kdim++) {
			
			    double lap_res_dctrl_ctrl_bd = 0.;
			      for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
				  lap_res_dctrl_ctrl_bd += gradSolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim] *
				                           phi_x_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin]][i_bdry * dim_offset_grad + jdim /*i_bdry + jdim*nve_bdry_max*/];
			      }//jdim

			
/*delta_state row */	   if (i_vol < nDofsV)      Res[kdim]                 [i_vol]  += - control_node_flag_iel_jface[kdim][i_vol] * penalty_ctrl * (sol_eldofs_Mat[SolPdeIndex[kdim + state_pos_begin]][i_vol] - sol_eldofs_Mat[SolPdeIndex[kdim + ctrl_pos_begin]][i_vol]);	    //u-g
/*delta_adjoint row */     if (i_vol < nDofsVadj)   Res[kdim + adj_pos_begin] [i_vol]  += 0.;	   
/*delta_control row */     if (i_vol < nDofsGctrl)  Res[kdim + ctrl_pos_begin][i_vol]  += - control_node_flag_iel_jface[kdim][i_vol] * weight_iqp_bdry * (
                                                                                          IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY * alpha * SolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * phi_bd_gss_fe[SolFEType_Mat[kdim +  ctrl_pos_begin]][i_bdry]
                                                                                        + IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY * beta * lap_res_dctrl_ctrl_bd
                                                                                        - IRe * grad_adj_dot_n_res[kdim]  * phi_bd_gss_fe[SolFEType_Mat[kdim +  ctrl_pos_begin]][i_bdry]
                                                                                        - /*(*sol->_Sol[SolIndex_Mat[theta_index]])(0)*/solTheta * phi_bd_gss_fe[SolFEType_Mat[kdim +  ctrl_pos_begin]][i_bdry] * normal[kdim]      //*sol->_Sol[SolIndex_Mat[theta_index]])(0) finds the global value from KKDof pos(72, 169,etc), SolVAReldof_theta gets the value in the boundary point which will be zero. Theta is just a const
                                                                                        );	    
		      }//kdim  

//============ Boundary-Boundary and Boundary-Volume Residuals - END ==================================================================================================

//============ Boundary-Boundary Jacobian i-BDRY/j-BDRY - BEGIN  ==================================================================================================
///@todo the error is in here
		      for(unsigned j_bdry = 0; j_bdry < nve_bdry_max; j_bdry ++) {
			  unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, jface, j_bdry);

				std::vector < double > lap_jac_dctrl_ctrl_bd(dim, 0.);
			     for (unsigned kdim = 0; kdim < dim; kdim++) {
				for (unsigned ldim = 0; ldim < dim_offset_grad /*space_dim*/; ldim++) {
                    lap_jac_dctrl_ctrl_bd[kdim] += phi_x_bd_gss_fe[SolFEType_Mat[ldim + ctrl_pos_begin]][i_bdry * dim_offset_grad + ldim] * 
                                                   phi_x_bd_gss_fe[SolFEType_Mat[ldim + ctrl_pos_begin]][j_bdry * dim_offset_grad + ldim];
                    }
                 }
                 
                 
			  for (unsigned  kdim = 0; kdim < dim; kdim++) {
                  
//DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
			    if(i_vol < nDofsV && j_vol < nDofsV && i_vol == j_vol)       		          Jac[kdim][kdim][i_vol * nDofsV + j_vol]	    += penalty_ctrl * control_node_flag_iel_jface[kdim][i_vol];  //u
			 
//BLOCK delta_state - control------------------------------------------------------------------------------------
			    if(i_vol < nDofsV && j_vol < nDofsGctrl && i_vol == j_vol) 	Jac[kdim][kdim + ctrl_pos_begin][i_vol * nDofsGctrl + j_vol]  += (-1.) * penalty_ctrl * control_node_flag_iel_jface[kdim][i_vol];  //-g

//DIAG BLOCK delta_control - control  --------------------------------------------------------------------------------------
			  if(i_vol < nDofsGctrl && j_vol < nDofsGctrl) {
				      Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i_vol * nDofsGctrl + j_vol] +=   control_node_flag_iel_jface[kdim][i_vol] * weight_iqp_bdry * (
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
			    grad_adj_dot_n_jac[ldim] = 0.;
			    for(unsigned d = 0; d < dim_offset_grad /*space_dim*/; d++) {
				grad_adj_dot_n_jac[ldim] += phi_x_vol_at_bdry_fe[SolFEType_Mat[ldim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + d] * normal[d];  //notice that the convention of the orders x y z is different from vol to bdry
			    }
		    }
			  //=============== grad dot n  =========================================    

			  for (unsigned kdim = 0; kdim < dim; kdim++) {
				Jac[kdim + ctrl_pos_begin][kdim + adj_pos_begin][i_vol*nDofsVadj + j] += control_node_flag_iel_jface[kdim][i_vol] * (-1.) * (weight_iqp_bdry  * phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin]][i_bdry] * IRe * grad_adj_dot_n_jac[kdim]);    		      
			  }
		} // end j loop for volume 
//============ Boundary-Volume Jacobian i-BDRY/j-VOL - END  ==================================================================================================

		    }//end i_bdry loop


		    
//============ Boundary Integral Constraint Residual - BEGIN ============================================================================================
	  for (unsigned  kdim = 0; kdim < dim; kdim++) {
// 		for(unsigned i=0; i < nDofsThetactrl; i ++) { avoid because it is an element dof
/*delta_theta row */ 	/* Res[theta_index][i]*/ Res_outer[0] +=  /*fake_theta_flag[i] **/ weight_iqp_bdry * SolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * normal[kdim] ;
// 		}  
	  }
		  
//============ Boundary Integral Constraint Residual - END ============================================================================================
		
//============ Boundary Integral Constraint Jacobian - BEGIN ============================================================================================
       		for(unsigned i_bdry = 0; i_bdry < nve_bdry_max; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
            
 		    for (unsigned  kdim = 0; kdim < dim; kdim++) { 
			  for(unsigned i =0; i < nDofsThetactrl; i ++) {
			    if(i_vol < nDofsGctrl) {
				double temp = weight_iqp_bdry * ( phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin]][i_bdry] * normal[kdim]);
//ROW_BLOCK delta_theta - control -- loop over i in the VOLUME (while j(/i_vol) is in the boundary) -------------------------------------------------------------------------------------------------------------
			      Jac[theta_index][ctrl_pos_begin + kdim][i*nDofsGctrl + i_vol]     += - temp; /*weight_iqp_bdry * ( phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin]][i_bdry] * normal[kdim])*/
//COLUMN_BLOCK delta_control - theta ---- loop over j in the VOLUME (while i(/i_vol) is in the boundary) ---------------------------------------------------------------------------------------------------
			      Jac[ctrl_pos_begin + kdim][theta_index][i_vol*nDofsThetactrl + i] += - control_node_flag_iel_jface[kdim][i_vol] /** phi_bd_gss_fe[SolFEType_Mat[theta_index]][i]*/*temp; /*weight_iqp_bdry * ( phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin]][i_bdry] * normal[kdim]);*/
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
    
 
 
 
// //======================= VolumeLoop with Integration (and Zero boundary control outside Gamma_c) - BEGIN =====================================================    
// 
// for(unsigned iqp = 0; iqp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); iqp++) {
// 	
// 	// *** get Jacobian and test function and test function derivatives ***
//        // *** get gauss point weight, test function and test function partial derivatives ***
//     elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), iqp, Jac_qp, JacI_qp, detJac_qp, space_dim);
//     weight_iqp = detJac_qp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[iqp];
//    
//      for(int fe=0; fe < NFE_FAMS; fe++) {
//     elem_all[qrule_i][ielGeom][fe]->shape_funcs_current_elem(iqp, JacI_qp,phi_gss_fe[fe],phi_x_gss_fe[fe], boost::none , space_dim);
//       }
//  
//  
// //geometry eval at gauss points - BEGIN ********************************
//  vector < double > coordX_gss(dim, 0.);
//  	for(unsigned k = 0; k <  dim; k++) {
// 	  for(unsigned i = 0; i < geom_element_iel.get_coords_at_dofs_3d()[k].size(); i++) { 
//          coordX_gss[k] += geom_element_iel.get_coords_at_dofs_3d()[k][i] * phi_gss_fe[ solType_coords ][i];
//       }
//     }
// //geometry eval at gauss points - END  ********************************
// 
// //begin unknowns eval at gauss points ********************************
// 	for(unsigned unk = 0; unk < n_unknowns; unk++) {
// 	    SolVAR_qp[unk] = 0.;
// 	    for(unsigned ivar2 = 0; ivar2 < dim_offset_grad /*space_dim*/; ivar2++) {
//            gradSolVAR_qp[unk][ivar2] = 0.; 
// 	    }
// 	  
// 	    for(unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[unk]; i++) {
//                         SolVAR_qp[unk] += phi_gss_fe[ SolFEType_Mat[unk] ][i]             * sol_eldofs_Mat[SolPdeIndex[unk]][i];
// 		for(unsigned ivar2 = 0; ivar2 < dim_offset_grad /*space_dim*/; ivar2++) {       
// 		    gradSolVAR_qp[unk][ivar2]  += phi_x_gss_fe[ SolFEType_Mat[unk] ][i * dim_offset_grad /*space_dim*/ + ivar2] * sol_eldofs_Mat[SolPdeIndex[unk]][i]; 
// 		}
// 	    }//ndofsunk
// 	  
// 	} //unk 
//  //end unknowns eval at gauss points ********************************
// 	
// //======= computation of RHS (force and desired velocity) using MMS - BEGIN =============================================== 
// //state values-------------------- //non-hom bdry
// vector <double>  exact_stateVel(dim);
//   mms_lid_driven::value_stateVel(coordX_gss, exact_stateVel);
// vector <double>  exact_lap_stateVel(dim);
//   mms_lid_driven::laplace_stateVel(coordX_gss, exact_lap_stateVel);
// vector < vector < double > > exact_grad_stateVel(dim);
// for (unsigned k = 0; k < dim; k++){ 
//     exact_grad_stateVel[k].resize(dim);
//     std::fill(exact_grad_stateVel[k].begin(), exact_grad_stateVel[k].end(), 0.);
// }
//   mms_lid_driven::gradient_stateVel(coordX_gss,exact_grad_stateVel);
// 
// //adjoint values--------------------//hom bdry
// vector <double>  exact_adjVel(dim);
//   mms_lid_driven::value_adjVel(coordX_gss, exact_adjVel);
// vector <double>  exact_lap_adjVel(dim);
//   mms_lid_driven::laplace_adjVel(coordX_gss, exact_lap_adjVel);
// vector < vector < double > > exact_grad_adjVel(dim);
// for (unsigned k = 0; k < dim; k++){ 
//     exact_grad_adjVel[k].resize(dim);
//     std::fill(exact_grad_adjVel[k].begin(), exact_grad_adjVel[k].end(), 0.);
// }
//   mms_lid_driven::gradient_adjVel(coordX_gss,exact_grad_adjVel);
// vector <double> exact_grad_adjPress(dim);
//   mms_lid_driven::gradient_adjPress(coordX_gss, exact_grad_adjPress);
// 
// //convection terms from delta_state-------------------------------------
// vector <double>  exact_conv_u_nabla_u(dim,0.);
// 
// for (unsigned k = 0; k < dim; k++){
//     for (unsigned i = 0; i < dim; i++){
//     exact_conv_u_nabla_u[k] += exact_grad_stateVel[k][i] * exact_stateVel[i] ; 
//     }
// }
// 
// //convection terms from delta_adjoint-------------------------
// vector <double>  exact_conv_u_nabla_uadj(dim,0.);
// vector <double>  exact_conv_nabla_uT_uadj(dim,0.);
// 
// for (unsigned k = 0; k < dim; k++){
//     for (unsigned i = 0; i < dim; i++){
//     exact_conv_u_nabla_uadj[k] += exact_grad_adjVel[k][i] * exact_stateVel[i] ; 
//     exact_conv_nabla_uT_uadj[k] += exact_grad_stateVel[i][k] * exact_adjVel[i];
//     }
// }
// 
// //force and desired velocity ---------------------------------------------
// vector <double> exactForce(dim,0.);
// vector <double> exactVel_d(dim,0.);
// for (unsigned k = 0; k < dim; k++){
//     exactForce[k] = - IRe * exact_lap_stateVel[k]  + advection_flag * exact_conv_u_nabla_u[k] + exact_grad_adjPress[k];
//     exactVel_d[k] =   exact_stateVel[k] + (1./cost_functional_coeff) * (IRe * exact_lap_adjVel[k] - exact_grad_adjPress[k]) 
//                     + (1./cost_functional_coeff) * advection_flag * (exact_conv_u_nabla_uadj[k] - exact_conv_nabla_uT_uadj[k]);
// }
// //======= computation of RHS (force and desired velocity) using MMS - END =============================================== 
//  
// //============ delta_state row - BEGIN ============================================================================================
// 
//   for (unsigned i = 0; i < nDofsV; i++) {
// // FIRST ROW
// 	for (unsigned  kdim = 0; kdim < dim; kdim++) { // velocity block row 
//                 double lap_res_du_u		= 0.; 
//                 double adv_res_uold_nablauold 	= 0.;
// 	      for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
// 		    lap_res_du_u 	       += gradSolVAR_qp[SolPdeIndex[kdim]][jdim] * phi_x_gss_fe[SolFEType_Mat[kdim]][i * dim_offset_grad /*space_dim*/ + jdim];
//           }
// 	      for (unsigned jdim = 0; jdim < dim; jdim++) {
// 		   adv_res_uold_nablauold  += SolVAR_qp[SolPdeIndex[jdim]]  * gradSolVAR_qp[SolPdeIndex[kdim]][jdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i];
// 	      }      
// 	      Res[kdim][i]   +=  (           
// #if exact_sol_flag == 0
//                                          + force[kdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i]
//  #endif                                      
//  #if exact_sol_flag == 1
//                                        + exactForce[kdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i]
//  #endif
//                                        - IRe*lap_res_du_u 
//                                        - advection_flag * adv_res_uold_nablauold 
//                                        + SolVAR_qp[SolPdeIndex[press_type_pos]] * phi_x_gss_fe[SolFEType_Mat[kdim]][i * dim_offset_grad /*space_dim*/ + kdim]) * weight_iqp; 
// 	}	    
// //DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
// 	for (unsigned j = 0; j < nDofsV; j++) {
// 		      vector < double > lap_jac_du_u(dim,0.);
// 		      vector < double > adv_uold_nablaunew(dim,0.);
// 	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
//             for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { 
// 		    lap_jac_du_u[kdim] += phi_x_gss_fe[SolFEType_Mat[kdim]][i * dim_offset_grad /*space_dim*/ + jdim]*phi_x_gss_fe[SolFEType_Mat[kdim]][j * dim_offset_grad/*space_dim*/ + jdim];
//             }
//           }
//         for (unsigned  kdim = 0; kdim < dim; kdim++) { 
//          for (unsigned  jdim = 0; jdim < dim; jdim++) { //diagonal blocks only
// 		    adv_uold_nablaunew[kdim] 	 += SolVAR_qp[SolPdeIndex[jdim]] * phi_x_gss_fe[ SolFEType_Mat[kdim] ][j * dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i];
//                 }  //jdim
// 	      }
// 	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
// 		Jac[kdim][kdim][i*nDofsV + j] += (   IRe * lap_jac_du_u[kdim] 
//                                             + advection_flag * adv_uold_nablaunew[kdim] 		// c(u_old, u_new, delta_lambda)
//                                             + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i]	 // c(u_new,u_old,delta_lambda) diagonal blocks  ..... unew_nablauold
//                                          ) * weight_iqp; 
//               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
// 		Jac[kdim][off_kdim][i*nDofsV + j] += (	advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[off_kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][off_kdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i]	// c(u_new,u_old,delta_lambda) off-diagonal blocks  ..... unew_nablauold
//                                              ) * weight_iqp;
// 	      }
// 	} //j_du_u loop
//      
// //BLOCK Pressure
//       for (unsigned j = 0; j < nDofsP; j++) {
// 	    for (unsigned  kdim = 0; kdim < dim; kdim++) {
// 	      Jac[kdim][press_type_pos][i*nDofsP + j] += -( phi_gss_fe[SolFEType_Mat[press_type_pos]][j] * phi_x_gss_fe[SolFEType_Mat[kdim]][i * dim_offset_grad /*space_dim*/ + kdim] ) * weight_iqp;
// 	    }
//       }//j_press loop
//    }//i_state loop
// 
// //DIV_state
//   for (unsigned i = 0; i < nDofsP; i++) {
// 		    double div_u_du_qp =0.;
//       for (unsigned  kdim = 0; kdim < dim; kdim++) {
// 	      div_u_du_qp += gradSolVAR_qp[SolPdeIndex[kdim]][kdim] ;
//       }
//       Res[press_type_pos][i]  +=  ( (div_u_du_qp) * phi_gss_fe[SolFEType_Mat[press_type_pos]][i] ) * weight_iqp;
//       for (unsigned j = 0; j < nDofsV; j++) {
// 	  for (unsigned  kdim = 0; kdim < dim; kdim++) {
// 	      Jac[press_type_pos][kdim][i*nDofsV + j] += -( phi_gss_fe[SolFEType_Mat[press_type_pos]][i] * phi_x_gss_fe[SolFEType_Mat[kdim]][j * dim_offset_grad /*space_dim*/ + kdim] ) * weight_iqp;
// 	  }
//       } //j loop
//    }//i_div_state
// //============ delta_state row - END ============================================================================================
// 
// 
//     
// //============ delta_adjoint row - BEGIN =============================================================================================
//   
//   for (unsigned i = 0; i < nDofsVadj; i++) {
// // SECOND ROW
//      for (unsigned kdim = 0; kdim < dim; kdim++) { 
// 		    double lap_res_dadj_adj 			= 0.;
// 		    double adv_res_phiadj_nablauold_uadjold 	= 0.;
// 		    double adv_res_uold_nablaphiadj_uadjold 	= 0.;
// 	   for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
// 		lap_res_dadj_adj 		             += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]][jdim]*phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim];
//        }
// 	   for (unsigned jdim = 0; jdim < dim; jdim++) {
// 		adv_res_phiadj_nablauold_uadjold     += phi_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i] * gradSolVAR_qp[SolPdeIndex[jdim]][kdim] 			* SolVAR_qp[SolPdeIndex[jdim + adj_pos_begin]];
// 		adv_res_uold_nablaphiadj_uadjold     += SolVAR_qp[SolPdeIndex[jdim]]		       * phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]];
// 	   }
// 	  Res[kdim + adj_pos_begin][i] += ( 
// #if exact_sol_flag == 0
//                             - cost_functional_coeff * target_flag * DesiredTargetVel()[kdim] 			      * phi_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i]
//  #endif                                      
//  #if exact_sol_flag == 1
//                             - cost_functional_coeff * target_flag * exactVel_d[kdim] 			      * phi_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i]
//  #endif
//                                         + cost_functional_coeff * target_flag * SolVAR_qp[SolPdeIndex[kdim]] * phi_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i]
//                                         - IRe*lap_res_dadj_adj
//                                         - advection_flag * adv_res_phiadj_nablauold_uadjold
//                                         - advection_flag * adv_res_uold_nablaphiadj_uadjold
//                                         + SolVAR_qp[SolPdeIndex[press_type_pos + adj_pos_begin]] * phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim]
//                                         ) * weight_iqp;
//       }
//       
// //BLOCK delta_adjoint - state------------------------------------------------------------------------------------------
//      for (unsigned j = 0; j < nDofsV; j++) {
// 	  for (unsigned kdim = 0; kdim < dim; kdim++) {
// 	      Jac[kdim + adj_pos_begin][kdim][i*nDofsV + j] += ( - cost_functional_coeff * target_flag * phi_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i] * phi_gss_fe[SolFEType_Mat[kdim]][j] 
//                                                              + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i]    * phi_x_gss_fe[ SolFEType_Mat[kdim] ][j*dim_offset_grad /*space_dim*/ + kdim] 		* SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]   //c(delta_u, u_new, lambda_old)  diagonal blocks  ......phiadj_nablaunew_uadjold 
//                                                              + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[kdim] ][j] 			* phi_x_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	 //c(u_new, delta_u, lambda_old) diagonal blocks  ......unew_nablaphiadj_uadjold
//                                                             ) * weight_iqp;
//               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
// 		Jac[kdim + adj_pos_begin][off_kdim][i*nDofsV + j] += (  advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i] * phi_x_gss_fe[ SolFEType_Mat[off_kdim] ][j*dim_offset_grad /*space_dim*/ + kdim]		      * SolVAR_qp[SolPdeIndex[off_kdim + adj_pos_begin]]   //c(delta_u, u_new, lambda_old)  off-diagonal blocks  ......phiadj_nablaunew_uadjold 
//                                                               + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[off_kdim] ][j] 		  * phi_x_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + off_kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	   //c(u_new, delta_u, lambda_old) off-diagonal blocks  ......unew_nablaphiadj_uadjold
//                                                              ) * weight_iqp;
// 	  }            
//      }//j_dadj_u loop
// 
// 
// //DIAG BLOCK delta_adjoint - adjoint---------------------------------------------------------------------------------
//      for (unsigned j = 0; j < nDofsVadj; j++) {
// 		    vector < double > lap_jac_dadj_adj(dim,0.);
// 		    vector < double > adv_uold_nablaphiadj_uadjnew(dim, 0.);
// 	  for (unsigned kdim = 0; kdim < dim; kdim++) {
//             for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { 
// 		  lap_jac_dadj_adj[kdim] += phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim] * phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + jdim];
//             }
//           }
//        for (unsigned  kdim = 0; kdim < dim; kdim++) { 
// 	   for (unsigned jdim = 0; jdim < dim; jdim++) { //diagonal blocks only
// 	     adv_uold_nablaphiadj_uadjnew[kdim]     += SolVAR_qp[SolPdeIndex[jdim]]  * phi_x_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][j] ;
// 	   }
// 	  }
// 	  for (unsigned kdim = 0; kdim < dim; kdim++) {
// 	      Jac[kdim + adj_pos_begin][kdim + adj_pos_begin][i*nDofsVadj + j] += (   IRe*lap_jac_dadj_adj[kdim] 
//                                                                                 + advection_flag * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim] * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][j]   //c(delta_u, u_old, lambda_new)  diagonal blocks  ......phiadj_nablauold_uadjnew  
//                                                                                 + advection_flag * adv_uold_nablaphiadj_uadjnew[kdim] 	//c(u_old, delta_u, lambda_new)
//                                                                                ) * weight_iqp;
//                unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
// 		  Jac[kdim + adj_pos_begin][off_kdim + adj_pos_begin][i*nDofsVadj + j] += ( advection_flag * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[off_kdim]][kdim] * phi_gss_fe[ SolFEType_Mat[off_kdim + adj_pos_begin] ][j]   //c(delta_u, u_old, lambda_new)  off-diagonal blocks  ......phiadj_nablauold_uadjnew   
//                                                                                   ) * weight_iqp;
// 	  }
//     } //j_dadj_adj loop
//       
// //BLOCK Pressure_adj
//     for (unsigned j = 0; j < nDofsPadj; j++) {
// 	  for (unsigned kdim = 0; kdim < dim; kdim++) {
// 	      Jac[kdim + adj_pos_begin][press_type_pos + adj_pos_begin][i*nDofsPadj + j] += -( phi_gss_fe[SolFEType_Mat[press_type_pos + adj_pos_begin]][j] * phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim] ) * weight_iqp;
// 	  }
//     }//j_press_adj loop
//   }//i_adj loop
// 
// //DIV_adj
//   for (unsigned i = 0; i < nDofsPadj; i++) {
// 		double div_adj_dadj_qp = 0.;
//       for (unsigned kdim = 0; kdim < dim; kdim++) {
// 	    div_adj_dadj_qp += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin ]][kdim] ;
//       }
//       Res[press_type_pos + adj_pos_begin][i] += ( (div_adj_dadj_qp) * phi_gss_fe[SolFEType_Mat[press_type_pos + adj_pos_begin]][i] ) * weight_iqp;
//       for (unsigned j = 0; j < nDofsVadj; j++) {
//         for (unsigned kdim = 0; kdim < dim; kdim++) {
//             Jac[press_type_pos + adj_pos_begin][kdim + adj_pos_begin][i*nDofsVadj + j] += - ( phi_gss_fe[SolFEType_Mat[press_type_pos + adj_pos_begin]][i] * phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + kdim] ) * weight_iqp;
//         }
//       }//j loop
//   }//i_div_adj
// 
// //============ delta_adjoint row - END =============================================================================================
// 
// //============ delta_control row - BEGIN  ==================================================================================================
// // delta_control
//     for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) {
//         
//          for (unsigned i = 0; i < nDofsGctrl; i++) {
//        Res[kdim + ctrl_pos_begin][i] += - penalty_outside_control_boundary * ( (1 - control_node_flag_iel_jface[kdim][i]) * (  sol_eldofs_Mat[SolPdeIndex[kdim + ctrl_pos_begin]][i] - 0.)  );              //enforce control zero outside the control boundary
// 
// 
// // //DIAG BLOCK delta_control - control--------------------------------------------------------------------------------------
//      for (unsigned j = 0; j < nDofsGctrl; j++) {
// 	    if (i == j) {
// 		Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsGctrl + j] += penalty_outside_control_boundary * (1 - control_node_flag_iel_jface[kdim][i]);              //enforce control zero outside the control boundary
//                   } //end i==j
//         }//j_dctrl_ctrl loop
//      }//i_ctrl loop
//   
//       }  //kdim
// 
// //============ delta_control row - END  ==================================================================================================
//  
//  
//       }  // end gauss point loop
//       
// //======================= VolumeLoop with Integration (and Zero boundary control outside Gamma_c) - END =====================================================    

    
    
// //======================= Loop without Integration For Fake Dofs related to Compatibility Condition - BEGIN =====================================================    
//     
//         //============ delta_theta - theta row ==================================================================================================
//   for (unsigned i = 0; i < nDofsThetactrl; i++) {
//              /* if ( fake_theta_flag[i] != 1 ) */             Res[ theta_index ][i]    = - (1 - fake_theta_flag[i]) * ( theta_value_outside_fake_element - sol_eldofs_Mat[SolPdeIndex[theta_index]][i]);  // Res_outer for the exact row (i.e. when fakeflag=1 , res =0(use Res_outer) and if not 1 this loop) and this is to take care of fake placement for the rest of dofs of theta values as 8
//      for (unsigned j = 0; j < nDofsThetactrl; j++) {
// 			         if(i == j)  Jac[ theta_index ][ theta_index ][i*nDofsThetactrl + j] = (1 - fake_theta_flag[i]) * 1.; //likewise Jac_outer (actually Jac itself works in the correct placement) for bdry integral and this is for rest of dofs
//              }//j_theta loop
//         }//i_theta loop
//    
//  //============ delta_theta row ==================================================================================================
// //======================= Loop without Integration For Fake Dofs related to Compatibility Condition - END =====================================================    



 //======================= From local to global - BEGIN =====================================================    

    //Sum the local matrices/vectors into the Global Matrix/Vector
    // FIRST ALL THE BLOCKS WITHOUT THETA ROW OR COLUMN 
    for(unsigned i_unk = 0; i_unk < n_unknowns-1; i_unk++) {
      RES->add_vector_blocked(Res[SolPdeIndex[i_unk]], L2G_dofmap_Mat[i_unk]);
        for(unsigned j_unk = 0; j_unk < n_unknowns-1; j_unk++) {
	  if(assembleMatrix) JAC->add_matrix_blocked( Jac[ SolPdeIndex[i_unk] ][ SolPdeIndex[j_unk] ], L2G_dofmap_Mat[i_unk], L2G_dofmap_Mat[j_unk]);
        }
    }
    
    // THEN THE BLOCKS WITH THETA ROW OR COLUMN 
	/*delta_theta-theta*/    JAC->add_matrix_blocked( Jac[ SolPdeIndex[n_unknowns-1] ][ SolPdeIndex[n_unknowns-1] ], L2G_dofmap_Mat[n_unknowns-1], L2G_dofmap_Mat[n_unknowns-1]);
	    
     if (does_iel_contain_a_bdry_control_face == 1) {
	      for (unsigned kdim = 0; kdim < dim; kdim++) {
                          /*delta_control*/       RES->add_vector_blocked(Res[SolPdeIndex[n_unknowns-2-kdim]], L2G_dofmap_Mat[n_unknowns-2-kdim]); 
		if(assembleMatrix) {
                          /*delta_theta-control*/ JAC->add_matrix_blocked( Jac[ SolPdeIndex[n_unknowns-1] ][ SolPdeIndex[n_unknowns-2-kdim] ], bdry_int_constr_pos_vec, L2G_dofmap_Mat[n_unknowns-2-kdim]);
                          /*delta_control-theta*/ JAC->add_matrix_blocked( Jac[ /*SolPdeIndex[n_unknowns-1] ][ SolPdeIndex[n_unknowns-2-kdim]*/SolPdeIndex[n_unknowns-2-kdim] ][ SolPdeIndex[n_unknowns-1] ], L2G_dofmap_Mat[n_unknowns-2-kdim], bdry_int_constr_pos_vec); 
		}
      }  //kdim
     }  //add control boundary element contributions
     
     
      if (does_iel_contain_a_bdry_control_face == 1) {
          /*delta_theta(bdry constr)*/         RES->add_vector_blocked(Res_outer, bdry_int_constr_pos_vec);
	  }
	  
     /* if (L2G_dofmap_Mat[n_unknowns-1][0] != bdry_int_constr_pos_vec[0]) */ /*delta_theta(fake)*/          RES->add_vector_blocked( Res[ SolPdeIndex[n_unknowns-1]],       L2G_dofmap_Mat[n_unknowns-1]);
	  

  
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
  
  
  JAC->close();
  RES->close();
  
  
    if (print_algebra_global) {
    assemble_jacobian< double, double >::print_global_jacobian(assembleMatrix, ml_prob, JAC, mlPdeSys->GetNonlinearIt());
//     assemble_jacobian< double, double >::print_global_residual(ml_prob, RES,  mlPdeSys->GetNonlinearIt());

    RES->close();
    std::ostringstream res_out; res_out << ml_prob.GetFilesHandler()->GetOutputPath() << "./" << "res_" << mlPdeSys->GetNonlinearIt()  << ".txt";
    pdeSys->print_with_structure_matlab_friendly(iproc, res_out.str().c_str(), RES);

  }

  
}





void ComputeIntegral(const MultiLevelProblem& ml_prob) {

   const NonLinearImplicitSystem & mlPdeSys   = ml_prob.get_system<NonLinearImplicitSystem> ("NSOpt");   
   const unsigned level = mlPdeSys.GetLevelToAssemble();
 

  const Mesh*          msh          	= ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
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
 
  std::vector<double> normal(dim_offset_grad /*space_dim*/, 0.);

  vector < vector < double > > coordX(dim);    // local coordinates
  vector< vector < double> > coordX_bd(dim);

  for (unsigned  k = 0; k < dim; k++) { 
        coordX[k].reserve(max_size);
        coordX_bd[k].reserve(max_size); 
  }
  
  double weight_iqp = 0.;
  double weight_iqp_bdry = 0.;
  
  
  //geometry *******************************

//STATE######################################################################
  //velocity *******************************
  vector < unsigned > solVIndex(dim);
  solVIndex[0] = ml_sol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = ml_sol->GetIndex("V");    // get the position of "V" in the ml_sol object

  if (dim == 3) solVIndex[2] = ml_sol->GetIndex("W");      // get the position of "V" in the ml_sol object

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
  

//CONTROL_@bdry######################################################################
  vector < unsigned > solVctrlIndex(dim);
  solVctrlIndex[0] = ml_sol->GetIndex("GX");    // get the position of "U" in the ml_sol object
  solVctrlIndex[1] = ml_sol->GetIndex("GY");    // get the position of "V" in the ml_sol object
  if (dim == 3) solVctrlIndex[2] = ml_sol->GetIndex("GZ");      // get the position of "V" in the ml_sol object

  unsigned solVctrlType = ml_sol->GetSolutionType(solVctrlIndex[0]);    // get the finite element type for "u"
  
  vector < vector < double > >  solVctrl(dim);    // local solution
  vector < double >   Vctrl_gss(dim, 0.);    //  solution
   
 for (unsigned  k = 0; k < dim; k++) {
    solVctrl[k].reserve(max_size);
  }

  
  vector <double> phiVctrl_gss_bd;  // local test function
  vector <double> phiVctrl_x_gss_bd; // local test function first order partial derivatives

  phiVctrl_gss_bd.reserve(max_size);
  phiVctrl_x_gss_bd.reserve(max_size * dim_offset_grad );
  
//CONTROL_@bdry######################################################################
  
//Theta value ######################################################################
   const unsigned solThetaIndex = ml_sol->GetIndex("THETA");
   const unsigned solThetaType = ml_sol->GetSolutionType(solThetaIndex);
   
//    double solTheta = (*sol->_Sol[solThetaIndex])(0)/*0.*/;
   //************** how to retrieve theta from proc0 ************************************* 
 double solTheta = get_theta_value(msh->n_processors(), sol, solThetaIndex);
//*************************************************** 
// 		     solTheta = (*sol->_Sol[solThetaIndex])(0);
//Theta value ######################################################################


// Vel_desired##################################################################
  vector <double> phiVdes_gss;  // local test function
  vector <double> phiVdes_x_gss; // local test function first order partial derivatives

  phiVdes_gss.reserve(max_size);
  phiVdes_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);

//   vector< vector < double > >  solVdes(dim);    // local solution
  vector <double>  solVdes(dim,0.);
  vector<double> Vdes_gss(dim, 0.);  
  
//  for (unsigned  k = 0; k < dim; k++) {
//     solVdes[k].reserve(max_size);
//   }
//   


// Vel_desired##################################################################

  
vector<double> integral(dim);

double  integral_target_alpha = 0.;

double	integral_beta   = 0.;
double	integral_gamma  = 0.;

double integral_g_dot_n = 0.;
  
//*************************************************** 
  //--- quadrature rules -------------------
  constexpr unsigned qrule_i = QRULE_I;
  
     std::vector < std::vector < double > >  JacI_qp(space_dim);
     std::vector < std::vector < double > >  Jac_qp(dim);
    for (unsigned d = 0; d < Jac_qp.size(); d++) {   Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp.size(); d++) { JacI_qp[d].resize(dim); }
    
    double detJac_qp;

     std::vector < std::vector < double > >  JacI_iqp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_iqp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_iqp_bdry.size(); d++) {   Jac_iqp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp_bdry.size(); d++) { JacI_iqp_bdry[d].resize(dim-1); }
    
    double detJac_iqp_bdry;
    
      //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all;
  ml_prob.get_all_abstract_fe_multiple(elem_all);
//*************************************************** 

  

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

   // geometry *****************************
      geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
      const short unsigned ielGeom = geom_element_iel.geom_type();
  // geometry end *****************************

// equation
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
//     unsigned nDofsVdes = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsVctrl = msh->GetElementDofNumber(iel, solVctrlType);    // number of solution element dofs
    
    unsigned nDofsThetactrl = msh->GetElementDofNumber(iel,solThetaType);
    
     
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

 
 //************ set control flag *********************
  int does_iel_contain_a_bdry_control_face = 0;
        does_iel_contain_a_bdry_control_face = ControlDomainFlag_bdry(geom_element_iel.get_elem_center_3d());
 //*************************************************** 

//========BoundaryLoop=====================================================================

  // Perform face loop over elements that contain some control face
  if (does_iel_contain_a_bdry_control_face == 1) {
	  
    double tau = 0.;
    vector<double> normal(dim_offset_grad /*space_dim*/,0);
	  
    for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {

       const unsigned ielGeom_bd = msh->GetElementFaceType(iel, jface);    
       
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();

	    // look for boundary faces
            const int bdry_index = el->GetFaceElementIndex(iel,jface);
            
	    if( bdry_index < 0) {
	   unsigned int face_in_rectangle_domain = -( msh->el->GetFaceElementIndex(iel,jface)+1);

	   if(  face_in_rectangle_domain == FACE_FOR_CONTROL) { //control face
	  
//=================================================== 
		
//========= initialize gauss quantities on the boundary ============================================
    vector < double >   Vctrl_bd_qp(dim, 0.);    //  solution@bdry
    vector < vector < double > > gradVctrl_bd_qp(dim);
      for (unsigned  k = 0; k < dim; k++) {
          gradVctrl_bd_qp[k].resize(dim_offset_grad /*space_dim*/);
          std::fill(gradVctrl_bd_qp[k].begin(), gradVctrl_bd_qp[k].end(), 0);
        }

//========= gauss_loop boundary===============================================================
	    for(unsigned iqp_bdry=0; iqp_bdry < ml_prob.GetQuadratureRule(ielGeom_bd).GetGaussPointsNumber(); iqp_bdry++) {
    elem_all[qrule_i][ielGeom_bd][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), iqp_bdry, Jac_iqp_bdry, JacI_iqp_bdry, detJac_iqp_bdry, space_dim);
	elem_all[qrule_i][ielGeom_bd][solType_coords]->compute_normal(Jac_iqp_bdry, normal);
    
    weight_iqp_bdry = detJac_iqp_bdry * ml_prob.GetQuadratureRule(ielGeom_bd).GetGaussWeightsPointer()[iqp_bdry];

    elem_all[qrule_i][ielGeom_bd][solVctrlType] ->shape_funcs_current_elem(iqp_bdry, JacI_iqp_bdry,phiVctrl_gss_bd,phiVctrl_x_gss_bd , boost::none , space_dim);

     
//========== compute gauss quantities on the boundary ===============================================
    for (unsigned  k = 0; k < dim; k++) {
	  Vctrl_bd_qp[k] = 0.;
	  for(unsigned ivar2 = 0; ivar2 < dim_offset_grad /*space_dim*/; ivar2++) { gradVctrl_bd_qp[k][ivar2] = 0.; }
	  
	  for (unsigned i = 0; i < nDofsVctrl; i++) {
                 const unsigned ndof_bdry = msh->GetElementFaceDofNumber(iel, jface, solVctrlType);

		   for(int i_bd = 0; i_bd < ndof_bdry; i_bd++) {
		       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bd);
		       Vctrl_bd_qp[k] += phiVctrl_gss_bd[i_bd] * solVctrl[k][i_vol];
		       for(unsigned ivar2 = 0; ivar2 < dim_offset_grad; ivar2++) {
			   gradVctrl_bd_qp[k][ivar2] += phiVctrl_x_gss_bd[i_bd * dim_offset_grad + ivar2 ] * solVctrl[k][i_vol]; 
		         }
		   }
      }
    }
 //end unknowns eval at gauss points ********************************
		  
//========== compute gauss quantities on the boundary ================================================
      for (unsigned  k = 0; k < dim; k++) {
	 integral_beta	+= ((Vctrl_bd_qp[k])*(Vctrl_bd_qp[k]) * weight_iqp_bdry);
	 integral_g_dot_n += Vctrl_bd_qp[k]*normal[k] * weight_iqp_bdry;
      }
      for (unsigned  k = 0; k < dim; k++) {
	for (unsigned  j = 0; j < dim; j++) {	
		integral_gamma	  += ((gradVctrl_bd_qp[k][j])*(gradVctrl_bd_qp[k][j]) * weight_iqp_bdry);
	}
      }


                }  //end iqp_bdryry loop
	  
             }    //end if control face
	 }  //end if boundary faces
      }  // loop over element faces //jface   
  } //end if control element flag

      // *** Gauss point loop ***
      for (unsigned iqp = 0; iqp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); iqp++) {
//STATE######## VolumeLoop #####################################################################	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), iqp, Jac_qp, JacI_qp, detJac_qp, space_dim);
    weight_iqp = detJac_qp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[iqp];
   
    elem_all[qrule_i][ielGeom][solVType]->shape_funcs_current_elem(iqp, JacI_qp, phiV_gss, phiV_x_gss, boost::none , space_dim);
    elem_all[qrule_i][ielGeom][solVType /*solVdes*/]->shape_funcs_current_elem(iqp, JacI_qp, phiVdes_gss, phiVdes_x_gss, boost::none , space_dim);

    
      for (unsigned  k = 0; k < dim; k++) {
	           V_gss[k] = 0.;
	           Vdes_gss[k] = 0.;
	    for (unsigned i = 0; i < nDofsV; i++) {
	   	   V_gss[k] += solV[k][i] * phiV_gss[i];
		   Vdes_gss[k] += solVdes[k]/*[i]*/ * phiVdes_gss[i];
	    }
	  }
	

      for (unsigned  k = 0; k < dim; k++) {
	      integral_target_alpha += (( target_flag ) *((V_gss[k]  - Vdes_gss[k]) * (V_gss[k]  - Vdes_gss[k])) * weight_iqp);
      }
      
      }// end gauss point loop
    } //end element loop  

       std::ostringstream filename_out; filename_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "Integral_computation"  << ".txt";

       std::ofstream intgr_fstream;
  if (paral::get_rank() == 0 ) {
      intgr_fstream.open(filename_out.str().c_str(),std::ios_base::app);
      intgr_fstream << " ***************************** Non Linear Iteration "<< mlPdeSys.GetNonlinearIt() << " *********************************** " <<  std::endl << std::endl;
      intgr_fstream << "The value of the target functional for " << "alpha " <<   std::setprecision(0) << std::scientific << cost_functional_coeff << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_target_alpha << std::endl;
      intgr_fstream << "The value of the L2 control for        " << "beta  " <<   std::setprecision(0) << std::scientific << alpha  << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_beta         << std::endl;
      intgr_fstream << "The value of the H1 control for        " << "gamma " <<   std::setprecision(0) << std::scientific << beta << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_gamma        << std::endl;
      intgr_fstream << "The value of the integral of g.n "<<    integral_g_dot_n << std::endl;
      intgr_fstream << "The value of the theta is                             " <<    std::setw(11) << std::setprecision(10) <<  solTheta << std::endl;
      intgr_fstream << "The value of the total integral is " << std::setw(11) << std::setprecision(10) <<  integral_target_alpha * cost_functional_coeff * 0.5  + integral_beta * alpha * 0.5 + integral_gamma * beta * 0.5 << std::endl;
      intgr_fstream <<  std::endl;
      intgr_fstream.close();  //you have to close to disassociate the file from the stream
}  
     
    return; 
	  
  
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
  const int n_vars_state = dim+1;
  const int n_unknowns = 3*n_vars_state; //(2.*dim)+1; //state , adjoint of velocity terms and one pressure term
  const int press_type_pos = dim;
  const int state_pos_begin = 0;
  const int adj_pos_begin   = dim+1;
  const int ctrl_pos_begin   = 2*(dim+1);
  const int theta_index = press_type_pos + ctrl_pos_begin;
  
  vector < std::string > Solname_Mat(n_unknowns);  // const char Solname_Mat[4][8] = {"U","V","W","P"};
  Solname_Mat              [state_pos_begin+0] =                "U";
  Solname_Mat              [state_pos_begin+1] =                "V";
  if (dim == 3) Solname_Mat[state_pos_begin+2] =                "W";
  Solname_Mat              [state_pos_begin + press_type_pos] = "P";
  
  Solname_Mat              [adj_pos_begin + 0] =              "UADJ";
  Solname_Mat              [adj_pos_begin + 1] =              "VADJ";
  if (dim == 3) Solname_Mat[adj_pos_begin + 2] =              "WADJ";
  Solname_Mat              [adj_pos_begin + press_type_pos] = "PADJ";

  Solname_Mat              [ctrl_pos_begin + 0] =              "GX";
  Solname_Mat              [ctrl_pos_begin + 1] =              "GY";
  if (dim == 3) Solname_Mat[ctrl_pos_begin + 2] =              "GZ";
  Solname_Mat              [ctrl_pos_begin + press_type_pos] = "THETA";
  
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
  double weight_iqp;
  
  
  //----------- dofs ------------------------------
  vector < vector < double > > sol_eldofs_Mat(n_unknowns);
  vector < vector < double > > gradsol_eldofs_Mat(n_unknowns);
  
  vector < vector < double > > SolVAR_coarser_prol_eldofs(n_unknowns);
  vector < vector < double > > gradSolVAR_coarser_prol_eldofs(n_unknowns);


  for(int k = 0; k < n_unknowns; k++) {
    sol_eldofs_Mat[k].reserve(max_size);
    gradsol_eldofs_Mat[k].reserve(max_size*dim); 
    
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
       sol_eldofs_Mat[k].resize(ndofs_unk);
       SolVAR_coarser_prol_eldofs[k].resize(ndofs_unk);
    for (unsigned i = 0; i < ndofs_unk; i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType_Mat[k]);    // global to global mapping between solution node and solution dof // via local to global solution node
       sol_eldofs_Mat[k][i] = (*sol->_Sol[SolIndex_Mat[k]])(solDof);      // global extraction and local storage for the solution
       SolVAR_coarser_prol_eldofs[k][i] = (*sol_coarser_prolongated->_Sol[SolIndex_Mat[k]])(solDof);      // global extraction and local storage for the solution
      }
    }
  //CTRL###################################################################

 
      // ********************** Gauss point loop *******************************
      for(unsigned i_qp = 0; i_qp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); i_qp++) {
	
 
      for(int fe=0; fe < NFE_FAMS; fe++) {
	msh->_finiteElement[ielGeom][fe]->Jacobian(coordX, i_qp, weight_iqp,phi_gss_fe[fe],phi_x_gss_fe[fe],phi_xx_gss_fe[fe]);
      }
         //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
  	msh->_finiteElement[ielGeom][BIQUADR_FE]->Jacobian(coordX,i_qp,weight_iqp,phi_gss_fe[BIQUADR_FE],phi_x_gss_fe[BIQUADR_FE],phi_xx_gss_fe[BIQUADR_FE]);

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
	    SolVAR_qp[unk] += phi_gss_fe[ SolFEType_Mat[unk] ][i] * sol_eldofs_Mat[unk][i];
	    SolVAR_coarser_prol_qp[unk] += phi_gss_fe[ SolFEType_Mat[unk] ][i] * SolVAR_coarser_prol_eldofs[unk][i];
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      gradSolVAR_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType_Mat[unk] ][i*dim+ivar2] * sol_eldofs_Mat[unk][i]; 
	      gradSolVAR_coarser_prol_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType_Mat[unk] ][i*dim+ivar2] * SolVAR_coarser_prol_eldofs[unk][i]; 
	    }
	  }
	  
	}  
 //end unknowns eval at gauss points ********************************


	for(unsigned unk = 0; unk < n_unknowns; unk++) {
        l2norm[unk] += ( SolVAR_qp[unk] - SolVAR_coarser_prol_qp[unk] ) * ( SolVAR_qp[unk] - SolVAR_coarser_prol_qp[unk] ) * weight_iqp ; 
        
     }
    
    
	for(unsigned unk = 0; unk < dim; unk++) {
        for(int j = 0; j < dim; j++){
        seminorm[unk] += (gradSolVAR_qp[unk][j] - gradSolVAR_coarser_prol_qp[unk][j] ) * ( gradSolVAR_qp[unk][j] - gradSolVAR_coarser_prol_qp[unk][j] ) * weight_iqp ;
        seminorm[unk + dim] += (gradSolVAR_qp[unk + adj_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + adj_pos_begin][j] ) * ( gradSolVAR_qp[unk + adj_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + adj_pos_begin][j] ) * weight_iqp ;
        seminorm[unk + 2*dim] += (gradSolVAR_qp[unk + ctrl_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j] ) * ( gradSolVAR_qp[unk + ctrl_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j] ) * weight_iqp ;
        
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
