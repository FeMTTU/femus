
#include <stdio.h>
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


//*********************** Sets Number of subdivisions in X and Y direction *****************************************

#define FACE_FOR_CONTROL  1
#include   "../nsopt_params.hpp"




#define exact_sol_flag 0 // 1 = if we want to use manufactured solution; 0 = if we use regular convention
#define compute_conv_flag 0 // 1 = if we want to compute the convergence and error ; 0 =  no error computation
#define no_of_ref 2     //mesh refinements

#define NO_OF_L2_NORMS 9   //U,V,P,UADJ,VADJ,PADJ,GX,GY,THETA
#define NO_OF_H1_NORMS 6    //U,V,UADJ,VADJ,GX, GY


using namespace femus;

 double penalty_outside_control_boundary = 1.e50;       // penalty for zero control outside Gamma_c
 double penalty_ctrl = 1.e10;         //penalty for u=q
 double theta_value_outside_fake_element = 0.;
  
 
bool SetBoundaryConditionOpt(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  //1: bottom  //2: right  //3: top  //4: left
  
  bool dirichlet = true;
   value = 0.;
 

                if (!strcmp(SolName, "GX"))       { if (facename == FACE_FOR_CONTROL) dirichlet = false; }
                if (!strcmp(SolName, "GY"))       { if (facename == FACE_FOR_CONTROL) dirichlet = false; }
                if (!strcmp(SolName, "GZ"))       { if (facename == FACE_FOR_CONTROL) dirichlet = false; }
                if (!strcmp(SolName, "THETA"))    { dirichlet = false; }
      
#if exact_sol_flag == 0
                if (!strcmp(SolName, "U"))       { if (facename == FACE_FOR_CONTROL) dirichlet = false; }
                if (!strcmp(SolName, "V"))       { if (facename == FACE_FOR_CONTROL) dirichlet = false; }
#endif
                if (!strcmp(SolName, "W"))       { if (facename == FACE_FOR_CONTROL) dirichlet = false; }
     
#if exact_sol_flag == 1
  //b.c. for manufactured lid driven cavity
  double pi = acos(-1.);
                if (!strcmp(SolName, "U"))       { if (facename == FACE_FOR_CONTROL) value =   sin(pi* x[0]) * sin(pi* x[0]) * cos(pi* x[1]) - sin(pi* x[0]) * sin(pi* x[0]); }
                if (!strcmp(SolName, "V"))       { if (facename == FACE_FOR_CONTROL) value = - sin(2. * pi * x[0]) * sin(pi* x[1]) + pi * x[1] * sin(2. * pi * x[0]); }
 #endif
               
  return dirichlet;

}

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
// //============== initial conditions =========
// double SetInitialCondition(const MultiLevelProblem * ml_prob, const std::vector <double> &x, const char SolName[]) {
//   
//   double value = 0.;
//   
//   if (x[1] < 1+ 1.e-5 && x[1] > 1 - 1.e-5 ) {
//                 if (!strcmp(SolName, "GX"))       { value = 1.; }
//                 if (!strcmp(SolName, "GY"))       { value = 1.; }
//   }
//   
//   return value;
// }
// //============== initial conditions =========



void AssembleNavierStokesOpt(MultiLevelProblem& ml_prob);    

void ComputeIntegral(const MultiLevelProblem& ml_prob);

double*  GetErrorNorm(const MultiLevelProblem& ml_prob, MultiLevelSolution* ml_sol, Solution* sol_coarser_prolongated);
// ||u_h - u_(h/2)||/||u_(h/2)-u_(h/4)|| = 2^alpha, alpha is order of conv 
//i.e. ||prol_(u_(i-1)) - u_(i)|| = err(i) => err(i-1)/err(i) = 2^alpha ,implemented as log(err(i)/err(i+1))/log2

void output_convergence_rate( double norm_i, double norm_ip1, std::string norm_name, unsigned maxNumberOfMeshes, int loop_i );


int main(int argc, char** args) {
  
  FemusInit mpinit(argc, args, MPI_COMM_WORLD); 	// init Petsc-MPI communicator
  
    // ======= Files ========================
  
  Files files; 
        files.CheckIODirectories();
        files.RedirectCout();
  
    // ======= Quad Rule ========================
    std::string fe_quad_rule("seventh");
    
  MultiLevelMesh mlMsh;			// define multilevel mesh
  MultiLevelMesh mlMsh_all_levels;
  double scalingFactor = 1.;		// read coarse level mesh and generate finers level meshes

   //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;
 // *** apparently needed by non-AD assemble only **********************
  // add fluid material
  Parameter parameter(Lref,Uref);
  
  // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  Fluid fluid(parameter,1,FLUID_DENSITY,"Newtonian");
  std::cout << "Fluid properties: " << std::endl;
  std::cout << fluid << std::endl;
  
// *************************
	
//   std::string input_file = "square_parametric.med";
//   std::string input_file = "square_4x5.med";
  std::string input_file = "Mesh_3_groups.med";
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_file;
  const std::string infile = mystream.str();
  
  //   MultiLevelMesh mlMsh;
 mlMsh.ReadCoarseMesh(infile.c_str(),fe_quad_rule.c_str(),Lref);
// //  mlMsh.RefineMesh(2, 2, NULL);
// //  mlMsh.EraseCoarseLevels(2-1);
#if compute_conv_flag == 1
 mlMsh_all_levels.ReadCoarseMesh(infile.c_str(),fe_quad_rule.c_str(),Lref);
#endif
//     mlMsh.GenerateCoarseBoxMesh(NSUB_X,NSUB_Y,0,0.,1.,0.,1.,0.,0.,QUAD9,fe_quad_rule.c_str());
//     mlMsh_all_levels.GenerateCoarseBoxMesh(NSUB_X,NSUB_Y,0,0.,1.,0.,1.,0.,0.,QUAD9,fe_quad_rule.c_str());
    
  unsigned dim = mlMsh.GetDimension();
  unsigned maxNumberOfMeshes;

  if (dim == 2) {
    maxNumberOfMeshes = no_of_ref;
  } else {
    maxNumberOfMeshes = 2;
  }


#if compute_conv_flag == 1
     double comp_conv[maxNumberOfMeshes][NO_OF_L2_NORMS+NO_OF_H1_NORMS];
 
  
        unsigned numberOfUniformLevels_finest = maxNumberOfMeshes;
        mlMsh_all_levels.RefineMesh(numberOfUniformLevels_finest, numberOfUniformLevels_finest, NULL);
//      mlMsh_all_levels.EraseCoarseLevels(numberOfUniformLevels - 2);  // need to keep at least two levels to send u_(i-1) projected(prolongated) into next refinement
        
        //store the fine solution  ==================
            MultiLevelSolution * ml_sol_all_levels;
            ml_sol_all_levels = new MultiLevelSolution (& mlMsh_all_levels);  //with the declaration outside and a "new" inside it persists outside the loop scopes
         // add variables to ml_sol_all_levels
        // state =====================  
            ml_sol_all_levels->AddSolution("U", LAGRANGE, SECOND);
            ml_sol_all_levels->AddSolution("V", LAGRANGE, SECOND);
            if (dim == 3) ml_sol_all_levels->AddSolution("W", LAGRANGE, SECOND);
            ml_sol_all_levels->AddSolution("P", LAGRANGE, FIRST);
        // adjoint =====================  
            ml_sol_all_levels->AddSolution("UADJ", LAGRANGE, SECOND);
            ml_sol_all_levels->AddSolution("VADJ", LAGRANGE, SECOND);
            if (dim == 3) ml_sol_all_levels->AddSolution("WADJ", LAGRANGE, SECOND);
            ml_sol_all_levels->AddSolution("PADJ", LAGRANGE, FIRST);
        // boundary condition =====================
            ml_sol_all_levels->AddSolution("GX", LAGRANGE, SECOND);
            ml_sol_all_levels->AddSolution("GY", LAGRANGE, SECOND);
            if (dim == 3) ml_sol_all_levels->AddSolution("GZ", LAGRANGE, SECOND);
            ml_sol_all_levels->AddSolution("THETA", DISCONTINUOUS_POLYNOMIAL, ZERO);
            ml_sol_all_levels->AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
            ml_sol_all_levels->AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
            
            ml_sol_all_levels->Initialize("All");
            ml_sol_all_levels->AttachSetBoundaryConditionFunction(SetBoundaryConditionOpt);
            ml_sol_all_levels->GenerateBdc("All");
#endif

         for (int i = 0; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

  unsigned numberOfUniformLevels = i + 1; 
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution ml_sol(&mlMsh);

  // add variables to ml_sol
  // state =====================  
  ml_sol.AddSolution("U", LAGRANGE, SECOND);
  ml_sol.AddSolution("V", LAGRANGE, SECOND);
  if (dim == 3) ml_sol.AddSolution("W", LAGRANGE, SECOND);
  ml_sol.AddSolution("P", LAGRANGE, FIRST);
  // adjoint =====================  
  ml_sol.AddSolution("UADJ", LAGRANGE, SECOND);
  ml_sol.AddSolution("VADJ", LAGRANGE, SECOND);
  if (dim == 3) ml_sol.AddSolution("WADJ", LAGRANGE, SECOND);
  ml_sol.AddSolution("PADJ", LAGRANGE, FIRST);
  // boundary condition =====================
  ml_sol.AddSolution("GX", LAGRANGE, SECOND);
  ml_sol.AddSolution("GY", LAGRANGE, SECOND);
  if (dim == 3) ml_sol.AddSolution("GZ", LAGRANGE, SECOND);
  ml_sol.AddSolution("THETA", DISCONTINUOUS_POLYNOMIAL, ZERO);
  // control ===================== 
  ml_sol.AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  ml_sol.AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  
   // ======= Problem  ==================
  MultiLevelProblem ml_prob(&ml_sol); 

  ml_prob.SetFilesHandler(&files);
  ml_prob.SetQuadratureRuleAllGeomElems(fe_quad_rule);
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;
  ml_prob.set_all_abstract_fe();
    
  
  // ======= Solution: Initial Conditions ==================
   ml_sol.Initialize("All");    // initialize all varaibles to zero
   ml_sol.Initialize("TargReg",     Solution_set_initial_conditions, & ml_prob);
   ml_sol.Initialize("ContReg",     Solution_set_initial_conditions, & ml_prob);

//   ml_sol.Initialize("GX", SetInitialCondition,&ml_prob);
//   ml_sol.Initialize("GY", SetInitialCondition,&ml_prob);
  
  
  // ======= Solution: Boundary Conditions ==================
  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionOpt);
  ml_sol.GenerateBdc("All");
  

 
  // add system OptBdryCtrl in ml_prob as a NonLinear Implicit System
  NonLinearImplicitSystem& system_opt    = ml_prob.add_system < NonLinearImplicitSystem > ("NSOpt");

  // ST ===================
  system_opt.AddSolutionToSystemPDE("U");
  system_opt.AddSolutionToSystemPDE("V");
  if (dim == 3) system_opt.AddSolutionToSystemPDE("W");
  system_opt.AddSolutionToSystemPDE("P");
//   ADJ ===================
  system_opt.AddSolutionToSystemPDE("UADJ");
  system_opt.AddSolutionToSystemPDE("VADJ");
  if (dim == 3) system_opt.AddSolutionToSystemPDE("WADJ");
  system_opt.AddSolutionToSystemPDE("PADJ");
  // BD ===================
  system_opt.AddSolutionToSystemPDE("GX");
  system_opt.AddSolutionToSystemPDE("GY");
  if (dim == 3)  system_opt.AddSolutionToSystemPDE("GZ");
  system_opt.AddSolutionToSystemPDE("THETA");
 
  
  // attach the assembling function to system
   system_opt.SetAssembleFunction(AssembleNavierStokesOpt);

   
  // initialize and solve the system
    system_opt.init();
    system_opt.ClearVariablesToBeSolved();
    system_opt.AddVariableToBeSolved("All");
  
    ml_sol.SetWriter(VTK);
    ml_sol.GetWriter()->SetDebugOutput(true);
    
    system_opt.SetDebugNonlinear(true);
    system_opt.SetDebugFunction(ComputeIntegral);
//   system_opt.SetMaxNumberOfNonLinearIterations(30);
//   system_opt.SetNonLinearConvergenceTolerance(1.e-15);
//     system_opt.SetMaxNumberOfLinearIterations(6);
//     system_opt.SetAbsoluteLinearConvergenceTolerance(1.e-14);
    system_opt.SetOuterSolver(PREONLY);
    system_opt.MGsolve();

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
       
 
  // print solutions
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
       output_convergence_rate(comp_conv[i][j], comp_conv[i + 1][j], norm_names_L2[j], maxNumberOfMeshes , i );
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
       output_convergence_rate(comp_conv[i][NO_OF_L2_NORMS + j], comp_conv[i + 1][NO_OF_L2_NORMS + j], norm_names_H1[j], maxNumberOfMeshes , i );
    }
  }
  std::cout << std::endl;
  std::cout << "=======================================================================" << std::endl;
#endif
  
  return 0;
}
 

void output_convergence_rate( double norm_i, double norm_ip1, std::string norm_name, unsigned maxNumberOfMeshes , int loop_i) {

    std::cout << loop_i + 1 << "\t\t" <<  std::setw(11) << std::setprecision(10) << norm_i << "\t\t\t\t" ;
  
    if (loop_i < maxNumberOfMeshes/*norm.size()*/ - 2) {
      std::cout << std::setprecision(3) << log( norm_i/ norm_ip1 ) / log(2.) << std::endl;
    }
  
}



//adj---------------------------------------------
void value_adjVel(const std::vector < double >& x, vector < double >& val_adjVel) {
  double pi = acos(-1.);
  val_adjVel[0] =   0.5 * sin(pi* x[0]) * sin(pi* x[0]) *  sin(2. * pi * x[1]); //u
  val_adjVel[1] = - 0.5 * sin(2. * pi * x[0]) * sin(pi* x[1]) * sin(pi* x[1]); //v
 };
 
double value_adjPress(const std::vector < double >& x) { 
  double pi = acos(-1.);
  return sin(2. * pi * x[0]) * sin(2. * pi * x[1]); //p
 };
 
 
void gradient_adjVel(const std::vector < double >& x, vector < vector < double > >& grad_adjVel) {
  double pi = acos(-1.);
  grad_adjVel[0][0]  =   0.5 * pi * sin(2. * pi * x[0]) * sin(2. * pi * x[1]); 
  grad_adjVel[0][1]  =   pi * sin(pi* x[0]) * sin(pi* x[0]) *  cos(2. * pi * x[1]);
  grad_adjVel[1][0]  = - pi * cos(2. * pi * x[0]) * sin(pi * x[1]) * sin(pi * x[1]); 
  grad_adjVel[1][1]  = - 0.5 * pi * sin(2. * pi * x[0]) * sin(2. * pi * x[1]);
 };

 void gradient_adjPress(const std::vector < double >& x, vector < double >& grad_adjPress) {
  double pi = acos(-1.);
  grad_adjPress[0]  =   2. * pi * cos(2. * pi * x[0]) * sin(2. * pi * x[1]); 
  grad_adjPress[1]  =   2. * pi * sin(2. * pi * x[0]) * cos(2. * pi * x[1]);
 };
 
 
void laplace_adjVel(const std::vector < double >& x, vector < double >& lap_adjVel) {
  double pi = acos(-1.);
  lap_adjVel[0] = pi * pi * cos(2. * pi * x[0]) * sin(2. * pi * x[1]) - 2. * pi * pi * sin(pi* x[0]) * sin(pi* x[0]) *  sin(2. * pi * x[1]);
  lap_adjVel[1] = 2. * pi * pi * sin(2. * pi * x[0]) * sin(pi* x[1]) * sin(pi* x[1]) - pi * pi * sin(2. * pi * x[0]) * cos(2. * pi * x[1]);
};
//adj---------------------------------------------

//state---------------------------------------------
void value_stateVel(const std::vector < double >& x, vector < double >& val_stateVel) {
  double pi = acos(-1.);
  val_stateVel[0] =   sin(pi* x[0]) * sin(pi* x[0]) * cos(pi* x[1]) - sin(pi* x[0]) * sin(pi* x[0]);
  val_stateVel[1] = - sin(2. * pi * x[0]) * sin(pi* x[1]) + pi * x[1] * sin(2. * pi * x[0]);
 };
 
 
void gradient_stateVel(const std::vector < double >& x, vector < vector < double > >& grad_stateVel) {
  double pi = acos(-1.);
  grad_stateVel[0][0]  =   pi * sin(2. * pi * x[0]) * cos(pi* x[1]) - pi * sin(2. * pi * x[0]);
  grad_stateVel[0][1]  = - pi * sin(pi* x[0]) * sin(pi* x[0]) *  sin(pi * x[1]); 
  grad_stateVel[1][0]  = - 2. * pi * cos(2. * pi * x[0]) * sin(pi* x[1]) + 2. * pi * pi * x[1] * cos(2. * pi * x[0]);   
  grad_stateVel[1][1]  = - pi * sin(2. * pi * x[0]) * cos(pi * x[1]) + pi * sin(2. * pi * x[0]); 
 };

  
void laplace_stateVel(const std::vector < double >& x, vector < double >& lap_stateVel) {
  double pi = acos(-1.);
  lap_stateVel[0] = - 2. * pi * pi * cos(2. * pi * x[0]) - 0.5 * pi * pi * cos(pi * x[1]) + 2.5 * pi * pi * cos(2. * pi* x[0]) * cos(pi* x[1]);
  lap_stateVel[1] = - 4. * pi * pi * pi * x[1] * sin(2. * pi * x[0]) + 5. * pi * pi * sin(2. * pi * x[0]) * sin(pi * x[1]);
};
//state---------------------------------------------


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
        

void AssembleNavierStokesOpt(MultiLevelProblem& ml_prob){
     
  //pointers
  NonLinearImplicitSystem& mlPdeSys  = ml_prob.get_system<NonLinearImplicitSystem>("NSOpt");
 const unsigned level = mlPdeSys.GetLevelToAssemble();

  bool assembleMatrix = mlPdeSys.GetAssembleMatrix(); 
   
  Solution*	 sol  	         = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearEquationSolver*  pdeSys	 = mlPdeSys._LinSolver[level];   
  const char* pdename            = mlPdeSys.name().c_str();
  
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;
  
  Mesh*		 msh    = ml_prob._ml_msh->GetLevel(level);
  elem*		 el	= msh->el;
  SparseMatrix*	 JAC	= pdeSys->_KK;
  NumericVector* RES 	= pdeSys->_RES;
  
  MatSetOption(static_cast< PetscMatrix* >(JAC)->mat(),MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    
  //data
  const unsigned dim 	= msh->GetDimension();
  unsigned dim2     = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned nel		= msh->GetNumberOfElements();
  unsigned igrid	= msh->GetLevel();
  unsigned iproc 	= msh->processor_id();

  const unsigned maxSize = static_cast< unsigned > (ceil(pow(3,dim)));

  //=============== Integration ========================================
  // quadratures ********************************
  double weight = 0.;
  double weight_bd = 0.;

   //=============== Geometry ========================================
   unsigned coordXType = 2; /*BIQUADR_FE*/// get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
   unsigned solType_coords = coordXType;  //FE_DOMAIN = 0; //we do linear FE this time // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
 
  CurrentElem < double > geom_element(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
  const unsigned int dim_offset_grad = /*dim*/  3  /*2*/    ;
 
  std::vector<double> normal(dim_offset_grad /*space_dim*/, 0.);
 
  // geometry *******************************************
  vector< vector < double> > coordX(dim);
  vector< vector < double> > coordX_bd(/*space_*/dim);
  for(int i=0;i<dim;i++) {   
       coordX[i].reserve(maxSize); 
       coordX_bd[i].reserve(maxSize); 
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
// //   phi_coords.reserve(maxSize);
// //   phi_coords_x.reserve(maxSize * space_dim);
// //   phi_coords_xx.reserve(maxSize * dim2);
// //   
// //   //boundary shape functions
// //   vector <double> phi_coords_bdry;  
// //   vector <double> phi_coords_x_bdry; 
// // 
// //   phi_coords_bdry.reserve(maxSize);
// //   phi_coords_x_bdry.reserve(maxSize * space_dim);
// //  //*************************************************** 
 
  // solution variables *******************************************
  
    std::vector<std::string> ctrl_name;
    ctrl_name.resize(dim);
    ctrl_name[0] = "GX";
    ctrl_name[1] = "GY";
     if (dim == 3)  ctrl_name[2] = "GZ";
 
  
  const int n_vars = dim+1;
  const int n_unknowns = 3*n_vars; //(2.*dim)+1; //state , adjoint of velocity terms and one pressure term
  const int vel_type_pos = 0;
  const int press_type_pos = dim;
  const int adj_vel_type_pos = vel_type_pos;
  const int state_pos_begin = 0;
  const int adj_pos_begin   = dim+1;
  const int ctrl_pos_begin   = 2*(dim+1);
  const int theta_index = press_type_pos + ctrl_pos_begin;

  vector < std::string > Solname(n_unknowns);  // const char Solname[4][8] = {"U","V","W","P"};
  Solname              [state_pos_begin+0] =                "U";
  Solname              [state_pos_begin+1] =                "V";
  if (dim == 3) Solname[state_pos_begin+2] =                "W";
  Solname              [state_pos_begin + press_type_pos] = "P";

  Solname              [adj_pos_begin + 0] =              "UADJ";
  Solname              [adj_pos_begin + 1] =              "VADJ";
  if (dim == 3) Solname[adj_pos_begin + 2] =              "WADJ";
  Solname              [adj_pos_begin + press_type_pos] = "PADJ";

  Solname              [ctrl_pos_begin + 0] =                  ctrl_name[0];
  Solname              [ctrl_pos_begin + 1] =                  ctrl_name[1];
  if (dim == 3) Solname[ctrl_pos_begin + 2] =                  ctrl_name[2];
  Solname              [theta_index ]       =               "THETA";
  
  vector < unsigned > SolPdeIndex(n_unknowns);
  vector < unsigned > SolIndex(n_unknowns);  
  vector < unsigned > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolPdeIndex[ivar]	= mlPdeSys.GetSolPdeIndex(Solname[ivar].c_str());
    SolIndex[ivar]	= ml_sol->GetIndex        (Solname[ivar].c_str());
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);
  }

  vector < double > Sol_n_el_dofs(n_unknowns);
  
  //==========================================================================================
  // velocity ************************************
  //-----------state------------------------------
  vector < vector < double > > phi_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_x_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_xx_gss_fe(NFE_FAMS);
  
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_gss_fe[fe].reserve(maxSize);
      phi_x_gss_fe[fe].reserve(maxSize*dim_offset_grad /*space_dim*/);
     phi_xx_gss_fe[fe].reserve(maxSize*dim2);
   }
   
   
  //boundary adjoint & ctrl shape functions  
  vector < vector < double > > phi_bd_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_x_bd_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_xx_bd_gss_fe_placeholder(NFE_FAMS);

  //bdry vol adj  evaluated at bdry points
   vector < vector < double > > phi_vol_at_bdry_fe(NFE_FAMS);
   vector < vector < double > > phi_x_vol_at_bdry_fe(NFE_FAMS);

    for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_bd_gss_fe[fe].reserve(maxSize);
      phi_x_bd_gss_fe[fe].reserve(maxSize*dim_offset_grad /*space_dim*/);
    phi_xx_bd_gss_fe_placeholder[fe].reserve(maxSize*dim2);
  //bdry vol adj  evaluated at bdry points
         phi_vol_at_bdry_fe[fe].reserve(maxSize);
       phi_x_vol_at_bdry_fe[fe].reserve(maxSize*dim_offset_grad /*space_dim*/);    
    }
   

  vector < vector < double > >  sol_adj_x_vol_at_bdry_gss(dim);
  vector < double > grad_dot_n_adj_res;
  vector < double > grad_adj_dot_n;
  for (int ldim =0; ldim < dim; ldim++) sol_adj_x_vol_at_bdry_gss[ldim].reserve(maxSize);
  grad_dot_n_adj_res.reserve(maxSize);
  grad_adj_dot_n.reserve(maxSize);
  
  //=================================================================================================
  
  
  // equation ***********************************
  vector < vector < int > > JACDof(n_unknowns); 
  vector < vector < double > > Res(n_unknowns);
  vector < vector < vector < double > > > Jac(n_unknowns);
  
  vector < vector < double > > Jac_outer(dim);
  vector < double > Res_outer(1);

  
  for(int i = 0; i < n_unknowns; i++) {     
    JACDof[i].reserve(maxSize);
       Res[i].reserve(maxSize);
  }
   
  if (assembleMatrix) {
    
    for(int i = 0; i < n_unknowns; i++) {
      Jac[i].resize(n_unknowns);    
      for(int j = 0; j < n_unknowns; j++) {
	Jac[i][j].reserve(maxSize*maxSize);	
      }
    }

         for(int i = 0; i < dim; i++) {  Jac_outer[i].reserve(maxSize); }
    
  }

  
  //----------- dofs ------------------------------
  vector < vector < double > > SolVAR_eldofs(n_unknowns); //sol_V,P_of_st,adj,ctrl
  vector < vector < double > > gradSolVAR_eldofs(n_unknowns);
  
  for(int k=0; k<n_unknowns; k++) {
    SolVAR_eldofs[k].reserve(maxSize);
    gradSolVAR_eldofs[k].reserve(maxSize*dim);    
  }

  //------------ at quadrature points ---------------------
    vector < double > SolVAR_qp(n_unknowns);   //sol_V,P_gss_of_st,adj,ctrl_ie@quadraturepoints
    vector < vector < double > > gradSolVAR_qp(n_unknowns);
    for(int k=0; k<n_unknowns; k++) {  gradSolVAR_qp[k].resize(dim_offset_grad /*space_dim*/);  }
      
      
  double IRe = ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();

  // Set to zero all the global structures
    RES->zero();
    if(assembleMatrix) JAC->zero();

 //*************************************************** 
     std::vector < std::vector < double > >  JacI_qp(space_dim);
     std::vector < std::vector < double > >  Jac_qp(dim);
    for (unsigned d = 0; d < Jac_qp.size(); d++) {   Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp.size(); d++) { JacI_qp[d].resize(dim); }
    
    double detJac_qp;

     std::vector < std::vector < double > >  JacI_qp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_qp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_qp_bdry.size(); d++) {   Jac_qp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp_bdry.size(); d++) { JacI_qp_bdry[d].resize(dim-1); }
    
    double detJac_qp_bdry;
    
    //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
//*************************************************** 
 
  
//************** how to retrieve theta from proc0 ************************************* 
 double solTheta = get_theta_value(msh->n_processors(), sol, SolIndex[theta_index]);
//*************************************************** 

  // ****************** element loop *******************
 
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

  // geometry *****************************
      geom_element.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
      const short unsigned ielGeom = geom_element.geom_type();
  // geometry end *****************************
  
  // equation *****************************
    unsigned nDofsV = msh->GetElementDofNumber(iel, SolFEType[vel_type_pos]);
    unsigned nDofsP = msh->GetElementDofNumber(iel, SolFEType[state_pos_begin + press_type_pos]);
    
    unsigned nDofsVadj = msh->GetElementDofNumber(iel,SolFEType[adj_pos_begin]);
    unsigned nDofsPadj = msh->GetElementDofNumber(iel,SolFEType[adj_pos_begin + press_type_pos]);

    unsigned nDofsGctrl = msh->GetElementDofNumber(iel,SolFEType[ctrl_pos_begin]);
    unsigned nDofsThetactrl = msh->GetElementDofNumber(iel,SolFEType[theta_index] );

    unsigned nDofsVP = dim * nDofsV + nDofsP; //same for state and adjoint
    unsigned nDofsVPctrl = dim * nDofsGctrl + nDofsThetactrl; //control
   
    unsigned nDofsVP_tot = 2*nDofsVP + (nDofsVPctrl);
  // equation end *****************************
  
  //***** set target domain flag ********************************** 
   geom_element.set_elem_center(iel, solType_coords);

   int target_flag = 0;
       target_flag = ElementTargetFlag(geom_element.get_elem_center()/*elem_center*/);
   //***************************************       
   
 //************ set control flag *********************
    int control_el_flag = 0;
        control_el_flag = ControlDomainFlag_bdry(geom_element.get_elem_center()/*elem_center*/);
    std::vector< std::vector<int> > control_node_flag(dim);
	    for(unsigned idim=0; idim < dim; idim++) {
	          control_node_flag[idim].resize(nDofsGctrl);
    std::fill(control_node_flag[idim].begin(), control_node_flag[idim].end(), 0);
	    }
 //*************************************************** 
  
   //STATE###################################################################  
    unsigned int fake_iel_flag = 0;
    unsigned int global_row_index_bdry_constr = pdeSys->KKoffset[SolPdeIndex[theta_index]][iproc];
  for (unsigned  k = 0; k < n_unknowns; k++) {
	unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);	//nDofs_V,P_of_st,adj,ctrl
	Sol_n_el_dofs[k]=ndofs_unk;
	SolVAR_eldofs[k].resize(ndofs_unk);	//sol_V,P_of_st,adj,ctrl
	JACDof[k].resize(ndofs_unk); 
    for (unsigned i = 0; i < ndofs_unk; i++) {
	      unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);    // global to global mapping between solution node and solution dof // via local to global solution node
	  SolVAR_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);      // global extraction and local storage for the solution
	         JACDof[k][i] = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
 
    if (k == SolPdeIndex[theta_index] && JACDof[k][i] == global_row_index_bdry_constr) {       fake_iel_flag = iel;  }
    }
  }
  //CTRL###################################################################
    
   //************ set fake theta flag: this flag tells me what degrees of freedom of the current element are fake for the theta variable *********************
    std::vector<int>  bdry_int_constr_pos_vec(1,global_row_index_bdry_constr); /*KKoffset[SolPdeIndex[PADJ]][iproc]*/
    std::vector<int> fake_theta_flag(nDofsThetactrl,0);
    for (unsigned i = 0; i < nDofsThetactrl; i++) {
      if ( JACDof[ SolPdeIndex[theta_index] ] [i] == bdry_int_constr_pos_vec[0]) { 	fake_theta_flag[i] = 1;       }
    }
    
   
 //************ end set fake theta flag *********************

// setting Jac and Res to zero ******************************* 
    for(int ivar=0; ivar<n_unknowns; ivar++) {
              Res[SolPdeIndex[ivar]].resize(Sol_n_el_dofs[ivar]);
      memset(&Res[SolPdeIndex[ivar]][0],0.,Sol_n_el_dofs[ivar]*sizeof(double));
    }
   
    for(int ivar=0; ivar<n_unknowns; ivar++) {
      for(int jvar=0; jvar<n_unknowns; jvar++) {
	    if(assembleMatrix){  //MISMATCH
                  Jac[ SolPdeIndex[ivar] ] [SolPdeIndex[jvar] ].resize(Sol_n_el_dofs[ivar]*Sol_n_el_dofs[jvar]);
		  memset(&Jac[ SolPdeIndex[ivar] ] [SolPdeIndex[jvar] ][0],0., Sol_n_el_dofs[ivar]*Sol_n_el_dofs[jvar]*sizeof(double));
           }
        }
     }
     
    for(int ivar = 0; ivar < dim; ivar++)     std::fill(Jac_outer[ivar].begin(), Jac_outer[ivar].end(), 0.); //did not use Jac_outer as Jac itself was placing the values as expected
    Res_outer[0] = 0.;
 // setting Jac and Res to zero ******************************* 

  
  
//========BoundaryLoop=====================================================================

  // Perform face loop over elements that contain some control face
  if (control_el_flag == 1) {
	  
      double tau = 0.;
      vector<double> normal(dim_offset_grad /*space_dim*/,0);
	       
	  // loop on faces of the current element

      for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {
          
       const unsigned ielGeom_bd = msh->GetElementFaceType(iel, jface);    
       const unsigned nve_bd = msh->GetElementFaceDofNumber(iel,jface,solType_coords);
       
       geom_element.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);
 
       geom_element.set_elem_center_bdry_3d();

	    // look for boundary faces
            const int bdry_index = el->GetFaceElementIndex(iel,jface);
            
	    if( bdry_index < 0) {
	      unsigned int face_in_rectangle_domain = -( msh->el->GetFaceElementIndex(iel,jface)+1);

// 	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"U",tau,face,0.) && tau!=0.){
          if(  face_in_rectangle_domain == FACE_FOR_CONTROL) { //control face
              
//=================================================== 
		   //we use the dirichlet flag to say: if dirichlet == true, we set 1 on the diagonal. if dirichlet == false, we put the boundary equation
		  std::vector<bool> dir_bool(dim);
		  for(unsigned idim=0; idim < dim; idim++) {
		      dir_bool[idim] = false; //ml_sol->GetBdcFunction()(geom_element.get_elem_center_bdry(),ctrl_name[idim].c_str(),tau,face_in_rectangle_domain,0.);
		  }
	  
	
//========= initialize gauss quantities on the boundary ============================================
		vector < double >                      SolVAR_bd_qp(n_unknowns);
		vector < vector < double > >       gradSolVAR_bd_qp(n_unknowns);
		for(int k=0; k<n_unknowns; k++) {  gradSolVAR_bd_qp[k].resize(dim_offset_grad /*space_dim*/);  }

//========= gauss_loop boundary===============================================================
		  for(unsigned ig_bd=0; ig_bd < ml_prob.GetQuadratureRule(ielGeom_bd).GetGaussPointsNumber(); ig_bd++) {
    
    elem_all[ielGeom_bd][solType_coords]->JacJacInv(geom_element.get_coords_at_dofs_bdry_3d(), ig_bd, Jac_qp_bdry, JacI_qp_bdry, detJac_qp_bdry, space_dim);
	elem_all[ielGeom_bd][solType_coords]->compute_normal(Jac_qp_bdry, normal);
    
    weight_bd = detJac_qp_bdry * ml_prob.GetQuadratureRule(ielGeom_bd).GetGaussWeightsPointer()[ig_bd];

    elem_all[ielGeom_bd][SolFEType[theta_index]] ->shape_funcs_current_elem(ig_bd, JacI_qp_bdry,phi_bd_gss_fe[SolFEType[theta_index]],phi_x_bd_gss_fe[SolFEType[theta_index]], phi_xx_bd_gss_fe_placeholder[SolFEType[theta_index]] , space_dim);
    elem_all[ielGeom_bd][SolFEType[ctrl_pos_begin]] ->shape_funcs_current_elem(ig_bd, JacI_qp_bdry,phi_bd_gss_fe[SolFEType[ctrl_pos_begin]],phi_x_bd_gss_fe[SolFEType[ctrl_pos_begin]], phi_xx_bd_gss_fe_placeholder[SolFEType[ctrl_pos_begin]] , space_dim);


    elem_all[ielGeom][solType_coords]->JacJacInv_vol_at_bdry_new(geom_element.get_coords_at_dofs_3d(), ig_bd, jface, Jac_qp/*not_needed_here*/, JacI_qp, detJac_qp/*not_needed_here*/, space_dim);
    elem_all[ielGeom][SolFEType[adj_pos_begin]]->shape_funcs_vol_at_bdry_current_elem(ig_bd, jface, JacI_qp, phi_vol_at_bdry_fe[SolFEType[adj_pos_begin]],phi_x_vol_at_bdry_fe[SolFEType[adj_pos_begin]], boost::none, space_dim);
     
		  
//========== compute gauss quantities on the boundary ===============================================
		    for (unsigned  kdim = 0; kdim < dim; kdim++) {
			    unsigned int ctrl_index = kdim + ctrl_pos_begin;
										  SolVAR_bd_qp[ SolPdeIndex[ctrl_index] ] = 0.;
			  for(unsigned ivar2=0; ivar2< dim_offset_grad /*space_dim*/; ivar2++) {  gradSolVAR_bd_qp[ SolPdeIndex[ctrl_index] ][ivar2] = 0.; }
	  
			  for(int i_bd = 0; i_bd < nve_bd; i_bd++) {
		                  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bd);
                                                                    SolVAR_bd_qp[SolPdeIndex[ctrl_index]]           += phi_bd_gss_fe  [ SolFEType[ctrl_index] ][i_bd]                  * SolVAR_eldofs[ SolPdeIndex[ ctrl_index ] ][i_vol];
			      for(unsigned ivar2=0; ivar2< dim_offset_grad /*space_dim*/; ivar2++) {  gradSolVAR_bd_qp[SolPdeIndex[ctrl_index]][ivar2] 	+= phi_x_bd_gss_fe[ SolFEType[ctrl_index] ][i_bd * dim_offset_grad /*space_dim*/ + ivar2 /*i_bd + ivar2 * nve_bd*/] * SolVAR_eldofs[ SolPdeIndex[ ctrl_index ] ][i_vol];     /* }*/
                        }
			  }
		    }//kdim
 //end unknowns eval at gauss points ********************************
 
//=============== grad dot n for residual ========================================= 
//     compute gauss quantities on the boundary through VOLUME interpolation
		for(unsigned ldim=0; ldim<dim; ldim++) {   sol_adj_x_vol_at_bdry_gss[ldim].resize(dim_offset_grad /*space_dim*/); }
		grad_dot_n_adj_res.resize(dim);
		grad_adj_dot_n.resize(dim);
		for(unsigned ldim=0; ldim<dim; ldim++) {   std::fill(sol_adj_x_vol_at_bdry_gss[ldim].begin(), sol_adj_x_vol_at_bdry_gss[ldim].end(), 0.);  }
		
		for(unsigned ldim=0; ldim<dim; ldim++) {  
		  for (int iv = 0; iv < nDofsVadj; iv++)  {
                     for (int d = 0; d < dim_offset_grad /*space_dim*/; d++) {
			   sol_adj_x_vol_at_bdry_gss[ldim][d] += SolVAR_eldofs[SolPdeIndex[ldim + adj_pos_begin]][iv] * phi_x_vol_at_bdry_fe[SolFEType[ldim + adj_pos_begin]][iv * dim_offset_grad /*space_dim*/ + d];//notice that the convention of the orders x y z is different from vol to bdry
                    }
		  }  
		      
		  grad_dot_n_adj_res[ldim] = 0.;
		  for(unsigned d=0; d < dim_offset_grad /*space_dim*/; d++) {
		      grad_dot_n_adj_res[ldim] += sol_adj_x_vol_at_bdry_gss[ldim][d]*normal[d];  
		  }
		}
		
//=============== grad dot n  for residual =========================================       
		  
//========== compute gauss quantities on the boundary ================================================

		
//============ Res _ Boundary Integral Constraint ============================================================================================
	  for (unsigned  kdim = 0; kdim < dim; kdim++) {
// 		for(unsigned i=0; i < nDofsThetactrl; i ++) { avoid because it is an element dof
/*delta_theta row */ 	/* Res[theta_index][i]*/ Res_outer[0] +=  /*fake_theta_flag[i] **/ weight_bd * SolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * normal[kdim] ;
// 		}  
	  }
		  
//============End of Res _ Boundary Integral Constraint ============================================================================================
		
  // *** phi_i loop ***
		for(unsigned i_bdry=0; i_bdry < nve_bd; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
                 
//=============== construct control node flag field on the go  =========================================    
	      /* (control_node_flag[i])       picks nodes on \Gamma_c
	         (1 - control_node_flag[i])   picks nodes on \Omega \setminus \Gamma_c
	       */
		    for(unsigned idim=0; idim<dim; idim++) {
			if (dir_bool[idim] == false) { 
// 		std::cout << " found boundary control nodes ==== " << std::endl;
			    for(unsigned k=0; k < control_node_flag[idim].size(); k++) {
				  control_node_flag[idim][i_vol] = 1;
			    }
			}
		    }
//=============== construct control node flag field on the go  =========================================    
		
		  
//Boundary Residuals  and Jacobians ==================	

		  
//============ Boundary Residuals============================================================================================
		  
		      for (unsigned  kdim = 0; kdim < dim; kdim++) {
			
			    double lap_res_dctrl_ctrl_bd = 0.;
			      for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
				  lap_res_dctrl_ctrl_bd += gradSolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim]*phi_x_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i_bdry * dim_offset_grad /*space_dim*/ + jdim /*i_bdry + jdim*nve_bd*/];
			      }//jdim

			
/*delta_state row */	    if(i_vol<nDofsV)     Res[kdim]                 [i_vol]  += - control_node_flag[kdim][i_vol] * penalty_ctrl * (SolVAR_eldofs[SolPdeIndex[kdim + state_pos_begin]][i_vol] - SolVAR_eldofs[SolPdeIndex[kdim + ctrl_pos_begin]][i_vol]);	    //u-g
/*delta_adjoint row */     if(i_vol<nDofsVadj)   Res[kdim + adj_pos_begin] [i_vol]  += 0.;	   
/*delta_control row */     if(i_vol<nDofsGctrl)  Res[kdim + ctrl_pos_begin][i_vol]  += - control_node_flag[kdim][i_vol] * weight_bd * (
                                                                                          beta_val* SolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * phi_bd_gss_fe[SolFEType[kdim +  ctrl_pos_begin]][i_bdry]
                                                                                        + gamma_val* lap_res_dctrl_ctrl_bd
                                                                                        - IRe * grad_dot_n_adj_res[kdim]  * phi_bd_gss_fe[SolFEType[kdim +  ctrl_pos_begin]][i_bdry]
                                                                                        - /*(*sol->_Sol[SolIndex[theta_index]])(0)*/solTheta * phi_bd_gss_fe[SolFEType[kdim +  ctrl_pos_begin]][i_bdry] * normal[kdim]      //*sol->_Sol[SolIndex[theta_index]])(0) finds the global value from KKDof pos(72, 169,etc), SolVAReldof_theta gets the value in the boundary point which will be zero. Theta is just a const
                                                                                        );	    
		      }//kdim  

//============ Boundary Residuals  ==================================================================================================


//============ Jac _ Boundary Integral Constraint ============================================================================================
		    for (unsigned  kdim = 0; kdim < dim; kdim++) { 
			  for(unsigned i =0; i < nDofsThetactrl; i ++) {
			    if(i_vol < nDofsGctrl) {
				double temp = weight_bd * ( phi_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i_bdry] * normal[kdim]);
//ROW_BLOCK delta_theta - control -- loop over i in the VOLUME (while j(/i_vol) is in the boundary) -------------------------------------------------------------------------------------------------------------
			      Jac[theta_index][ctrl_pos_begin + kdim][i*nDofsGctrl + i_vol]     += - temp; /*weight_bd * ( phi_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i_bdry] * normal[kdim])*/
//COLUMN_BLOCK delta_control - theta ---- loop over j in the VOLUME (while i(/i_vol) is in the boundary) ---------------------------------------------------------------------------------------------------
			      Jac[ctrl_pos_begin + kdim][theta_index][i_vol*nDofsThetactrl + i] += - control_node_flag[kdim][i_vol] /** phi_bd_gss_fe[SolFEType[theta_index]][i]*/*temp; /*weight_bd * ( phi_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i_bdry] * normal[kdim]);*/
			    }//endif
			  }// i 
		    }//kdim
//============ End of Jac _ Boundary Integral Constraint ============================================================================================



//============  ==================================================================================================
		  
		      for(unsigned j_bdry=0; j_bdry < nve_bd; j_bdry ++) {
			  unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, jface, j_bdry);

			  for (unsigned  kdim = 0; kdim < dim; kdim++) { 
//DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
			    if(i_vol < nDofsV && j_vol < nDofsV && i_vol == j_vol)       		          Jac[kdim][kdim][i_vol*nDofsV + j_vol]	    += penalty_ctrl * control_node_flag[kdim][i_vol];  //u
			 
//BLOCK delta_state - control------------------------------------------------------------------------------------
			    if(i_vol < nDofsV && j_vol < nDofsGctrl && i_vol == j_vol) 	Jac[kdim][kdim + ctrl_pos_begin][i_vol*nDofsGctrl + j_vol]  += (-1.) * penalty_ctrl *control_node_flag[kdim][i_vol];  //-g
			  
			  }//kdim		      


//DIAG BLOCK delta_control - control  --------------------------------------------------------------------------------------
			  if(i_vol < nDofsGctrl && j_vol < nDofsGctrl) {
				vector < double > lap_jac_dctrl_ctrl_bd(dim, 0.);
			     for (unsigned kdim = 0; kdim < dim; kdim++) {
				for (unsigned ldim = 0; ldim < dim_offset_grad /*space_dim*/; ldim++) {
                    lap_jac_dctrl_ctrl_bd[kdim] += phi_x_bd_gss_fe[SolFEType[ldim + ctrl_pos_begin]][i_bdry * dim_offset_grad /*space_dim*/ + ldim /*i_bdry + ldim*nve_bd*/]*phi_x_bd_gss_fe[SolFEType[ldim + ctrl_pos_begin]][j_bdry * dim_offset_grad /*space_dim*/ + ldim /*j_bdry + ldim*nve_bd*/];
                    }
                 }

			      for (unsigned kdim = 0; kdim < dim; kdim++) {
				      Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i_vol*nDofsGctrl + j_vol] +=   control_node_flag[kdim][i_vol] * weight_bd *(
                                                                                                       beta_val * phi_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin] ][i_bdry] * phi_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin] ][j_bdry]
                                                                                                    + gamma_val * lap_jac_dctrl_ctrl_bd[kdim]
                                                                                                    );
			      }//kdim
			  }//endif
                   
                       }//end j_bdry loop
		    
//BLOCK delta_control - adjoint------------------------------------------------------------------------------------------------
//===================loop over j in the VOLUME (while i is in the boundary)	      
		for(unsigned j=0; j < nDofsVadj; j ++) {
			//=============== grad dot n  =========================================    
		    for(unsigned ldim=0; ldim<dim; ldim++) {
			    grad_adj_dot_n[ldim] = 0.;
			    for(unsigned d=0; d<dim_offset_grad /*space_dim*/; d++) {
				grad_adj_dot_n[ldim] += phi_x_vol_at_bdry_fe[SolFEType[ldim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + d]*normal[d];  //notice that the convention of the orders x y z is different from vol to bdry
			    }
		    }
			
			  //=============== grad dot n  =========================================    

			  for (unsigned kdim = 0; kdim < dim; kdim++) {
				Jac[kdim + ctrl_pos_begin][kdim + adj_pos_begin][i_vol*nDofsVadj + j] += control_node_flag[kdim][i_vol] * (-1.) * (weight_bd  * phi_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i_bdry]* IRe * grad_adj_dot_n[kdim]);    		      
			  }
		} // end j loop for volume 
		  
		    }//end i_bdry loop

                }  //end ig_bdry loop
	  
             }    //end if control face
	 }  //end if boundary faces
      }  // loop over element faces //jface   
  } //end if control element flag

//End Boundary Residuals  and Jacobians ==================	
    
    
    
    
//======================= Loop without Integration =====================================================    
    
        //============ delta_theta - theta row ==================================================================================================
  for (unsigned i = 0; i < nDofsThetactrl; i++) {
             /* if ( fake_theta_flag[i] != 1 ) */             Res[ theta_index ][i]    = - (1 - fake_theta_flag[i]) * ( theta_value_outside_fake_element - SolVAR_eldofs[SolPdeIndex[theta_index]][i]);  // Res_outer for the exact row (i.e. when fakeflag=1 , res =0(use Res_outer) and if not 1 this loop) and this is to take care of fake placement for the rest of dofs of theta values as 8
     for (unsigned j = 0; j < nDofsThetactrl; j++) {
			         if(i==j)  Jac[ theta_index ][ theta_index ][i*nDofsThetactrl + j] = (1 - fake_theta_flag[i]) * 1.; //likewise Jac_outer (actually Jac itself works in the correct placement) for bdry integral and this is for rest of dofs
             }//j_theta loop
        }//i_theta loop
   
 //============ delta_theta row ==================================================================================================
 //======================= Loop without Integration =====================================================    

 
 
 
//======================= VolumeLoop with Integration (and fake boundary) =====================================================    
// ********************** Gauss point loop *******************************
 for(unsigned ig=0;ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {
	
	// *** get Jacobian and test function and test function derivatives ***
       // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[ielGeom][solType_coords]->JacJacInv(geom_element.get_coords_at_dofs_3d(), ig, Jac_qp, JacI_qp, detJac_qp, space_dim);
    weight = detJac_qp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[ig];
   
     for(int fe=0; fe < NFE_FAMS; fe++) {
    elem_all[ielGeom][fe]->shape_funcs_current_elem(ig, JacI_qp,phi_gss_fe[fe],phi_x_gss_fe[fe],phi_xx_gss_fe[fe] , space_dim);
      }
 
 
 vector < double > coordX_gss(dim, 0.);
 	for(unsigned k = 0; k <  dim; k++) {
	  for(unsigned i = 0; i < Sol_n_el_dofs[k]; i++) {
         coordX_gss[k] += coordX[k][i] * phi_gss_fe[ SolFEType[k] ][i];
      }
    }
//begin unknowns eval at gauss points ********************************
	for(unsigned unk = 0; unk < /*n_vars*/ n_unknowns; unk++) {
	    SolVAR_qp[unk] = 0.;
	    for(unsigned ivar2=0; ivar2 < dim_offset_grad /*space_dim*/; ivar2++){ 
		gradSolVAR_qp[unk][ivar2] = 0.; 
	    }
	  
	    for(unsigned i = 0; i < Sol_n_el_dofs[unk]; i++) {
                        SolVAR_qp[unk] += phi_gss_fe[ SolFEType[unk] ][i]             * SolVAR_eldofs[SolPdeIndex[unk]][i];
		for(unsigned ivar2 = 0; ivar2 < dim_offset_grad /*space_dim*/; ivar2++) {       
		    gradSolVAR_qp[unk][ivar2]  += phi_x_gss_fe[ SolFEType[unk] ][i * dim_offset_grad /*space_dim*/ + ivar2] * SolVAR_eldofs[SolPdeIndex[unk]][i]; 
		}
	    }//ndofsunk
	  
	} //unk 
 //end unknowns eval at gauss points ********************************
	
//computation of RHS (force and desired velocity) using MMS=============================================== 
//state values-------------------- //non-hom bdry
vector <double>  exact_stateVel(dim);
value_stateVel(coordX_gss, exact_stateVel);
vector <double>  exact_lap_stateVel(dim);
laplace_stateVel(coordX_gss, exact_lap_stateVel);
vector < vector < double > > exact_grad_stateVel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_stateVel[k].resize(dim);
    std::fill(exact_grad_stateVel[k].begin(), exact_grad_stateVel[k].end(), 0.);
}
gradient_stateVel(coordX_gss,exact_grad_stateVel);

//adjoint values--------------------//hom bdry
vector <double>  exact_adjVel(dim);
value_adjVel(coordX_gss, exact_adjVel);
vector <double>  exact_lap_adjVel(dim);
laplace_adjVel(coordX_gss, exact_lap_adjVel);
vector < vector < double > > exact_grad_adjVel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_adjVel[k].resize(dim);
    std::fill(exact_grad_adjVel[k].begin(), exact_grad_adjVel[k].end(), 0.);
}
gradient_adjVel(coordX_gss,exact_grad_adjVel);
vector <double> exact_grad_adjPress(dim);
gradient_adjPress(coordX_gss, exact_grad_adjPress);

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
    exactVel_d[k] =   exact_stateVel[k] + (1./alpha_val) * (IRe * exact_lap_adjVel[k] - exact_grad_adjPress[k]) 
                    + (1./alpha_val) * advection_flag * (exact_conv_u_nabla_uadj[k] - exact_conv_nabla_uT_uadj[k]);
}
//computation of RHS (force and desired velocity) using MMS=============================================== 
 
//============ delta_state row ============================================================================================

  for (unsigned i = 0; i < nDofsV; i++) {
// FIRST ROW
	for (unsigned  kdim = 0; kdim < dim; kdim++) { // velocity block row 
                double lap_res_du_u		= 0.; 
                double adv_res_uold_nablauold 	= 0.;
	      for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		    lap_res_du_u 	       += gradSolVAR_qp[SolPdeIndex[kdim]][jdim] * phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad /*space_dim*/ + jdim];
          }
	      for (unsigned jdim = 0; jdim < dim; jdim++) {
		   adv_res_uold_nablauold  += SolVAR_qp[SolPdeIndex[jdim]]  * gradSolVAR_qp[SolPdeIndex[kdim]][jdim] * phi_gss_fe[ SolFEType[kdim] ][i];
	      }      
	      Res[kdim][i]   +=  (           
#if exact_sol_flag == 0
                                         + force[kdim] * phi_gss_fe[ SolFEType[kdim] ][i]
 #endif                                      
 #if exact_sol_flag == 1
                                       + exactForce[kdim] * phi_gss_fe[ SolFEType[kdim] ][i]
 #endif
                                       - IRe*lap_res_du_u 
                                       - advection_flag * adv_res_uold_nablauold 
                                       + SolVAR_qp[SolPdeIndex[press_type_pos]] * phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad /*space_dim*/ + kdim]) * weight; 
	}	    
//DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
	for (unsigned j = 0; j < nDofsV; j++) {
		      vector < double > lap_jac_du_u(dim,0.);
		      vector < double > adv_uold_nablaunew(dim,0.);
	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { 
		    lap_jac_du_u[kdim] += phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad /*space_dim*/ + jdim]*phi_x_gss_fe[SolFEType[kdim]][j * dim_offset_grad/*space_dim*/ + jdim];
            }
          }
        for (unsigned  kdim = 0; kdim < dim; kdim++) { 
         for (unsigned  jdim = 0; jdim < dim; jdim++) { //diagonal blocks only
		    adv_uold_nablaunew[kdim] 	 += SolVAR_qp[SolPdeIndex[jdim]] * phi_x_gss_fe[ SolFEType[kdim] ][j * dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim] ][i];
                }  //jdim
	      }
	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
		Jac[kdim][kdim][i*nDofsV + j] += (   IRe * lap_jac_du_u[kdim] 
                                            + advection_flag * adv_uold_nablaunew[kdim] 		// c(u_old, u_new, delta_lambda)
                                            + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim] * phi_gss_fe[ SolFEType[kdim] ][i]	 // c(u_new,u_old,delta_lambda) diagonal blocks  ..... unew_nablauold
                                         ) * weight; 
              unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
		Jac[kdim][off_kdim][i*nDofsV + j] += (	advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][off_kdim] * phi_gss_fe[ SolFEType[kdim] ][i]	// c(u_new,u_old,delta_lambda) off-diagonal blocks  ..... unew_nablauold
                                             ) * weight;
	      }
	} //j_du_u loop
     
//BLOCK Pressure
      for (unsigned j = 0; j < nDofsP; j++) {
	    for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim][press_type_pos][i*nDofsP + j] += -( phi_gss_fe[SolFEType[press_type_pos]][j] * phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad /*space_dim*/ + kdim] ) * weight;
	    }
      }//j_press loop
   }//i_state loop

//DIV_state
  for (unsigned i = 0; i < nDofsP; i++) {
		    double div_u_du_qp =0.;
      for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      div_u_du_qp += gradSolVAR_qp[SolPdeIndex[kdim]][kdim] ;
      }
      Res[press_type_pos][i]  +=  ( (div_u_du_qp) * phi_gss_fe[SolFEType[press_type_pos]][i] ) * weight;
      for (unsigned j = 0; j < nDofsV; j++) {
	  for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      Jac[press_type_pos][kdim][i*nDofsV + j] += -( phi_gss_fe[SolFEType[press_type_pos]][i] * phi_x_gss_fe[SolFEType[kdim]][j * dim_offset_grad /*space_dim*/ + kdim] ) * weight;
	  }
      } //j loop
   }//i_div_state
    //============ delta_state row ============================================================================================


    
//============ delta_adjoint row =============================================================================================
  
  for (unsigned i = 0; i < nDofsVadj; i++) {
// SECOND ROW
     for (unsigned kdim = 0; kdim < dim; kdim++) { 
		    double lap_res_dadj_adj 			= 0.;
		    double adv_res_phiadj_nablauold_uadjold 	= 0.;
		    double adv_res_uold_nablaphiadj_uadjold 	= 0.;
	   for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		lap_res_dadj_adj 		             += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]][jdim]*phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim];
       }
	   for (unsigned jdim = 0; jdim < dim; jdim++) {
		adv_res_phiadj_nablauold_uadjold     += phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i] * gradSolVAR_qp[SolPdeIndex[jdim]][kdim] 			* SolVAR_qp[SolPdeIndex[jdim + adj_pos_begin]];
		adv_res_uold_nablaphiadj_uadjold     += SolVAR_qp[SolPdeIndex[jdim]]		       * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]];
	   }
	  Res[kdim + adj_pos_begin][i] += ( 
#if exact_sol_flag == 0
                            - alpha_val * target_flag * Vel_desired[kdim] 			      * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]
 #endif                                      
 #if exact_sol_flag == 1
                            - alpha_val * target_flag * exactVel_d[kdim] 			      * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]
 #endif
                                        + alpha_val * target_flag * SolVAR_qp[SolPdeIndex[kdim]] * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]
                                        - IRe*lap_res_dadj_adj
                                        - advection_flag * adv_res_phiadj_nablauold_uadjold
                                        - advection_flag * adv_res_uold_nablaphiadj_uadjold
                                        + SolVAR_qp[SolPdeIndex[press_type_pos + adj_pos_begin]] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim]
                                        ) * weight;
      }
      
//BLOCK delta_adjoint - state------------------------------------------------------------------------------------------
     for (unsigned j = 0; j < nDofsV; j++) {
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + adj_pos_begin][kdim][i*nDofsV + j] += ( - alpha_val * target_flag * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i] * phi_gss_fe[SolFEType[kdim]][j] 
                                                             + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i]    * phi_x_gss_fe[ SolFEType[kdim] ][j*dim_offset_grad /*space_dim*/ + kdim] 		* SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]   //c(delta_u, u_new, lambda_old)  diagonal blocks  ......phiadj_nablaunew_uadjold 
                                                             + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim] ][j] 			* phi_x_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	 //c(u_new, delta_u, lambda_old) diagonal blocks  ......unew_nablaphiadj_uadjold
                                                            ) * weight;
              unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
		Jac[kdim + adj_pos_begin][off_kdim][i*nDofsV + j] += (  advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i] * phi_x_gss_fe[ SolFEType[off_kdim] ][j*dim_offset_grad /*space_dim*/ + kdim]		      * SolVAR_qp[SolPdeIndex[off_kdim + adj_pos_begin]]   //c(delta_u, u_new, lambda_old)  off-diagonal blocks  ......phiadj_nablaunew_uadjold 
                                                              + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim] ][j] 		  * phi_x_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + off_kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	   //c(u_new, delta_u, lambda_old) off-diagonal blocks  ......unew_nablaphiadj_uadjold
                                                             ) * weight;
	  }            
     }//j_dadj_u loop


//DIAG BLOCK delta_adjoint - adjoint---------------------------------------------------------------------------------
     for (unsigned j = 0; j < nDofsVadj; j++) {
		    vector < double > lap_jac_dadj_adj(dim,0.);
		    vector < double > adv_uold_nablaphiadj_uadjnew(dim, 0.);
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { 
		  lap_jac_dadj_adj[kdim] += phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + jdim];
            }
          }
       for (unsigned  kdim = 0; kdim < dim; kdim++) { 
	   for (unsigned jdim = 0; jdim < dim; jdim++) { //diagonal blocks only
	     adv_uold_nablaphiadj_uadjnew[kdim]     += SolVAR_qp[SolPdeIndex[jdim]]  * phi_x_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j] ;
	   }
	  }
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + adj_pos_begin][kdim + adj_pos_begin][i*nDofsVadj + j] += (   IRe*lap_jac_dadj_adj[kdim] 
                                                                                + advection_flag * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim] * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j]   //c(delta_u, u_old, lambda_new)  diagonal blocks  ......phiadj_nablauold_uadjnew  
                                                                                + advection_flag * adv_uold_nablaphiadj_uadjnew[kdim] 	//c(u_old, delta_u, lambda_new)
                                                                               ) * weight;
               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
		  Jac[kdim + adj_pos_begin][off_kdim + adj_pos_begin][i*nDofsVadj + j] += ( advection_flag * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[off_kdim]][kdim] * phi_gss_fe[ SolFEType[off_kdim + adj_pos_begin] ][j]   //c(delta_u, u_old, lambda_new)  off-diagonal blocks  ......phiadj_nablauold_uadjnew   
                                                                                  ) * weight;
	  }
    } //j_dadj_adj loop
      
//BLOCK Pressure_adj
    for (unsigned j = 0; j < nDofsPadj; j++) {
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + adj_pos_begin][press_type_pos + adj_pos_begin][i*nDofsPadj + j] += -( phi_gss_fe[SolFEType[press_type_pos + adj_pos_begin]][j] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim] ) * weight;
	  }
    }//j_press_adj loop
  }//i_adj loop

//DIV_adj
  for (unsigned i = 0; i < nDofsPadj; i++) {
		double div_adj_dadj_qp = 0.;
      for (unsigned kdim = 0; kdim < dim; kdim++) {
	    div_adj_dadj_qp += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin ]][kdim] ;
      }
      Res[press_type_pos + adj_pos_begin][i] += ( (div_adj_dadj_qp) * phi_gss_fe[SolFEType[press_type_pos + adj_pos_begin]][i] ) * weight;
      for (unsigned j = 0; j < nDofsVadj; j++) {
        for (unsigned kdim = 0; kdim < dim; kdim++) {
            Jac[press_type_pos + adj_pos_begin][kdim + adj_pos_begin][i*nDofsVadj + j] += - ( phi_gss_fe[SolFEType[press_type_pos + adj_pos_begin]][i] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + kdim] ) * weight;
        }
      }//j loop
  }//i_div_adj

      //============ delta_adjoint row =============================================================================================

      //============ delta_control row ==================================================================================================
// delta_control
    for (unsigned kdim = 0; kdim < dim; kdim++) { 
         for (unsigned i = 0; i < nDofsGctrl; i++) {
       Res[kdim + ctrl_pos_begin][i] += - penalty_outside_control_boundary * ( (1 - control_node_flag[kdim][i]) * (  SolVAR_eldofs[SolPdeIndex[kdim + ctrl_pos_begin]][i] - 0.)  );              //enforce control zero outside the control boundary


// //DIAG BLOCK delta_control - control--------------------------------------------------------------------------------------
     for (unsigned j = 0; j < nDofsGctrl; j++) {
	    if (i==j) {
		Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsGctrl + j] += penalty_outside_control_boundary *(1 - control_node_flag[kdim][i]);              //enforce control zero outside the control boundary
                  } //end i==j
      }//j_dctrl_ctrl loop
  }//i_ctrl loop
      }  //kdim

 //============ delta_control row ==================================================================================================
 
 
      }  // end gauss point loop
      

    
    
// //  //--------------------PRINTING------------------------------------------------------------------------------------
// //  // Add the local Matrix/Vector into the global Matrix/Vector
// //   std::cout << " -------------------------Element = " << iel << " ----------------------Jac -------------------------- " << std::endl;      
// //   for(unsigned i_unk=0; i_unk < n_unknowns; i_unk++) {
// //     unsigned ndofs_unk_i = msh->GetElementDofNumber(iel, SolFEType[i_unk]);	//nDofs_V,P_of_st,adj,ctrl//Sol_n_el_dofs[k]
// //     std::cout << " ======== Row ==== " << i_unk << " Unk_i ===================================== " << std::endl;      
// //     for(unsigned j_unk=0; j_unk < n_unknowns; j_unk++) {
// // 	unsigned ndofs_unk_j = msh->GetElementDofNumber(iel, SolFEType[j_unk]);	//nDofs_V,P_of_st,adj,ctrl
// // 	std::cout << " ======= Column ==== " << j_unk << " Unk_j ================== " << std::endl;      
// // 	for (unsigned i = 0; i < ndofs_unk_i; i++) {
// // 	      std::cout << " " << std::setfill(' ') << std::setw(10) << Res[SolPdeIndex[i_unk]][ i ];
// // 	      std::cout << std::endl;
// // // // //             for (unsigned j = 0; j < ndofs_unk_j; j++) {
// // // // // 	    std::cout << " " << std::setfill(' ') << std::setw(10) << Jac[ SolPdeIndex[i_unk] ][ SolPdeIndex[j_unk] ][ i*ndofs_unk_i + j ];
// // // // // 	    }  //j end
// // // // // 	    std::cout << std::endl;
// // 	 } //i end
// //      } //j_unk end
// //    } //i_unk end
// //  //--------------------PRINTING------------------------------------------------------------------------------------

      //***************************************************************************************************************

    //Sum the local matrices/vectors into the Global Matrix/Vector
    // FIRST ALL THE BLOCKS WITHOUT THETA ROW OR COLUMN 
    for(unsigned i_unk=0; i_unk < n_unknowns-1; i_unk++) {
      RES->add_vector_blocked(Res[SolPdeIndex[i_unk]],JACDof[i_unk]);
        for(unsigned j_unk=0; j_unk < n_unknowns-1; j_unk++) {
	  if(assembleMatrix) JAC->add_matrix_blocked( Jac[ SolPdeIndex[i_unk] ][ SolPdeIndex[j_unk] ], JACDof[i_unk], JACDof[j_unk]);
        }
    }
    
    // THEN THE BLOCKS WITH THETA ROW OR COLUMN 
	/*delta_theta-theta*/    JAC->add_matrix_blocked( Jac[ SolPdeIndex[n_unknowns-1] ][ SolPdeIndex[n_unknowns-1] ], JACDof[n_unknowns-1], JACDof[n_unknowns-1]);
	    
     if (control_el_flag == 1) {
	      for (unsigned kdim = 0; kdim < dim; kdim++) {
                          /*delta_control*/       RES->add_vector_blocked(Res[SolPdeIndex[n_unknowns-2-kdim]],JACDof[n_unknowns-2-kdim]); 
		if(assembleMatrix) {
                          /*delta_theta-control*/ JAC->add_matrix_blocked( Jac[ SolPdeIndex[n_unknowns-1] ][ SolPdeIndex[n_unknowns-2-kdim] ], bdry_int_constr_pos_vec, JACDof[n_unknowns-2-kdim]);
                          /*delta_control-theta*/ JAC->add_matrix_blocked( Jac[ /*SolPdeIndex[n_unknowns-1] ][ SolPdeIndex[n_unknowns-2-kdim]*/SolPdeIndex[n_unknowns-2-kdim] ][ SolPdeIndex[n_unknowns-1] ], JACDof[n_unknowns-2-kdim], bdry_int_constr_pos_vec); 
		}
	      }  //kdim
     }  //add control boundary element contributions
     
     
          if (control_el_flag == 1) {
          /*delta_theta(bdry constr)*/         RES->add_vector_blocked(Res_outer,bdry_int_constr_pos_vec);
	  }
	  
     /* if (JACDof[n_unknowns-1][0] != bdry_int_constr_pos_vec[0]) */ /*delta_theta(fake)*/          RES->add_vector_blocked( Res[ SolPdeIndex[n_unknowns-1]],       JACDof[n_unknowns-1]);
	  
   //--------------------------------------------------------------------------------------------------------  
  } //end list of elements loop for each subdomain
  
  
  JAC->close();
  RES->close();
  
//     std::ostringstream mat_out; mat_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/matrix_non_ad" << mlPdeSys.GetNonlinearIt()  << ".txt";
//   JAC->print_matlab(mat_out.str(),"ascii");
//     std::ostringstream res_out; res_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/res_non_ad_" << mlPdeSys.GetNonlinearIt()  << ".txt";
//     std::filebuf res_fb;
//    res_fb.open (res_out.str().c_str(),std::ios::out);
//     std::ostream  res_file_stream(&res_fb);
//   RES->print(res_file_stream);
 
//   JAC->print();
//   RES->print();
  // ***************** END ASSEMBLY *******************
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
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //geometry *******************************
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  unsigned solType_coords = coordXType;  //FE_DOMAIN = 0; //we do linear FE this time // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
 
  CurrentElem < double > geom_element(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
  const unsigned int dim_offset_grad = /*dim*/  3  /*2*/    ;
 
  std::vector<double> normal(dim_offset_grad /*space_dim*/, 0.);

  vector < vector < double > > coordX(dim);    // local coordinates
  vector< vector < double> > coordX_bd(dim);

  for (unsigned  k = 0; k < dim; k++) { 
        coordX[k].reserve(maxSize);
        coordX_bd[k].reserve(maxSize); 
  }
  
  double weight = 0.;
  double weight_bd = 0.;
  
  
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
    solV[k].reserve(maxSize);
  }

  
  vector <double> phiV_gss;  // local test function
  vector <double> phiV_x_gss; // local test function first order partial derivatives
  vector <double> phiV_xx_gss; // local test function second order partial derivatives

  phiV_gss.reserve(maxSize);
  phiV_x_gss.reserve(maxSize * dim_offset_grad /*space_dim*/);
  phiV_xx_gss.reserve(maxSize * dim2);
  
  
  //velocity *******************************
   

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
    solVctrl[k].reserve(maxSize);
  }

  
  vector <double> phiVctrl_gss_bd;  // local test function
  vector <double> phiVctrl_x_gss_bd; // local test function first order partial derivatives
  vector <double> phiVctrl_xx_gss_bd; // local test function second order partial derivatives

  phiVctrl_gss_bd.reserve(maxSize);
  phiVctrl_x_gss_bd.reserve(maxSize * dim_offset_grad /*space_dim*/);
  phiVctrl_xx_gss_bd.reserve(maxSize * dim2);
  
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
  vector <double> phiVdes_xx_gss; // local test function second order partial derivatives

  phiVdes_gss.reserve(maxSize);
  phiVdes_x_gss.reserve(maxSize * dim_offset_grad /*space_dim*/);
  phiVdes_xx_gss.reserve(maxSize * dim2);

//   vector< vector < double > >  solVdes(dim);    // local solution
  vector <double>  solVdes(dim,0.);
  vector<double> Vdes_gss(dim, 0.);  
  
//  for (unsigned  k = 0; k < dim; k++) {
//     solVdes[k].reserve(maxSize);
//   }
//   
//   double* Vdes_gss [3] = Vel_desired/*= 0.*/;


// Vel_desired##################################################################

  
vector<double> integral(dim);

double  integral_target_alpha = 0.;

double	integral_beta   = 0.;
double	integral_gamma  = 0.;

double integral_g_dot_n = 0.;
  
//*************************************************** 
     std::vector < std::vector < double > >  JacI_qp(space_dim);
     std::vector < std::vector < double > >  Jac_qp(dim);
    for (unsigned d = 0; d < Jac_qp.size(); d++) {   Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp.size(); d++) { JacI_qp[d].resize(dim); }
    
    double detJac_qp;

     std::vector < std::vector < double > >  JacI_qp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_qp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_qp_bdry.size(); d++) {   Jac_qp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp_bdry.size(); d++) { JacI_qp_bdry[d].resize(dim-1); }
    
    double detJac_qp_bdry;
    
    //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
//*************************************************** 

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

   // geometry *****************************
      geom_element.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
      const short unsigned ielGeom = geom_element.geom_type();
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
  geom_element.set_elem_center(iel, solType_coords);

   int target_flag = 0;
   target_flag = ElementTargetFlag(geom_element.get_elem_center());
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
        solVdes[k]/*[i]*/ = Vel_desired[k] /*(*sol->_Sol[solVIndex[k]])(solVdesDof)*/;      // global extraction and local storage for the solution
     }
//     }
 //DESIRED VEL###################################################################

 
 //************ set control flag *********************
  int control_el_flag = 0;
        control_el_flag = ControlDomainFlag_bdry(geom_element.get_elem_center());
  std::vector< std::vector<int> > control_node_flag(dim);
	    for(unsigned idim=0; idim<dim; idim++) {
	      control_node_flag[idim].resize(nDofsVctrl);
   /*if (control_el_flag == 0)*/ std::fill(control_node_flag[idim].begin(), control_node_flag[idim].end(), 0);
	    }
 //*************************************************** 

//========BoundaryLoop=====================================================================

  // Perform face loop over elements that contain some control face
  if (control_el_flag == 1) {
	  
    double tau=0.;
    vector<double> normal(dim_offset_grad /*space_dim*/,0);
	  
    for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {

       const unsigned ielGeom_bd = msh->GetElementFaceType(iel, jface);    
       const unsigned nve_bd = msh->GetElementFaceDofNumber(iel,jface,solType_coords);
       
       geom_element.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);
 
       geom_element.set_elem_center_bdry_3d();

	    // look for boundary faces
            const int bdry_index = el->GetFaceElementIndex(iel,jface);
            
	    if( bdry_index < 0) {
	   unsigned int face_in_rectangle_domain = -( msh->el->GetFaceElementIndex(iel,jface)+1);

	   if(  face_in_rectangle_domain == FACE_FOR_CONTROL) { //control face
// //=================================================== 
// 		//we use the dirichlet flag to say: if dirichlet = true, we set 1 on the diagonal. if dirichlet = false, we put the boundary equation
// 	    std::vector<bool> dir_bool; dir_bool.resize(dim);
// 	    for(unsigned idim=0; idim<dim; idim++) {
// 		dir_bool[idim] = ml_sol->GetBdcFunction()(xyz_bdc,ctrl_name[idim].c_str(),tau,face,0.);
// 	    }
	  
//=================================================== 
		
//========= initialize gauss quantities on the boundary ============================================
    vector < double >   Vctrl_bd_qp(dim, 0.);    //  solution@bdry
    vector < vector < double > > gradVctrl_bd_qp(dim);
      for (unsigned  k = 0; k < dim; k++) {
          gradVctrl_bd_qp[k].resize(dim_offset_grad /*space_dim*/);
          std::fill(gradVctrl_bd_qp[k].begin(), gradVctrl_bd_qp[k].end(), 0);
        }

//========= gauss_loop boundary===============================================================
	    for(unsigned ig_bd=0; ig_bd < ml_prob.GetQuadratureRule(ielGeom_bd).GetGaussPointsNumber(); ig_bd++) {
    elem_all[ielGeom_bd][solType_coords]->JacJacInv(geom_element.get_coords_at_dofs_bdry_3d(), ig_bd, Jac_qp_bdry, JacI_qp_bdry, detJac_qp_bdry, space_dim);
	elem_all[ielGeom_bd][solType_coords]->compute_normal(Jac_qp_bdry, normal);
    
    weight_bd = detJac_qp_bdry * ml_prob.GetQuadratureRule(ielGeom_bd).GetGaussWeightsPointer()[ig_bd];

    elem_all[ielGeom_bd][solVctrlType] ->shape_funcs_current_elem(ig_bd, JacI_qp_bdry,phiVctrl_gss_bd,phiVctrl_x_gss_bd , phiVctrl_xx_gss_bd , space_dim);

     
//========== compute gauss quantities on the boundary ===============================================
    for (unsigned  k = 0; k < dim; k++) {
	  Vctrl_bd_qp[k] = 0.;
	  for(unsigned ivar2=0; ivar2<dim_offset_grad /*space_dim*/; ivar2++) { gradVctrl_bd_qp[k][ivar2] = 0.; }
	  
	  for (unsigned i = 0; i < nDofsVctrl; i++) {
		   for(int i_bd = 0; i_bd < nve_bd; i_bd++) {
		       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bd);
		       Vctrl_bd_qp[k] += phiVctrl_gss_bd[i_bd] * solVctrl[k][i_vol];
		       for(unsigned ivar2=0; ivar2<dim_offset_grad /*space_dim*/; ivar2++) {
			   gradVctrl_bd_qp[k][ivar2] += phiVctrl_x_gss_bd[i_bd * dim_offset_grad /*space_dim*/ + ivar2 /*i_bd + ivar2 * nve_bd*/] * solVctrl[k][i_vol]; 
		       }
		   }
	}
    }
 //end unknowns eval at gauss points ********************************
		  
//========== compute gauss quantities on the boundary ================================================
      for (unsigned  k = 0; k < dim; k++) {
	 integral_beta	+= ((Vctrl_bd_qp[k])*(Vctrl_bd_qp[k])*weight_bd);
	 integral_g_dot_n += Vctrl_bd_qp[k]*normal[k]*weight_bd;
      }
      for (unsigned  k = 0; k < dim; k++) {
	for (unsigned  j = 0; j < dim; j++) {	
		integral_gamma	  += ((gradVctrl_bd_qp[k][j])*(gradVctrl_bd_qp[k][j])*weight_bd);
	}
      }


                }  //end ig_bdry loop
	  
             }    //end if control face
	 }  //end if boundary faces
      }  // loop over element faces //jface   
  } //end if control element flag

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {
//STATE######## VolumeLoop #####################################################################	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[ielGeom][solType_coords]->JacJacInv(geom_element.get_coords_at_dofs_3d(), ig, Jac_qp, JacI_qp, detJac_qp, space_dim);
    weight = detJac_qp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[ig];
   
    elem_all[ielGeom][solVType]->shape_funcs_current_elem(ig, JacI_qp, phiV_gss, phiV_x_gss, phiV_xx_gss , space_dim);
    elem_all[ielGeom][solVType /*solVdes*/]->shape_funcs_current_elem(ig, JacI_qp, phiVdes_gss, phiVdes_x_gss, phiVdes_xx_gss , space_dim);

    
      for (unsigned  k = 0; k < dim; k++) {
	           V_gss[k] = 0.;
	           Vdes_gss[k] = 0.;
	    for (unsigned i = 0; i < nDofsV; i++) {
	   	   V_gss[k] += solV[k][i] * phiV_gss[i];
		   Vdes_gss[k] += solVdes[k]/*[i]*/ * phiVdes_gss[i];
	    }
	  }
	

      for (unsigned  k = 0; k < dim; k++) {
	      integral_target_alpha+=(( target_flag ) *((V_gss[k]  - Vdes_gss[k]) * (V_gss[k]  - Vdes_gss[k]))*weight);
      }
      
      }// end gauss point loop
    } //end element loop  

       std::ostringstream filename_out; filename_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "Integral_computation"  << ".txt";

       std::ofstream intgr_fstream;
  if (paral::get_rank() == 0 ) {
      intgr_fstream.open(filename_out.str().c_str(),std::ios_base::app);
      intgr_fstream << " ***************************** Non Linear Iteration "<< mlPdeSys.GetNonlinearIt() << " *********************************** " <<  std::endl << std::endl;
      intgr_fstream << "The value of the target functional for " << "alpha " <<   std::setprecision(0) << std::scientific << alpha_val << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_target_alpha << std::endl;
      intgr_fstream << "The value of the L2 control for        " << "beta  " <<   std::setprecision(0) << std::scientific << beta_val  << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_beta         << std::endl;
      intgr_fstream << "The value of the H1 control for        " << "gamma " <<   std::setprecision(0) << std::scientific << gamma_val << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_gamma        << std::endl;
      intgr_fstream << "The value of the integral of g.n "<<    integral_g_dot_n << std::endl;
      intgr_fstream << "The value of the theta is                             " <<    std::setw(11) << std::setprecision(10) <<  solTheta << std::endl;
      intgr_fstream << "The value of the total integral is " << std::setw(11) << std::setprecision(10) <<  integral_target_alpha * alpha_val*0.5  + integral_beta *beta_val*0.5 + integral_gamma *gamma_val*0.5 << std::endl;
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
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //geometry *******************************
  vector < vector < double > > coordX(dim);    // local coordinates

  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) { 
    coordX[k].reserve(maxSize);
  }
   
  //geometry *******************************

 // solution variables *******************************************
  const int n_vars = dim+1;
  const int n_unknowns = 3*n_vars; //(2.*dim)+1; //state , adjoint of velocity terms and one pressure term
  const int vel_type_pos = 0;
  const int press_type_pos = dim;
  const int adj_vel_type_pos = vel_type_pos;
  const int state_pos_begin = 0;
  const int adj_pos_begin   = dim+1;
  const int ctrl_pos_begin   = 2*(dim+1);
  const int theta_index = press_type_pos + ctrl_pos_begin;
  
  vector < std::string > Solname(n_unknowns);  // const char Solname[4][8] = {"U","V","W","P"};
  Solname              [state_pos_begin+0] =                "U";
  Solname              [state_pos_begin+1] =                "V";
  if (dim == 3) Solname[state_pos_begin+2] =                "W";
  Solname              [state_pos_begin + press_type_pos] = "P";
  
  Solname              [adj_pos_begin + 0] =              "UADJ";
  Solname              [adj_pos_begin + 1] =              "VADJ";
  if (dim == 3) Solname[adj_pos_begin + 2] =              "WADJ";
  Solname              [adj_pos_begin + press_type_pos] = "PADJ";

  Solname              [ctrl_pos_begin + 0] =              "GX";
  Solname              [ctrl_pos_begin + 1] =              "GY";
  if (dim == 3) Solname[ctrl_pos_begin + 2] =              "GZ";
  Solname              [ctrl_pos_begin + press_type_pos] = "THETA";
  
  vector < unsigned > SolIndex(n_unknowns);  
  vector < unsigned > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolIndex[ivar]	= ml_sol->GetIndex        (Solname[ivar].c_str());
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);
  }

  vector < double > Sol_n_el_dofs(n_unknowns);
  
  //==========================================================================================
  // velocity ************************************
  vector < vector < double > > phi_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_x_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_xx_gss_fe(NFE_FAMS);
 
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_gss_fe[fe].reserve(maxSize);
      phi_x_gss_fe[fe].reserve(maxSize*dim);
     phi_xx_gss_fe[fe].reserve(maxSize*(3*(dim-1)));
   }
  
  //=================================================================================================
  
  // quadratures ********************************
  double weight;
  
  
  //----------- dofs ------------------------------
  vector < vector < double > > SolVAR_eldofs(n_unknowns);
  vector < vector < double > > gradSolVAR_eldofs(n_unknowns);
  
  vector < vector < double > > SolVAR_coarser_prol_eldofs(n_unknowns);
  vector < vector < double > > gradSolVAR_coarser_prol_eldofs(n_unknowns);


  for(int k = 0; k < n_unknowns; k++) {
    SolVAR_eldofs[k].reserve(maxSize);
    gradSolVAR_eldofs[k].reserve(maxSize*dim); 
    
    SolVAR_coarser_prol_eldofs[k].reserve(maxSize);
    gradSolVAR_coarser_prol_eldofs[k].reserve(maxSize*dim);    
  }

  //------------ at quadrature points ---------------------
  vector < double > SolVAR_qp(n_unknowns);
  vector < double > SolVAR_coarser_prol_qp(n_unknowns);
  vector < vector < double > > gradSolVAR_qp(n_unknowns);
  vector < vector < double > > gradSolVAR_coarser_prol_qp(n_unknowns);
  for(int k = 0; k < n_unknowns; k++) {
      gradSolVAR_qp[k].reserve(maxSize);  
      gradSolVAR_coarser_prol_qp[k].reserve(maxSize);  
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

   unsigned nDofsGctrl = msh->GetElementDofNumber(iel,SolFEType[ctrl_pos_begin]);
    unsigned nDofsThetactrl = msh->GetElementDofNumber(iel,SolFEType[theta_index] );

    unsigned nDofsVP = dim * nDofsV + nDofsP; //same for state and adjoint
    unsigned nDofsVPctrl = dim * nDofsGctrl + nDofsThetactrl; //control
   
    unsigned nDofsVP_tot = 2*nDofsVP + (nDofsVPctrl);
  // equation end *****************************


   //STATE###################################################################  
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
	Sol_n_el_dofs[k]=ndofs_unk;
       SolVAR_eldofs[k].resize(ndofs_unk);
       SolVAR_coarser_prol_eldofs[k].resize(ndofs_unk);
    for (unsigned i = 0; i < ndofs_unk; i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);    // global to global mapping between solution node and solution dof // via local to global solution node
       SolVAR_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);      // global extraction and local storage for the solution
       SolVAR_coarser_prol_eldofs[k][i] = (*sol_coarser_prolongated->_Sol[SolIndex[k]])(solDof);      // global extraction and local storage for the solution
      }
    }
  //CTRL###################################################################

 
      // ********************** Gauss point loop *******************************
      for(unsigned ig=0;ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {
	
 
      for(int fe=0; fe < NFE_FAMS; fe++) {
	msh->_finiteElement[ielGeom][fe]->Jacobian(coordX,ig,weight,phi_gss_fe[fe],phi_x_gss_fe[fe],phi_xx_gss_fe[fe]);
      }
         //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
  	msh->_finiteElement[ielGeom][BIQUADR_FE]->Jacobian(coordX,ig,weight,phi_gss_fe[BIQUADR_FE],phi_x_gss_fe[BIQUADR_FE],phi_xx_gss_fe[BIQUADR_FE]);

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
	    SolVAR_qp[unk] += phi_gss_fe[ SolFEType[unk] ][i] * SolVAR_eldofs[unk][i];
	    SolVAR_coarser_prol_qp[unk] += phi_gss_fe[ SolFEType[unk] ][i] * SolVAR_coarser_prol_eldofs[unk][i];
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      gradSolVAR_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType[unk] ][i*dim+ivar2] * SolVAR_eldofs[unk][i]; 
	      gradSolVAR_coarser_prol_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType[unk] ][i*dim+ivar2] * SolVAR_coarser_prol_eldofs[unk][i]; 
	    }
	  }
	  
	}  
 //end unknowns eval at gauss points ********************************


	for(unsigned unk = 0; unk < n_unknowns; unk++) {
        l2norm[unk] += ( SolVAR_qp[unk] - SolVAR_coarser_prol_qp[unk] ) * ( SolVAR_qp[unk] - SolVAR_coarser_prol_qp[unk] ) * weight ; 
        
     }
    
    
	for(unsigned unk = 0; unk < dim; unk++) {
        for(int j = 0; j < dim; j++){
        seminorm[unk] += (gradSolVAR_qp[unk][j] - gradSolVAR_coarser_prol_qp[unk][j] ) * ( gradSolVAR_qp[unk][j] - gradSolVAR_coarser_prol_qp[unk][j] ) * weight ;
        seminorm[unk + dim] += (gradSolVAR_qp[unk + adj_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + adj_pos_begin][j] ) * ( gradSolVAR_qp[unk + adj_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + adj_pos_begin][j] ) * weight ;
        seminorm[unk + 2*dim] += (gradSolVAR_qp[unk + ctrl_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j] ) * ( gradSolVAR_qp[unk + ctrl_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j] ) * weight ;
        
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
