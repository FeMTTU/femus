#include "MultiLevelProblem.hpp"
#include "TransientSystem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"
#include "FemusInit.hpp"
#include "SparseMatrix.hpp"
#include "VTKWriter.hpp"
#include "FElemTypeEnum.hpp"
#include "Files.hpp"
#include "../include/FSIassembly.hpp"
#include "../include/IncompressibleFSIAssembly.hpp"


using std::cout;
using std::endl;
using namespace femus;

// User defined functions
//int AssembleMatrixResFSI(NonLinearMultiLevelProblem &nl_td_ml_prob, unsigned level, const unsigned &gridn, const unsigned &ipde, const bool &assemble_matrix);

double SetVariableTimeStep(const double time);

bool SetBoundaryCondition(const double &x, const double &y, const double &z,const char name[], 
			  double &value, const int FaceName, const double = 0.);

bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber,const int &level);

//------------------------------------------------------------------------------------------------------------------

int main(int argc,char **args) {
    
  bool Vanka=0, Gmres=0, Asm=0;
  if(argc >= 2) {
    if( !strcmp("vanka",args[1])) 	Vanka=1;
    else if( !strcmp("gmres",args[1])) 	Gmres=1;
    else if( !strcmp("asm",args[1])) 	Asm=1;
    
    if(Vanka+Gmres+Asm==0) {
      cout << "wrong input arguments!" << endl;
      exit(0);
    }
  }
  else {
    cout << "No input argument set default smoother = Asm" << endl;
    Asm=1;
  }
  
   
  /// Init Petsc-MPI communicator
  FemusInit mpinit(argc,args,MPI_COMM_WORLD);

  Files files; 
        files.CheckIODirectories();
        files.RedirectCout();

  unsigned short nm,nr;
  std::cout<<"#MULTIGRID levels? (>=1) \n";
  //std::cin>>nm;
  nm=4;

  std::cout<<"#MAX_REFINEMENT levels? (>=0) \n";
  //std::cin>>nr;
  nr=0;
  int tmp=nm;
  nm+=nr;
  nr=tmp;

  char *infile = new char [50];
  
  //input file
  sprintf(infile,"./input/mesh.comsolbenchmark.neu");

  const double Lref = 1.;
  double Uref = 1.;
  double rhof = 1000.;
  double muf = 0.001;
  double rhos = 7850;
  double ni = 0.5;
  double E = 200000;
  
//   sprintf(infile,"./input/fsifirst.neu");
// 
//   double Lref = 1.;
//   double Uref = 1.;
//   double rhof = 1000.;
//   double muf = 1.;
//   double rhos = 1000;
//   double ni = 0.4;
//   double E = 1400000;
  
  MultiLevelMesh ml_msh(nm,nr,infile,"fifth",Lref,SetRefinementFlag);
  
  MultiLevelSolution ml_sol(&ml_msh);
  
  //Start System Variables
  ml_sol.AddSolution("DX",LAGRANGE,SECOND,1);
  ml_sol.AddSolution("DY",LAGRANGE,SECOND,1);
  ml_sol.AssociatePropertyToSolution("DX","Displacement"); // Add this line
  ml_sol.AssociatePropertyToSolution("DY","Displacement"); // Add this line 
  ml_sol.AddSolution("U",LAGRANGE,SECOND,1);
  ml_sol.AddSolution("V",LAGRANGE,SECOND,1);
  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  ml_sol.AddSolution("P",DISCONTINOUS_POLYNOMIAL,FIRST,1);
  ml_sol.AssociatePropertyToSolution("P","Pressure"); // Add this line

  //Initialize (update Init(...) function)
  ml_sol.Initialize("All");

  //Set Boundary (update Dirichlet(...) function)
  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  ml_sol.GenerateBdc("DX","Steady");
  ml_sol.GenerateBdc("DY","Steady");
  ml_sol.GenerateBdc("U","Steady");
  ml_sol.GenerateBdc("V","Steady");
  ml_sol.GenerateBdc("P","Steady");
  
  MultiLevelProblem ml_prob(&ml_sol);

  Parameter par(Lref,Uref);
  
  // Generate Solid Object  
  //Solid solid(par,E,ni,rhos,"Linear_elastic");
  //Solid solid(par,E,ni,rhos,"Neo-Hookean");
  //Solid solid(par,E,ni,rhos,"Neo-Hookean");
  Solid solid(par,E,ni,rhos,"Mooney-Rivlin");
  //Solid solid(par,E,ni,rhos,"Neo-Hookean-BW-Penalty");
  
  
  cout << "Solid properties: " << endl;
  cout << solid << endl;
  
  // Generate Fluid Object
  Fluid fluid(par,muf,rhof,"Newtonian");
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;

  // Add fluid object
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;
  
  // Add Solid Object
  ml_prob.parameters.set<Solid>("Solid") = solid;

  ml_msh.MarkStructureNode();
   
  //create systems
  // add the system FSI to the MultiLevel problem
  MonolithicFSINonLinearImplicitSystem & system = ml_prob.add_system<MonolithicFSINonLinearImplicitSystem> ("Fluid-Structure-Interaction");
  system.AddSolutionToSystemPDE("DX");
  system.AddSolutionToSystemPDE("DY");
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  system.AddSolutionToSystemPDE("P");
  
//   bool sparsity_pattern_matrix[5][5]={{1, 0, 1, 0, 0},
// 				      {0, 1, 0, 1, 0},
// 				      {1, 1, 1, 1, 1},
// 				      {1, 1, 1, 1, 1},
// 				      {1, 1, 0, 0, 1}};
//   vector < bool > sparsity_pattern (sparsity_pattern_matrix[0],sparsity_pattern_matrix[0]+25*sizeof(bool));
//   system.SetSparsityPattern(sparsity_pattern);  
   
  // System Fluid-Structure-Interaction
  system.SetAssembleFunction(IncompressibleFSIAssemblyAD_DD);  
  
  system.SetMaxNumberOfLinearIterations(2);
  system.SetMaxNumberOfNonLinearIterations(10);
  system.SetMgType(F_CYCLE);
  system.SetAbsoluteConvergenceTolerance(1.e-10);
  system.SetNonLinearConvergenceTolerance(1.e-05);
  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);
   
  //Set Smoother Options
  if(Gmres) 		system.SetMgSmoother(GMRES_SMOOTHER);
  else if(Asm) 		system.SetMgSmoother(ASM_SMOOTHER);
  else if(Vanka)	system.SetMgSmoother(VANKA_SMOOTHER);
  
  // init all the systems
  system.init();
  
  
  system.SetSolverFineGrids(GMRES);
  system.SetPreconditionerFineGrids(ILU_PRECOND); 
  system.SetTolerances(1.e-12,1.e-20,1.e+50,10);
  
  
//   system.SetSolverFineGrids(GMRES);
//   system.SetPreconditionerFineGrids(LU_PRECOND); 
//   system.SetTolerances(1.e-12,1.e-20,1.e+50,20);
  
  system.ClearVariablesToBeSolved();
  //system.AddVariableToBeSolved("All");
  system.AddVariableToBeSolved("DX");
  system.AddVariableToBeSolved("DY");
  system.AddVariableToBeSolved("U");
  system.AddVariableToBeSolved("V");
  system.AddVariableToBeSolved("P");
  
  //for Vanka and ASM smoothers
  system.SetNumberOfSchurVariables(1);
  system.SetElementBlockNumber(3);   
  //for Gmres smoother
  //system.SetDirichletBCsHandling(PENALTY); 
  system.SetDirichletBCsHandling(ELIMINATION);   
   
  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  VTKWriter vtkio(ml_sol);
  vtkio.SetMovingMesh(mov_vars);
  
  // Solving Fluid-Structure-Interaction system
  std::cout << std::endl;
  std::cout << " *********** Fluid-Structure-Interaction ************  " << std::endl;
  system.solve();
   
  //print solution 
  std::vector<std::string> print_vars;
  print_vars.push_back("DX");
  print_vars.push_back("DY");
  print_vars.push_back("U");
  print_vars.push_back("V");
  print_vars.push_back("P");
      
  vtkio.write_system_solutions(files.GetOutputPath(),"biquadratic",print_vars);

  // Destroy all the new systems
  ml_prob.clear();
   
  delete [] infile;
  return 0;
}











// int main(int argc,char **args) {
//   
//   /// Init Petsc-MPI communicator
//   FemusInit mpinit(argc,args,MPI_COMM_WORLD);
// 
//   unsigned short nm,nr;
//   std::cout<<"#MULTIGRID levels? (>=1) \n";
//   //std::cin>>nm;
//   nm=1;
// 
//   std::cout<<"#MAX_REFINEMENT levels? (>=0) \n";
//   //std::cin>>nr;
//   nr=0;
//   int tmp=nm;
//   nm+=nr;
//   nr=tmp;
// 
//   char *infile = new char [50];
// 
//   //input file
//   sprintf(infile,"./input/mesh.comsolbenchmark.neu");
// 
//   const double Lref = 1.;
//   double Uref = 1.;
//   double rhof = 1000.;
//   double muf = 0.001;
//   double rhos = 7850;
//   double ni = 0.3;
//   double E = 200000;
//   
//   MultiLevelMesh ml_msh(nm,nr,infile,"fifth",Lref,SetRefinementFlag);
//   
//   MultiLevelSolution ml_sol(&ml_msh);
//   
//    //Start System Variables
//   ml_sol.AddSolution("DX",LAGRANGE,SECOND,2);
//   ml_sol.AddSolution("DY",LAGRANGE,SECOND,2);
//   ml_sol.AssociatePropertyToSolution("DX","Displacement"); // Add this line
//   ml_sol.AssociatePropertyToSolution("DY","Displacement"); // Add this line 
//   ml_sol.AddSolution("U",LAGRANGE,SECOND,2);
//   ml_sol.AddSolution("V",LAGRANGE,SECOND,2);
//   ml_sol.AddSolution("AX",LAGRANGE,SECOND,1,0);
//   ml_sol.AddSolution("AY",LAGRANGE,SECOND,1,0);
//   // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
//   ml_sol.AddSolution("P",DISCONTINOUS_POLYNOMIAL,FIRST,1);
//   ml_sol.AssociatePropertyToSolution("P","Pressure"); // Add this line
// 
//   //Initialize (update Init(...) function)
//   ml_sol.Initialize("All");
// 
//   //Set Boundary (update Dirichlet(...) function)
//   ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
//   ml_sol.GenerateBdc("DX","Steady");
//   ml_sol.GenerateBdc("DY","Steady");
//   ml_sol.GenerateBdc("U","Time_dependent");
//   ml_sol.GenerateBdc("V","Steady");
//   ml_sol.GenerateBdc("AX","Steady");
//   ml_sol.GenerateBdc("AY","Steady");
//   ml_sol.GenerateBdc("P","Steady");
//   
//   
//   
//   
//   MultiLevelProblem ml_probl(&ml_msh, &ml_sol);
//   
// 
//   Parameter par(Lref,Uref);
//   
//   // Generate Solid Object
//   Solid solid(par,E,ni,rhos,"Neo-Hookean");
//   cout << "Solid properties: " << endl;
//   cout << solid << endl;
//   
//   //Generate Fluid Object
//   Fluid fluid(par,muf,rhof,"Newtonian");
//   cout << "Fluid properties: " << endl;
//   cout << fluid << endl;
// 
//   // Add fluid object
//   ml_probl.parameters.set<Fluid>("Fluid") = fluid;
//   
//   // Add Solid Object
//   ml_probl.parameters.set<Solid>("Solid") = solid;
// 
//  
//   
// //   std::vector<std::string> mov_vars;
// //   mov_vars.push_back("DX");
// //   mov_vars.push_back("DY");
// //   ml_probl.SetMovingMesh(mov_vars);
//   ml_msh.MarkStructureNode();
//   
//   //create systems
//   // add the system FSI to the MultiLevel problem
//   TransientMonolithicFSINonlinearImplicitSystem & system = ml_probl.add_system<TransientMonolithicFSINonlinearImplicitSystem> ("Fluid-Structure-Interaction",VANKA_SMOOTHER);
//   system.AddSolutionToSystemPDE("DX");
//   system.AddSolutionToSystemPDE("DY");
//   system.AddSolutionToSystemPDE("U");
//   system.AddSolutionToSystemPDE("V");
//   system.AddSolutionToSystemPDE("P");
//   
//   // init all the systems
//   system.init();
//  
//   // System Navier-Stokes
//   system.SetAssembleFunction(AssembleMatrixResFSI);  
//   system.SetMaxNumberOfLinearIterations(1);
//   system.SetAbsoluteConvergenceTolerance(1.e-8);  
//   system.SetMaxNumberOfNonLinearIterations(5);  
//   system.SetNonLinearConvergenceTolerance(1.e-4);
//   system.SetDirichletBCsHandling(PENALTY);
//   
//   system.SetMgType(F_CYCLE);
//   system.SetNumberPreSmoothingStep(0);
//   system.SetNumberPostSmoothingStep(3);
//   //system.SetMgSmoother(VANKA_SMOOTHER);
//   system.ClearVariablesToBeSolved();
//   system.AddVariableToBeSolved("All");
//   /*system.AddVariableToVankaIndex("DY");
//   system.AddVariableToVankaIndex("U");
//   system.AddVariableToVankaIndex("V");
//   system.AddVariableToVankaIndex("P");*/
//   system.SetSolverFineGrids(GMRES);
//   system.SetPreconditionerFineGrids(MLU_PRECOND); // Use MLU_PRECOND (MumpsLU) instead LU_PRECOND (pestcLu), it's more robust!!!
//   system.SetNumberOfSchurVariables(1);
//   system.SetTolerances(1.e-12,1.e-20,1.e+50,1);
//   //system.SetSchurTolerances(1.e-12,1.e-20,1.e+50,4);
//   system.SetElementBlockNumber(4);                
// 
//   // time loop parameter
//   system.SetIntervalTime(0.005);
//   system.AttachGetTimeIntervalFunction(SetVariableTimeStep);
//   const unsigned int n_timesteps = 30;
//   const unsigned int write_interval = 1;
//   
//   std::vector<std::string> mov_vars;
//   mov_vars.push_back("DX");
//   mov_vars.push_back("DY");
// //   ml_probl.SetMovingMesh(mov_vars);
//   VTKWriter vtkio(ml_sol);
//   vtkio.SetMovingMesh(mov_vars);
//   
//   for (unsigned time_step = 0; time_step < n_timesteps; time_step++) {
//    
//     // Solving Fluid-Structure-Interaction system
//     std::cout << std::endl;
//     std::cout << " *********** Fluid-Structure-Interaction ************  " << std::endl;
//     system.solve();
//    
//     //The update of the acceleration must be done before the update of the other variables
//     system.NewmarkAccUpdate();
//     
//     //update Solution
//     system.UpdateSolution();
// 
//     // print solution
//     if ( !(time_step%write_interval) ) {
//         
//       //print solution 
//       std::vector<std::string> print_vars;
//       print_vars.push_back("DX");
//       print_vars.push_back("DY");
//       print_vars.push_back("U");
//       print_vars.push_back("V");
//       print_vars.push_back("P");
//       
// //       ml_probl.printsol_vtu_inline("biquadratic",print_vars,time_step);
//       vtkio.write_system_solutions("biquadratic",print_vars,time_step);
//     }
//   
//   } //end loop timestep
//   
// 
//   // Destroy all the new systems
//   ml_probl.clear();
//    
//   delete[] infile;
//   return 0;
// }

//---------------------------------------------------------------------------------------------------------------------

bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber,const int &level) {
  bool refine=0;

  //refinemenet based on Elemen Group Number
  if (ElemGroupNumber==5) refine=1;
  if (ElemGroupNumber==6) refine=1;

  return refine;

}

//---------------------------------------------------------------------------------------------------------------------

double SetVariableTimeStep(const double time) {
  if (time < 0.75) {
    return 0.005;
  }
  else {
    return 0.25;
  }
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryCondition(const double &x, const double &y, const double &z,const char name[], double &value, const int FaceName, const double time) {
  bool test=1; //Dirichlet
  value=0.;
  //   cout << "Time bdc : " <<  time << endl;
  if (!strcmp(name,"U")) {
    if (1==FaceName) { //inflow
      test=1;
      //comsol Benchmark
      //value = (0.05*time*time)/(sqrt( (0.04 - time*time)*(0.04 - time*time) + (0.1*time)*(0.1*time) ))*y*(0.0001-y)*4.*100000000;
      value = 0.05*y*(0.0001-y)*4.*100000000;
    }
    else if (2==FaceName ) {  //outflow
      test=0;
      //    test=1;
      value=0.;
    }
    else if (3==FaceName ) {  // no-slip fluid wall
      test=1;
      value=0.;
    }
    else if (4==FaceName ) {  // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if (!strcmp(name,"V")) {
    if (1==FaceName) {          //inflow
      test=1;
      value=0.;
    }
    else if (2==FaceName ) {    //outflow
      test=0;
      //    test=1;
      value=0.;
    }
    else if (3==FaceName ) {    // no-slip fluid wall
      test=1;
      value=0;
    }
    else if (4==FaceName ) {    // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if (!strcmp(name,"W")) {
    if (1==FaceName) {
      test=1;
      value=0.;
    }
    else if (2==FaceName ) {
      test=1;
      value=0.;
    }
    else if (3==FaceName ) {
      test=1;
      value=0.;
    }
    else if (4==FaceName ) {
      test=1;
      value=0.;
    }
  }
  else if (!strcmp(name,"P")) {
    if (1==FaceName) {
      test=0;
      value=0.;
    }
    else if (2==FaceName ) {
      test=0;
      value=0.;
    }
    else if (3==FaceName ) {
      test=0;
      value=0.;
    }
    else if (4==FaceName ) {
      test=0;
      value=0.;
    }
  }
  else if (!strcmp(name,"DX")) {
    if (1==FaceName) {       //inflow
      test=1;
      value=0.;
    }
    else if (2==FaceName ) { //outflow
      test=1;
      value=0.;
    }
    else if (3==FaceName ) { // no-slip Top fluid wall
      test=0;
      value=0;
    }
    else if (4==FaceName ) { // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if (!strcmp(name,"DY")) {
    if (1==FaceName) {       //inflow
      test=0;
      value=0.;
    }
    else if (2==FaceName ) { //outflow
      test=0;
      value=0.;
    }
    else if (3==FaceName ) { // no-slip fluid wall
      test=1;
      value=0.;
    }
    else if (4==FaceName ) { // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if (!strcmp(name,"DZ")) {
    if (1==FaceName) {       //inflow
      test=1;
      value=0.;
    }
    else if (2==FaceName ) { //outflow
      test=1;
      value=0.;
    }
    else if (3==FaceName ) { // no-slip fluid wall
      test=1;
      value=0.;
    }
    else if (4==FaceName ) { // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if (!strcmp(name,"AX")) {
    test=0;
    value=0;
  }
  else if (!strcmp(name,"AY")) {
    test=0;
    value=0;
  }
  else if (!strcmp(name,"AZ")) {
    test=0;
    value=0;
  }
  return test;
}


/*
//--------------------------------------------------------------------------------------------------------------------

void AssembleMatrixResFSI(NonLinearMultiLevelProblem &nl_td_ml_prob2, unsigned level, const unsigned &gridn, const unsigned &ipde, const bool &assemble_matrix) {

  clock_t AssemblyTime=0;
  clock_t start_time, end_time;
  PetscErrorCode ierr;
  
  //Static conversion from the Base class (NonLinearMultiLevelProblem) to the Derived class (NonLinearMultiLevelProblemTimeLoop)
  NonLinearTimeDependentMultiLevelProblem& nl_td_ml_prob = static_cast<NonLinearTimeDependentMultiLevelProblem&>(nl_td_ml_prob2);
 
  const char* pdename= nl_td_ml_prob.GetThisPdeName(ipde);
  
  
  //pointers and references
  Solution*		mysolution	=	nl_td_ml_prob._solution[level];
  LinearSolver*		mylsyspde	=	nl_td_ml_prob._LinSolver[ipde][level];
  mesh*			mymsh		=	nl_td_ml_prob._msh[level];
  elem*			myel		=	mymsh->el;
  SparseMatrix*		myKK		=	mylsyspde->_KK;
  NumericVector*	myRES		=	mylsyspde->_RES;
  vector <int>&		myKKIndex	=	mylsyspde->KKIndex;
  
  
  const unsigned dim = mymsh->GetDimension();
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  
  // local objects
  vector<double> SolVAR(3*dim+1);
  vector<double> SolOldVAR(3*dim+1);
  vector<vector<double> > GradSolVAR(2*dim);
  vector<vector<double> > GradSolOldVAR(2*dim);
  for(int i=0;i<2*dim;i++){
    GradSolVAR[i].resize(dim);
    GradSolOldVAR[i].resize(dim);
  }
  vector<vector<double> > GradSolhatVAR(dim);
  vector<vector<double> > GradSolOldhatVAR(dim);
  for(int i=0;i<dim;i++){
    GradSolhatVAR[i].resize(dim);
    GradSolOldhatVAR[i].resize(dim);
  }
    
  vector <int> metis_node1;
  vector <int> metis_node2;
  vector <bool> solidmark;
  
  vector <double > phi;
  vector <double > phi_hat;
  vector <double > phi_old;
  
  vector <double> gradphi;
  vector <double> gradphi_hat;
  vector <double> gradphi_old;
  
  metis_node1.reserve(max_size);
  metis_node2.reserve(max_size);
  solidmark.reserve(max_size);
  phi.reserve(max_size);
  phi_hat.reserve(max_size);
  phi_old.reserve(max_size);
  gradphi.reserve(max_size*dim);
  gradphi_hat.reserve(max_size*dim);
  gradphi_old.reserve(max_size*dim);
  
  const double * phi1 = NULL;
    
  double Weight=0.;
  double Weight_nojac=0.;
  double Weight_hat=0.;
  double Weight_old=0.;
  
  vector <vector < double> > vx(dim);
  vector <vector < double> > vx_hat(dim);
  vector <vector < double> > vx_old(dim);
  
  for(int i=0;i<dim;i++){
    vx[i].reserve(max_size);
    vx_hat[i].reserve(max_size);
    vx_old[i].reserve(max_size);
  }
   
  vector< vector< double > > Rhs(2*dim+1);
  vector< vector< vector< double > > > B(2*dim+1); 
  for(int i=0;i<2*dim+1;i++){
    B[i].resize(2*dim+1);
  }
  
  vector< vector< int > > dofsVAR(2*dim+1); 
  
  // algorithm parameters
  double eps_pen = 1.e40;
  bool   newton = 0;
  bool   penalty = mylsyspde->GetStabilization();
  
  // ------------------------------------------------------------------------
  // Physical parameters
  double rhof = nl_td_ml_prob._fluid->get_density();
  double rhos = nl_td_ml_prob._solid->get_density();
  double mu_lame = nl_td_ml_prob._solid->get_lame_shear_modulus();
  double lambda_lame = nl_td_ml_prob._solid->get_lame_lambda();
  double mus= mu_lame/rhof;
  double IRe = nl_td_ml_prob._fluid->get_IReynolds_number();
  double lambda= lambda_lame / rhof;
  double betafsi= rhos / rhof;
  double betans=1.;
  int    solid_model= nl_td_ml_prob._solid->get_physical_model();

  //physical quantity
  double Jnp1_hat;
  double Jn_hat;
  double I_bleft;
  double I_e;
  double Cauchy[3][3];
  double Cauchy_old[3][3];
  double b_left[3][3];
  double e[3][3];
  double e_old[3][3];
  double tg_stiff_matrix[3][3];
  
  const double Id2th[3][3]= {
    { 1.,0.,0.},
    { 0.,1.,0.},
    { 0.,0.,1.}
  };
  //initialization C tensor: Saint-Venaint Kirchoff model
  double C_mat[3][3][3][3];
  for (int I=0; I<3; ++I) {
    for (int J=0; J<3; ++J) {
      for (int K=0; K<3; ++K) {
        for (int L=0; L<3; ++L) {
          C_mat[I][J][K][L] = 2.*mus*Id2th[I][K]*Id2th[J][L];
        }
      }
    }
  }

  // ale map
  double _lambda_map=0.;
  double _mu_ale[3] = {1.,1.,1.};
 
  // gravity
  double _gravity[3]={0.,0.,0.};

  // -----------------------------------------------------------------

  // time discretization algorithm paramters
  double dt = nl_td_ml_prob.GetTimeStep();
  const double gamma = 0.5;
  const double gammaratio = (1.-gamma)/gamma;
  double _theta_ns=1.0;

  // space discretization parameters
  unsigned order_ind2 = nl_td_ml_prob.SolType[nl_td_ml_prob.GetIndex("U")];  
  unsigned end_ind2   = mymsh->GetEndIndex(order_ind2);

  unsigned order_ind1 = nl_td_ml_prob.SolType[nl_td_ml_prob.GetIndex("P")];  
  unsigned end_ind1   = mymsh->GetEndIndex(order_ind1);

  // mesh and procs
  unsigned nel    = mymsh->GetElementNumber();
  unsigned igrid  = mymsh->GetGridNumber();
  unsigned iproc  = mymsh->GetProcID();
  unsigned nprocs = mymsh->GetNumProcs();

  //----------------------------------------------------------------------------------
  //variable-name handling
  const char varname[10][3] = {"DX","DY","DZ","U","V","W","P","AX","AY","AZ"};
  const char coordname[3][2] = {"X","Y","Z"};
  vector <unsigned> indexVAR(2*dim+1);
  vector <unsigned> indCOORD(dim);
  vector <unsigned> indVAR(3*dim+1);  
  vector <unsigned> SolType(3*dim+1);  
  
  for(unsigned ivar=0; ivar<dim; ivar++) {
    indCOORD[ivar]=nl_td_ml_prob.GetIndex(&coordname[ivar][0]);
    indVAR[ivar]=nl_td_ml_prob.GetIndex(&varname[ivar][0]);
    indVAR[ivar+dim]=nl_td_ml_prob.GetIndex(&varname[ivar+3][0]);
    indVAR[ivar+2*dim+1]=nl_td_ml_prob.GetIndex(&varname[ivar+7][0]);
    
    SolType[ivar]=nl_td_ml_prob.GetSolType(&varname[ivar][0]);
    SolType[ivar+dim]=nl_td_ml_prob.GetSolType(&varname[ivar+3][0]);
    SolType[ivar+2*dim+1]=nl_td_ml_prob.GetSolType(&varname[ivar+7][0]);
    
    indexVAR[ivar]=nl_td_ml_prob.GetSolPdeIndex("FSI",&varname[ivar][0]);
    indexVAR[ivar+dim]=nl_td_ml_prob.GetSolPdeIndex("FSI",&varname[ivar+3][0]);
  }
  indexVAR[2*dim]=nl_td_ml_prob.GetSolPdeIndex("FSI",&varname[6][0]);
  indVAR[2*dim]=nl_td_ml_prob.GetIndex(&varname[6][0]);
  SolType[2*dim]=nl_td_ml_prob.GetSolType(&varname[6][0]);
  //----------------------------------------------------------------------------------
  
  start_time=clock();
    
  myKK->zero();
  
  /// *** element loop ***
  for(int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

    unsigned kel        = mymsh->IS_Mts2Gmt_elem[iel]; 
    short unsigned kelt = myel->GetElementType(kel);
    unsigned nve        = myel->GetElementDofNumber(kel,end_ind2);
    unsigned nve1       = myel->GetElementDofNumber(kel,end_ind1);
    int flag_mat        = myel->GetElementMaterial(kel);

    //*******************************************************************************************************
    
    //initialization of everything is in common fluid and solid
    
    //Rhs
    for(int i=0; i<2*dim; i++) {
      dofsVAR[i].resize(nve);	
      Rhs[indexVAR[i]].resize(nve);
      memset(&Rhs[indexVAR[i]][0],0,nve*sizeof(double));
    }
    dofsVAR[2*dim].resize(nve1);
    Rhs[indexVAR[2*dim]].resize(nve1);
    memset(&Rhs[indexVAR[2*dim]][0],0,nve1*sizeof(double));
      
    //Kinematic relation (solid) and ALE Map (fluid)
    for(int i=0; i<dim; i++) {
      B[indexVAR[i]][indexVAR[i]].resize(nve*nve);
      memset(&B[indexVAR[i]][indexVAR[i]][0],0,nve*nve*sizeof(double));
    }

    //Stiffness Matrix (solid)
    for(int i=0; i<dim; i++) {
      for(int j=0; j<dim; j++) {
	B[indexVAR[dim+i]][indexVAR[j]].resize(nve*nve);
	memset(&B[indexVAR[dim+i]][indexVAR[j]][0],0,nve*nve*sizeof(double));
      }
    }
      
    //Mass Matrix (solid and fluid) and Diffusion Matrix (fluid only) 
    for(int i=0; i<dim; i++) {
      B[indexVAR[dim+i]][indexVAR[dim+i]].resize(nve*nve);
      memset(&B[indexVAR[dim+i]][indexVAR[dim+i]][0],0,nve*nve*sizeof(double));
    }

    //Pressure gradient (fluid and solid)
    for(int i=0; i<dim; i++) {
      B[indexVAR[dim+i]][indexVAR[2*dim]].resize(nve*nve1);
      memset(&B[indexVAR[dim+i]][indexVAR[2*dim]][0],0,nve*nve1*sizeof(double));
    }

    // Pressure Mass Matrix
    B[indexVAR[2*dim]][indexVAR[2*dim]].resize(nve1*nve1);
    memset(&B[indexVAR[2*dim]][indexVAR[2*dim]][0],0,nve1*nve1*sizeof(double));
    
    if (flag_mat==2) { //initialization for fluid only
      // Fluid Continuity Matrix: divergence of the velocity
      for(int i=0; i<dim; i++) {
	B[indexVAR[2*dim]][indexVAR[dim+i]].resize(nve1*nve);
	memset(&B[indexVAR[2*dim]][indexVAR[dim+i]][0],0,nve1*nve*sizeof(double));
      }
    }
    else{ // initialization for solid only
      // Kinematic relation
      for(int i=0; i<dim; i++) {
	B[indexVAR[i]][indexVAR[dim+i]].resize(nve*nve);
	memset(&B[indexVAR[i]][indexVAR[dim+i]][0],0,nve*nve*sizeof(double));
      }
      // Solid Continuity Matrix: divergence of the displacemnet
      for(int i=0; i<dim; i++) {
	B[indexVAR[2*dim]][indexVAR[i]].resize(nve1*nve);
	memset(&B[indexVAR[2*dim]][indexVAR[i]][0],0,nve1*nve*sizeof(double));
      }
    }
    
    // ----------------------------------------------------------------------------------------
    // coordinates, displacement, velocity dofs
    
    metis_node2.resize(nve);
    metis_node1.resize(nve1);
    solidmark.resize(nve);
    phi.resize(nve);
    phi_hat.resize(nve);
    phi_old.resize(nve);
    gradphi.resize(nve*dim);
    gradphi_hat.resize(nve*dim);
    gradphi_old.resize(nve*dim);
        
    for(int i=0;i<dim;i++){
      vx[i].resize(nve);
      vx_old[i].resize(nve);
      vx_hat[i].resize(nve);
    }
    
    for (unsigned i=0;i<nve;i++) {
      // gambit nodes
      unsigned inode=myel->GetElementVertexIndex(kel,i)-1u;
      // dof metis
      unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
      metis_node2[i]=inode_Metis;
      
      //unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
      // flag to know if the node "inode" lays on the fluid-solid interface
      solidmark[i]=myel->GetNodeRegion(inode); // to check
      for(int j=0; j<dim; j++) {
	//Updated coordinates (Moving frame)
        vx[j][i]= (*mysolution->_Sol[indCOORD[j]])(inode_Metis) + (*mysolution->_Sol[indVAR[j]])(inode_Metis);
	//Old coordinates (Moving frame)
        vx_old[j][i]= (*mysolution->_Sol[indCOORD[j]])(inode_Metis) + (*mysolution->_SolOld[indVAR[j]])(inode_Metis);
	//Fixed coordinates (Reference frame)
	vx_hat[j][i]= (*mysolution->_Sol[indCOORD[j]])(inode_Metis);  
	// displacement dofs
	dofsVAR[j][i]= mylsyspde->GetKKDof(indVAR[j],indexVAR[j],inode); 
	// velocity dofs
	dofsVAR[j+dim][i]= mylsyspde->GetKKDof(indVAR[j+dim],indexVAR[j+dim],inode);   
      }
    }

    // pressure dofs
    for (unsigned i=0;i<nve1;i++) {
      unsigned inode=(order_ind1<3)?(myel->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
      metis_node1[i]=mymsh->GetMetisDof(inode,SolType[2*dim]);
      dofsVAR[2*dim][i]=mylsyspde->GetKKDof(indVAR[2*dim],indexVAR[2*dim],inode);
    }
    // ----------------------------------------------------------------------------------------
       
    if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {
      /// *** Gauss point loop ***
      for (unsigned ig=0;ig < nl_td_ml_prob.type_elem[kelt][order_ind2]->GetGaussPointNumber(); ig++) {

	// *** get Jacobian and test function and test function derivatives in the moving frame***
	(nl_td_ml_prob.type_elem[kelt][order_ind2]->*(nl_td_ml_prob.type_elem[kelt][order_ind2])->Jacobian_ptr)(vx,ig,Weight,phi,gradphi);
	(nl_td_ml_prob.type_elem[kelt][order_ind2]->*(nl_td_ml_prob.type_elem[kelt][order_ind2])->Jacobian_ptr)(vx_old,ig,Weight_old,phi_old,gradphi_old);
	(nl_td_ml_prob.type_elem[kelt][order_ind2]->*(nl_td_ml_prob.type_elem[kelt][order_ind2])->Jacobian_ptr)(vx_hat,ig,Weight_hat,phi_hat,gradphi_hat);
	phi1=nl_td_ml_prob.type_elem[kelt][order_ind1]->GetPhi(ig);
	if (flag_mat==2) Weight_nojac = nl_td_ml_prob.type_elem[kelt][order_ind2]->GetGaussWeight(ig);

	// ---------------------------------------------------------------------------
	// displacement and velocity
	for(int i=0; i<2*dim; i++){
	  SolVAR[i]=0.;
	  SolOldVAR[i]=0.;
	  for(int j=0; j<dim; j++) {
	    GradSolVAR[i][j]=0.;
	    GradSolOldVAR[i][j]=0.;
	    if(i<dim){
	      GradSolhatVAR[i][j]=0.;
	      GradSolOldhatVAR[i][j]=0.;
	    }
	  }
	    
	  for (unsigned inode=0; inode<nve; inode++) {
	    unsigned sol_dof = metis_node2[inode];
	      
	    double soli = (*mysolution->_Sol[indVAR[i]])(sol_dof);
	    SolVAR[i]+=phi[inode]*soli;
	      
	    double soli_old = (*mysolution->_SolOld[indVAR[i]])(sol_dof);
	    SolOldVAR[i]+=phi[inode]*soli_old;
	      
	    for(int j=0; j<dim; j++) {
	      GradSolVAR[i][j]+=gradphi[inode*dim+j]*soli;
	      GradSolOldVAR[i][j]+=gradphi[inode*dim+j]*soli_old;
	      if(i<dim){ 
		GradSolhatVAR[i][j]   +=gradphi_hat[inode*dim+j]*soli;
		GradSolOldhatVAR[i][j]+=gradphi_hat[inode*dim+j]*soli_old;
	      }
	    }	      
	  }
	}  
		  
	// pressure
	SolVAR[2*dim]=0.;
	for (unsigned inode=0; inode<nve1; inode++) {
	  double soli = (*mysolution->_Sol[indVAR[2*dim]])(metis_node1[inode]);
	  SolVAR[2*dim]+=phi1[inode]*soli;
	}
  
	// acceleration (solid only)  
	if(flag_mat!=2){
	  for(int i=2*dim+1; i<3*dim+1; i++){
	    SolVAR[i]=0.;
	    for (unsigned inode=0; inode<nve; inode++) {
	      double soli = (*mysolution->_Sol[indVAR[i]])(metis_node2[inode]);
	      SolVAR[i]+=phi[i]*soli;  
	    }
	  }
	}
  
	// ---------------------------------------------------------------------------
	 	  
	//BEGIN FLUID ASSEMBLY ============
	  
	if(flag_mat==2){
	  
          //divergence of the velocity
	  double div_vel=0.;
	  double div_w=0.;
	  for(int i=0; i<dim; i++) {
	    div_vel+=GradSolVAR[dim+i][i];
	    div_w+=(GradSolVAR[i][i]-GradSolOldVAR[i][i])*(1./dt);
	  }
	  
          { // Laplace operator + adection operator + Mass operator

            const double *gradfi=&gradphi[0];
            const double *fi=&phi[0];

            // *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim,fi++) {

              //BEGIN RESIDUALS A block ===========================
 	      double LapmapVAR[3] = {0., 0., 0.};
 	      for(int idim=0; idim<dim; idim++) {
 		for(int idim2=0; idim2<dim; idim2++) {
 	          LapmapVAR[idim] += dt*( _mu_ale[idim2]*gradphi_hat[i*dim+idim2]*GradSolhatVAR[idim][idim2] );
 		}
 	      }
      
	      // Residual ALE equations
 	      for(int idim=0; idim<dim; idim++) {
 	        Rhs[indexVAR[idim]][i]+=(!solidmark[i])*(-LapmapVAR[idim]*Weight_nojac);
 	      }
	      
	      //-------------------------------------------------------------------------------
	      
	      double LapvelVAR[3]={0.,0.,0.};
	      double AdvaleVAR[3]={0.,0.,0.};
	      for(int idim=0.; idim<dim; idim++) {
		for(int idim2=0.; idim2<dim; idim2++) {
		  LapvelVAR[idim]+=gradphi[i*dim+idim2]*GradSolVAR[dim+idim][idim2];
		  AdvaleVAR[idim]+=((SolVAR[dim+idim2]*dt - (SolVAR[idim2]-SolOldVAR[idim2]))*GradSolVAR[dim+idim][idim2])*phi[i];;
		}
	      }
	      // Residual Momentum equations 
	      for(int idim=0; idim<dim; idim++) {
		Rhs[indexVAR[dim+idim]][i]+= (
					      -phi[i]*betans*SolVAR[dim+idim]*Weight
					      +phi[i]*betans*SolOldVAR[dim+idim]*Weight_old
					      -AdvaleVAR[idim]*0.5*(Weight+Weight_old)
					      -0.5*dt*div_vel*SolVAR[dim+idim]*phi[i]*0.5*(Weight+Weight_old)
					      +dt*div_w*SolVAR[dim+idim]*phi[i]*0.5*(Weight+Weight_old)
					      -dt*IRe*LapvelVAR[idim]*Weight
					      +dt*SolVAR[2*dim]*gradphi[i*dim+idim]*Weight
					      );
	      }
	      //END RESIDUALS A block ===========================
	      
              const double *gradfj=&gradphi[0];
              const double *fj=&phi[0];
              //  *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim,fj++) {

                //Laplacian
                double Lap=0.;
		for(int idim=0; idim<dim; idim++) {
		  Lap+=(*(gradfi+idim))*(*(gradfj+idim));
		}
                double LapXweight = Lap*Weight;

                //advection term I
		double Adv1=0.;
		for(int idim=0; idim<dim; idim++) {
		  Adv1+= ((SolVAR[dim+idim]*dt - (SolVAR[idim]-SolOldVAR[idim]))*(*(gradfj+idim))*(*(fi)))*0.5*(Weight+Weight_old);
		}
				
		double div_stab = 0.5*div_vel*((*fi))*((*fj))*0.5*(Weight+Weight_old);
		
		double div_ale = div_w*((*fi))*((*fj))*0.5*(Weight+Weight_old);

                double Mass = ((*fi))*((*fj))*Weight;

		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[dim+idim]][i*nve+j] += dt*IRe*LapXweight + Adv1 + dt*div_stab - dt*div_ale + betans*Mass;
		}
		
                for(int idim=0; idim<dim; idim++) {
		  for(int idim2=0; idim2<dim; idim2++) {
		    B[indexVAR[dim+idim]][indexVAR[idim2]][i*nve+j] += betans*SolVAR[dim+idim]*((*fi))*(*(gradfi+idim2))*Weight;
		  }
		}

		double Lap_ale=0.;
		for(int idim=0; idim<dim; idim++) {
		  Lap_ale+=_mu_ale[idim]*(*(gradfi+idim))*(*(gradfj+idim));
		}
		  
		// Laplacian ALE map
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[0+idim]][indexVAR[idim]][i*nve+j] += (!solidmark[i])*dt*Lap_ale*Weight_nojac;  
		}

              } // end phi_j loop
            } // end phi loop
          } // end A + Bt

          { //Gradient of Pressure operator

            const double *gradfi=&gradphi[0];
            const double *fi=phi1;

            /// *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim,fi++) {

              const double *fj=phi1;
              /// *** phi_j loop ***
              for (unsigned j=0; j<nve1; j++,fj++) {
                for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[2*dim]][i*nve1+j] -= dt*((*(gradfi+idim))*(*fj))*Weight;
		}
              } // end phi_j loop
            } // end phi_i loop
          } // End Bt
        
          { // Divergence of the Velocity operator
            const double *fi=phi1;
            // *** phi_i loop ***
            for (unsigned i=0; i<nve1; i++,fi++) {

              //BEGIN RESIDUALS B block ===========================
              Rhs[indexVAR[2*dim]][i] += -(-((*fi))*div_vel)*Weight;
              //END RESIDUALS  B block ===========================

              const double *gradfj=&gradphi[0];
              // *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim) {
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[2*dim]][indexVAR[idim+dim]][i*nve+j] -= ((*fi)*(*(gradfj+idim)))*Weight;
		}
              }
            }
          }
	}   
	//END FLUID ASSEMBLY ============
	//*******************************************************************************************************
	//BEGIN SOLID ASSEMBLY ============
	  
	else{
	  
	  //------------------------------------------------------------------------------------------------------------
          if (solid_model==0) {
	    
	    //computation of the stress tensor
	    for(int i=0;i<dim;i++){
	      for(int j=0;j<dim;j++){
		e[i][j]=0.5*(GradSolhatVAR[i][j]+GradSolhatVAR[j][i]);
		e_old[i][j]=0.5*(GradSolOldhatVAR[i][j]+GradSolOldhatVAR[j][i]);
	      }
	    }
	    I_e=0;
	    for(int i=0;i<dim;i++){
	      I_e += e[i][i];
	    }
	    
	    for (int i=0; i<dim; i++) {
              for (int j=0; j<dim; j++) {
                //incompressible
                Cauchy[i][j] = 2*mus*e[i][j];
		Cauchy_old[i][j]=2*mus*e_old[i][j];
              }
            }
	    
	    
//             //computation of the stress tensor
//             e[0][0] = GradSolhatVAR[0][0];
//             e[0][1] = 0.5*(GradSolhatVAR[0][1] + GradSolhatVAR[1][0]);
//             e[0][2] = 0.;
//             e[0][2] = 0.*0.5*(GradSolhatVAR[0][2] + GradSolhatVAR[2][0]);
// 
//             e[1][0] = 0.5*(GradSolhatVAR[1][0] + GradSolhatVAR[0][1]);
//             e[1][1] = GradSolhatVAR[1][1];
//             e[1][2] = 0.;
//             e[1][2] = 0.*0.5*(GradSolhatVAR[1][2] + GradSolhatVAR[2][1]);
// 
//             e[2][0] = 0.*0.5*(GradSolhatVAR[0][2] + GradSolhatVAR[2][0]);
//             e[2][1] = 0.*0.5*(GradSolhatVAR[1][2] + GradSolhatVAR[2][1]);
//             e[2][2] = 0.*GradSolhatVAR[2][2];
// 
//             I_e = e[0][0] + e[1][1] + e[2][2];
// 
// 
//             Cauchy stress tensor
//             for (int irow=0; irow<3; ++irow) {
//               for (int jcol=0; jcol<3; ++jcol) {
//                 //compressible
// 		// 	   Cauchy[irow][jcol] = _lambda*I_e*Id2th[irow][jcol] + 2*_mus*e[irow][jcol];
//                 //incompressible
//                 Cauchy[irow][jcol] = 2*mus*e[irow][jcol];
//               }
//             }
// 
//             //computation of the older stress tensor
//             e_old[0][0] = GradSolOldhatVAR[0][0];
//             e_old[0][1] = 0.5*(GradSolOldhatVAR[0][1] + GradSolOldhatVAR[1][0]);
//             e_old[0][2] = 0.;
//             e_old[0][2] = 0.*0.5*(GradSolOldhatVAR[0][2] + GradSolOldhatVAR[2][0]);
// 
//             e_old[1][0] = 0.5*(GradSolOldhatVAR[1][0] + GradSolOldhatVAR[0][1]);
//             e_old[1][1] = GradSolOldhatVAR[1][1];
//             e_old[1][2] = 0.;
//             e_old[1][2] = 0.*0.5*(GradSolOldhatVAR[1][2] + GradSolOldhatVAR[2][1]);
// 
//             e_old[2][0] = 0.*0.5*(GradSolOldhatVAR[0][2] + GradSolOldhatVAR[2][0]);
//             e_old[2][1] = 0.*0.5*(GradSolOldhatVAR[1][2] + GradSolOldhatVAR[2][1]);
//             e_old[2][2] = 0.*GradSolOldhatVAR[2][2];
// 
//             //Iold_e = e_old[0][0] + e_old[1][1] + e_old[2][2];
// 
//             // Cauchy stress tensor
//             for (int irow=0; irow<3; ++irow) {
//               for (int jcol=0; jcol<3; ++jcol) {
//                 //compressible
// 		// 	   Cauchy[irow][jcol] = _lambda*I_e*Id2th[irow][jcol] + 2*_mus*e[irow][jcol];
//                 //incompressible
//                 Cauchy_old[irow][jcol] = 2*mus*e_old[irow][jcol];
//               }
//		}
          }

          else if (solid_model==1) {
	    double F[3][3]={{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
	    for(int i=0;i<dim;i++){
	      for(int j=0;j<dim;j++){
		F[i][j]+=GradSolhatVAR[i][j];
	      }
	    }
	    
	    Jnp1_hat =  F[0][0]*F[1][1]*F[2][2] + F[0][1]*F[1][2]*F[2][0] + F[0][2]*F[1][0]*F[2][1]
		      - F[2][0]*F[1][1]*F[0][2] - F[2][1]*F[1][2]*F[0][0] - F[2][2]*F[1][0]*F[0][1];

	    // computation of the the three deformation tensor b
	    for (int I=0; I<3; ++I) {
	      for (int J=0; J<3; ++J) {
		b_left[I][J]=0.;
		for (int K=0; K<3; ++K) {
		  //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
		  b_left[I][J] += F[I][K]*F[J][K];
		}
		Cauchy[I][J] = (mus/Jnp1_hat)*(b_left[I][J] - Id2th[I][J]);
	      }
	    }

	    I_bleft = b_left[0][0] + b_left[1][1] + b_left[2][2];
	    
	    //compressible case
	    //             for (int ii=0; ii<3; ++ii) {
	    //               for (int jj=0; jj<3; ++jj) {
	    //                 //Cauchy stress tensor
	    //                 Cauchy[ii][jj] = (_mus/Jnp1_hat)*(b_left[ii][jj] - Id2th[ii][jj]) + (_lambda/Jnp1_hat)*log(Jnp1_hat)*Id2th[ii][jj];
	    //                 for (int k=0; k<3; ++k) {
	    //                   for (int l=0; l<3; ++l) {
	    //                   }
	    //                 }
	    //               }
	    //             }
	    
	    //for the incompressible(nearly incompressible) case
	    for (int ii=0; ii<3; ++ii) {
	      for (int jj=0; jj<3; ++jj) {
		for (int kk=0; kk<3; ++kk) {
		  for (int ll=0; ll<3; ++ll) {
		    C_mat[ii][jj][kk][ll] = 2.*mus*pow(Jnp1_hat,-1.6666666666666)*(
										  0.333333333333*I_bleft*Id2th[ii][kk]*Id2th[jj][ll]              //1/3*I_c*i
										  // 	                        +0.111111111111*I_C*Id2th[i][j]*Id2th[k][l]             //1/9*I_b*IxI
										  // 				-0.333333333333*b_left[i][j]*Id2th[k][l]                //-1/3*b*I
										  // 				-0.333333333333*Id2th[i][j]*b_left[k][l]                //-1/3*b*I
										  )
		      -SolVAR[2*dim]*(Id2th[ii][jj]*Id2th[kk][ll]-2.*Id2th[ii][kk]*Id2th[jj][ll] );  // -p(IxI-2i)
		  }
		}
	      }
	    }
	    
	    //Old deformation gradient
	    double F_old[3][3]={{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
	    for(int i=0;i<dim;i++){
	      for(int j=0;j<dim;j++){
		F_old[i][j]+=GradSolhatVAR[i][j];
	      }
	    }
            
            Jn_hat =  F_old[0][0]*F_old[1][1]*F_old[2][2] + F_old[0][1]*F_old[1][2]*F_old[2][0] + F_old[0][2]*F_old[1][0]*F_old[2][1]
		    - F_old[2][0]*F_old[1][1]*F_old[0][2] - F_old[2][1]*F_old[1][2]*F_old[0][0] - F_old[2][2]*F_old[1][0]*F_old[0][1] ;

            // computation of the the three deformation tensor b
            for (int I=0; I<3; ++I) {
              for (int J=0; J<3; ++J) {
                b_left[I][J]=0.;
                for (int k=0; k<3; ++k) {
                  //left Cauchy-Green deformation tensor or F_oldinger tensor (b = F_old*F_old^T)
                  b_left[I][J] += F_old[I][k]*F_old[J][k];
                }
                Cauchy_old[I][J] = (mus/Jn_hat)*(b_left[I][J] - Id2th[I][J]);
              }
            }
          }
          //----------------------------------------------------------------------------------------------------------------------------

          /////////////

          ///Mass + Stiffness operator
          {
            const double *gradfi=&gradphi[0];
            const double *fi=&phi[0];

            /// *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim,fi++) {


              //BEGIN RESIDUALS A + Bt block ===========================
	      
	      // Residual ALE equations
 	      for(int idim=0; idim<dim; idim++) {
	        Rhs[indexVAR[idim]][i] += (
					   -phi[i]*eps_pen*(
							    +SolVAR[idim] - SolOldVAR[idim]
							    -dt*gamma*SolVAR[dim+idim]
							    -dt*(1.-gamma)*SolOldVAR[dim+idim]
							    )
					   )*Weight_hat;
              }
              
              double CauchyDIR[3]={0.,0.,0.};
	      for(int idim=0.; idim<dim; idim++) {
		for(int idim2=0.; idim2<dim; idim2++) {
		  CauchyDIR[idim]+= gradphi[i*dim+idim2]*Cauchy[idim][idim2];
		}
	      }

              // Residual Momentum equations
              for(int idim=0; idim<dim; idim++) {
	        Rhs[indexVAR[dim+idim]][i] += (
					       phi[i]*dt*_gravity[idim]*Weight_hat
					       -phi[i]*betafsi*(SolVAR[dim+idim] - SolOldVAR[dim+idim])*Weight_hat*(1./gamma)
					       +phi[i]*betafsi*SolVAR[2*dim+1+idim]*Weight_hat*gammaratio*dt
					       -dt*CauchyDIR[idim]*Weight
					       +dt*SolVAR[2*dim]*gradphi[i*dim+idim]*Weight
					       );

              }
              
              //---------------------------------------------------------------------------------------------------------------------------------

              //END RESIDUALS A + Bt block ===========================

              const double *gradfj=&gradphi[0];
              const double *fj=&phi[0];
              // *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim,fj++) {

                /// Mass term
                // (v_n+1,psi)
                //gamma = 1 Backward Euler
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[dim+idim]][i*nve+j] += betafsi*(*(fi))*(*(fj))*Weight_hat*(1./gamma);
		}

                //Da collaudare il 3D
                for (int icount=0; icount<dim; ++icount) {
                  for (int jcount=0; jcount<dim; ++jcount) {
                    tg_stiff_matrix[icount][jcount] = 0.;
                    for (int kcount=0; kcount<dim; ++kcount) {
                      for (int lcount=0; lcount<dim; ++lcount) {
                        tg_stiff_matrix[icount][jcount] += (*(gradfi+kcount))*0.25*(
										    C_mat[icount][kcount][jcount][lcount]+C_mat[icount][kcount][lcount][jcount]
										    +C_mat[kcount][icount][jcount][lcount]+C_mat[kcount][icount][lcount][jcount]
										    )*(*(gradfj+lcount));
                      }
                    }
                  }
                }
                
                //geometric tangent stiffness matrix
                double geom_tg_stiff_matrx = 0.;
                for(int kcount=0; kcount<dim; ++kcount) {
                  for(int lcount=0; lcount<dim; ++lcount) {
                    geom_tg_stiff_matrx += (*(gradfi+kcount))*Cauchy[kcount][lcount]*(*(gradfj+lcount));
                  }
                }

                /// Stiffness operator -- Elasticity equation (Linear or not)
                for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[0+idim]][i*nve+j] += dt*geom_tg_stiff_matrx*Weight;
 		  for(int idim2=0; idim2<dim; idim2++) {
 		    B[indexVAR[dim+idim]][indexVAR[0+idim2]][i*nve+j] += dt*tg_stiff_matrix[0+idim][0+idim2]*Weight;
 		  }
 		}

                /// Kinematic equation v = du/dt
                //   -_theta*dt*(v_n+1,eta)
		for(int idim=0; idim<dim; idim++) {
		  //   -_theta*dt*(v_n+1,eta)
		  B[indexVAR[0+idim]][indexVAR[dim+idim]][i*nve+j] -= eps_pen*gamma*dt*(*(fi))*(*(fj))*Weight_hat;
		  // (u_n+1,eta)
		  B[indexVAR[0+idim]][indexVAR[0+idim]][i*nve+j] += eps_pen*(*(fi))*(*(fj))*Weight_hat;
		}
	      }
            }
          }
          ////////////
          { ///Gradient of Pressure
            const double *gradfi=&gradphi[0];
            // *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim) {
              const double *fj=phi1;
              // *** phi_j loop ***
              for (unsigned j=0; j<nve1; j++,fj++) {
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[2*dim]][i*nve1+j] -= dt*((*(gradfi+idim))*(*fj))*Weight;
		}
              }
            }
          }
          ////////////
          { ///Divergence of the Displacement
            const double *fi=phi1;
            // *** phi_i loop ***
            for (unsigned i=0; i<nve1; i++,fi++) {

              //BEGIN RESIDUALS B block ===========================

              if (solid_model==0) {
                Rhs[indexVAR[2*dim]][i] += -(-((*fi))*(I_e + (1./lambda)*SolVAR[2*dim] ) )*Weight_hat;
              }
              else if (solid_model==1) {
                Rhs[indexVAR[2*dim]][i] += -(-((*fi))*( log(Jnp1_hat)/Jnp1_hat + (1./lambda)*SolVAR[2*dim] ) )*Weight_hat;
              }

              //END RESIDUALS B block ===========================

              const double *gradfj=&gradphi[0];
              // *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim) {
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[2*dim]][indexVAR[idim]][i*nve+j] -= ((*fi)*(*(gradfj+idim)))*Weight;
		}
              }
            }
          }
          //  /////////////
          {  ///Pressure Mass term
            const double *fi=phi1;
            // *** phi_i loop ***
            for (unsigned i=0; i<nve1; i++,fi++) {
              const double *fj=phi1;
              // *** phi_j loop ***
              for (unsigned j=0; j<nve1; j++,fj++) {
                B[indexVAR[2*dim]][indexVAR[2*dim]][i*nve1+j] -= (1./lambda)*((*fi)*(*fj))*Weight_hat;
              }
            }
          }  //end pressure mass term
	  //---------------------------------------------------------------------------------------------------------------------------------
	}
	  
	//END SOLID ASSEMBLY ============
      }
    }

    //BEGIN local to global assembly 
    // ALE mapping
    for(int i=0; i<dim; i++) {
      myRES->add_vector_blocked(Rhs[indexVAR[i]],dofsVAR[i]);
      myKK ->add_matrix_blocked(B[indexVAR[i]][indexVAR[i]],dofsVAR[i],dofsVAR[i]);  
      if(flag_mat!=2){ //Solid only
	myKK->add_matrix_blocked(B[indexVAR[i]][indexVAR[i+dim]],dofsVAR[i],dofsVAR[i+dim]);
      }
    }
    
    // Momentum equation
    for(int i=0; i<dim; i++) {
      myRES->add_vector_blocked(Rhs[indexVAR[dim+i]],dofsVAR[dim+i]);
      myKK->add_matrix_blocked(B[indexVAR[dim+i]][indexVAR[dim+i]],dofsVAR[dim+i],dofsVAR[dim+i]);
      myKK->add_matrix_blocked(B[indexVAR[dim+i]][indexVAR[2*dim]],dofsVAR[dim+i],dofsVAR[2*dim]);
      for(int j=0; j<dim; j++) {
	myKK->add_matrix_blocked(B[indexVAR[dim+i]][indexVAR[j]],dofsVAR[dim+i],dofsVAR[j]);
      }
    }
   
    //P-continuity equation
    myRES->add_vector_blocked(Rhs[indexVAR[2*dim]],dofsVAR[2*dim]);
    for(int i=0; i<dim; i++) {
      if(flag_mat==2){ //Fluid only
	myKK->add_matrix_blocked(B[indexVAR[2*dim]][indexVAR[dim+i]],dofsVAR[2*dim],dofsVAR[dim+i]);
      }
      else{ //Solid only
	myKK->add_matrix_blocked(B[indexVAR[2*dim]][indexVAR[i]],dofsVAR[2*dim],dofsVAR[i]);  
      } 
    }
    myKK->add_matrix_blocked(B[indexVAR[2*dim]][indexVAR[2*dim]],dofsVAR[2*dim],dofsVAR[2*dim]);
    //END local to global assembly
   
  } //end list of elements loop


  //----------------------------------------------------------------------------------------------------------
  myKK->close();
  
  myRES->close();
  
  // *************************************
  end_time=clock();
  AssemblyTime+=(end_time-start_time);
  // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

  //----------------------------------------------------------------------------------------------

}*/



//===============================================================================================================
//====================================== TIME DEPENDENT, OTHER VERSION ==========================================
//===============================================================================================================





/*
void AssembleMatrixResFSI(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assemble_matrix) {
    
  clock_t AssemblyTime=0;
  clock_t start_time, end_time;
  
  //pointers and references
  MultiLevelSolution*	 ml_sol	                      = ml_prob._ml_sol;
  Solution*	 mysolution  	                      = ml_sol->GetSolutionLevel(level);
  MonolithicFSINonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<MonolithicFSINonLinearImplicitSystem>("Fluid-Structure-Interaction");
  LinearEquationSolver*  myLinEqSolver	              = my_nnlin_impl_sys._LinSolver[level];   
  mesh		*mymsh		=  ml_prob._ml_msh->GetLevel(level);
  elem		*myel		=  mymsh->el;
  SparseMatrix	*myKK		=  myLinEqSolver->_KK;
  NumericVector *myRES		=  myLinEqSolver->_RES;
  vector <int>	&myKKIndex	=  myLinEqSolver->KKIndex;
    
  const unsigned dim = mymsh->GetDimension();
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  
  // local objects
  vector<double> SolVAR(2*dim+1);
  vector<vector<double> > GradSolVAR(2*dim);
  for(int i=0;i<2*dim;i++){
    GradSolVAR[i].resize(dim);
  }
  vector<vector<double> > GradSolhatVAR(dim);
  for(int i=0;i<dim;i++){
    GradSolhatVAR[i].resize(dim);
  }
    
  vector <int> metis_node1;
  vector <int> metis_node2;
  vector <bool> solidmark;
  
  vector <double > phi;
  vector <double > phi_hat;
  
  vector <double> gradphi;
  vector <double> gradphi_hat;
  
  metis_node1.reserve(max_size);
  metis_node2.reserve(max_size);
  solidmark.reserve(max_size);
  phi.reserve(max_size);
  phi_hat.reserve(max_size);
  gradphi.reserve(max_size*dim);
  gradphi_hat.reserve(max_size*dim);
  
  const double *phi1;
    
  double Weight=0.;
  double Weight_nojac=0.;
  double Weight_hat=0.;
  
  vector <vector < double> > vx(dim);
  vector <vector < double> > vx_hat(dim);
  
  for(int i=0;i<dim;i++){
    vx[i].reserve(max_size);
    vx_hat[i].reserve(max_size);
  }
   
  vector< vector< double > > Rhs(2*dim+1);
  vector< vector< vector< double > > > B(2*dim+1); 
  for(int i=0;i<2*dim+1;i++){
    B[i].resize(2*dim+1);
  }
  vector< vector< int > > dofsVAR(2*dim+1); 
  
  // ------------------------------------------------------------------------
  // Physical parameters
  double rhof	 	= ml_prob.parameters.get<Fluid>("Fluid").get_density();            
  double rhos 		= ml_prob.parameters.get<Solid>("Solid").get_density();            
  double mu_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_shear_modulus(); 
  double lambda_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_lambda();        
  double mus		= mu_lame/rhof;
  double IRe 		= ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();   
  double lambda		= lambda_lame / rhof;
  double betafsi	= rhos / rhof;
  double betans		= 1.;
  int    solid_model	= ml_prob.parameters.get<Solid>("Solid").get_physical_model();     

  //physical quantity
  double Jnp1_hat;
  double Jn_hat;
  double I_bleft;
  double I_e;
  double Cauchy[3][3];
  double Cauchy_old[3][3]; 
  double tg_stiff_matrix[3][3];
  //initialization C tensor: Saint-Venaint Kirchoff model : solid_model==1;
  const double Id2th[3][3]= {{ 1.,0.,0.}, { 0.,1.,0.}, { 0.,0.,1.}};
  double C_mat[3][3][3][3];
  for (int I=0; I<3; ++I) {
    for (int J=0; J<3; ++J) {
      for (int K=0; K<3; ++K) {
        for (int L=0; L<3; ++L) {
          C_mat[I][J][K][L] = 2.*mus*Id2th[I][K]*Id2th[J][L];
        }
      }
    }
  }

  // ale map
  double _lambda_map=0.;
  double _mu_ale[3] = {1.,1.,1.};
 
  // gravity
  double _gravity[3]={0.,-1.,0.};
  
  // newton algorithm
  bool nwtn_alg = true;

  // -----------------------------------------------------------------
  // space discretization parameters
  unsigned order_ind2 = ml_sol->GetSolutionType(ml_sol->GetIndex("U"));  
  unsigned end_ind2   = mymsh->GetEndIndex(order_ind2);

  unsigned order_ind1 = ml_sol->GetSolutionType(ml_sol->GetIndex("P"));  
  unsigned end_ind1   = mymsh->GetEndIndex(order_ind1);

  // mesh and procs
  unsigned nel    = mymsh->GetElementNumber();
  unsigned igrid  = mymsh->GetGridNumber();
  unsigned iproc  = mymsh->processor_id();

  //----------------------------------------------------------------------------------
  //variable-name handling
  const char varname[7][3] = {"DX","DY","DZ","U","V","W","P"};
  vector <unsigned> indexVAR(2*dim+1);
  vector <unsigned> indVAR(2*dim+1);  
  vector <unsigned> SolType(2*dim+1);  
  
  for(unsigned ivar=0; ivar<dim; ivar++) {
    indVAR[ivar]=ml_sol->GetIndex(&varname[ivar][0]);
    indVAR[ivar+dim]=ml_sol->GetIndex(&varname[ivar+3][0]);
    SolType[ivar]=ml_sol->GetSolutionType(&varname[ivar][0]);
    SolType[ivar+dim]=ml_sol->GetSolutionType(&varname[ivar+3][0]);
    indexVAR[ivar]=my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);
    indexVAR[ivar+dim]=my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar+3][0]);
  }
  indexVAR[2*dim]=my_nnlin_impl_sys.GetSolPdeIndex(&varname[6][0]);
  indVAR[2*dim]=ml_sol->GetIndex(&varname[6][0]);
  SolType[2*dim]=ml_sol->GetSolutionType(&varname[6][0]);
  //----------------------------------------------------------------------------------
  
  start_time=clock();
    
  myKK->zero();
  
  /// *** element loop ***
  for(int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

    unsigned kel        = mymsh->IS_Mts2Gmt_elem[iel]; 
    short unsigned kelt = myel->GetElementType(kel);
    unsigned nve        = myel->GetElementDofNumber(kel,end_ind2);
    unsigned nve1       = myel->GetElementDofNumber(kel,end_ind1);
    int flag_mat        = myel->GetElementMaterial(kel);

    //*******************************************************************************************************
    
    //initialization of everything is in common fluid and solid
    
    //Rhs
    for(int i=0; i<2*dim; i++) {
      dofsVAR[i].resize(nve);	
      Rhs[indexVAR[i]].resize(nve);
      memset(&Rhs[indexVAR[i]][0],0,nve*sizeof(double));
    }
    dofsVAR[2*dim].resize(nve1);
    Rhs[indexVAR[2*dim]].resize(nve1);
    memset(&Rhs[indexVAR[2*dim]][0],0,nve1*sizeof(double));
      
    //Kinematic relation (solid) and ALE Map (fluid)
    for(int i=0; i<dim; i++) {
      B[indexVAR[i]][indexVAR[i]].resize(nve*nve);
      memset(&B[indexVAR[i]][indexVAR[i]][0],0,nve*nve*sizeof(double));
    }

    //Stiffness Matrix (solid)
    for(int i=0; i<dim; i++) {
      for(int j=0; j<dim; j++) {
	B[indexVAR[dim+i]][indexVAR[j]].resize(nve*nve);
	memset(&B[indexVAR[dim+i]][indexVAR[j]][0],0,nve*nve*sizeof(double));
	
// 	B[indexVAR[dim+i]][indexVAR[dim+j]].resize(nve*nve);
// 	memset(&B[indexVAR[dim+i]][indexVAR[dim+j]][0],0,nve*nve*sizeof(double));
      }
    }
      
    //Mass Matrix (solid and fluid) and Diffusion Matrix (fluid only) 
    for(int i=0; i<dim; i++) {
      B[indexVAR[dim+i]][indexVAR[dim+i]].resize(nve*nve);
      memset(&B[indexVAR[dim+i]][indexVAR[dim+i]][0],0,nve*nve*sizeof(double));
      if(nwtn_alg== true) {
	for(int idim2=1; idim2<dim; idim2++) {
	  B[indexVAR[dim+i]][indexVAR[dim+(i+idim2)%dim]].resize(nve*nve);
	  memset(&B[indexVAR[dim+i]][indexVAR[dim+(i+idim2)%dim]][0],0,nve*nve*sizeof(double));
	}
      }
    }

    //Pressure gradient (fluid and solid)
    for(int i=0; i<dim; i++) {
      B[indexVAR[dim+i]][indexVAR[2*dim]].resize(nve*nve1);
      memset(&B[indexVAR[dim+i]][indexVAR[2*dim]][0],0,nve*nve1*sizeof(double));
    }

    // Pressure Mass Matrix
    B[indexVAR[2*dim]][indexVAR[2*dim]].resize(nve1*nve1);
    memset(&B[indexVAR[2*dim]][indexVAR[2*dim]][0],0,nve1*nve1*sizeof(double));
    
    if (flag_mat==2) { //initialization for fluid only
      // Fluid Continuity Matrix: divergence of the velocity
      for(int i=0; i<dim; i++) {
	B[indexVAR[2*dim]][indexVAR[dim+i]].resize(nve1*nve);
	memset(&B[indexVAR[2*dim]][indexVAR[dim+i]][0],0,nve1*nve*sizeof(double));
      }
    }
    else{ // initialization for solid only
      // Kinematic relation
      for(int i=0; i<dim; i++) {
	B[indexVAR[i]][indexVAR[dim+i]].resize(nve*nve);
	memset(&B[indexVAR[i]][indexVAR[dim+i]][0],0,nve*nve*sizeof(double));
      }
      // Solid Continuity Matrix: divergence of the displacemnet
      for(int i=0; i<dim; i++) {
	B[indexVAR[2*dim]][indexVAR[i]].resize(nve1*nve);
	memset(&B[indexVAR[2*dim]][indexVAR[i]][0],0,nve1*nve*sizeof(double));
      }
    }
    
    // ----------------------------------------------------------------------------------------
    // coordinates, displacement, velocity dofs
    
    metis_node2.resize(nve);
    metis_node1.resize(nve1);
    solidmark.resize(nve);
    phi.resize(nve);
    phi_hat.resize(nve);
    gradphi.resize(nve*dim);
    gradphi_hat.resize(nve*dim);
        
    for(int i=0;i<dim;i++){
      vx[i].resize(nve);
      vx_hat[i].resize(nve);
    }
    
    for (unsigned i=0;i<nve;i++) {
      // gambit nodes
      unsigned inode=myel->GetElementVertexIndex(kel,i)-1u;
      // dof metis
      unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
      metis_node2[i]=inode_Metis;
      
      //unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
      // flag to know if the node "inode" lays on the fluid-solid interface
      solidmark[i]=myel->GetNodeRegion(inode); // to check
      for(int j=0; j<dim; j++) {
	//Updated coordinates (Moving frame)
        vx[j][i]= (*mymsh->_coordinate->_Sol[j])(inode_Metis) + (*mysolution->_Sol[indVAR[j]])(inode_Metis);
	//Fixed coordinates (Reference frame)
	vx_hat[j][i]= (*mymsh->_coordinate->_Sol[j])(inode_Metis);  
	// displacement dofs
	dofsVAR[j][i]= myLinEqSolver->GetKKDof(indVAR[j],indexVAR[j],inode); 
	// velocity dofs
	dofsVAR[j+dim][i]= myLinEqSolver->GetKKDof(indVAR[j+dim],indexVAR[j+dim],inode);   
      }
    }

    // pressure dofs
    for (unsigned i=0;i<nve1;i++) {
      unsigned inode=(order_ind1<3)?(myel->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
      metis_node1[i]=mymsh->GetMetisDof(inode,SolType[2*dim]);
      dofsVAR[2*dim][i]=myLinEqSolver->GetKKDof(indVAR[2*dim],indexVAR[2*dim],inode);
    }
    // ----------------------------------------------------------------------------------------
       
    if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {
      /// *** Gauss point loop ***
      for (unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->GetGaussPointNumber(); ig++) {

	// *** get Jacobian and test function and test function derivatives in the moving frame***
	(ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->*(ml_prob._ml_msh->_finiteElement[kelt][order_ind2])->Jacobian_ptr)(vx,ig,Weight,phi,gradphi);
	(ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->*(ml_prob._ml_msh->_finiteElement[kelt][order_ind2])->Jacobian_ptr)(vx_hat,ig,Weight_hat,phi_hat,gradphi_hat);
	phi1=ml_prob._ml_msh->_finiteElement[kelt][order_ind1]->GetPhi(ig);
	if (flag_mat==2) Weight_nojac = ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->GetGaussWeight(ig);

	// ---------------------------------------------------------------------------
	// displacement and velocity
	for(int i=0; i<2*dim; i++){
	  SolVAR[i]=0.;
	  for(int j=0; j<dim; j++) {
	    GradSolVAR[i][j]=0.;
	    if(i<dim){
	      GradSolhatVAR[i][j]=0.;
	    }
	  }
	    
	  for (unsigned inode=0; inode<nve; inode++) {
	    unsigned sol_dof = metis_node2[inode];
	      
	    double soli = (*mysolution->_Sol[indVAR[i]])(sol_dof);
	    SolVAR[i]+=phi[inode]*soli;
	      
	    for(int j=0; j<dim; j++) {
	      GradSolVAR[i][j]+=gradphi[inode*dim+j]*soli;
	      if(i<dim){ 
		GradSolhatVAR[i][j]   +=gradphi_hat[inode*dim+j]*soli;
	      }
	    }	      
	  }
	}  
		  
	// pressure
	SolVAR[2*dim]=0.;
	for (unsigned inode=0; inode<nve1; inode++) {
	  double soli = (*mysolution->_Sol[indVAR[2*dim]])(metis_node1[inode]);
	  SolVAR[2*dim]+=phi1[inode]*soli;
	}
  
	// ---------------------------------------------------------------------------
	 	  
	//BEGIN FLUID ASSEMBLY ============
	  
	if(flag_mat==2){
	  
          { // Laplace operator + adection operator + Mass operator

            const double *gradfi=&gradphi[0];
            const double *fi=&phi[0];

            // *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim,fi++) {

              //BEGIN RESIDUALS A + Bt block ===========================
 	      
	      // begin redidual Laplacian ALE map
	      double LapmapVAR[3] = {0., 0., 0.};
 	      for(int idim=0; idim<dim; idim++) {
 		for(int jdim=0; jdim<dim; jdim++) {
 	          LapmapVAR[idim] += (_mu_ale[jdim]*GradSolVAR[idim][jdim]*gradphi[i*dim+jdim]);
 		}
 	      }
      	      for(int idim=0; idim<dim; idim++) {
 	        Rhs[indexVAR[idim]][i]+=(!solidmark[i])*(-LapmapVAR[idim]*Weight_nojac);
 	      }
	      // end redidual Laplacian ALE map

	      // begin redidual Navier-Stokes    
	      double LapvelVAR[3]={0.,0.,0.};
	      double AdvaleVAR[3]={0.,0.,0.};
	      for(int idim=0.; idim<dim; idim++) {
		for(int jdim=0.; jdim<dim; jdim++) {
		  LapvelVAR[idim]+=GradSolVAR[dim+idim][jdim]*gradphi[i*dim+jdim];
		  AdvaleVAR[idim]+=SolVAR[dim+jdim]*GradSolVAR[dim+idim][jdim]*phi[i];
		}
	      }
	      // Residual Momentum equations 
	      for(int idim=0; idim<dim; idim++) {
		Rhs[indexVAR[dim+idim]][i]+= (-AdvaleVAR[idim]			 // advection term	
					      -IRe*LapvelVAR[idim]		 // viscous dissipation
					      +SolVAR[2*dim]*gradphi[i*dim+idim] // pressure gradient
					      )*Weight;
	      // End redidual Navier-Stokes    
	      }
	      //END RESIDUALS A + Bt block ===========================
	      
	      //BEGIN A block ===========================
	      const double *gradfj=&gradphi[0];
              const double *fj=&phi[0];
              //  *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim,fj++) {

		// begin Laplacian ALE map
		double Lap_ale=0.;
		for(int jdim=0; jdim<dim; jdim++) {
		  Lap_ale+=_mu_ale[jdim]*(*(gradfi+jdim))*(*(gradfj+jdim));
		}	
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[idim]][indexVAR[idim]][i*nve+j] += (!solidmark[i])*Lap_ale*Weight_nojac;  
		}
		// end Laplacian ALE map
		
                //Laplacian
                double Lap=0.;
		for(int jdim=0; jdim<dim; jdim++) {
		  Lap+=(*(gradfi+jdim))*(*(gradfj+jdim))*Weight;
		}
               
                //advection term I
		double Adv1 = 0.;
		for(int jdim=0; jdim<dim; jdim++) {
		  Adv1+= SolVAR[dim+jdim]*(*(gradfj+jdim))*(*(fi))*Weight;
		}
		
		//advection term II
		double Adv2 = ((*fi))*((*fj))*Weight;
				
                for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[dim+idim]][i*nve+j] += (IRe*Lap + Adv1);
		  // Advection term II
                  if(nwtn_alg== true) {
		    B[indexVAR[dim+idim]][indexVAR[dim+idim]][i*nve+j]  += Adv2*GradSolVAR[dim+idim][idim];
		    for(unsigned jdim=1; jdim<dim; jdim++) {
		      B[indexVAR[dim+idim]][indexVAR[dim + (idim+jdim)%dim]][i*nve+j] += Adv2*GradSolVAR[dim+idim][(idim+jdim)%dim];
		    }
		  }
	        }
              } // end phi_j loop
              //END A block ===========================
            } // end phi loop
          } // end A 
	 
	 
	  //BEGIN Bt block ===========================
          { //Gradient of Pressure operator

            const double *gradfi=&gradphi[0];
            const double *fi=phi1;
            /// *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim,fi++) {
              const double *fj=phi1;
              /// *** phi_j loop ***
              for (unsigned j=0; j<nve1; j++,fj++) {
                for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[2*dim]][i*nve1+j] -= ((*(gradfi+idim))*(*fj))*Weight;
		}
              } // end phi_j loop
            } // end phi_i loop
          } // End Bt
          //END Bt block ===========================
          
          //BEGIN B block ===========================
          {  	    
	    //divergence of the velocity
	    double div_vel=0.;
	    for(int i=0; i<dim; i++) {
	      div_vel+=GradSolVAR[dim+i][i];
	    }
	    // Divergence of the Velocity operator
            const double *fi=phi1;
            // *** phi_i loop ***
            for (unsigned i=0; i<nve1; i++,fi++) {

              //BEGIN RESIDUALS B block ===========================
              Rhs[indexVAR[2*dim]][i] += -(-((*fi))*div_vel)*Weight;
              //END RESIDUALS  B block ===========================

              const double *gradfj=&gradphi[0];
              // *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim) {
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[2*dim]][indexVAR[idim+dim]][i*nve+j] -= ((*fi)*(*(gradfj+idim)))*Weight;
		}
              }
            }
          }
          //END B block ===========================
	}   
	//END FLUID ASSEMBLY ============
	//*******************************************************************************************************
	//BEGIN SOLID ASSEMBLY ============
	  
	else{
	  //------------------------------------------------------------------------------------------------------------
          if (solid_model==0) {
	    double e[3][3];
	    //computation of the stress tensor
	    for(int i=0;i<dim;i++){
	      for(int j=0;j<dim;j++){
		e[i][j]=0.5*(GradSolhatVAR[i][j]+GradSolhatVAR[j][i]);
	      }
	    }
	    I_e=0;
	    for(int i=0;i<dim;i++){
	      I_e += e[i][i];
	    }
	    
	    for (int i=0; i<dim; i++) {
              for (int j=0; j<dim; j++) {
                //incompressible
                Cauchy[i][j] = 2*mus*e[i][j];
	      }
            }
          }
	  
          else if (solid_model==1) {
	    double F[3][3]={{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
	    double b_left[3][3];      
	    
	    for(int i=0;i<dim;i++){
	      for(int j=0;j<dim;j++){
		F[i][j]+=GradSolhatVAR[i][j];
	      }
	    }
	   	    	    	    
	    Jnp1_hat =  F[0][0]*F[1][1]*F[2][2] + F[0][1]*F[1][2]*F[2][0] + F[0][2]*F[1][0]*F[2][1]
		      - F[2][0]*F[1][1]*F[0][2] - F[2][1]*F[1][2]*F[0][0] - F[2][2]*F[1][0]*F[0][1];		      
	   
	    // computation of the the three deformation tensor b
	    for (int I=0; I<3; ++I) {
	      for (int J=0; J<3; ++J) {
		b_left[I][J]=0.;
		for (int K=0; K<3; ++K) {
		  //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
		  b_left[I][J] += F[I][K]*F[J][K];
		}
		Cauchy[I][J] = (mus/Jnp1_hat)*(b_left[I][J] - Id2th[I][J]);
	      }
	    }

	    
	    
	    I_bleft = b_left[0][0] + b_left[1][1] + b_left[2][2];
	    
	    //compressible case
	    //             for (int ii=0; ii<3; ++ii) {
	    //               for (int jj=0; jj<3; ++jj) {
	    //                 //Cauchy stress tensor
	    //                 Cauchy[ii][jj] = (_mus/Jnp1_hat)*(b_left[ii][jj] - Id2th[ii][jj]) + (_lambda/Jnp1_hat)*log(Jnp1_hat)*Id2th[ii][jj];
	    //                 for (int k=0; k<3; ++k) {
	    //                   for (int l=0; l<3; ++l) {
	    //                   }
	    //                 }
	    //               }
	    //             }
	    
	    //for the incompressible(nearly incompressible) case
	    for (int ii=0; ii<3; ++ii) {
	      for (int jj=0; jj<3; ++jj) {
		for (int kk=0; kk<3; ++kk) {
		  for (int ll=0; ll<3; ++ll) {
		    C_mat[ii][jj][kk][ll] = 2.*mus*pow(Jnp1_hat,-1.6666666666666)*(
										  0.333333333333*I_bleft*Id2th[ii][kk]*Id2th[jj][ll]              //1/3*I_c*i
										  // 	                        +0.111111111111*I_C*Id2th[i][j]*Id2th[k][l]             //1/9*I_b*IxI
										  // 				-0.333333333333*b_left[i][j]*Id2th[k][l]                //-1/3*b*I
										  // 				-0.333333333333*Id2th[i][j]*b_left[k][l]                //-1/3*b*I
										  )
		      -SolVAR[2*dim]*(Id2th[ii][jj]*Id2th[kk][ll]-2.*Id2th[ii][kk]*Id2th[jj][ll] );  // -p(IxI-2i)
		  }
		}
	      }
	    }
	  }
          //----------------------------------------------------------------------------------------------------------------------------

          /////////////

          ///Mass + Stiffness operator
          {
            const double *gradfi=&gradphi[0];
            const double *fi=&phi[0];

            /// *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim,fi++) {


              //BEGIN RESIDUALS A + Bt block ===========================
	      
	      // Residual ALE equations
 	      for(int idim=0; idim<dim; idim++) {
	        Rhs[indexVAR[idim]][i] += (-phi[i]*(-SolVAR[dim+idim] ))*Weight_hat;
              }
              
              double CauchyDIR[3]={0.,0.,0.};
	      for(int idim=0.; idim<dim; idim++) {
		for(int jdim=0.; jdim<dim; jdim++) {
		  CauchyDIR[idim]+= gradphi[i*dim+jdim]*Cauchy[idim][jdim];
		}
	      }

              // Residual Momentum equations
              for(int idim=0; idim<dim; idim++) {
	        Rhs[indexVAR[dim+idim]][i] += (
					       phi[i]*_gravity[idim]*Weight_hat
					       -CauchyDIR[idim]*Weight
					       +SolVAR[2*dim]*gradphi[i*dim+idim]*Weight
					       );

              }
              
              //---------------------------------------------------------------------------------------------------------------------------------

              //END RESIDUALS A + Bt block ===========================

              const double *gradfj=&gradphi[0];
              const double *fj=&phi[0];
              // *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim,fj++) {

		// tangent stiffness matrix
                for (int idim=0; idim<dim; ++idim) {
                  for (int jdim=0; jdim<dim; ++jdim) {
                    tg_stiff_matrix[idim][jdim] = 0.;
                    for (int kdim=0; kdim<dim; ++kdim) {
                      for (int ldim=0; ldim<dim; ++ldim) {
                        tg_stiff_matrix[idim][jdim] += (*(gradfi+kdim))*0.25*(
									       C_mat[idim][kdim][jdim][ldim]+C_mat[idim][kdim][ldim][jdim]
									      +C_mat[kdim][idim][jdim][ldim]+C_mat[kdim][idim][ldim][jdim]
									      )*(*(gradfj+ldim));
                      }
                    }
                  }
                }
                
                //geometric tangent stiffness matrix
                double geom_tg_stiff_matrix = 0.;
                for(int idim=0; idim<dim; ++idim) {
                  for(int jdim=0; jdim<dim; ++jdim) {
                    geom_tg_stiff_matrix += (*(gradfi+idim))*Cauchy[idim][jdim]*(*(gradfj+jdim));
                  }
                }

                /// Stiffness operator -- Elasticity equation (Linear or not)
                for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[idim]][i*nve+j] += geom_tg_stiff_matrix*Weight;
 		  for(int jdim=0; jdim<dim; jdim++) {
 		    B[indexVAR[dim+idim]][indexVAR[jdim]][i*nve+j] += tg_stiff_matrix[idim][jdim]*Weight;
 		  }
 		}

                /// Kinematic equation v = du/dt --> In the steady state we write \deltau^n+1 - \deltav^n+1 = v - 0
                //   
		for(int idim=0; idim<dim; idim++) {
		  // -(v_n+1,eta)
		  B[indexVAR[0+idim]][indexVAR[dim+idim]][i*nve+j] -= (*(fi))*(*(fj))*Weight_hat;
		  //  (u_n+1,eta)
		  B[indexVAR[0+idim]][indexVAR[0+idim]][i*nve+j] += (*(fi))*(*(fj))*Weight_hat;
		}
	      }
            }
          }
          ////////////
          { ///Gradient of Pressure
            const double *gradfi=&gradphi[0];
            // *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim) {
              const double *fj=phi1;
              // *** phi_j loop ***
              for (unsigned j=0; j<nve1; j++,fj++) {
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[2*dim]][i*nve1+j] -= ((*(gradfi+idim))*(*fj))*Weight;
		}
              }
            }
          }
          ////////////
          { ///Divergence of the Displacement
            const double *fi=phi1;
            // *** phi_i loop ***
            for (unsigned i=0; i<nve1; i++,fi++) {

              //BEGIN RESIDUALS B block ===========================

              if (solid_model==0) {
                Rhs[indexVAR[2*dim]][i] += -(-((*fi))*(I_e + (1./lambda)*SolVAR[2*dim] ) )*Weight_hat;
              }
              else if (solid_model==1) {
                Rhs[indexVAR[2*dim]][i] += -(-((*fi))*( log(Jnp1_hat)/Jnp1_hat + (1./lambda)*SolVAR[2*dim] ) )*Weight_hat;
              }

              //END RESIDUALS B block ===========================

              const double *gradfj=&gradphi[0];
              // *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim) {
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[2*dim]][indexVAR[idim]][i*nve+j] -= ((*fi)*(*(gradfj+idim)))*Weight;
		}
              }
            }
          }
          //  /////////////
          {  ///Pressure Mass term
            const double *fi=phi1;
            // *** phi_i loop ***
            for (unsigned i=0; i<nve1; i++,fi++) {
              const double *fj=phi1;
              // *** phi_j loop ***
              for (unsigned j=0; j<nve1; j++,fj++) {
                B[indexVAR[2*dim]][indexVAR[2*dim]][i*nve1+j] -= (1./lambda)*((*fi)*(*fj))*Weight_hat;
              }
            }
          }  //end pressure mass term
	  //---------------------------------------------------------------------------------------------------------------------------------
	}  
	//END SOLID ASSEMBLY ============
      }
    }

    //BEGIN local to global assembly 
    // ALE mapping
    for(int i=0; i<dim; i++) {
      myRES->add_vector_blocked(Rhs[indexVAR[i]],dofsVAR[i]);
      myKK ->add_matrix_blocked(B[indexVAR[i]][indexVAR[i]],dofsVAR[i],dofsVAR[i]);  
      if(flag_mat!=2){ //Solid only
	myKK->add_matrix_blocked(B[indexVAR[i]][indexVAR[i+dim]],dofsVAR[i],dofsVAR[i+dim]);
      }
    }
    
    // Momentum equation
    for(int i=0; i<dim; i++) {
      myRES->add_vector_blocked(Rhs[indexVAR[dim+i]],dofsVAR[dim+i]);
      myKK->add_matrix_blocked(B[indexVAR[dim+i]][indexVAR[dim+i]],dofsVAR[dim+i],dofsVAR[dim+i]);
      if(nwtn_alg== true){
	for(unsigned idim2=1; idim2<dim; idim2++) {
	  myKK->add_matrix_blocked(B[indexVAR[dim+i]][indexVAR[dim+(i+idim2)%dim]],dofsVAR[dim+i],dofsVAR[dim+(i+idim2)%dim]);  
	}
      }
      myKK->add_matrix_blocked(B[indexVAR[dim+i]][indexVAR[2*dim]],dofsVAR[dim+i],dofsVAR[2*dim]);
      for(int j=0; j<dim; j++) {
	myKK->add_matrix_blocked(B[indexVAR[dim+i]][indexVAR[j]],dofsVAR[dim+i],dofsVAR[j]);
      }
    }
   
    //P-continuity equation
    myRES->add_vector_blocked(Rhs[indexVAR[2*dim]],dofsVAR[2*dim]);
    for(int i=0; i<dim; i++) {
      if(flag_mat==2){ //Fluid only
	myKK->add_matrix_blocked(B[indexVAR[2*dim]][indexVAR[dim+i]],dofsVAR[2*dim],dofsVAR[dim+i]);
      }
      else{ //Solid only
	myKK->add_matrix_blocked(B[indexVAR[2*dim]][indexVAR[i]],dofsVAR[2*dim],dofsVAR[i]);  
      } 
    }
    myKK->add_matrix_blocked(B[indexVAR[2*dim]][indexVAR[2*dim]],dofsVAR[2*dim],dofsVAR[2*dim]);
    //END local to global assembly
   
  } //end list of elements loop

  // close residual vector and matrix

  myKK->close();
  myRES->close();
  
  // *************************************
  end_time=clock();
  AssemblyTime+=(end_time-start_time);
  // ***************** END ASSEMBLY RESIDUAL + MATRIX ********************/

//}


