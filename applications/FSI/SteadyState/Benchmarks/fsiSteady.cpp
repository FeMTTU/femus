#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"
#include "FemusInit.hpp"
#include "SparseMatrix.hpp"
#include "FElemTypeEnum.hpp"
#include "Files.hpp"
#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "../include/FSISteadyStateAssembly.hpp"

using namespace std;
using namespace femus;

bool SetBoundaryConditionTurek_2D_FSI_and_solid(const std::vector < double >& x,const char name[],
						double &value, const int FaceName, const double = 0.);
bool SetBoundaryConditionBathe_2D_FSI(const std::vector < double >& x,const char name[],
				      double &value, const int FaceName, const double = 0.);
bool SetBoundaryConditionBathe_3D_FSI_and_fluid(const std::vector < double >& x,const char name[],
						double &value, const int facename, const double time);

bool SetBoundaryConditionBathe_3D_solid(const std::vector < double >& x,const char name[],
					double &value, const int facename, const double time);
bool SetBoundaryConditionComsol_2D_FSI(const std::vector < double >& x,const char name[],
				       double &value, const int FaceName, const double = 0.);

double InitalValueU(const std::vector < double >& x);

bool SetRefinementFlag(const std::vector < double >& x, const int &ElemGroupNumber,const int &level);

//------------------------------------------------------------------------------------------------------------------

int main(int argc,char **args) {

  // ******* Init Petsc-MPI communicator *******
  FemusInit mpinit(argc,args,MPI_COMM_WORLD);

  // process options 
  int dimension=2;
  char infile[256] = "";
  size_t len_infile_name = 256;
  int simulation = 1;
  double Lref=1., Uref=1., rhof=1., muf=1., rhos=1., ni=0., E=1.;
  int numofmeshlevels = 1;
  int numofrefinements = 0;
  std::string gauss_integration_order = "fifth";
  
  // ******* reading input parameters *******
  PetscOptionsBegin(PETSC_COMM_WORLD, "", "FSI steady problem options", "Unstructured mesh");
  
  PetscOptionsInt("-dim", "The dimension of the problem", "fsiSteady.cpp", dimension, &dimension, NULL);
  printf(" dim: %i\n", dimension);
  
  PetscOptionsInt("-sim", "The type of the Simulation", "fsiSteady.cpp", simulation, &simulation, NULL);
  printf(" sim: %i\n", simulation);
  
  PetscOptionsInt("-nlevel", "The number of mesh levels", "fsiSteady.cpp", numofmeshlevels , &numofmeshlevels, NULL);
  printf(" nlevel: %i\n", numofmeshlevels);
  
  PetscOptionsInt("-nrefinement", "The number of refinements", "fsiSteady.cpp", numofrefinements , &numofrefinements, NULL);
  printf(" nrefinement: %i\n", numofrefinements);
  
  PetscOptionsString("-input", "The name of the input file", "fsiSteady.cpp", "./mesh.neu", infile, len_infile_name, NULL);
  printf(" input: %s\n", infile);
  
  PetscOptionsReal("-rhof", "The density of the fluid", "fsiSteady.cpp", rhof, &rhof, NULL);
  printf(" rhof: %f\n", rhof);
  
  PetscOptionsReal("-rhof", "The density of the solid", "fsiSteady.cpp", rhos, &rhos, NULL);
  printf(" rhos: %f\n", rhos);
  
  PetscOptionsReal("-E", "The young module of the solid", "fsiSteady.cpp", E, &E, NULL);
  printf(" E: %f\n", E);
  
  PetscOptionsReal("-muf", "The viscosity of the fluid", "fsiSteady.cpp", muf, &muf, NULL);
  printf(" muf: %f\n", muf);
  
  PetscOptionsReal("-ni", "The Poisson coefficient of the Solid", "fsiSteady.cpp", ni, &ni, NULL);
  printf(" ni: %f\n", ni);
  
  printf("\n");
  
  PetscOptionsEnd();
  
  
  // ******* Init multilevel mesh from mesh.neu file *******
  MultiLevelMesh ml_msh(numofrefinements, numofrefinements, infile, gauss_integration_order.c_str(), Lref, 
			SetRefinementFlag);

  ml_msh.EraseCoarseLevels(numofrefinements - numofmeshlevels);

  ml_msh.PrintInfo();

  // mark Solid nodes
  ml_msh.MarkStructureNode();
  
  
  // ******* Fluid and Solid Parameters *******
  Parameter par(Lref,Uref);

  // Generate Solid Object
  Solid solid;
  solid = Solid(par,E,ni,rhos,"Mooney-Rivlin");

  cout << "Solid properties: " << endl;
  cout << solid << endl;

  // Generate Fluid Object
  Fluid fluid(par,muf,rhof,"Newtonian");
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;
  

  // ******* Init multilevel solution ******
  MultiLevelSolution ml_sol(&ml_msh);

  // ******* Add solution variables to multilevel solution and pair them *******
  ml_sol.AddSolution("DX",LAGRANGE,SECOND,1);
  ml_sol.AddSolution("DY",LAGRANGE,SECOND,1);
  if (dimension==3) ml_sol.AddSolution("DZ",LAGRANGE,SECOND,1);

  ml_sol.AddSolution("U",LAGRANGE,SECOND,1);
  ml_sol.AddSolution("V",LAGRANGE,SECOND,1);
  if (dimension==3) ml_sol.AddSolution("W",LAGRANGE,SECOND,1);

  // Pair each velocity variable with the corresponding displacement variable
  ml_sol.PairSolution("U","DX"); // Add this line
  ml_sol.PairSolution("V","DY"); // Add this line
  if (dimension==3) ml_sol.PairSolution("W","DZ"); // Add this line

  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  ml_sol.AddSolution("P",DISCONTINOUS_POLYNOMIAL,FIRST,1);
  ml_sol.AssociatePropertyToSolution("P","Pressure",false); // Add this line

  // ******* Initialize solution *******
  ml_sol.Initialize("All");
  if (1 == simulation )
    ml_sol.Initialize("U",InitalValueU);

  // ******* Set boundary functions *******
  if(1==simulation || 2==simulation)
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionTurek_2D_FSI_and_solid);
  else if( 3==simulation)
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionBathe_2D_FSI);
  else if (4==simulation || 6==simulation)
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionBathe_3D_FSI_and_fluid);
  else if (5==simulation)
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionBathe_3D_solid);
  else if (7 == simulation)
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionComsol_2D_FSI);

  // ******* Set boundary conditions *******
  ml_sol.GenerateBdc("DX","Steady");
  ml_sol.GenerateBdc("DY","Steady");
  if (dimension==3) ml_sol.GenerateBdc("DZ","Steady");
  ml_sol.GenerateBdc("U","Steady");
  ml_sol.GenerateBdc("V","Steady");
  if (dimension==3) ml_sol.GenerateBdc("W","Steady");
  ml_sol.GenerateBdc("P","Steady");

  
  // ******* Define the FSI Multilevel Problem *******
  MultiLevelProblem ml_prob(&ml_sol);
  // Add fluid object
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;
  // Add Solid Object
  ml_prob.parameters.set<Solid>("Solid") = solid;


  // ******* Add FSI system to the MultiLevel problem *******
  MonolithicFSINonLinearImplicitSystem & system = ml_prob.add_system<MonolithicFSINonLinearImplicitSystem> ("Fluid-Structure-Interaction");
  system.AddSolutionToSystemPDE("DX");
  system.AddSolutionToSystemPDE("DY");
  if (dimension==3) system.AddSolutionToSystemPDE("DZ");
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  if (dimension==3) system.AddSolutionToSystemPDE("W");
  system.AddSolutionToSystemPDE("P");

  // ******* System Fluid-Structure-Interaction Assembly *******
  system.SetAssembleFunction(FSISteadyStateAssembly);


  // Solver settings
  // ******* set Non-linear MG-Solver *******
  system.SetMaxNumberOfLinearIterations(3);
  system.SetLinearConvergenceTolerance(1.e-10);
  
  system.SetMaxNumberOfNonLinearIterations(10);
  system.SetNonLinearConvergenceTolerance(1.e-9);
 
  // Set type of multigrid preconditioner (V-cycle, F-cycle)
  system.SetMgType(F_CYCLE);
  system.SetNumberPreSmoothingStep(0);
  system.SetNumberPostSmoothingStep(2);

  // ******* Set Smoother *******
  // Set Preconditioner of the smoother (name to be changed)
  system.SetMgSmoother(ASM_SMOOTHER);
  
  // System init
  system.init();
  
  // Set the preconditioner for each ASM block
  system.SetPreconditionerFineGrids(ILU_PRECOND);
  // Set block size for the ASM smoother
  system.SetElementBlockNumber(2);
  // Set Solver of the smoother (name to be changed)
  system.SetSolverFineGrids(RICHARDSON);
  
  // set the tolerances for the GMRES outer solver
  system.SetTolerances(1.e-12,1.e-20,1.e+50,5);
  
  // ******* Add variables to be solved *******
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");

  // ******* Set the last (1) variables in system (i.e. P) to be a schur variable *******
  system.SetNumberOfSchurVariables(1);


  // ******* Solve *******
  std::cout << std::endl;
  
  std::cout << " *********** Solving... ************  " << std::endl;
  
  system.MGsolve();

  // ******* Print solution *******
  ml_sol.SetWriter(VTK);

  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  mov_vars.push_back("DZ");
  ml_sol.GetWriter()->SetMovingMesh(mov_vars);

  std::vector<std::string> print_vars;
  print_vars.push_back("All");

  ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR,"biquadratic",print_vars);

  // ******* Clear all systems *******
  ml_prob.clear();
  
  return 0;
}


bool SetRefinementFlag(const std::vector < double >& x, const int &elemgroupnumber,const int &level) {
  bool refine=0;

  //refinemenet based on elemen group number
  if (elemgroupnumber==5) refine=1;
  if (elemgroupnumber==6) refine=1;
  if (elemgroupnumber==7 && level<5) refine=1;

  return refine;

}

double InitalValueU(const std::vector < double >& x) {
  double xc = 0.2;
  double yc = 0.2;
  double r = 0.05;
  double r2 = r*r;
  double xMxc2 = (x[0]-xc)*(x[0]-xc);
  double OMxc2 = (0.-xc)*(0.-xc);
  double yMyc2 = (x[1]-yc)*(x[1]-yc);

  double H = 0.41;
  double L = 2.5;
  double um = 0.2;
  return (xMxc2+yMyc2-r2)/(OMxc2+yMyc2-r2)*(1.5*um*4.0/0.1681*x[1]*(H-x[1]))*exp(-L*x[0]);
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionTurek_2D_FSI_and_solid(const std::vector < double >& x,const char name[], double &value, const int facename, const double time) {
  bool test=1; //dirichlet
  value=0.;
  if(!strcmp(name,"U")) {
    if(1==facename){   //inflow
      test=1;
      double um = 0.2;
      value=1.5*um*4.0/0.1681*x[1]*(0.41-x[1]);
    }
    else if(2==facename ){  //outflow
      test=0;
      //    test=1;
      value=0.;
    }
    else if(3==facename ){  // no-slip fluid wall
      test=1;
      value=0.;
    }
    else if(4==facename ){  // no-slip solid wall
      test=1;
      value=0.;
    }
    else if(6==facename ){   // beam case zero stress
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"V")){
    if(1==facename){            //inflow
      test=1;
      value=0.;
    }
    else if(2==facename ){      //outflow
      test=0;
      //    test=1;
      value=0.;
    }
    else if(3==facename ){      // no-slip fluid wall
      test=1;
      value=0;
    }
    else if(4==facename ){      // no-slip solid wall
      test=1;
      value=0.;
    }
    else if(6==facename ){   // beam case zero stress
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"P")){
    if(1==facename){
      test=0;
      value=0.;
    }
    else if(2==facename ){
      test=0;
      value=0.;
    }
    else if(3==facename ){
      test=0;
      value=0.;
    }
    else if(4==facename ){
      test=0;
      value=0.;
    }
    else if(6==facename ){   // beam case zero stress
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"DX")){
    if(1==facename){         //inflow
      test=1;
      value=0.;
    }
    else if(2==facename ){   //outflow
      test=1;
      value=0.;
    }
    else if(3==facename ){   // no-slip fluid wall
      test=0; //0
      value=0.;
    }
    else if(4==facename ){   // no-slip solid wall
      test=1;
      value=0.;
    }
    else if(6==facename ){   // beam case zero stress
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"DY")){
    if(1==facename){         //inflow
      test=0; // 0
      value=0.;
    }
    else if(2==facename ){   //outflow
      test=0; // 0
      value=0.;
    }
    else if(3==facename ){   // no-slip fluid wall
      test=1;
      value=0.;
    }
    else if(4==facename ){   // no-slip solid wall
      test=1;
      value=0.;
    }
    else if(6==facename ){   // beam case zero stress
      test=0;
      value=0.;
    }
  }
  return test;
}

bool SetBoundaryConditionBathe_2D_FSI(const std::vector < double >& x,const char name[], double &value, const int facename, const double time) {
  bool test=1; //dirichlet
  value=0.;
  if(!strcmp(name,"U")) {
    if(1==facename){   //top
      test=0;
      value=0;
    }
    else if(2==facename ){  //top side
      test=0;
      value=0.;
    }
    else if(3==facename ){  //top bottom
      test=1;
      value=0;
    }
    else if(4==facename ){  //solid side
      test=1;
      value=0;
    }
    else if(5==facename ){  //bottom side
      test=1;
      value=0;
    }
    else if(6==facename ){  //bottom
      test=0;
      value=200000;
    }
  }
  else if(!strcmp(name,"V")){
    if(1==facename){   //top
      test=0;
      value=0;
    }
    else if(2==facename ){  //top side
      test=0;
      value=0.;
    }
    else if(3==facename ){  //top bottom
      test=1;
      value=0;
    }
    else if(4==facename ){  //solid side
      test=1;
      value=0;
    }
    else if(5==facename ){  //bottom side
      test=1;
      value=0;
    }
    else if(6==facename ){  //bottom
      test=0;
      value=0;
    }
  }
  else if(!strcmp(name,"P")){
    if(facename==facename){
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"DX")){
    if(1==facename){   //top
      test=0;
      value=0;
    }
    else if(2==facename ){  //top side
      test=1;
      value=0.;
    }
    else if(3==facename ){  //top bottom
      test=0;
      value=0;
    }
    else if(4==facename ){  //solid side
      test=1;
      value=0;
    }
    else if(5==facename ){  //bottom side
      test=1;
      value=0;
    }
    else if(6==facename ){  //bottom
      test=0;
      value=0;
    }
  }
  else if(!strcmp(name,"DY")){
    if(1==facename){   //top
      test=1;
      value=0;
    }
    else if(2==facename ){  //top side
      test=0;
      value=0.;
    }
    else if(3==facename ){  //top bottom
      test=1;
      value=0;
    }
    else if(4==facename ){  //solid side
      test=1;
      value=0;
    }
    else if(5==facename ){  //bottom side
      test=1;
      value=0;
    }
    else if(6==facename ){  //bottom
      test=1;
      value=0;
    }
  }

  return test;
}



bool SetBoundaryConditionBathe_3D_FSI_and_fluid(const std::vector < double >& x,const char name[], double &value, const int facename, const double time) {
  bool test=1; //dirichlet
  value=0.;

  if(!strcmp(name,"U")) {
    if(1==facename){   //inflow
      //test=1;
      //double r=sqrt(y*y+z*z);
      //value=1000*(0.05-r)*(0.05+r);
      test=0;
      value=15*1.5*1000;
    }
    else if(2==facename){  //outflow
      test=0;
      value=13*1.5*1000;
    }
    else if(3==facename || 4==facename ){  // clamped solid
      test=1;
      value=0.;
    }
    else if(5==facename ){  // free solid
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"V")){
    if(1==facename){   //inflow
      test=0;
      value=0;
    }
    else if(2==facename){  //outflow
      test=0;
      value=0;
    }
    else if(3==facename || 4==facename ){  // clamped solid
      test=1;
      value=0.;
    }
    else if(5==facename ){  // free solid
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"W")){
    if(1==facename){   //inflow
      test=0;
      value=0;
    }
    else if(2==facename){  //outflow
      test=0;
      value=0;
    }
    else if(3==facename || 4==facename ){  // clamped solid
      test=1;
      value=0.;
    }
    else if(5==facename ){  // free solid
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"P")){
    if(1==facename){   //outflow
      test=0;
      value=0;
    }
    else if(2==facename){  //inflow
      test=0;
      value=0;
    }
    else if(3==facename || 4==facename ){  // clamped solid
      test=0;
      value=0.;
    }
    else if(5==facename ){  // free solid
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"DX")){
    if(1==facename){   //outflow
      test=1;
      value=0;
    }
    else if(2==facename){  //inflow
      test=1;
      value=0;
    }
    else if(3==facename || 4==facename ){  // clamped solid
      test=1;
      value=0.;
    }
    else if(5==facename ){  // free solid
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"DY")){
    if(1==facename){   //outflow
      test=0;
      value=0;
    }
    else if(2==facename){  //inflow
      test=0;
      value=0;
    }
    else if(3==facename || 4==facename ){  // clamped solid
      test=1;
      value=0.;
    }
    else if(5==facename ){  // free solid
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"DZ")){
    if(1==facename){   //outflow
      test=0;
      value=0;
    }
    else if(2==facename){  //inflow
      test=0;
      value=0;
    }
    else if(3==facename || 4==facename ){  // clamped solid
      test=1;
      value=0.;
    }
    else if(5==facename ){  // free solid
      test=0;
      value=0.;
    }
  }

  return test;
}

bool SetBoundaryConditionBathe_3D_solid(const std::vector < double >& x,const char name[], double &value, const int facename, const double time) {
  bool test=1; //dirichlet
  value=0.;

  if(!strcmp(name,"U")) {
    if(2==facename){  //stress
      test=0;
      value=15*1.5*1000;
    }
    else if(3==facename || 4==facename ){  // clamped solid
      test=1;
      value=0.;
    }
    else if(5==facename ){  // free solid
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"V")){
    if(2==facename){  //stress
      test=0;
      value=0;
    }
    else if(3==facename || 4==facename ){  // clamped solid
      test=1;
      value=0.;
    }
    else if(5==facename ){  // free solid
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"W")){
    if(2==facename){  //stress
      test=0;
      value=0;
    }
    else if(3==facename || 4==facename ){  // clamped solid
      test=1;
      value=0.;
    }
    else if(5==facename ){  // free solid
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"P")){
    if(2==facename){  //stress
      test=0;
      value=0;
    }
    else if(3==facename || 4==facename ){  // clamped solid
      test=0;
      value=0.;
    }
    else if(5==facename ){  // free solid
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"DX")){
    if(2==facename){  //stress
      test=0;
      value=0;
    }
    else if(3==facename || 4==facename ){  // clamped solid
      test=1;
      value=0.;
    }
    else if(5==facename ){  // free solid
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"DY")){
    if(2==facename){  //stress
      test=0;
      value=0;
    }
    else if(3==facename || 4==facename ){  // clamped solid
      test=1;
      value=0.;
    }
    else if(5==facename ){  // free solid
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"DZ")){
    if(2==facename){  //stress
      test=0;
      value=0;
    }
    else if(3==facename || 4==facename ){  // clamped solid
      test=1;
      value=0.;
    }
    else if(5==facename ){  // free solid
      test=0;
      value=0.;
    }
  }

  return test;
}



//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionComsol_2D_FSI(const std::vector < double >& x,const char name[], double &value, const int FaceName, const double time) {
  bool test=1; //Dirichlet
  value=0.;
  //   cout << "Time bdc : " <<  time << endl;
  if (!strcmp(name,"U")) {
    if (1==FaceName) { //inflow
      test=1;
      //comsol Benchmark
      //value = (0.05*time*time)/(sqrt( (0.04 - time*time)*(0.04 - time*time) + (0.1*time)*(0.1*time) ))*y*(0.0001-y)*4.*100000000;
      value = 0.05*x[1]*(0.0001-x[1])*4.*100000000;
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
  return test;
}

