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
#include "CurrentElem.hpp"
#include "../include/FSISteadyStateAssembly.hpp"

double scale=1000.;

using namespace std;
using namespace femus;

bool SetBoundaryConditionTurek_2D_FSI_and_solid(const vector< double > &x, const char name[],
						double &value, const int FaceName, const double = 0.);
bool SetBoundaryConditionBathe_2D_FSI(const vector< double > &x, const char name[],
				      double &value, const int FaceName, const double = 0.);
bool SetBoundaryConditionBathe_3D_FSI_and_fluid(const vector< double > &x, const char name[],
						double &value, const int facename, const double time);

bool SetBoundaryConditionBathe_3D_solid(const vector< double > &x, const char name[],
					double &value, const int facename, const double time);
bool SetBoundaryConditionComsol_2D_FSI(const vector< double > &x, const char name[],
				       double &value, const int FaceName, const double = 0.);


bool SetRefinementFlag(const vector< double > &x, const int &ElemGroupNumber,const int &level);

//------------------------------------------------------------------------------------------------------------------

int main(int argc,char **args) {

  // ******* Init Petsc-MPI communicator *******
  FemusInit mpinit(argc,args,MPI_COMM_WORLD);

  //Files files;
  //files.CheckIODirectories();
  //files.RedirectCout();

  // ******* Extract the problem dimension and simulation identifier based on the inline input *******
  unsigned simulation;
  bool dimension2D;

  if(argc >= 2) {
    if( !strcmp("turek_2D_FSI",args[1])) {   /** FSI steady-state Turek benchmark */
      simulation=1;
      dimension2D=1;
    }
    else if( !strcmp("turek_2D_solid",args[1])) {    /** Solid Turek beam benchmark test. Beware: activate gravity in assembly */
      simulation=2;
      dimension2D=1;
    }
    else if( !strcmp("bathe_2D_FSI",args[1])){  /** Bathe 2D membrane benchmark */
      simulation=3;
      dimension2D=1;
    }
    else if( !strcmp("bathe_3D_FSI",args[1])){  /** Bathe 3D cylinder FSI benchmark */
      simulation=4;
      dimension2D=0;
    }
    else if( !strcmp("bathe_3D_solid",args[1])) { /** Bathe 3D solid, for debugging */
      simulation=5;
      dimension2D=0;
    }
    else if( !strcmp("bathe_3D_fluid",args[1])) { /** Bathe 3D fluid, for debugging */
      simulation=6;
      dimension2D=0;
    }
    else if( !strcmp("comsol_2D_FSI",args[1])) { /** Comsol 2D vertical beam benchmark */
      simulation=7;
      dimension2D=1;
    }

    else{
      cout << "wrong input arguments!\n";
      cout << "please specify the simulation you want to run, options are\n";
      cout << "turek_2D_FSI \n turek_2D_solid \n bathe_2D_FSI \n bathe_3D_FSI \n bathe_3D_solid \n bathe_3D_fluid \n comsol_2D_FSI \n";
      abort();
    }
  }
  else {
    cout << "no input arguments!\n";
    cout << "please specify the simulation you want to run, options are\n";
    cout << "turek_2D_FSI \n turek_2D_solid \n bathe_2D_FSI \n bathe_3D_FSI \n bathe_3D_solid \n bathe_3D_fluid \n comsol_2D_FSI \n";
    abort();
  }

  // ******* Extract the preconditioner type based on the inline input *******
  bool Gmres=0, Asm=0;
  if(argc >= 3) {
    if( !strcmp("gmres",args[2])) 	Gmres=1;
    else if( !strcmp("asm",args[2])) 	Asm=1;

    if(Gmres+Asm==0) {
      cout << "wrong input arguments!" << endl;
      abort();
    }
  }
  else {
    cout << "No input argument set default smoother = Asm" << endl;
    Asm=1;
  }


  // ******* Extract the mesh.neu file name based on the simulation identifier *******
  std::string infile;

  if(1 == simulation){
    infile = "./input/turek.neu";
  }
  else if(2 == simulation){
    infile = "./input/beam.neu";
  }
  else if(3 == simulation){
    infile = "./input/drum.neu";
  }
  else if(4 == simulation) {
    infile = "./input/bathe_FSI.neu";
  }
  else if(5 == simulation){
    infile = "./input/bathe_shell.neu";
  }
  else if(6 == simulation){
    infile = "./input/bathe_cylinder.neu";
  }
  else if(7 == simulation){
    infile = "./input/comsolbenchmark.neu";
  }


  // ******* Set physics parameters *******
  double Lref, Uref, rhof, muf, rhos, ni, E;

  Lref = 1.;
  Uref = 1.;

  if(simulation<3){ //turek 2D
    rhof = 1000.;
    muf = 1.;
    rhos = 1000;
    ni = 0.5;
    E = 1400000;
  }
  else if(simulation==3){ //bathe 2D
    rhof = 1000;
    muf = 0.04;
    rhos = 800;
    ni = 0.5;
    E = 140000000;
  }
  else if(simulation<7){ //bathe 3D
    rhof = 100.;
    muf = 1.;
    rhos = 800;
    ni = 0.5;
    E = 1800000;
  }
  else if(simulation==7){ //comsol
    rhof = 1000.;
    muf = 0.001;
    rhos = 7850;
    ni = 0.5;
    E = 200000;
  }

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

  // ******* Init multilevel mesh from mesh.neu file *******
  unsigned short numberOfUniformRefinedMeshes, numberOfAMRLevels;

  if(simulation < 3)
    numberOfUniformRefinedMeshes=3;
  else if(simulation == 3 || simulation == 7)
    numberOfUniformRefinedMeshes=4;
  else if(simulation < 7)
    numberOfUniformRefinedMeshes=2;

  numberOfAMRLevels = 0;

  MultiLevelMesh ml_msh(numberOfUniformRefinedMeshes, numberOfUniformRefinedMeshes + numberOfAMRLevels,
			infile.c_str(),"fifth",Lref,SetRefinementFlag);

  // ******* Init multilevel solution ******
  MultiLevelSolution ml_sol(&ml_msh);

  // ******* Add solution variables to multilevel solution and pair them *******
  ml_sol.AddSolution("DX",LAGRANGE,SECOND,1);
  ml_sol.AddSolution("DY",LAGRANGE,SECOND,1);
  if (!dimension2D) ml_sol.AddSolution("DZ",LAGRANGE,SECOND,1);

  ml_sol.AddSolution("U",LAGRANGE,SECOND,1);
  ml_sol.AddSolution("V",LAGRANGE,SECOND,1);
  if (!dimension2D) ml_sol.AddSolution("W",LAGRANGE,SECOND,1);

  // Pair each velocity varible with the corresponding displacement variable
  ml_sol.PairSolution("U","DX"); // Add this line
  ml_sol.PairSolution("V","DY"); // Add this line
  if (!dimension2D) ml_sol.PairSolution("W","DZ"); // Add this line

  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  ml_sol.AddSolution("P",DISCONTINUOUS_POLYNOMIAL,FIRST,1);
  ml_sol.AssociatePropertyToSolution("P","Pressure",false); // Add this line

  // ******* Initialize solution *******
  ml_sol.Initialize("All");

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
  if (!dimension2D) ml_sol.GenerateBdc("DZ","Steady");
  ml_sol.GenerateBdc("U","Steady");
  ml_sol.GenerateBdc("V","Steady");
  if (!dimension2D) ml_sol.GenerateBdc("W","Steady");
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
  if (!dimension2D) system.AddSolutionToSystemPDE("DZ");
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  if (!dimension2D) system.AddSolutionToSystemPDE("W");
  system.AddSolutionToSystemPDE("P");

  // ******* System Fluid-Structure-Interaction Assembly *******
  system.SetAssembleFunction(FSISteadyStateAssembly);

  // ******* set MG-Solver *******
  system.SetMgType(F_CYCLE);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-10);
  system.SetNonLinearConvergenceTolerance(1.e-10);
  if( simulation == 7 )
    system.SetNonLinearConvergenceTolerance(1.e-5);
  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);
  if( simulation < 3 || simulation == 7 ) {
    system.SetMaxNumberOfLinearIterations(2);
    system.SetMaxNumberOfNonLinearIterations(10);
  }
  else {
    system.SetMaxNumberOfLinearIterations(8);
    system.SetMaxNumberOfNonLinearIterations(15);
  }

  // ******* Set Preconditioner *******
  if(Gmres) 		system.SetLinearEquationSolverType(FEMuS_DEFAULT);
  else if(Asm) 		system.SetLinearEquationSolverType(FEMuS_ASM);
  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids(GMRES);
  if( simulation < 3 || simulation > 5 )
    system.SetPreconditionerFineGrids(ILU_PRECOND);
  else
    system.SetPreconditionerFineGrids(MLU_PRECOND);
  system.SetTolerances(1.e-12,1.e-20,1.e+50,20);

  // ******* Add variables to be solved *******
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");

  // ******* Set the last (1) variables in system (i.e. P) to be a schur variable *******
  system.SetNumberOfSchurVariables(1);

  // ******* Set block size for the ASM smoothers *******
  if(simulation < 3){
    system.SetElementBlockNumber(2);
  }
  else if(simulation < 7 ){
    system.SetElementBlockNumberFluid(2);
    system.SetElementBlockSolidAll();
  }
  else if(simulation == 7 ){
    system.SetElementBlockNumber(3);
  }

  // ******* For Gmres Preconditioner only *******
  system.SetDirichletBCsHandling(ELIMINATION);

  // ******* Solve *******
  std::cout << std::endl;
  std::cout << " *********** Fluid-Structure-Interaction ************  " << std::endl;
  system.SetOuterSolver(PREONLY);
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


  ml_sol.GetWriter()->SetDebugOutput( true );
  //ml_sol.GetWriter()->ParallelWrite(DEFAULT_OUTPUTDIR,"biquadratic",print_vars);
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

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionTurek_2D_FSI_and_solid(const vector< double > &xx,const char name[], double &value, const int facename, const double time) {
  double x=xx[0];
  double y=xx[1];
  double z=xx[2];
  bool test=1; //dirichlet
  value=0.;
  if(!strcmp(name,"U")) {
    if(1==facename){   //inflow
      test=1;
      double um = 0.2;
      value=1.5*um*4.0/0.1681*y*(0.41-y);
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

bool SetBoundaryConditionBathe_2D_FSI(const vector< double > &xx,const char name[], double &value, const int facename, const double time) {
  double x=xx[0];
  double y=xx[1];
  double z=xx[2];

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



bool SetBoundaryConditionBathe_3D_FSI_and_fluid(const vector< double > &xx,const char name[], double &value, const int facename, const double time) {
  double x=xx[0];
  double y=xx[1];
  double z=xx[2];


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

bool SetBoundaryConditionBathe_3D_solid(const vector< double > &xx,const char name[], double &value, const int facename, const double time) {
  double x=xx[0];
  double y=xx[1];
  double z=xx[2];
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

bool SetBoundaryConditionComsol_2D_FSI(const vector< double > &xx,const char name[], double &value, const int FaceName, const double time) {
  double x=xx[0];
  double y=xx[1];
  double z=xx[2];
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
  return test;
}

