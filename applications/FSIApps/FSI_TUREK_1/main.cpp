#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"
#include "FemTTUInit.hpp"
#include "SparseMatrix.hpp"
#include "VTKWriter.hpp"
#include "FElemTypeEnum.hpp"
#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "../include/IncompressibleFSIAssembly.hpp"

using std::cout;
using std::endl;
using namespace femus;

void AssembleMatrixResFSI(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix);

bool SetBoundaryConditionTurek(const double &x, const double &y, const double &z,const char name[], 
		double &value, const int FaceName, const double = 0.);
bool SetBoundaryConditionBatheCylinder(const double &x, const double &y, const double &z,const char name[], 
				       double &value, const int facename, const double time);

bool SetBoundaryConditionBatheShell(const double &x, const double &y, const double &z,const char name[], 
				       double &value, const int facename, const double time);
 
bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber,const int &level);

//------------------------------------------------------------------------------------------------------------------

int main(int argc,char **args) {
    
  unsigned simulation;
  bool dimension2D;
  
  if(argc >= 2) {
    if( !strcmp("turek",args[1])) {
      simulation=1; 
      dimension2D=1;  
    }
    else if( !strcmp("beam",args[1])) {
      simulation=2; 
      dimension2D=1;  
    }
    else if( !strcmp("bathe_FSI",args[1])){
      simulation=3; 
      dimension2D=0;  
    }
    else if( !strcmp("bathe_shell",args[1])) {
      simulation=4; 
      dimension2D=0;  
    }
    else if( !strcmp("bathe_cylinder",args[1])) {
      simulation=5; 
      dimension2D=0;  
    }
    else{    
      cout << "wrong input arguments!\n";
      cout << "please specify the simulation you want to run, options are\n";
      cout << "turek\nbeam\nbathe_cylinder\nbathe_shell\n";
      exit(0);
    }
  }
  else {   
    cout << "wrong input arguments!\n";
    cout << "please specify the simulation you want to run, options are\n";
    cout << "turek\nbeam\nbathe_FSI\nbathe_shell\nbathe_cylinder\n";
    exit(0);
  }
   
  
  bool Vanka=0, Gmres=0, Asm=0;
  if(argc >= 3) {
    if( !strcmp("vanka",args[2])) 	Vanka=1;
    else if( !strcmp("gmres",args[2])) 	Gmres=1;
    else if( !strcmp("asm",args[2])) 	Asm=1;
    
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
  FemTTUInit mpinit(argc,args,MPI_COMM_WORLD);

  unsigned short nm,nr;
  std::cout<<"#MULTIGRID levels? (>=1) \n";
  //std::cin>>nm;
  if(simulation<3)
    nm=3;
  else if(simulation<6)
    nm=2;

  std::cout<<"#MAX_REFINEMENT levels? (>=0) \n";
  //std::cin>>nr;
  nr=0;
  int tmp=nm;
  nm+=nr;
  nr=tmp;

  char *infile = new char [50];
  
  if(1 == simulation){
    sprintf(infile,"./input/turek.neu");
  }
  else if(2 == simulation){
    sprintf(infile,"./input/beam.neu");
  }  
  else if(3 == simulation) {
    sprintf(infile,"./input/bathe_FSI.neu");
  }
  else if(4==simulation){
    sprintf(infile,"./input/bathe_shell.neu");
  }
  else if(5==simulation){
    sprintf(infile,"./input/bathe_cylinder.neu");
  }
  
   double Lref, Uref, rhof, muf, rhos, ni, E; 
  
  
  if(simulation<3){
    Lref = 1.;	
    Uref = 1.;
    rhof = 1000.;
    muf = 1.;
    rhos = 1000;
    ni = 0.5;
    E = 1400000; 
  }
  else if(simulation<6){
    Lref = 1.;
    Uref = 1.;
    rhof = 100.;
    muf = 1.;
    rhos = 800;
    ni = 0.5;
    E = 1800000;
  }
  
  Parameter par(Lref,Uref);
  
  // Generate Solid Object
  Solid solid;
  if(simulation<3){
    //solid = Solid(par,E,ni,rhos,"Neo-Hookean-MassPenalty");
    solid = Solid(par,E,ni,rhos,"Mooney-Rivlin-MassPenalty"); 
    //solid = Solid(par,E,ni,rhos,"Neo-Hookean-BW-Penalty");
  }
  else if(simulation < 6){	
    //solid = Solid(par,E,ni,rhos,"Neo-Hookean");
    //Solid solid(par,E,ni,rhos,"Neo-Hookean-BW");
    //solid = Solid(par,E,ni,rhos,"Neo-Hookean-BW-Penalty");
    //solid = Solid(par,E,ni,rhos,"Neo-Hookean-AB-Penalty"); //Allan Bower
    solid = Solid(par,E,ni,rhos,"Mooney-Rivlin"); 
  }
  
  //Solid solid(par,E,ni,rhos,"Linear_elastic");
  //Solid solid(par,E,ni,rhos,"Neo-Hookean");
  //Solid solid(par,E,ni,rhos,"Neo-Hookean-BW");
  //Solid solid(par,E,ni,rhos,"Neo-Hookean-BW-Penalty");
  
  cout << "Solid properties: " << endl;
  cout << solid << endl;
  
  // Generate Fluid Object
  Fluid fluid(par,muf,rhof,"Newtonian");
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;

  MultiLevelMesh ml_msh(nm,nr,infile,"fifth",Lref,SetRefinementFlag);
  
  MultiLevelSolution ml_sol(&ml_msh);
  
  //Start System Variables
  ml_sol.AddSolution("DX",LAGRANGE,SECOND,1);
  ml_sol.AddSolution("DY",LAGRANGE,SECOND,1);
  if (!dimension2D) ml_sol.AddSolution("DZ",LAGRANGE,SECOND,1);
  ml_sol.AssociatePropertyToSolution("DX","Displacement"); // Add this line
  ml_sol.AssociatePropertyToSolution("DY","Displacement"); // Add this line 
  if (!dimension2D) ml_sol.AssociatePropertyToSolution("DZ","Displacement"); // Add this line 
  ml_sol.AddSolution("U",LAGRANGE,SECOND,1);
  ml_sol.AddSolution("V",LAGRANGE,SECOND,1);
  if (!dimension2D) ml_sol.AddSolution("W",LAGRANGE,SECOND,1);
  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  ml_sol.AddSolution("P",DISCONTINOUS_POLYNOMIAL,FIRST,1);
  //ml_sol.AddSolution("P",LAGRANGE,FIRST,1);
  ml_sol.AssociatePropertyToSolution("P","Pressure"); // Add this line

  //Initialize (update Init(...) function)
  ml_sol.Initialize("All");

  //Set Boundary (update Dirichlet(...) function)
  if(1==simulation || 2==simulation)
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionTurek);
  else if (3==simulation || 5==simulation)
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionBatheCylinder);
  else if (4==simulation)
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionBatheShell);

  ml_sol.GenerateBdc("DX","Steady");
  ml_sol.GenerateBdc("DY","Steady");
  if (!dimension2D) ml_sol.GenerateBdc("DZ","Steady");
  ml_sol.GenerateBdc("U","Steady");
  ml_sol.GenerateBdc("V","Steady");
  if (!dimension2D) ml_sol.GenerateBdc("W","Steady");
  ml_sol.GenerateBdc("P","Steady");
  
  MultiLevelProblem ml_prob(&ml_msh,&ml_sol);
  // Add fluid object
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;
  // Add Solid Object
  ml_prob.parameters.set<Solid>("Solid") = solid;
  // mark Solid nodes
  ml_msh.MarkStructureNode();
 
  //create systems
  // add the system FSI to the MultiLevel problem
  MonolithicFSINonLinearImplicitSystem & system = ml_prob.add_system<MonolithicFSINonLinearImplicitSystem> ("Fluid-Structure-Interaction");
  system.AddSolutionToSytemPDE("DX");
  system.AddSolutionToSytemPDE("DY");
  if (!dimension2D) system.AddSolutionToSytemPDE("DZ");
  system.AddSolutionToSytemPDE("U");
  system.AddSolutionToSytemPDE("V");
  if (!dimension2D) system.AddSolutionToSytemPDE("W");
  system.AddSolutionToSytemPDE("P");
  
  if(dimension2D){
    bool sparsity_pattern_matrix[5][5]={{1, 0, 1, 0, 0},
					{0, 1, 0, 1, 0},
					{1, 1, 1, 1, 1},
					{1, 1, 1, 1, 1},
					{1, 1, 0, 0, 1}};
    vector < bool > sparsity_pattern (sparsity_pattern_matrix[0],sparsity_pattern_matrix[0]+25*sizeof(bool));
    system.SetSparsityPattern(sparsity_pattern);  
  }
  else{				   
    bool sparsity_pattern_matrix[7][7]={{1, 0, 0, 1, 0, 0, 0},
					{0, 1, 0, 0, 1, 0, 0},
					{0, 0, 1, 0, 0, 1, 0},
					{1, 1, 1, 1, 1, 1, 1},
					{1, 1, 1, 1, 1, 1, 1},
					{1, 1, 1, 1, 1, 1, 1},
					{1, 1, 1, 0, 0, 0, 1}};
    vector < bool > sparsity_pattern (sparsity_pattern_matrix[0],sparsity_pattern_matrix[0]+49*sizeof(bool));
    system.SetSparsityPattern(sparsity_pattern);  
  }
   
  // System Fluid-Structure-Interaction
  system.AttachAssembleFunction(IncompressibleFSIAssemblyAD);  
  //system.AttachAssembleFunction(AssembleMatrixResFSI);  
  
  if(simulation < 3){
    system.SetMaxNumberOfLinearIterations(2);
    system.SetMaxNumberOfNonLinearIterations(10);
  }
  else if(simulation < 6){	
    system.SetMaxNumberOfLinearIterations(8);
    system.SetMaxNumberOfNonLinearIterations(15); 
  }
  system.SetMgType(F_CYCLE);
  system.SetAbsoluteConvergenceTolerance(1.e-10);
  system.SetNonLinearConvergenceTolerance(1.e-10);
  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);
   
  //Set Smoother Options
  if(Gmres) 		system.SetMgSmoother(GMRES_SMOOTHER);
  else if(Asm) 		system.SetMgSmoother(ASM_SMOOTHER);
  else if(Vanka)	system.SetMgSmoother(VANKA_SMOOTHER);
  
  // init all the systems
  system.init();
  
  system.SetSolverFineGrids(GMRES);
  if(1==simulation || 2==simulation || 5==simulation)
    system.SetPreconditionerFineGrids(ILU_PRECOND); 
  else
    system.SetPreconditionerFineGrids(MLU_PRECOND); 
 
  system.SetTolerances(1.e-12,1.e-20,1.e+50,20);
 
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");
  
  system.AddVariableToBeSolved("DX");
  system.AddVariableToBeSolved("DY");
  if (!dimension2D)  system.AddVariableToBeSolved("DZ");
  system.AddVariableToBeSolved("U");
  system.AddVariableToBeSolved("V");
  if (!dimension2D)  system.AddVariableToBeSolved("W");
  system.AddVariableToBeSolved("P");
    
  //for Vanka and ASM smoothers
  system.SetNumberOfSchurVariables(1);
  system.SetElementBlockNumber(2);   
  //for Gmres smoother
  //system.SetDirichletBCsHandling(PENALTY); 
  system.SetDirichletBCsHandling(ELIMINATION);   
   
  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  mov_vars.push_back("DZ");
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
  if (!dimension2D) print_vars.push_back("DZ");
  print_vars.push_back("U");
  print_vars.push_back("V");
  if (!dimension2D) print_vars.push_back("W");
  print_vars.push_back("P");
      
  vtkio.write_system_solutions("biquadratic",print_vars);

  // Destroy all the new systems
  ml_prob.clear();
   
  delete [] infile;
  return 0;
}


bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &elemgroupnumber,const int &level) {
  bool refine=0;

  //refinemenet based on elemen group number
  if (elemgroupnumber==5) refine=1;
  if (elemgroupnumber==6) refine=1;
  if (elemgroupnumber==7 && level<5) refine=1;

  return refine;

}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionTurek(const double &x, const double &y, const double &z,const char name[], double &value, const int facename, const double time) {
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
  else if(!strcmp(name,"W")){
    if(1==facename){
      test=1;
      value=0.;
    }  
    else if(2==facename ){  
      test=1;
      value=0.;
    }
    else if(3==facename ){  
      test=1;
      value=0.;
    }
    else if(4==facename ){  
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
  else if(!strcmp(name,"DZ")){
    if(1==facename){         //inflow
      test=1;
      value=0.;
    }  
    else if(2==facename ){   //outflow
     test=1;
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



bool SetBoundaryConditionBatheCylinder(const double &x, const double &y, const double &z,const char name[], double &value, const int facename, const double time) {
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

bool SetBoundaryConditionBatheShell(const double &x, const double &y, const double &z,const char name[], double &value, const int facename, const double time) {
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


void AssembleMatrixResFSI(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix) {
    
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
      for (unsigned ig=0;ig < ml_prob._ml_msh->_type_elem[kelt][order_ind2]->GetGaussPointNumber(); ig++) {

	// *** get Jacobian and test function and test function derivatives in the moving frame***
	(ml_prob._ml_msh->_type_elem[kelt][order_ind2]->*(ml_prob._ml_msh->_type_elem[kelt][order_ind2])->Jacobian_ptr)(vx,ig,Weight,phi,gradphi);
	(ml_prob._ml_msh->_type_elem[kelt][order_ind2]->*(ml_prob._ml_msh->_type_elem[kelt][order_ind2])->Jacobian_ptr)(vx_hat,ig,Weight_hat,phi_hat,gradphi_hat);
	phi1=ml_prob._ml_msh->_type_elem[kelt][order_ind1]->GetPhi(ig);
	if (flag_mat==2) Weight_nojac = ml_prob._ml_msh->_type_elem[kelt][order_ind2]->GetGaussWeight(ig);

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
  // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}
