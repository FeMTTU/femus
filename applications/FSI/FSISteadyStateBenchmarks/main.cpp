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
#include "../include/IncompressibleFSIAssembly.hpp"

double scale=1000.;

using namespace std;
using namespace femus;

void AssembleMatrixResFSI(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assemble_matrix);

bool SetBoundaryConditionTurek(const double &x, const double &y, const double &z,const char name[], 
		double &value, const int FaceName, const double = 0.);
bool SetBoundaryConditionDrum(const double &x, const double &y, const double &z,const char name[], 
		double &value, const int FaceName, const double = 0.);
bool SetBoundaryConditionBatheCylinder(const double &x, const double &y, const double &z,const char name[], 
				       double &value, const int facename, const double time);

bool SetBoundaryConditionBatheShell(const double &x, const double &y, const double &z,const char name[], 
				       double &value, const int facename, const double time);
bool SetBoundaryConditionComsol(const double &x, const double &y, const double &z,const char name[], 
		double &value, const int FaceName, const double = 0.);

 
bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber,const int &level);

//------------------------------------------------------------------------------------------------------------------

int main(int argc,char **args) {
    
  /// Init Petsc-MPI communicator
  FemusInit mpinit(argc,args,MPI_COMM_WORLD);

  Files files; 
  files.CheckIODirectories();
  files.RedirectCout();
  
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
  
   
  
	
  unsigned short nm,nr;
  std::cout<<"#MULTIGRID levels? (>=1) \n";
  //std::cin>>nm;
  if(simulation < 3)
    nm=3;
  else if(simulation == 3 || simulation == 7)
    nm=4;
  else if(simulation < 7)
    nm=2;


  std::cout<<"#MAX_REFINEMENT levels? (>=0) \n";
  //std::cin>>nr;
  nr=0;
  int tmp=nm;
  nm+=nr;
  nr=tmp;

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
  
   double Lref, Uref, rhof, muf, rhos, ni, E; 
  
    Lref = 1.;	
    Uref = 1.;
 
  if(simulation<3){ //turek FSI
    rhof = 1000.;
    muf = 1.;
    rhos = 1000;
    ni = 0.5;
    E = 1400000; 
  }
  else if(simulation==3){ //bathe membrane
    rhof = 1000;
    muf = 0.04;
    rhos = 800;
    ni = 0.5;
    E = 140000000; 
  }
  else if(simulation<7){ //bathe cylinder
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
  if(simulation<3){
    solid = Solid(par,E,ni,rhos,"Mooney-Rivlin-MassPenalty"); 
    //solid = Solid(par,E,ni,rhos,"Neo-Hookean-BW-MassPenalty");
    //solid = Solid(par,E,ni,rhos,"Neo-Hookean-BW-Penalty");
  }
  else if(simulation==3 || simulation == 7){
    solid = Solid(par,E,ni,rhos,"Mooney-Rivlin"); 
    //solid = Solid(par,E,ni,rhos,"Neo-Hookean-BW-MassPenalty");
    //solid = Solid(par,E,ni,rhos,"Neo-Hookean-BW-Penalty");
  }
  else if(simulation < 7){	
    solid = Solid(par,E,ni,rhos,"Mooney-Rivlin"); 
    //solid = Solid(par,E,ni,rhos,"Neo-Hookean");
    //Solid solid(par,E,ni,rhos,"Neo-Hookean-BW");
    //solid = Solid(par,E,ni,rhos,"Neo-Hookean-BW-Penalty");
    //solid = Solid(par,E,ni,rhos,"Neo-Hookean-AB-Penalty"); //Allan Bower
   
  }
  
  
  cout << "Solid properties: " << endl;
  cout << solid << endl;
  
  // Generate Fluid Object
  Fluid fluid(par,muf,rhof,"Newtonian");
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;

  MultiLevelMesh ml_msh(nm,nr,infile.c_str(),"fifth",Lref,SetRefinementFlag);
  
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
  else if( 3==simulation)
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionDrum);
  else if (4==simulation || 6==simulation)
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionBatheCylinder);
  else if (5==simulation)
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionBatheShell);
  else if (7 == simulation)
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionComsol);

  ml_sol.GenerateBdc("DX","Steady");
  ml_sol.GenerateBdc("DY","Steady");
  if (!dimension2D) ml_sol.GenerateBdc("DZ","Steady");
  ml_sol.GenerateBdc("U","Steady");
  ml_sol.GenerateBdc("V","Steady");
  if (!dimension2D) ml_sol.GenerateBdc("W","Steady");
  ml_sol.GenerateBdc("P","Steady");
  
  MultiLevelProblem ml_prob(&ml_sol);
  // Add fluid object
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;
  // Add Solid Object
  ml_prob.parameters.set<Solid>("Solid") = solid;
  // mark Solid nodes
  ml_msh.MarkStructureNode();
 
  //create systems
  // add the system FSI to the MultiLevel problem
  MonolithicFSINonLinearImplicitSystem & system = ml_prob.add_system<MonolithicFSINonLinearImplicitSystem> ("Fluid-Structure-Interaction");
  system.AddSolutionToSystemPDE("DX");
  system.AddSolutionToSystemPDE("DY");
  if (!dimension2D) system.AddSolutionToSystemPDE("DZ");
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  if (!dimension2D) system.AddSolutionToSystemPDE("W");
  system.AddSolutionToSystemPDE("P");
  
//   if(dimension2D){
//     bool sparsity_pattern_matrix[5][5]={{1, 0, 1, 0, 0},
// 					{0, 1, 0, 1, 0},
// 					{1, 1, 1, 1, 1},
// 					{1, 1, 1, 1, 1},
// 					{1, 1, 0, 0, 1}};
//     vector < bool > sparsity_pattern (sparsity_pattern_matrix[0],sparsity_pattern_matrix[0]+25*sizeof(bool));
//     system.SetSparsityPattern(sparsity_pattern);  
//   }
//   else{				   
//     bool sparsity_pattern_matrix[7][7]={{1, 0, 0, 1, 0, 0, 0},
// 					{0, 1, 0, 0, 1, 0, 0},
// 					{0, 0, 1, 0, 0, 1, 0},
// 					{1, 1, 1, 1, 1, 1, 1},
// 					{1, 1, 1, 1, 1, 1, 1},
// 					{1, 1, 1, 1, 1, 1, 1},
// 					{1, 1, 1, 0, 0, 0, 1}};
//     vector < bool > sparsity_pattern (sparsity_pattern_matrix[0],sparsity_pattern_matrix[0]+49*sizeof(bool));
//     system.SetSparsityPattern(sparsity_pattern);  
//   }
   
  // System Fluid-Structure-Interaction
  system.SetAssembleFunction(IncompressibleFSIAssemblyAD_DD);  
  //system.SetAssembleFunction(AssembleMatrixResFSI);  
  
  system.SetMgType(F_CYCLE);
  system.SetAbsoluteConvergenceTolerance(1.e-10);
  system.SetNonLinearConvergenceTolerance(1.e-10);
  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);
  if(simulation < 3 || simulation == 7){
    system.SetMaxNumberOfLinearIterations(2);
    system.SetMaxNumberOfNonLinearIterations(10);
  }
  else if(simulation < 7){	
    system.SetMaxNumberOfLinearIterations(8);
    system.SetMaxNumberOfNonLinearIterations(15); 
  }
  
  if(simulation == 7)  system.SetNonLinearConvergenceTolerance(1.e-5);

   
  //Set Smoother Options
  if(Gmres) 		system.SetMgSmoother(GMRES_SMOOTHER);
  else if(Asm) 		system.SetMgSmoother(ASM_SMOOTHER);
  else if(Vanka)	system.SetMgSmoother(VANKA_SMOOTHER);
  
  // init all the systems
  system.init();
  
  system.SetSolverFineGrids(GMRES);
  if(3 >= simulation || 6 == simulation || 7 == simulation )
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
  if(simulation < 3){
    system.SetElementBlockNumber(2);
  }
  else if(simulation ==3 || !dimension2D){
    //system.SetElementBlockNumber("All");
    //system.SetElementBlockNumber(2);
    system.SetElementBlockNumberFluid(2);
    //system.SetElementBlockNumberSolid(2);
    //system.SetElementBlockFluidAll();
    system.SetElementBlockSolidAll();
  }
  else if(simulation == 7 ){
     system.SetElementBlockNumber(3);   
  }
  
  //for Gmres smoother
  //system.SetDirichletBCsHandling(PENALTY); 
  system.SetDirichletBCsHandling(ELIMINATION);   
   
  ml_sol.SetWriter(GMV);

  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  mov_vars.push_back("DZ");
  ml_sol.GetWriter()->SetMovingMesh(mov_vars);
  
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
      
  ml_sol.GetWriter()->write_system_solutions(files.GetOutputPath(),"biquadratic",print_vars);

  // Destroy all the new systems
  ml_prob.clear();
   
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

bool SetBoundaryConditionDrum(const double &x, const double &y, const double &z,const char name[], double &value, const int facename, const double time) {
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



//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionComsol(const double &x, const double &y, const double &z,const char name[], double &value, const int FaceName, const double time) {
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

