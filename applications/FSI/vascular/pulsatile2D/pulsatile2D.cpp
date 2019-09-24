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
#include "TransientSystem.hpp"
#include "VTKWriter.hpp"
#include "../include/FSITimeDependentAssemblySupgNonConservativeTwoPressures.hpp"
#include <cmath>
double scale = 1000.;

using namespace std;
using namespace femus;

double SetVariableTimeStep(const double time);

bool SetBoundaryConditionTurek2D(const std::vector < double >& x, const char name[],
                                 double &value, const int facename, const double time);

bool SetBoundaryConditionThrombus2D(const std::vector < double >& x, const char name[],
                                    double &value, const int facename, const double time);

bool SetBoundaryConditionAorticBifurcation(const std::vector < double >& x, const char name[],
    double &value, const int facename, const double time);

bool SetBoundaryConditionThrombus2DPorous(const std::vector < double >& x, const char name[],
    double &value, const int facename, const double time);

bool SetBoundaryConditionVeinValve(const std::vector < double >& x, const char name[],
                                   double &value, const int facename, const double time);

void GetSolutionNorm(MultiLevelSolution& mlSol, const unsigned & group, std::vector <double> &data);

void PrintConvergenceInfo(char *stdOutfile, const unsigned &numberOfUniformRefinedMeshes, const int &nprocs);
//------------------------------------------------------------------------------------------------------------------

int main(int argc, char **args)
{

  
  
  // ******* Init Petsc-MPI communicator *******
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
  clock_t start_time = clock();

  valve = false;
  twoPressure = true;
  
  unsigned simulation = 0;

  if(argc >= 2) {
    if(!strcmp("0", args[1])) {       /** FSI Turek2D no stent */
      simulation = 0;
    }
    else if(!strcmp("1", args[1])) {       /** FSI Turek porous */
      simulation = 1;
    }
    else if(!strcmp("2", args[1])) {       /** FSI Turek stents 60 micron */
      simulation = 2;
    }
    else if(!strcmp("3", args[1])) {       /** FSI Turek 11 stents 60 micron */
      simulation = 3;
    }
    else if(!strcmp("4", args[1])) {       /** FSI AAA thrombus 2D */
      simulation = 4;
    }
    else if(!strcmp("5", args[1])) {       /** FSI Aortic Bifurcation */
      simulation = 5;
    }
    else if(!strcmp("6", args[1])) {       /** FSI AAA thrombus 2D Porous */
      simulation = 6;
    }
    else if(!strcmp("7", args[1])) {       /** FSI Vein Valve */
      simulation = 7;
    }
    else{
      simulation = 5;
    }
  }

  //Files files;
  //files.CheckIODirectories();
  //files.RedirectCout();

  // ******* Extract the problem dimension and simulation identifier based on the inline input *******


  // ******* Extract the mesh.neu file name based on the simulation identifier *******
  std::string infile;
  if(simulation == 0) {
    infile = "./../input/steady&pulsatile/2D/Turek_moreElem.neu";
  }
  else if(simulation == 1) {
    infile = "./../input/steady&pulsatile/2D/Turek_porous_60micron.neu";
  }
  else if(simulation == 2) {
    infile = "./../input/steady&pulsatile/2D/Turek_stents_60micron.neu";
  }
  else if(simulation == 3) {
    infile = "./../input/steady&pulsatile/2D/Turek_11stents_60micron.neu";
  }
  else if(simulation == 4) {
    infile = "./../input/steady&pulsatile/2D/AAA_thrombus_2D.neu";
  }
  else if(simulation == 5) {
    infile = "./../input/steady&pulsatile/2D/aortic_bifurcation.neu";
  }
  else if(simulation == 6) {
    infile = "./../input/steady&pulsatile/2D/AAA_thrombus_2D_porous.neu";
  }
  else if(simulation == 7) {
    //infile = "./../input/valve/valveOld/vein_valve_closed.neu";
    //infile = "./../input/valve/valveOld/vein_valve_thiner.neu";
    infile = "./../input/valve/valveOld/vein_valve_modifiedFluid.neu";
    //infile = "./../input/valve/valveOld/vein_valve_new.neu";
  }

  // ******* Set physics parameters *******
  double Lref, Uref, rhof, muf, rhos, ni, E, E1;

  Lref = 1.;
  Uref = 1.;


  if(simulation == 7) {
    rhof = 1060.;
    muf = 2.2 * 1.0e-3;
    rhos = 960;
    ni = 0.5;
    //E = 3.3 * 1.0e6; //vein young modulus
    E = 4.3874951 * 1.0e12;
    E1 = 15 * 1.0e6; //leaflet young modulus
  }
  else {
    rhof = 1035.;
    muf = 3.5 * 1.0e-3; //wrong=3.38*1.0e-4*rhof, note:3.38*1.0e-6*rhof=3.5*1.0e-3
    rhos = 1120;
    ni = 0.5;
    E = 100 * 1.e6; //simulation=5
    //E = 1. * 1.0e6; //simulation=0,1,2,3
    E1 = 50000;
  }

  Parameter par(Lref, Uref);

  // Generate Solid Object
  Solid solid;
  solid = Solid(par, E, ni, rhos, "Mooney-Rivlin");

  Solid solid1;
  solid1 = Solid(par, E1, ni, rhos, "Mooney-Rivlin");

  cout << "Solid properties: " << endl;
  cout << solid << endl;

  // Generate Fluid Object
  Fluid fluid(par, muf, rhof, "Newtonian");
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;

  // ******* Init multilevel mesh from mesh.neu file *******
  unsigned short numberOfUniformRefinedMeshes, numberOfAMRLevels;

  numberOfUniformRefinedMeshes = 3;
  
  numberOfAMRLevels = 0;

  std::cout << 0 << std::endl;

  MultiLevelMesh ml_msh(numberOfUniformRefinedMeshes + numberOfAMRLevels, numberOfUniformRefinedMeshes,
                        infile.c_str(), "fifth", Lref, NULL);

  //ml_msh.EraseCoarseLevels(numberOfUniformRefinedMeshes - 2);

  ml_msh.PrintInfo();

  // mark Solid nodes

  //ml_msh.MarkStructureNode();

  // ******* Init multilevel solution ******
  MultiLevelSolution ml_sol(&ml_msh);

  // ******* Add solution variables to multilevel solution and pair them *******
  ml_sol.AddSolution("DX", LAGRANGE, SECOND, 2);
  ml_sol.AddSolution("DY", LAGRANGE, SECOND, 2);
 
  ml_sol.AddSolution("U", LAGRANGE, SECOND, 2);
  ml_sol.AddSolution("V", LAGRANGE, SECOND, 2);

  // Pair each velocity variable with the corresponding displacement variable
  ml_sol.PairSolution("U", "DX");    // Add this line
  ml_sol.PairSolution("V", "DY");    // Add this line

  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  ml_sol.AddSolution("PS", DISCONTINUOUS_POLYNOMIAL, FIRST, 2);
  ml_sol.AssociatePropertyToSolution("PS", "Pressure", false);    // Add this line
  
  ml_sol.AddSolution("PF", DISCONTINUOUS_POLYNOMIAL, FIRST, 2);
  ml_sol.AssociatePropertyToSolution("PF", "Pressure", false);    // Add this line

  ml_sol.AddSolution("lmbd", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

  ml_sol.AddSolution("Um", LAGRANGE, SECOND, 0, false);
  ml_sol.AddSolution("Vm", LAGRANGE, SECOND, 0, false);
  
//   ml_sol.AddSolution("AX", LAGRANGE, SECOND, 2);
//   ml_sol.AddSolution("AY", LAGRANGE, SECOND, 2);
    
  // ******* Initialize solution *******
  ml_sol.Initialize("All");

  if(simulation == 0 || simulation == 1 || simulation == 2 || simulation == 3) {
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionTurek2D);
  }
  else if(simulation == 4 || simulation == 6) {
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionThrombus2D);
  }
  else if(simulation == 5) {
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionAorticBifurcation);
  }
  else if(simulation == 7) {
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionVeinValve);
  }

  // ******* Set boundary conditions *******
  ml_sol.GenerateBdc("DX", "Steady");
  ml_sol.GenerateBdc("DY", "Steady");
  ml_sol.GenerateBdc("PF", "Steady");

  if(simulation == 4 || simulation == 5 || simulation == 6) {
    ml_sol.GenerateBdc("U", "Steady");
    ml_sol.GenerateBdc("V", "Time_dependent");
  }
  else if(simulation == 7) {
    ml_sol.GenerateBdc("U", "Steady");
    ml_sol.GenerateBdc("V", "Steady");
  }
  else {
    ml_sol.GenerateBdc("U", "Time_dependent");
    ml_sol.GenerateBdc("V", "Steady");
  }
  if(simulation == 7) {
    ml_sol.GenerateBdc("PS", "Steady");
  }
  else {
    ml_sol.GenerateBdc("PS", "Time_dependent");
  }

//   for(unsigned level = 0; level < numberOfUniformRefinedMeshes; level++ ){
//     SetLambda(ml_sol, level , SECOND, ELASTICITY);
//   }

  // ******* Define the FSI Multilevel Problem *******

  MultiLevelProblem ml_prob(&ml_sol);
  // Add fluid object
  ml_prob.parameters.set<Fluid> ("Fluid") = fluid;
  // Add Solid Object
  ml_prob.parameters.set<Solid> ("Solid") = solid;
  ml_prob.parameters.set<Solid> ("Solid1") = solid1;

  // ******* Add FSI system to the MultiLevel problem *******
  TransientMonolithicFSINonlinearImplicitSystem & system = ml_prob.add_system<TransientMonolithicFSINonlinearImplicitSystem> ("Fluid-Structure-Interaction");
  system.AddSolutionToSystemPDE("DX");
  system.AddSolutionToSystemPDE("DY");

  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");

  system.AddSolutionToSystemPDE("PS");
  
  if(twoPressure) system.AddSolutionToSystemPDE("PF");

  // ******* System Fluid-Structure-Interaction Assembly *******
  system.SetAssembleFunction(FSITimeDependentAssemblySupgNew2);

  // ******* set MG-Solver *******
  system.SetMgType(F_CYCLE);
  
  system.SetNonLinearConvergenceTolerance(1.e-7);
  //system.SetResidualUpdateConvergenceTolerance ( 1.e-15 );

  if(twoPressure){
    system.SetMaxNumberOfNonLinearIterations(20); //20
    //system.SetMaxNumberOfResidualUpdatesForNonlinearIteration ( 4 );

    system.SetMaxNumberOfLinearIterations(1);
    system.SetAbsoluteLinearConvergenceTolerance(1.e-50);

    system.SetNumberPreSmoothingStep(1);
    system.SetNumberPostSmoothingStep(1);
  }
  else{
    system.SetMaxNumberOfNonLinearIterations(5); //10
    //system.SetMaxNumberOfResidualUpdatesForNonlinearIteration ( 4 );
  
    system.SetMaxNumberOfLinearIterations(3); //6
    system.SetAbsoluteLinearConvergenceTolerance(1.e-13);

    system.SetNumberPreSmoothingStep(0);
    system.SetNumberPostSmoothingStep(2);
  }

  // ******* Set Preconditioner *******

  system.SetLinearEquationSolverType(FEMuS_ASM);

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids(RICHARDSON);
  //system.SetSolverFineGrids(GMRES);
  if(twoPressure) system.SetRichardsonScaleFactor(0.4);

  if(twoPressure)
    system.SetPreconditionerFineGrids(MLU_PRECOND);
  else
    system.SetPreconditionerFineGrids(ILU_PRECOND);
  
  if(twoPressure)
    system.SetTolerances(1.e-10, 1.e-9, 1.e+50, 40, 40);
  else
    system.SetTolerances(1.e-12, 1.e-20, 1.e+50, 20, 10);

  // ******* Add variables to be solved *******
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");

  // ******* Set the last (1) variables in system (i.e. P) to be a schur variable *******
  if(twoPressure)
    system.SetNumberOfSchurVariables(2);
  else 
    system.SetNumberOfSchurVariables(1);

  // ******* Set block size for the ASM smoothers *******
  if(twoPressure)
    system.SetElementBlockNumber(4);
  else
    system.SetElementBlockNumber(2);

  // ******* Print solution *******
  ml_sol.SetWriter(VTK);


  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  //mov_vars.push_back("DZ");
  ml_sol.GetWriter()->SetMovingMesh(mov_vars);

  std::vector<std::string> print_vars;
  print_vars.push_back("All");

  ml_sol.GetWriter()->SetDebugOutput(true);
  ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);


  // ******* Solve *******
  std::cout << std::endl;
  std::cout << " *********** Fluid-Structure-Interaction ************  " << std::endl;

  // time loop parameter
  system.AttachGetTimeIntervalFunction ( SetVariableTimeStep );
  const unsigned int n_timesteps = 128;

  std::vector < std::vector <double> > data(n_timesteps);

  system.ResetComputationalTime();
  
  for(unsigned time_step = 0; time_step < n_timesteps; time_step++) {
    for(unsigned level = 0; level < numberOfUniformRefinedMeshes; level++) {
      SetLambdaNew(ml_sol, level , SECOND, ELASTICITY);
    }
    data[time_step].resize(5);
    if(time_step > 0)
      system.SetMgType(V_CYCLE);
    system.CopySolutionToOldSolution();

    system.MGsolve();
    system.PrintComputationalTime();
    
    StoreMeshVelocity(ml_prob);
    
    //data[time_step][0] = time_step / 16.;
    //data[time_step][0] = time_step / 20.;
    //data[time_step][0] = time_step / 32.;
    //data[time_step][0] = time_step / ( 64 * 1.4 );
//     if(simulation == 0 || simulation == 1 || simulation == 2 || simulation == 3) {
//       GetSolutionNorm(ml_sol, 9, data[time_step]);
//     }
//     else if(simulation == 4) {    //AAA_thrombus, 15=thrombus
//       GetSolutionNorm(ml_sol, 7, data[time_step]);
//     }
//     else if(simulation == 6) {    //AAA_thrombus_porous, 15=thrombus
//       GetSolutionNorm(ml_sol, 7, data[time_step]);
//     }
    ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, time_step + 1);
  }


  int  iproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
//   if(iproc == 0) {
//     std::ofstream outf;
//     if(simulation == 0) {
//       outf.open("DataPrint_Turek.txt");
//     }
//     else if(simulation == 1) {
//       outf.open("DataPrint_TurekPorous.txt");
//     }
//     else if(simulation == 2) {
//       outf.open("DataPrint_TurekStents.txt");
//     }
//     else if(simulation == 3) {
//       outf.open("DataPrint_Turek11Stents.txt");
//     }
//     else if(simulation == 4) {
//       outf.open("DataPrint_AAA_thrombus_2D.txt");
//     }
//     else if(simulation == 6) {
//       outf.open("DataPrint_AAA_thrombus_2D_porous.txt");
//     }


//     if(!outf) {
//       std::cout << "Error in opening file DataPrint.txt";
//       return 1;
//     }
//     for(unsigned k = 0; k < n_timesteps; k++) {
//       outf << data[k][0] << "\t" << data[k][1] << "\t" << data[k][2] << "\t" << data[k][3] << "\t" << data[k][4] << std::endl;
//     }
//     outf.close();
//   }

 

  // ******* Clear all systems *******
  ml_prob.clear();
  std::cout << " TOTAL TIME:\t" << \
          static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl;
	  
//   int  nprocs;	    
//   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
//   if(iproc == 0){
//     char stdOutputName[100];
//     sprintf(stdOutputName, "stdoutput_level%d_nprocs%d_turek2D.txt",numberOfUniformRefinedMeshes, nprocs);
//     PrintConvergenceInfo(stdOutputName, numberOfUniformRefinedMeshes, nprocs);
//   } 
	  
  return 0;
}

//---------------------------------------------------------------------------------------------------------------------

double SetVariableTimeStep(const double time)
{
  //double dt = 1. / ( 64 * 1.4 );
  double dt = 1./32;
//   double shiftedTime = time - floor(time);
//   if( time > 1 && shiftedTime >= 0.125 && shiftedTime < 0.25){
//     dt = 1./256;
//   }
//   std::cout << " Shifted Time = " << shiftedTime << " dt = " << dt<< std::endl; 
  
//   double PI = acos(-1.);
//   double dt = 1./32.*(1.25 + 0.75 * sin(2.*PI*time - 3./4. * PI));
//   std::cout << "time = " << time <<  " dt * 32 = " << dt * 32<<std::endl;
  //double dt = 1. / 20;
  //double dt = 1./16.;
  //double dt = 1.0e6;

//   if( turek_FSI == 2 ){
//     if ( time < 9 ) dt = 0.05;
//     else dt = 0.025;
//   }
//   else if ( turek_FSI == 3 ){
//     //if	    ( time < 5. ) dt = 0.1;
//     //else
//     if ( time < 6. ) dt = 0.01;
//     else             dt = 0.01;
//   }
//   else if ( simulation == 3 ) dt=0.001;
//   else if ( simulation == 4 ) dt=0.1;
//   else if ( simulation == 5 ) dt=0.1;
//   else if ( simulation == 6 ) dt=0.1;
//   else if ( simulation == 7 ) dt=0.001;
//   else{
//     std::cout << "Warning this simulation case has not been considered yet for the time dependent case"<<std::endl;
//     abort();
//   }
  return dt;
}


//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionTurek2D(const std::vector < double >& x, const char name[], double &value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;


//   std::ifstream inf;
//   inf.open ( "./input/womersleyProfile_velMax65cms.txt" );
//   if ( !inf ) {
//     std::cout << "velocity file ./input/womersleyProfile_velMax65cms.txt can not be opened\n";
//     exit ( 0 );
//   }
//
//   std::vector<double> vel ( 64 );
//
//   for ( unsigned i = 0; i < 64; i++ ) {
//     inf >> vel[i];
//   }
//   inf.close();
//
//   double period = 1. / 1.4;
//   double dt = period / 64;
//
//   double time1 = time - floor ( time / period ) * period;
//
//   unsigned j = static_cast < unsigned > ( floor ( time1 / dt ) );

  //git pstd::cout<< name << " " << time <<" "<< j <<" "<<  vel[j] << std::endl;

  double PI = acos(-1.);
  double ramp = (time < 1) ? sin(PI / 2 * time) : 1.;
  //double ramp = ( time < period ) ? sin ( PI / 2 * time / period ) : 1.;

  if(!strcmp(name, "U")) {
    if(1 == facename) {
      //value = 0.05 * (x[1] * 1000 - 6) * ( x[1] * 1000 - 8); //inflow
      value = 0.05 * (x[1] * 1000 - 6) * (x[1] * 1000 - 8) * (1. + 0.75 * sin(2. * PI * time)) * ramp; //inflow
      //value = ( x[1] * 1000 - 6 ) * ( x[1] * 1000 - 8 ) * vel[j] * ramp; //inflow
    }
    else if(2 == facename || 5 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if(!strcmp(name, "V")) {
    //if ( 2 == facename ) {
    //if ( 5 == facename ) {
    if(2 == facename || 5 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if(!strcmp(name, "PS")) {
    test = 0;
    value = 0.;
  }
  else if(!strcmp(name, "PF")) {
    test = 0;
    value = 0.;
  }
  else if(!strcmp(name, "DX")) {
    //if(2 == facename || 4 == facename || 5 == facename || 6 == facename) {
    if(5 == facename || 6 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "DY")) {
    //if(1 == facename || 3 == facename || 5 == facename || 6 == facename) {
    if(5 == facename || 6 == facename) {
      test = 0;
      value = 0;
    }
  }

  return test;

}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionThrombus2D(const std::vector < double >& x, const char name[], double &value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos(-1.);

  double ramp = (time < 1) ? sin(PI / 2 * time) : 1.;

  if(!strcmp(name, "V")) {
    if(1 == facename) {
      double r2 = (x[0] * 100.) * (x[0] * 100.);
      //value = -0.01/.9 * (.9 - r2); //inflow
      value = -0.05 / .81 * (.81 - r2) * (1. + 0.75 * sin(2.*PI * time)) * ramp;        //inflow
    }
    if(2 == facename || 5 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if(!strcmp(name, "U")) {
    if(2 == facename) {
      test = 0;
      value = (10000 + 2500 * sin(2 * PI * time)) * ramp;;
    }
    else if(5 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "PS")) {
    test = 0;
    value = 0.;
  }
  else if(!strcmp(name, "PF")) {
    test = 0;
    value = 0.;
  }
  else if(!strcmp(name, "DX")) {
    if(5 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "DY")) {
    if(5 == facename) {
      test = 0;
      value = 0;
    }
  }

  return test;
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionAorticBifurcation(const std::vector < double >& x, const char name[], double &value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos(-1.);

  double ramp = (time < 1) ? sin(PI / 2 * time) : 1.;

  if(!strcmp(name, "V")) {
    if(1 == facename) {
      double r2 = (x[0] * 100.) * (x[0] * 100.);
      //value = -0.01/.9 * (.9 - r2); //inflow
      value = -0.1 / .81 * (.81 - r2) * (1. + 0.25 * sin(2.*PI * time)) * ramp;        //inflow
    }
    if(2 == facename || 3 == facename || 7 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if(!strcmp(name, "U")) {
    if(2 == facename || 3 == facename || 7 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "PS")) {
    test = 0;
    value = 0.;
//     if ( 2 == facename || 3 == facename ) {
//       //value = ( 12500 + 2500 * sin ( 2 * PI * time ) ) * ramp;
//       value = 50 * ramp;
//     }
  }
  else if(!strcmp(name, "DX")) {
    if(7 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "DY")) {
    if(7 == facename) {
      test = 0;
      value = 0;
    }
  }

  return test;
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionVeinValve(const std::vector < double >& x, const char name[], double &value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos(-1.);
  double ramp = (time < 1) ? sin(PI / 2 * time) : 1.;

  if(!strcmp(name, "V")) {
    if(1 == facename || 2 == facename || 6 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "PS")) {
    test = 0;
    value = 0.;
    if(1 == facename) {
      //value = -1;
      //value = ( /*2.5*/ + 2.5 * sin ( 2 * PI * time ) ) * ramp;
      //value = ( 5 + 3 * sin ( 2 * PI * time ) ) * ramp; //+ 4.5
      //value = ( 6 + 3 * sin ( 2 * PI * time ) ) * ramp; //+ 4.5
      //value = ( 12 + 9 * sin ( 2 * PI * time ) ) * ramp; //runna
      //value = ( 24 + 21 * sin ( 2 * PI * time ) ) * ramp; //runna
      value = (0 + 2 * sin(2 * PI * time)) * ramp;      //+ 4.5
    }
    else if(2 == facename) {
      //value = 1;
      //value = ( /*2.5*/ - 2.5 * sin ( 2 * PI * time ) ) * ramp;
      //value = ( 4 - 1 * sin ( 2 * PI * time ) ) * ramp; //- 4.5
      //value = ( 5 - 3 * sin ( 2 * PI * time ) ) * ramp; //non runna
      value = (0 - 2 * sin(2 * PI * time)) * ramp;      //- 4.5
    }
  }
  else if(!strcmp(name, "DX")) {
    if(5 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "DY")) {
    if(5 == facename || 6 == facename) {
      test = 0;
      value = 0;
    }
  }

  return test;

}

//---------------------------------------------------------------------------------------------------------------------

void GetSolutionNorm(MultiLevelSolution& mlSol, const unsigned & group, std::vector <double> &data)
{

  int  iproc, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  NumericVector* p2;
  NumericVector* v2;
  NumericVector* vol;
  NumericVector* vol0;
  p2 = NumericVector::build().release();
  v2 = NumericVector::build().release();
  vol = NumericVector::build().release();
  vol0 = NumericVector::build().release();

  if(nprocs == 1) {
    p2->init(nprocs, 1, false, SERIAL);
    v2->init(nprocs, 1, false, SERIAL);
    vol->init(nprocs, 1, false, SERIAL);
    vol0->init(nprocs, 1, false, SERIAL);
  }
  else {
    p2->init(nprocs, 1, false, PARALLEL);
    v2->init(nprocs, 1, false, PARALLEL);
    vol->init(nprocs, 1, false, PARALLEL);
    vol0->init(nprocs, 1, false, PARALLEL);
  }

  p2->zero();
  v2->zero();
  vol->zero();
  vol0->zero();

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* solution  = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);


  const unsigned dim = msh->GetDimension();


  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

  vector< double > solP;
  vector< vector < double> >  solV(dim);
  vector< vector < double> > x0(dim);
  vector< vector < double> > x(dim);

  solP.reserve(max_size);
  for(unsigned d = 0; d < dim; d++) {
    solV[d].reserve(max_size);
    x0[d].reserve(max_size);
    x[d].reserve(max_size);
  }
  double weight;
  double weight0;

  vector <double> phiV;
  vector <double> gradphiV;
  vector <double> nablaphiV;

  double *phiP;

  phiV.reserve(max_size);
  gradphiV.reserve(max_size * dim);
  nablaphiV.reserve(max_size * (3 * (dim - 1) + !(dim - 1)));

  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol.GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol.GetIndex("V");    // get the position of "V" in the ml_sol object
  if(dim == 3) solVIndex[2] = mlSol.GetIndex("W");       // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol.GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  vector < unsigned > solDIndex(dim);
  solDIndex[0] = mlSol.GetIndex("DX");    // get the position of "U" in the ml_sol object
  solDIndex[1] = mlSol.GetIndex("DY");    // get the position of "V" in the ml_sol object
  if(dim == 3) solDIndex[2] = mlSol.GetIndex("DZ");       // get the position of "V" in the ml_sol object

  unsigned solDType = mlSol.GetSolutionType(solDIndex[0]);

  unsigned solPIndex;
  solPIndex = mlSol.GetIndex("PS");
  unsigned solPType = mlSol.GetSolutionType(solPIndex);

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    if(msh->GetElementGroup(iel) == group) {
      short unsigned ielt = msh->GetElementType(iel);
      unsigned ndofV = msh->GetElementDofNumber(iel, solVType);
      unsigned ndofP = msh->GetElementDofNumber(iel, solPType);
      unsigned ndofD = msh->GetElementDofNumber(iel, solDType);
      // resize

      phiV.resize(ndofV);
      gradphiV.resize(ndofV * dim);
      nablaphiV.resize(ndofV * (3 * (dim - 1) + !(dim - 1)));

      solP.resize(ndofP);
      for(int d = 0; d < dim; d++) {
        solV[d].resize(ndofV);
        x0[d].resize(ndofD);
        x[d].resize(ndofD);
      }
      // get local to global mappings
      for(unsigned i = 0; i < ndofD; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solDType);
        for(unsigned d = 0; d < dim; d++) {
          x0[d][i] = (*msh->_topology->_Sol[d])(idof);

          x[d][i] = (*msh->_topology->_Sol[d])(idof) +
                    (*solution->_Sol[solDIndex[d]])(idof);
        }
      }

      for(unsigned i = 0; i < ndofV; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof
        for(unsigned  d = 0; d < dim; d++) {
          solV[d][i] = (*solution->_Sol[solVIndex[d]])(idof);      // global extraction and local storage for the solution
        }
      }



      for(unsigned i = 0; i < ndofP; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solPType);
        solP[i] = (*solution->_Sol[solPIndex])(idof);
      }


      for(unsigned ig = 0; ig < mlSol._mlMesh->_finiteElement[ielt][solVType]->GetGaussPointNumber(); ig++) {
        // *** get Jacobian and test function and test function derivatives ***
        msh->_finiteElement[ielt][solVType]->Jacobian(x0, ig, weight0, phiV, gradphiV, nablaphiV);
        msh->_finiteElement[ielt][solVType]->Jacobian(x, ig, weight, phiV, gradphiV, nablaphiV);
        phiP = msh->_finiteElement[ielt][solPType]->GetPhi(ig);

        vol0->add(iproc, weight0);
        vol->add(iproc, weight);

        std::vector < double> SolV2(dim, 0.);
        for(unsigned i = 0; i < ndofV; i++) {
          for(unsigned d = 0; d < dim; d++) {
            SolV2[d] += solV[d][i] * phiV[i];
          }
        }

        double V2 = 0.;
        for(unsigned d = 0; d < dim; d++) {
          V2 += SolV2[d] * SolV2[d];
        }
        v2->add(iproc, V2 * weight);

        double P2 = 0;
        for(unsigned i = 0; i < ndofP; i++) {
          P2 += solP[i] * phiP[i];
        }
        P2 *= P2;
        p2->add(iproc, P2 * weight);
      }
    }
  }

  p2->close();
  v2->close();
  vol0->close();
  vol->close();

  double p2_l2 = p2->l1_norm();
  double v2_l2 = v2->l1_norm();
  double VOL0 = vol0->l1_norm();
  double VOL = vol->l1_norm();

  std::cout.precision(14);
  std::scientific;
  std::cout << " vol0 = " << VOL0 << std::endl;
  std::cout << " vol = " << VOL << std::endl;
  std::cout << " (vol-vol0)/vol0 = " << (VOL - VOL0) / VOL0 << std::endl;
  std::cout << " p_l2 norm / vol = " << sqrt(p2_l2 / VOL)  << std::endl;
  std::cout << " v_l2 norm / vol = " << sqrt(v2_l2 / VOL)  << std::endl;

  data[1] = (VOL - VOL0) / VOL0;
  data[2] = VOL;
  data[3] = sqrt(p2_l2 / VOL);
  data[4] = sqrt(v2_l2 / VOL);

  delete p2;
  delete v2;
  delete vol;

}

//---------------------------------------------------------------------------------------------------------------------

void PrintConvergenceInfo(char *stdOutfile, const unsigned &level, const int &nprocs){

  std::cout<<"END_COMPUTATION\n"<<std::flush;

  std::ifstream inf;
  inf.open(stdOutfile);
  if (!inf) {
    std::cout<<"Redirected standard output file not found\n";
    std::cout<<"add option -std_output std_out_filename > std_out_filename\n";
    return;
  }

  std::ofstream outf;
  char outFileName[100];
  sprintf(outFileName, "turek2D_convergence_level%d_nprocs%d.txt",level, nprocs);

  outf.open(outFileName, std::ofstream::app);
  outf << std::endl << std::endl;
  outf << "Number_of_refinements="<<level<<std::endl;
  outf << "Simulation_Time,Nonlinear_Iteration,resid_norm0,resid_normN,N,convergence";

  std::string str1;
  inf >> str1;
  double simulationTime=0.;
  while (str1.compare("END_COMPUTATION") != 0) {

    if (str1.compare("Simulation") == 0){
      inf >> str1;
      if (str1.compare("Time:") == 0){
        inf >> simulationTime;
      }
    }
    else if (str1.compare("Nonlinear") == 0) {
      inf >> str1;
      if (str1.compare("iteration") == 0) {
        inf >> str1;
        outf << std::endl << simulationTime<<","<<str1;
      }
    }
    else if (str1.compare("KSP") == 0){
      inf >> str1;
      if (str1.compare("preconditioned") == 0){
        inf >> str1;
        if (str1.compare("resid") == 0){
          inf >> str1;
          if (str1.compare("norm") == 0){
            double norm0 = 1.;
            double normN = 1.;
            unsigned counter = 0;
            inf >> norm0;
            outf <<","<< norm0;
            for (unsigned i = 0; i < 11; i++){
              inf >> str1;
            }
            while(str1.compare("norm") == 0){
              inf >> normN;
              counter++;
              for (unsigned i = 0; i < 11; i++){
                inf >> str1;
              }
            }
            outf <<","<< normN;
            if(counter != 0){
              outf << "," <<counter<< "," << pow(normN/norm0,1./counter);
            }
            else{
              outf << "Invalid solver, set -outer_ksp_solver \"gmres\"";
            }
          }
        }
      }
    }
    inf >> str1;
  }

  outf.close();
  inf.close();

}

