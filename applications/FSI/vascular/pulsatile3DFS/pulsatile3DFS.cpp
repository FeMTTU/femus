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
//#include "../include/FSITimeDependentAssemblySupg_OLD.hpp"
//#include "../../include/FSISteadyStateAssembly.hpp"
#include <cmath>

double scale = 1000.;

using namespace std;
using namespace femus;

double SetVariableTimeStep ( const double time );

bool SetBoundaryConditionOmino ( const std::vector < double > & x, const char name[],
                                 double & value, const int facename, const double time );

bool SetBoundaryConditionAorta ( const std::vector < double > & x, const char name[],
                                 double & value, const int facename, const double time );

bool SetBoundaryCondition ( const std::vector < double > & x, const char name[],
                            double & value, const int facename, const double time );

bool SetBoundaryConditionTurek ( const std::vector < double > & x, const char name[],
                                 double & value, const int facename, const double time );

bool SetBoundaryConditionPorous ( const std::vector < double > & x, const char name[],
                                  double & value, const int facename, const double time );

bool SetBoundaryConditionThrombus ( const std::vector < double > & x, const char name[],
                                    double & value, const int facename, const double time );

bool SetBoundaryConditionOminoPorous ( const std::vector < double > & x, const char name[],
                                       double & value, const int facename, const double time );

bool SetBoundaryConditionTubo ( const std::vector < double > & x, const char name[],
                                double & value, const int facename, const double time );

void GetSolutionNorm ( MultiLevelSolution& mlSol, const unsigned & group, std::vector <double> &data );

void PrintConvergenceInfo(char *stdOutfile, const unsigned &numberOfUniformRefinedMeshes, const int &nprocs);
//------------------------------------------------------------------------------------------------------------------

int main ( int argc, char ** args )
{

  // ******* Init Petsc-MPI communicator *******
  FemusInit mpinit ( argc, args, MPI_COMM_WORLD );
  
  clock_t start_time = clock();
  
  valve = false;
  twoPressure = true;

  unsigned simulation = 0;

  if ( argc >= 2 ) {
    if ( !strcmp ( "0", args[1] ) ) { /** FSI Turek3D no stent */
      simulation = 0;
    }
    else if ( !strcmp ( "1", args[1] ) ) {   /** FSI Omino no stent */
      simulation = 1;
    }
    else if ( !strcmp ( "2", args[1] ) ) {   /** FSI Thoracic Aortic Aneurysm */
      simulation = 2;
    }
    else if ( !strcmp ( "3", args[1] ) ) {   /** FSI Abdominal Aortic Aneurysm */
      simulation = 3;
    }
    else if ( !strcmp ( "4", args[1] ) ) {   /** FSI Turek 3D porous */
      simulation = 4;
    }
    else if ( !strcmp ( "5", args[1] ) ) {   /** FSI tubo 3D */
      simulation = 5;
    }
    else{
      simulation = 0;
    }
  }

  bool dimension2D = false;

  //Files files;
  //files.CheckIODirectories();
  //files.RedirectCout();

  // ******* Extract the problem dimension and simulation identifier based on the inline input *******


  // ******* Extract the mesh.neu file name based on the simulation identifier *******
  std::string infile;
  if ( simulation == 0 ) {
    infile = "./../input/steady&pulsatile/3D/Turek_3D_F2.neu"; //Turek_3D_D.neu non ha gruppo 10
  }
  else if ( simulation == 1 ) {
    infile = "./../input/steady&pulsatile/3D/aneurysm_omino.neu";
  }
  else if ( simulation == 2 ) {
    infile = "./../input/steady&pulsatile/3D/aneurisma_aorta.neu";
  }
  else if ( simulation == 3 ) {
    infile = "./../input/steady&pulsatile/3D/AAA_thrombus.neu";
  }
  else if ( simulation == 4 ) {
    infile = "./../input/steady&pulsatile/3D/Turek_3D_porous.neu";
  }
  else if ( simulation == 5 ) {
    infile = "./../input/steady&pulsatile/3D/tubo3D.neu";
  }

  //std::string infile = "./input/turek_porous_scaled.neu";
  //std::string infile = "./input/turek_porous_omino.neu";
  //std::string infile = "./input/AAA_thrombus.neu";

  // ******* Set physics parameters *******
  double Lref, Uref, rhof, muf, rhos, ni, E, E1;

  Lref = 1.;
  Uref = 1.;

  rhof = 1035.;
  muf = 3.5 * 1.0e-3; //3.38 * 1.0e-6 * rhof;
  rhos = 1120;
  ni = 0.5;
  E = 1. * 1.0e6;
  //E = 12000; //E=6000;
  E1 = 3000;

  // Maximum aneurysm_omino deformation (velocity = 0.1)
//   rhof = 1035.;
//   muf = 3.38 * 1.0e-6 * rhof;
//   rhos = 1120;
//   ni = 0.5;
//   E = 6000;

  // Maximum Turek_3D_D deformation (velocity = 0.2)
//   rhof = 1035.;
//   muf = 3.38 * 1.0e-6 * rhof;
//   rhos = 1120;
//   ni = 0.5;
//   E = 6000;

  //E = 60000 * 1.e1;

  Parameter par ( Lref, Uref );

  // Generate Solid Object
  Solid solid;
  solid = Solid ( par, E, ni, rhos, "Mooney-Rivlin" );

  Solid solid1;
  solid1 = Solid ( par, E1, ni, rhos, "Mooney-Rivlin" );

  cout << "Solid properties: " << endl;
  cout << solid << endl;

  // Generate Fluid Object
  Fluid fluid ( par, muf, rhof, "Newtonian" );
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;

  // ******* Init multilevel mesh from mesh.neu file *******
  unsigned short numberOfUniformRefinedMeshes, numberOfAMRLevels;

  numberOfUniformRefinedMeshes = 3;
  numberOfAMRLevels = 0;

  std::cout << 0 << std::endl;

  MultiLevelMesh ml_msh ( numberOfUniformRefinedMeshes + numberOfAMRLevels, numberOfUniformRefinedMeshes,
                          infile.c_str(), "fifth", Lref, NULL );
  
  unsigned dim = ml_msh.GetDimension();

  //ml_msh.EraseCoarseLevels(numberOfUniformRefinedMeshes - 1);

  ml_msh.PrintInfo();

  // mark Solid nodes

  //ml_msh.MarkStructureNode();

  // ******* Init multilevel solution ******
  MultiLevelSolution ml_sol ( &ml_msh );

  // ******* Add solution variables to multilevel solution and pair them *******
  ml_sol.AddSolution ( "DX", LAGRANGE, SECOND, 2 );
  ml_sol.AddSolution ( "DY", LAGRANGE, SECOND, 2 );
  if ( !dimension2D ) ml_sol.AddSolution ( "DZ", LAGRANGE, SECOND, 2 );

  ml_sol.AddSolution ( "U", LAGRANGE, SECOND, 2 );
  ml_sol.AddSolution ( "V", LAGRANGE, SECOND, 2 );
  if ( !dimension2D ) ml_sol.AddSolution ( "W", LAGRANGE, SECOND, 2 );

  // Pair each velocity variable with the corresponding displacement variable
  ml_sol.PairSolution ( "U", "DX" ); // Add this line
  ml_sol.PairSolution ( "V", "DY" ); // Add this line
  if ( !dimension2D ) ml_sol.PairSolution ( "W", "DZ" ); // Add this line

  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  ml_sol.AddSolution ( "PS", DISCONTINUOUS_POLYNOMIAL, FIRST, 2 );
  ml_sol.AssociatePropertyToSolution ( "PS", "Pressure", false ); // Add this line
  
  ml_sol.AddSolution("PF", DISCONTINUOUS_POLYNOMIAL, FIRST, 2);
  ml_sol.AssociatePropertyToSolution("PF", "Pressure", false);    // Add this line

  ml_sol.AddSolution ( "lmbd", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false );
  
  ml_sol.AddSolution("Um", LAGRANGE, SECOND, 0, false);
  ml_sol.AddSolution("Vm", LAGRANGE, SECOND, 0, false);
  if ( !dimension2D ) ml_sol.AddSolution("Wm", LAGRANGE, SECOND, 0, false);
  
//   ml_sol.AddSolution("AX", LAGRANGE, SECOND, 2);
//   ml_sol.AddSolution("AY", LAGRANGE, SECOND, 2);
//   if(dim == 3) ml_sol.AddSolution("AZ", LAGRANGE, SECOND, 2);

  // ******* Initialize solution *******
  ml_sol.Initialize ( "All" );

  if ( simulation == 0 || simulation == 4 ) {
    ml_sol.AttachSetBoundaryConditionFunction ( SetBoundaryConditionTurek );
  }
  else if ( simulation == 1 ) {
    ml_sol.AttachSetBoundaryConditionFunction ( SetBoundaryConditionOmino );
  }
  else if ( simulation == 2 ) {
    ml_sol.AttachSetBoundaryConditionFunction ( SetBoundaryConditionAorta );
  }
  else if ( simulation == 3 ) {
    ml_sol.AttachSetBoundaryConditionFunction ( SetBoundaryConditionThrombus );
  }
  else if ( simulation == 5 ) {
    ml_sol.AttachSetBoundaryConditionFunction ( SetBoundaryConditionTubo );
  }


  //ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionPorous);
  //ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionOminoPorous);


  // ******* Set boundary conditions *******
  ml_sol.GenerateBdc ( "DX", "Steady" );
  ml_sol.GenerateBdc ( "DY", "Steady" );
  if ( !dimension2D ) ml_sol.GenerateBdc ( "DZ", "Steady" );

  if ( simulation == 0 || simulation == 4 || simulation == 5 ) {
    ml_sol.GenerateBdc ( "U", "Time_dependent" );
    ml_sol.GenerateBdc ( "V", "Steady" );
  }
  else {
    ml_sol.GenerateBdc ( "U", "Steady" );
    ml_sol.GenerateBdc ( "V", "Time_dependent" );
  }

  if ( !dimension2D ) ml_sol.GenerateBdc ( "W", "Steady" );
  ml_sol.GenerateBdc ( "PS", "Steady" );
  ml_sol.GenerateBdc ("PF", "Steady");

//   for(unsigned level = 0; level < numberOfUniformRefinedMeshes; level++ ){
//     SetLambda(ml_sol, level , SECOND, ELASTICITY);
//   }

  // ******* Define the FSI Multilevel Problem *******

  MultiLevelProblem ml_prob ( &ml_sol );
  // Add fluid object
  ml_prob.parameters.set<Fluid> ( "Fluid" ) = fluid;
  // Add Solid Object
  ml_prob.parameters.set<Solid> ( "Solid" ) = solid;
  ml_prob.parameters.set<Solid> ( "Solid1" ) = solid1;

  // ******* Add FSI system to the MultiLevel problem *******
  TransientMonolithicFSINonlinearImplicitSystem & system = ml_prob.add_system<TransientMonolithicFSINonlinearImplicitSystem> ( "Fluid-Structure-Interaction" );
  system.AddSolutionToSystemPDE ( "DX" );
  system.AddSolutionToSystemPDE ( "DY" );
  if ( !dimension2D ) system.AddSolutionToSystemPDE ( "DZ" );
  system.AddSolutionToSystemPDE ( "U" );
  system.AddSolutionToSystemPDE ( "V" );
  if ( !dimension2D ) system.AddSolutionToSystemPDE ( "W" );
  system.AddSolutionToSystemPDE ( "PS" );
  if (twoPressure) system.AddSolutionToSystemPDE("PF");
  
  //BEGIN buid fieldSplitTree (only for FieldSplitPreconditioner)
  std::vector < unsigned > fieldVelPf(dim + 1);
  fieldVelPf[0] = system.GetSolPdeIndex("U");
  fieldVelPf[1] = system.GetSolPdeIndex("V");
  fieldVelPf[2] = system.GetSolPdeIndex("W");
  fieldVelPf[dim] = system.GetSolPdeIndex("PF");

  std::vector < unsigned > solutionTypeVelPf(dim + 1);
  solutionTypeVelPf[0] = ml_sol.GetSolutionType("U");
  solutionTypeVelPf[1] = ml_sol.GetSolutionType("V");
  solutionTypeVelPf[2] = ml_sol.GetSolutionType("W");
  solutionTypeVelPf[dim] = ml_sol.GetSolutionType("PF");

  FieldSplitTree VelPf(PREONLY, ASM_PRECOND, fieldVelPf, solutionTypeVelPf, "VelPf");
  VelPf.SetAsmStandard(false);
  VelPf.SetAsmBlockSizeSolid(10);
  //if(dim == 2){
   // VelPf.SetAsmBlockSizeFluid(4);
  //}
  //else{
    VelPf.SetAsmBlockSizeFluid(3);
  //}
    
  VelPf.SetAsmBlockPreconditionerSolid(ILU_PRECOND);
  VelPf.SetAsmBlockPreconditionerFluid(MLU_PRECOND);
  VelPf.SetAsmNumeberOfSchurVariables(1);

  
  std::vector < unsigned > fieldDispPs(dim + 1);
  fieldDispPs[0] = system.GetSolPdeIndex("DX");
  fieldDispPs[1] = system.GetSolPdeIndex("DY");
  fieldDispPs[2] = system.GetSolPdeIndex("DZ");
  fieldDispPs[dim] = system.GetSolPdeIndex("PS");

  std::vector < unsigned > solutionTypeDispPs(dim + 1);
  solutionTypeDispPs[0] = ml_sol.GetSolutionType("DX");
  solutionTypeDispPs[1] = ml_sol.GetSolutionType("DY");
  solutionTypeDispPs[2] = ml_sol.GetSolutionType("DZ");
  solutionTypeDispPs[dim] = ml_sol.GetSolutionType("PS");

  FieldSplitTree DispPs(PREONLY, ASM_PRECOND, fieldDispPs, solutionTypeDispPs, "DispPs");
  DispPs.SetAsmStandard(false);
  //if(dim == 2){
    //DispPs.SetAsmBlockSize(4);
  //}
  //else{
    DispPs.SetAsmBlockSize(3);
  //}
  DispPs.SetAsmBlockPreconditionerSolid(MLU_PRECOND);
  DispPs.SetAsmBlockPreconditionerFluid(MLU_PRECOND);
  DispPs.SetAsmNumeberOfSchurVariables(1);
  
  std::vector < FieldSplitTree *> FS2;
  FS2.reserve(2);

  FS2.push_back(&DispPs);  //displacement first
  FS2.push_back(&VelPf);  // velocity second

  FieldSplitTree FSI(RICHARDSON, FIELDSPLIT_PRECOND, FS2, "FSI");
  FSI.SetRichardsonScaleFactor(.4);

  //END buid fieldSplitTree
  
  system.SetLinearEquationSolverType(FEMuS_FIELDSPLIT);   // Field-Split preconditioned
  
  system.SetOuterSolver(LGMRES); //system.SetOuterKSPSolver("lgmres");

  // ******* System Fluid-Structure-Interaction Assembly *******
  system.SetAssembleFunction ( FSITimeDependentAssemblySupgNew2 );

  // ******* set MG-Solver *******
  system.SetMgType ( F_CYCLE );  
  
  system.SetNonLinearConvergenceTolerance(1.e-7);
  system.SetMaxNumberOfNonLinearIterations(20);
  system.SetMaxNumberOfLinearIterations(1);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-50);
 
  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);
  
//   old parameters   
//   system.SetNonLinearConvergenceTolerance(1.e-9);
//   system.SetResidualUpdateConvergenceTolerance(1.e-15);
//   system.SetMaxNumberOfNonLinearIterations(4);
//   system.SetMaxNumberOfResidualUpdatesForNonlinearIteration(1);
//   
//   system.SetNumberPreSmoothingStep ( 0 );
//   system.SetNumberPostSmoothingStep ( 2 );

  // ******* Set Preconditioner *******

  //system.SetLinearEquationSolverType ( FEMuS_ASM );

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids ( RICHARDSON );
  system.SetRichardsonScaleFactor(.4);

  //system.SetPreconditionerFineGrids ( ILU_PRECOND );
  
  system.SetFieldSplitTree(&FSI);

  //system.SetTolerances ( 1.e-12, 1.e-20, 1.e+50, 20, 10 );
  
  //if(dim==2){
    //system.SetTolerances(1.e-10, 1.e-8, 1.e+50, 40, 40);
  //}
  //else{
    system.SetTolerances(1.e-10, 1.e-12, 1.e+50, 40, 40);
  //}

  // ******* Add variables to be solved *******
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved ( "All" );

  // ******* Set the last (1) variables in system (i.e. P) to be a schur variable *******
  //system.SetNumberOfSchurVariables ( 1 );

  // ******* Set block size for the ASM smoothers *******
  //system.SetElementBlockNumber ( 2 );

  // ******* For Gmres Preconditioner only *******
  //system.SetDirichletBCsHandling(ELIMINATION);


  // ******* Print solution *******
  ml_sol.SetWriter ( VTK );


  std::vector<std::string> mov_vars;
  mov_vars.push_back ( "DX" );
  mov_vars.push_back ( "DY" );
  mov_vars.push_back ( "DZ" );
  ml_sol.GetWriter()->SetMovingMesh ( mov_vars );

  std::vector<std::string> print_vars;
  print_vars.push_back ( "All" );

  ml_sol.GetWriter()->SetDebugOutput ( true );
  ml_sol.GetWriter()->Write ( DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0 );


  // ******* Solve *******
  std::cout << std::endl;
  std::cout << " *********** Fluid-Structure-Interaction ************  " << std::endl;

  // time loop parameter
  system.AttachGetTimeIntervalFunction ( SetVariableTimeStep );
  const unsigned int n_timesteps = 128;

  std::vector < std::vector <double> > data ( n_timesteps );
  
  system.ResetComputationalTime();

  for ( unsigned time_step = 0; time_step < n_timesteps; time_step++ ) {
    for ( unsigned level = 0; level < numberOfUniformRefinedMeshes; level++ ) {
      SetLambdaNew ( ml_sol, level , SECOND, ELASTICITY );
    }
    data[time_step].resize ( 5 );
    if ( time_step > 0 )
      system.SetMgType ( V_CYCLE );
    system.CopySolutionToOldSolution();
    system.MGsolve();
    system.PrintComputationalTime();
    
    StoreMeshVelocity(ml_prob);
    //data[time_step][0] = time_step/16.;
    //data[time_step][0] = time_step / 32.;
    //data[time_step][0] = time_step / (64*1.4);
//     if ( simulation == 0 || simulation == 4 ) {
//       GetSolutionNorm ( ml_sol, 9, data[time_step] );
//     }
//     else if ( simulation == 2 ) {   //aneurisma_aorta
//       GetSolutionNorm ( ml_sol, 14, data[time_step] );
//     }
//     else if ( simulation == 3 ) {   //AAA_thrombus, 15=thrombus
//       GetSolutionNorm ( ml_sol, 7, data[time_step] );
//     }
//     ml_sol.GetWriter()->Write ( DEFAULT_OUTPUTDIR, "biquadratic", print_vars, time_step + 1 );
  }

  int  iproc;
  MPI_Comm_rank ( MPI_COMM_WORLD, &iproc );
//   if ( iproc == 0 ) {
//     std::ofstream outf;
//     if ( simulation == 0 ) {
//       outf.open ( "DataPrint_Turek_3D.txt" );
//     }
//     else if ( simulation == 2 ) {
//       outf.open ( "DataPrint_aorta.txt" );
//     }
//     else if ( simulation == 3 ) {
//       outf.open ( "DataPrint_AAA_thrombus_3D.txt" );
//     }
//     else if ( simulation == 4 ) {
//       outf.open ( "DataPrint_Turek_3D_Porous.txt" );
//     }
//     else {
//       outf.open ( "DataPrint.txt" );
//     }


//     if ( !outf ) {
//       std::cout << "Error in opening file DataPrint.txt";
//       return 1;
//     }
//     for ( unsigned k = 0; k < n_timesteps; k++ ) {
//       outf << data[k][0] << "\t" << data[k][1] << "\t" << data[k][2] << "\t" << data[k][3] << "\t" << data[k][4] << std::endl;
//     }
//     outf.close();
//   }

  
  
  // ******* Clear all systems *******
  ml_prob.clear();
  std::cout << " TOTAL TIME:\t" << \
          static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl;
	  
  int  nprocs;	    
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if(iproc == 0){
    char stdOutputName[100];
    sprintf(stdOutputName, "stdoutput_level%d_nprocs%d_turek3DFS.txt",numberOfUniformRefinedMeshes, nprocs);
    PrintConvergenceInfo(stdOutputName, numberOfUniformRefinedMeshes, nprocs);
  } 	  	  
	  
  return 0;
}

//---------------------------------------------------------------------------------------------------------------------

double SetVariableTimeStep ( const double time )
{
  //double dt = 1./16.;
  double dt = 1. / 32.;
  //double dt = 1./(64*1.4);
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

bool SetBoundaryCondition ( const std::vector < double > & x, const char name[], double & value, const int facename, const double time )
{
  bool test = 1; //dirichlet
  value = 0.;

  if ( !strcmp ( name, "U" ) || !strcmp ( name, "V" ) || !strcmp ( name, "W" ) ) {

    if ( 1 == facename ) {
      test = 0;
      value = 10.;
    }
    else if ( 2 == facename ) {
      test = 0;
      value = 10.;

    }
    else if ( 3 == facename ) {
      test = 0;
      value = 20.;
    }
  }

  else if ( !strcmp ( name, "PS" ) ) {
    test = 0;
    value = 0.;
  }

  else if ( !strcmp ( name, "DX" ) || !strcmp ( name, "DY" ) || !strcmp ( name, "DZ" ) ) {
    if ( 7 == facename ) {
      test = 0;
      value = 1.;
    }
  }

  return test;
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionTurek ( const std::vector < double > & x, const char name[], double & value, const int facename, const double time )
{
  bool test = 1; //dirichlet
  value = 0.;

//   std::ifstream inf;
//   inf.open ( "./../input/steady&pulsatile/3D/womersleyProfile_velMax40cms.txt" );
//   if ( !inf ) {
//     std::cout << "velocity file ./../input/steady&pulsatile/3D/womersleyProfile_velMax40cms.txt can not be opened\n";
//     exit ( 0 );
//   }
//   std::ifstream inf2;
//   inf2.open ( "./../input/steady&pulsatile/3D/OutflowResistence64_R0p001_f84.txt" );
//   if ( !inf2 ) {
//     std::cout << "pressure file ./../input/steady&pulsatile/3D/OutflowResistence64_R0p001_f84.txt can not be opened\n";
//     exit ( 0 );
//   }
// 
//   std::vector<double> vel ( 64 );
//   std::vector<double> pressure ( 64 );
// 
//   for ( unsigned i = 0; i < 64; i++ ) {
//     inf >> vel[i];
//     inf2 >> pressure[i];
//   }
//   inf.close();
//   inf2.close();
// 
//   double period = 1. / 1.4;
//   double dt = period / 64;
// 
//   double time1 = time - floor ( time / period ) * period;
// 
//   unsigned j = static_cast < unsigned > ( floor ( time1 / dt ) );

  double PI = acos ( -1. );
  double ramp = (time < 1) ? sin(PI / 2 * time) : 1.;
  //double ramp = ( time < period ) ? sin ( PI / 2 * time / period ) : 1.;

  if ( !strcmp ( name, "U" ) ) {
    if ( 1 == facename ) {
      double r2 = ( ( x[1] * 1000. ) - 7. ) * ( ( x[1] * 1000. ) - 7. ) + ( x[2] * 1000. ) * ( x[2] * 1000. );
      //value = -0.3 * (1. - r2); //inflow
      value = -0.3 * (1. - r2) * (1. + 0.75 * sin(2.*PI * time)) * ramp; //inflow
      //value = - ( 1. - r2 ) * vel[j] * ramp; //inflow
      //std::cout << value << " " << time << " " << ramp << std::endl;
      //value=25;
    }
    else if ( 2 == facename || 5 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "V" ) || !strcmp ( name, "W" ) ) {
    if ( 2 == facename || 5 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "PS" ) ) {
    test = 0;
    value = 0.;
//     if ( 2 == facename ) {
//       value = pressure[j] * ramp;
//       //value = 5000 * ramp;
//       //value = (5000 + 1000 * sin(2 * PI * time)) * ramp;
//     }
  }
  else if ( !strcmp ( name, "PF" ) ) {
    test = 0;
    value = 0.;
  }
  else if ( !strcmp ( name, "DX" ) ) {
    if ( 5 == facename || 6 == facename ) {
      test = 0;
      value = 0;
    }
  }
  else if ( !strcmp ( name, "DY" ) ) {
    if ( 5 == facename || 6 == facename ) {
      test = 0;
      value = 0;
    }
  }
  else if ( !strcmp ( name, "DZ" ) ) {
    if ( 5 == facename || 6 == facename ) {
      test = 0;
      value = 0;
    }
  }

  return test;
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionPorous ( const std::vector < double > & x, const char name[], double & value, const int facename, const double time )
{
  bool test = 1; //dirichlet
  value = 0.;

  if ( !strcmp ( name, "U" )  || !strcmp ( name, "W" ) ) {
    if ( 2 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "V" ) ) {
    if ( 1 == facename ) {
      double r2 = ( x[0] * 1000. ) * ( x[0] * 1000. ) + ( x[2] * 1000. ) * ( x[2] * 1000. );
      value = 0.25 * ( 1. - r2 ); //inflow
    }
    else if ( 2 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "PS" ) ) {
    test = 0;
    value = 0.;
  }
  else if ( !strcmp ( name, "DX" ) ) {
    if ( 1 == facename || 2 == facename || 5 == facename ) {
      test = 0;
      value = 0;
    }
  }
  else if ( !strcmp ( name, "DY" ) ) {
    if ( 5 == facename ) {
      test = 0;
      value = 0;
    }
  }
  else if ( !strcmp ( name, "DZ" ) ) {
    if ( 1 == facename || 2 == facename || 5 == facename ) {
      test = 0;
      value = 0;
    }
  }

  return test;
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionOminoPorous ( const std::vector < double > & x, const char name[], double & value, const int facename, const double time )
{
  bool test = 1; //dirichlet
  value = 0.;

  if ( !strcmp ( name, "U" )  || !strcmp ( name, "V" ) ) {
    if ( 1 == facename || 2 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "W" ) ) {
    if ( 6 == facename ) {
      double r2 = ( x[0] / .000375 ) * ( x[0] / .000375 ) + ( x[1] / .000375 ) * ( x[1] / .000375 );
      value = 1.0 * ( 1. - r2 ); //inflow
    }
    else if ( 1 == facename || 2 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "PS" ) ) {
    test = 0;
    value = 0.;
  }
  else if ( !strcmp ( name, "DX" ) ) {
    if ( 1 == facename || 2 == facename || 6 == facename || 5 == facename ) {
      test = 0;
      value = 0;
    }
  }
  else if ( !strcmp ( name, "DY" ) ) {
    if ( 5 == facename || 6 == facename ) {
      test = 0;
      value = 0;
    }
  }
  else if ( !strcmp ( name, "DZ" ) ) {
    if ( 1 == facename || 2 == facename || 5 == facename ) {
      test = 0;
      value = 0;
    }
  }

  return test;
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionOmino ( const std::vector < double > & x, const char name[], double & value, const int facename, const double time )
{
  bool test = 1; //dirichlet
  value = 0.;
  double PI = acos ( -1. );
  double ramp = ( time < 1 ) ? sin ( PI / 2 * time ) : 1.;

  if ( !strcmp ( name, "V" ) ) {
    if ( 3 == facename ) {
      double r2 = ( ( x[2] * 1000. ) + 0.403 ) * ( ( x[2] * 1000. ) + 0.403 ) + ( ( x[0] * 1000. ) + 0.589 ) * ( ( x[0] * 1000. ) + 0.589 );
      value = 0.1 * ( 1. - r2 ) * ( 1. + 0.75 * sin ( 2.*PI * time ) ) * ramp; //inflow
      //value = 0.1;
    }
    else if ( 1 == facename || 2 == facename || 7 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "U" ) ) {
    if ( 1 == facename || 2 == facename ) {
      test = 0;
      //value = 10 * ramp;
      value = ( 10000 + 2500 * sin ( 2 * PI * time ) ) * ramp;
    }
    else if ( 7 == facename ) {
      test = 0;
      value = 0;
    }
  }
  else if ( !strcmp ( name, "W" ) ) {
    if ( 1 == facename || 2 == facename || 7 == facename ) {
      test = 0;
      value = 0;
    }
  }
  else if ( !strcmp ( name, "PS" ) ) {
    if ( 1 == facename || 2 == facename ) {
      test = 0;
      value = 0;
    }
  }
  else if ( !strcmp ( name, "DX" ) || !strcmp ( name, "DY" ) || !strcmp ( name, "DZ" ) ) {
    if ( 7 == facename ) {
      test = 0;
      value = 0.;
    }
  }

  return test;
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionAorta ( const std::vector < double > & x, const char name[], double & value, const int facename, const double time )
{
  bool test = 1; //dirichlet
  value = 0.;
  double PI = acos ( -1. );
  double ramp = ( time < 1 ) ? sin ( PI / 2 * time ) : 1.;

  if ( !strcmp ( name, "V" ) ) {
    if ( 1 == facename || 2 == facename || 3 == facename || 4 == facename || 11 == facename ) {
      test = 0;
      value = 0;
    }
    else if ( 5 == facename ) {
//       test = 0;
//       value = 0.5 *1.e-01; //Pressure value/rhof
      double r2 = ( ( x[0] + 0.075563 ) / 0.0104 ) * ( ( x[0] + 0.075563 ) / 0.0104 ) + ( x[2] / 0.0104 ) * ( x[2] / 0.0104 );
      value = 0.02 * ( 1. - r2 ) * ( 1. + 0.75 * sin ( 2.*PI * time ) ) * ramp; //inflow
    }
  }

  else if ( !strcmp ( name, "U" ) ) {
    if ( 1 == facename ) { // || 5 == facename) {
      test = 0;
      value = ( 10000 + 2500 * sin ( 2 * PI * time ) ) * ramp;
    }
    else if ( 2 == facename || 3 == facename || 4 == facename ) {   // || 5 == facename) {
      test = 0;
      value = ( 10500 + 2500 * sin ( 2 * PI * time ) ) * ramp;
    }
    else if ( 11 == facename ) {
      test = 0;
      value = 0.;
    }
  }

  else if ( !strcmp ( name, "W" ) ) {
    if ( 1 == facename || 2 == facename || 3 == facename || 4 == facename || 11 == facename ) { // || 5 == facename) {
      test = 0;
      value = 0;
    }
  }


  else if ( !strcmp ( name, "PS" ) ) {
    test = 0;
    value = 0.;
  }

  else if ( !strcmp ( name, "DX" ) || !strcmp ( name, "DY" ) || !strcmp ( name, "DZ" ) ) {
    if ( 2 == facename || 3 == facename || 4 == facename ||
         7 == facename || 9 == facename || 8 == facename || 11 == facename ) {
      test = 0;
      value = 0;
    }
  }

  return test;
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionThrombus ( const std::vector < double > & x, const char name[], double & value, const int facename, const double time )
{
  bool test = 1; //dirichlet
  value = 0.;
  double PI = acos ( -1. );

  double ramp = ( time < 1 ) ? sin ( PI / 2 * time ) : 1.;
  if ( !strcmp ( name, "V" ) ) {
    if ( 1 == facename ) {
      double r2 = ( x[0] * 100. ) * ( x[0] * 100. ) + ( x[2] * 100. ) * ( x[2] * 100. );
      value = -0.01 / .81 * ( .81 - r2 ) * ( 1. + 0.75 * sin ( 2.*PI * time ) ) * ramp; //inflow
    }
    else if ( 2 == facename || 5 == facename ) {
      test = 0;
      value = 0;
    }
  }
  else if ( !strcmp ( name, "U" ) ) {
    if ( 2 == facename ) {
      test = 0;
      value = ( 10000 + 2500 * sin ( 2 * PI * time ) ) * ramp;;
    }
    else if ( 5 == facename ) {
      test = 0;
      value = 0;
    }
  }
  else if ( !strcmp ( name, "W" ) ) {
    if ( 2 == facename || 5 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "PS" ) ) {
    test = 0;
    value = 0.;
  }
  else if ( !strcmp ( name, "PF" ) ) {
    test = 0;
    value = 0.;
  }
  else if ( !strcmp ( name, "DX" ) ) {
    if ( 5 == facename ) {
      test = 0;
      value = 0;
    }
  }
  else if ( !strcmp ( name, "DY" ) ) {
    if ( 5 == facename ) {
      test = 0;
      value = 0;
    }
  }
  else if ( !strcmp ( name, "DZ" ) ) {
    if ( 5 == facename ) {
      test = 0;
      value = 0;
    }
  }

  return test;
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionTubo ( const std::vector < double > & x, const char name[], double & value, const int facename, const double time )
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos ( -1. );

  if ( !strcmp ( name, "U" ) ) {
    double ramp = ( time < 1 ) ? sin ( PI / 2 * time ) : 1.;
    if ( 2 == facename ) {
      double r2 = ((x[1] - 0.0196 ) * ( x[1] - 0.0196 ) + ( x[2] * x[2] ) ) / (0.0035 * 0.0035);
      value = 0.04 * ( 1. - r2 ) * ( 1. + 0.75 * sin ( 2.*PI * time ) ) * ramp; //inflow
      //std::cout << value << " " << time << " " << ramp << std::endl;
      //value=25;
    }
    else if ( 1 == facename ) {
      test = 0;
      //value = 11335 * ramp;
      value = ( 10000 + 2500 * sin ( 2 * PI * time ) ) * ramp;
    }
    else if ( 5 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "V" ) || !strcmp ( name, "W" ) ) {
    if ( 1 == facename || 5 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "PS" ) ) {
    test = 0;
    value = 0.;
  }
  else if ( !strcmp ( name, "DX" ) || !strcmp ( name, "DY" ) || !strcmp ( name, "DZ" ) ) {
    if ( 5 == facename ) {
      test = 0;
      value = 0;
    }
  }

  return test;
}

//---------------------------------------------------------------------------------------------------------------------

void GetSolutionNorm ( MultiLevelSolution& mlSol, const unsigned & group, std::vector <double> &data )
{

  int  iproc, nprocs;
  MPI_Comm_rank ( MPI_COMM_WORLD, &iproc );
  MPI_Comm_size ( MPI_COMM_WORLD, &nprocs );

  NumericVector* p2;
  NumericVector* v2;
  NumericVector* vol;
  NumericVector* vol0;
  p2 = NumericVector::build().release();
  v2 = NumericVector::build().release();
  vol = NumericVector::build().release();
  vol0 = NumericVector::build().release();

  if ( nprocs == 1 ) {
    p2->init ( nprocs, 1, false, SERIAL );
    v2->init ( nprocs, 1, false, SERIAL );
    vol->init ( nprocs, 1, false, SERIAL );
    vol0->init ( nprocs, 1, false, SERIAL );
  }
  else {
    p2->init ( nprocs, 1, false, PARALLEL );
    v2->init ( nprocs, 1, false, PARALLEL );
    vol->init ( nprocs, 1, false, PARALLEL );
    vol0->init ( nprocs, 1, false, PARALLEL );
  }

  p2->zero();
  v2->zero();
  vol->zero();
  vol0->zero();

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* solution  = mlSol.GetSolutionLevel ( level );
  Mesh* msh = mlSol._mlMesh->GetLevel ( level );


  const unsigned dim = msh->GetDimension();


  const unsigned max_size = static_cast< unsigned > ( ceil ( pow ( 3, dim ) ) );

  vector< double > solP;
  vector< vector < double> >  solV ( dim );
  vector< vector < double> > x0 ( dim );
  vector< vector < double> > x ( dim );

  solP.reserve ( max_size );
  for ( unsigned d = 0; d < dim; d++ ) {
    solV[d].reserve ( max_size );
    x0[d].reserve ( max_size );
    x[d].reserve ( max_size );
  }
  double weight;
  double weight0;

  vector <double> phiV;
  vector <double> gradphiV;
  vector <double> nablaphiV;

  double *phiP;

  phiV.reserve ( max_size );
  gradphiV.reserve ( max_size * dim );
  nablaphiV.reserve ( max_size * ( 3 * ( dim - 1 ) + ! ( dim - 1 ) ) );

  vector < unsigned > solVIndex ( dim );
  solVIndex[0] = mlSol.GetIndex ( "U" ); // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol.GetIndex ( "V" ); // get the position of "V" in the ml_sol object
  if ( dim == 3 ) solVIndex[2] = mlSol.GetIndex ( "W" ); // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol.GetSolutionType ( solVIndex[0] ); // get the finite element type for "u"

  vector < unsigned > solDIndex ( dim );
  solDIndex[0] = mlSol.GetIndex ( "DX" ); // get the position of "U" in the ml_sol object
  solDIndex[1] = mlSol.GetIndex ( "DY" ); // get the position of "V" in the ml_sol object
  if ( dim == 3 ) solDIndex[2] = mlSol.GetIndex ( "DZ" ); // get the position of "V" in the ml_sol object

  unsigned solDType = mlSol.GetSolutionType ( solDIndex[0] );

  unsigned solPIndex;
  solPIndex = mlSol.GetIndex ( "PS" );
  unsigned solPType = mlSol.GetSolutionType ( solPIndex );

  for ( int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++ ) {
    if ( msh->GetElementGroup ( iel ) == group ) {
      short unsigned ielt = msh->GetElementType ( iel );
      unsigned ndofV = msh->GetElementDofNumber ( iel, solVType );
      unsigned ndofP = msh->GetElementDofNumber ( iel, solPType );
      unsigned ndofD = msh->GetElementDofNumber ( iel, solDType );
      // resize

      phiV.resize ( ndofV );
      gradphiV.resize ( ndofV * dim );
      nablaphiV.resize ( ndofV * ( 3 * ( dim - 1 ) + ! ( dim - 1 ) ) );

      solP.resize ( ndofP );
      for ( int d = 0; d < dim; d++ ) {
        solV[d].resize ( ndofV );
        x0[d].resize ( ndofD );
        x[d].resize ( ndofD );
      }
      // get local to global mappings
      for ( unsigned i = 0; i < ndofD; i++ ) {
        unsigned idof = msh->GetSolutionDof ( i, iel, solDType );
        for ( unsigned d = 0; d < dim; d++ ) {
          x0[d][i] = ( *msh->_topology->_Sol[d] ) ( idof );

          x[d][i] = ( *msh->_topology->_Sol[d] ) ( idof ) +
                    ( *solution->_Sol[solDIndex[d]] ) ( idof );
        }
      }

      for ( unsigned i = 0; i < ndofV; i++ ) {
        unsigned idof = msh->GetSolutionDof ( i, iel, solVType ); // global to global mapping between solution node and solution dof
        for ( unsigned  d = 0; d < dim; d++ ) {
          solV[d][i] = ( *solution->_Sol[solVIndex[d]] ) ( idof ); // global extraction and local storage for the solution
        }
      }



      for ( unsigned i = 0; i < ndofP; i++ ) {
        unsigned idof = msh->GetSolutionDof ( i, iel, solPType );
        solP[i] = ( *solution->_Sol[solPIndex] ) ( idof );
      }


      for ( unsigned ig = 0; ig < mlSol._mlMesh->_finiteElement[ielt][solVType]->GetGaussPointNumber(); ig++ ) {
        // *** get Jacobian and test function and test function derivatives ***
        msh->_finiteElement[ielt][solVType]->Jacobian ( x0, ig, weight0, phiV, gradphiV, nablaphiV );
        msh->_finiteElement[ielt][solVType]->Jacobian ( x, ig, weight, phiV, gradphiV, nablaphiV );
        phiP = msh->_finiteElement[ielt][solPType]->GetPhi ( ig );

        vol0->add ( iproc, weight0 );
        vol->add ( iproc, weight );

        std::vector < double> SolV2 ( dim, 0. );
        for ( unsigned i = 0; i < ndofV; i++ ) {
          for ( unsigned d = 0; d < dim; d++ ) {
            SolV2[d] += solV[d][i] * phiV[i];
          }
        }

        double V2 = 0.;
        for ( unsigned d = 0; d < dim; d++ ) {
          V2 += SolV2[d] * SolV2[d];
        }
        v2->add ( iproc, V2 * weight );

        double P2 = 0;
        for ( unsigned i = 0; i < ndofP; i++ ) {
          P2 += solP[i] * phiP[i];
        }
        P2 *= P2;
        p2->add ( iproc, P2 * weight );
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

  std::cout.precision ( 14 );
  std::scientific;
  std::cout << " vol0 = " << VOL0 << std::endl;
  std::cout << " vol = " << VOL << std::endl;
  std::cout << " (vol-vol0)/vol0 = " << ( VOL - VOL0 ) / VOL0 << std::endl;
  std::cout << " p_l2 norm / vol = " << sqrt(p2_l2 / VOL)  << std::endl;
  std::cout << " v_l2 norm / vol = " << sqrt(v2_l2 / VOL)  << std::endl;

  data[1] = ( VOL - VOL0 ) / VOL0;
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
  sprintf(outFileName, "turek3D_convergence_level%d_nprocs%dFS.txt",level, nprocs);

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
