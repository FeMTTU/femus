
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "XDMFWriter.hpp"
#include "NonLinearImplicitSystem.hpp"


#include "03_navier_stokes.hpp"



using namespace femus;

// bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int faceIndex, const double time) {
//   bool dirichlet = true; //dirichlet
//   value = 0.;
//
//   if (!strcmp(SolName, "T")) {
//     if (faceIndex == 2) {
//       value = 1.;
//     } else if (faceIndex == 3) {
//       dirichlet = false; //Neumann
//     }
//   } else if (!strcmp(SolName, "P")) {
//     dirichlet = false;
//   }
//
//   return dirichlet;
// }

unsigned dim;

bool SetBoundaryCondition(const MultiLevelProblem * mlProb,
                          const std::vector < double >& x, 
                          const char name[],
                          double& value,
                          const int FaceName, 
                          const double time ) {
  bool test = 1; //Dirichlet
  value = 0.;
  //   cout << "Time bdc : " <<  time << endl;
  if( !strcmp( name, "U" ) ) {
    if( 1 == FaceName ) { //inflow
      test = 1;
      if( dim == 2 ) {
        double um2D = 0.3;
        value = um2D * ( 4.0 / ( 0.1681 ) ) * x[1] * ( 0.41 - x[1] );
      }
      else {
        double um3D = 0.45;
        value = ( um3D * ( 16.0 / ( 0.02825761 ) ) * x[1] * x[2] * ( 0.41 - x[1] ) * ( 0.41 - x[2] ) );
      }
    }
    else if( 2 == FaceName ) { //outflow
      test = 0;
      //    test=1;
      value = 0.;
    }
    else if( 3 == FaceName ) { // no-slip fluid wall
      test = 1;
      value = 0.;
    }
    else if( 4 == FaceName ) { // no-slip solid wall
      test = 1;
      value = 0.;
    }
  }
  else if( !strcmp( name, "V" ) ) {
    if( 1 == FaceName ) {       //inflow
      test = 1;
      value = 0.;
    }
    else if( 2 == FaceName ) {  //outflow
      test = 0;
      //    test=1;
      value = 0.;
    }
    else if( 3 == FaceName ) {  // no-slip fluid wall
      test = 1;
      value = 0.;
    }
    else if( 4 == FaceName ) {  // no-slip solid wall
      test = 1;
      value = 0.;
    }
  }
  else if( !strcmp( name, "W" ) ) {
    if( 1 == FaceName ) {      //inflow
      test = 1;
      value = 0.;
    }
    else if( 2 == FaceName ) { //outflow
      //test=1;
      test = 0;
      value = 0.;
    }
    else if( 3 == FaceName ) { // no-slip fluid wall
      test = 1;
      value = 0.;
    }
    else if( 4 == FaceName ) { // no-slip fluid wall
      test = 1;
      value = 0.;
    }
  }
  else if( !strcmp( name, "P" ) ) {
    if( 1 == FaceName ) {
      test = 0;
      value = 0.;
    }
    else if( 2 == FaceName ) {
      test = 0;
      value = 0.;
    }
    else if( 3 == FaceName ) {
      test = 0;
      value = 0.;
    }
    else if( 4 == FaceName ) {
      test = 0;
      value = 0.;
    }
  }
  return test;
}

// //------------------------------------------------------------------------------------------------------------
unsigned numberOfUniformLevels;

bool SetRefinementFlag( const std::vector < double >& x, const int& elemgroupnumber, const int& level ) {

  bool refine = 0;

  //------------------ 3D --------------------------//
  //if (elemgroupnumber == 6 && level < 2) refine = 1;

  if( elemgroupnumber == 7 && level < numberOfUniformLevels ) refine = 1;

  if( elemgroupnumber == 8 && level < numberOfUniformLevels + 1 ) refine = 1;

  //------------------ 2D --------------------------//
  //if (elemgroupnumber == 7 && level < 2) refine = 1;

  //if (elemgroupnumber == 6 && level < 3) refine = 1;

  //if (elemgroupnumber == 5 && level < 4) refine = 1;

//   if (elemgroupnumber==6 && level<1) refine=1;
//   if (elemgroupnumber==7 && level<2) refine=1;
//   if (elemgroupnumber==8 && level<3) refine=1;

  return refine;

}

// //------------------------------------------------------------------------------------------------------------



int main( int argc, char** args ) {

  // init Petsc-MPI communicator
  FemusInit mpinit( argc, args, MPI_COMM_WORLD );

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  // mlMsh.ReadCoarseMesh("./input/cylinder2Dnew.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh( "./input/cylinder3Dnew.neu", "seventh", scalingFactor );
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  dim = mlMsh.GetDimension();

  numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 2;
  mlMsh.RefineMesh( numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag );

  MultiLevelSolution mlSol( &mlMsh );

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb( &mlSol );


  // add variables to mlSol
  mlSol.AddSolution( "U", LAGRANGE, SECOND );
  mlSol.AddSolution( "V", LAGRANGE, SECOND );
  if( dim == 3 ) mlSol.AddSolution( "W", LAGRANGE, SECOND );

  //mlSol.AddSolution("P", LAGRANGE, FIRST);
  mlSol.AddSolution( "P",  DISCONTINUOUS_POLYNOMIAL, FIRST );

  mlSol.AssociatePropertyToSolution( "P", "Pressure", false );
  mlSol.Initialize( "All" );

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction( SetBoundaryCondition );

  mlSol.GenerateBdc( "All", "Steady", & mlProb );

  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ( "NS" );

  system.AddSolutionToSystemPDE( "U" );
  system.AddSolutionToSystemPDE( "V" );
  if( dim == 3 ) system.AddSolutionToSystemPDE( "W" );

  system.AddSolutionToSystemPDE( "P" );

  system.SetLinearEquationSolverType(FEMuS_DEFAULT);
  // // // system.SetLinearEquationSolverType( FEMuS_ASM ); /// Additive Swartz Method  ///@todo does not work
  
  // attach the assembling function to system
  system.SetAssembleFunction( femus::AssembleNavierStokes_AD /*AssembleIncompressibleNavierStokes */);

  system.SetMaxNumberOfNonLinearIterations(10);
  system.SetNonLinearConvergenceTolerance(1.e-8);
  system.SetMaxNumberOfResidualUpdatesForNonlinearIteration(10);
  system.SetResidualUpdateConvergenceTolerance(1.e-15);

  system.SetMgType( F_CYCLE );

  system.SetNumberPreSmoothingStep( 1 );
  system.SetNumberPostSmoothingStep( 1 );
  // initialize and solve the system
  system.init();

  system.SetSolverFineGrids(GMRES);
  system.SetTolerances( 1.e-20, 1.e-20, 1.e+50, 50, 10 );
  
//   system.SetSolverFineGrids(RICHARDSON);
//   system.SetRichardsonScaleFactor(.5);
//   system.SetTolerances( 1.e-20, 1.e-20, 1.e+50, 50, 10 );
  
  system.SetPreconditionerFineGrids( ILU_PRECOND );

  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved( "All" );
  system.SetNumberOfSchurVariables( 1 );
  system.SetElementBlockNumber( 2 );
  //system.UseSamePreconditioner();
  
  system.SetOuterSolver(PREONLY);
  system.MGsolve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back( "All" );

  VTKWriter vtkIO( &mlSol );
  vtkIO.SetDebugOutput( true );
  vtkIO.Write( Files::_application_output_directory, "biquadratic", variablesToBePrinted );

  GMVWriter gmvIO( &mlSol );
  gmvIO.SetDebugOutput( true );
  gmvIO.Write( Files::_application_output_directory, "biquadratic", variablesToBePrinted );

  XDMFWriter xdmfIO( &mlSol );
  xdmfIO.SetDebugOutput( true );
  xdmfIO.Write( Files::_application_output_directory, "biquadratic", variablesToBePrinted );



  return 0;
}



