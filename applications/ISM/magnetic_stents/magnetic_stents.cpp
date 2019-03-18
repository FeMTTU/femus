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
#include "GMVWriter.hpp"
#include "XDMFWriter.hpp"
#include "Line.hpp"
#include "adept.h"
#include "../../FSI/vascular/include/FSITimeDependentAssemblySupgNonConservativeTwoPressures.hpp"
#include <cmath>


#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

using boost::math::ellint_1;
using boost::math::ellint_2;

double scale = 1000.;

using namespace std;
using namespace femus;

double SetVariableTimeStep ( const double time );

bool SetBoundaryConditionAorticBifurcation ( const std::vector < double >& x, const char name[],
    double &value, const int facename, const double time );

bool SetBoundaryConditionTubo3D ( const std::vector < double >& x, const char name[],
                                  double &value, const int facename, const double time );

bool SetBoundaryConditionCarotidBifurcation ( const std::vector < double >& x, const char name[],
    double &value, const int facename, const double time );

bool SetBoundaryConditionMagneticStents ( const std::vector < double >& x, const char name[],
    double &value, const int facename, const double time );

void UpdateMeshCoordinates ( MultiLevelMesh &mlMesh, MultiLevelSolution& mlSol );

void MagneticForceWire ( const std::vector <double> & xMarker, std::vector <double> &Fm, const unsigned &material );
void MagneticForceSC ( const std::vector <double> & xMarker, std::vector <double> &Fm, const unsigned &material );
void MagneticForceStents ( const std::vector <double> & xMarker, std::vector <double> &Fm, const unsigned &material );
//------------------------------------------------------------------------------------------------------------------

unsigned partSim;
unsigned configuration;

int main ( int argc, char **args )
{

  // ******* Init Petsc-MPI communicator *******
  FemusInit mpinit ( argc, args, MPI_COMM_WORLD );

  twoPressure = true;
  valve = false;

  unsigned simulation = 0;

  if ( argc >= 2 ) {
    if ( !strcmp ( "0", args[1] ) ) { /** FSI Aortic Bifurcation 2D*/
      simulation = 0;
    }
    else if ( !strcmp ( "1", args[1] ) ) { /** FSI Tubo 3D */
      simulation = 1;
    }
    else if ( !strcmp ( "2", args[1] ) ) { /** FSI Carotid Bifurcation 3D*/
      simulation = 2;
    }
    else if ( !strcmp ( "3", args[1] ) ) { /** FSI Carotid Bifurcation 3D*/
      simulation = 3;
    }
  }

  //Files files;
  //files.CheckIODirectories();
  //files.RedirectCout();

  // ******* Extract the problem dimension and simulation identifier based on the inline input *******


  // ******* Extract the mesh.neu file name based on the simulation identifier *******
  std::string infile;
  bool dimension2D = true;

  if ( simulation == 0 ) {
    infile = "./input/aortic_bifurcation.neu";
  }
  else if ( simulation == 1 ) {
    infile = "./input/tubo3D.neu";
    //infile = "./input/tubo3D_thick.neu"; //used only once, not important
    dimension2D = false;
  }
  else if ( simulation == 2 ) {
    infile = "./input/carotid_bifurcation_3D.neu";
    dimension2D = false;
  }
  else if ( simulation == 3 ) {
    infile = "./input/magnetic_stents_solidWires.neu";
  }


  // ******* Set physics parameters *******
  double Lref, Uref, rhof, muf, rhos, ni, E, E1;

  Lref = 1.;
  Uref = 1.;

  rhof = 1035.;
  muf = 3.5 * 1.0e-3; //wrong=3.38*1.0e-4*rhof, note:3.38*1.0e-6*rhof=3.5*1.0e-3
  rhos = 1120;
  ni = 0.5;
  if ( simulation == 0 ) { //aortic_bifurcation
    E = 100 * 1.e6;
  }
  else if ( simulation == 1 ) { //tubo3D
    //E = 1000;
    E = 1. * 1.e12 ;
  }
  else if ( simulation == 2 ) { //carotide
    //E = 1.e6 * 1.e6; //CFD case
    E = 1 * 1.e6; // FSI with E = 1 MPa
    //E = 0.5 * 1.e6; // FSI with E = 0.5 MPa
  }
  else if ( simulation == 3 ) { //magnetic_stents
    E = 1. * 1.e6;
  }
  E1 = 1000;


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

  numberOfUniformRefinedMeshes = 1;
  numberOfAMRLevels = 0;

  std::cout << 0 << std::endl;

  MultiLevelMesh ml_msh ( numberOfUniformRefinedMeshes + numberOfAMRLevels, numberOfUniformRefinedMeshes,
                          infile.c_str(), "fifth", Lref, NULL );

//   MultiLevelMesh ml_msh1(numberOfUniformRefinedMeshes + numberOfAMRLevels, numberOfUniformRefinedMeshes,
//                          infile.c_str(), "fifth", Lref, NULL);

//   ml_msh.EraseCoarseLevels(numberOfUniformRefinedMeshes - 1);
//   numberOfUniformRefinedMeshes = 1;


  ml_msh.PrintInfo();

  // mark Solid nodes
  // ml_msh.MarkStructureNode();

  // ******* Init multilevel solution ******
  MultiLevelSolution ml_sol ( &ml_msh );

  ml_sol.SetIfFSI ( true );

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

  ml_sol.AddSolution ( "PF", DISCONTINUOUS_POLYNOMIAL, FIRST, 2 );
  ml_sol.AssociatePropertyToSolution ( "PF", "Pressure", false ); // Add this line

  ml_sol.AddSolution ( "lmbd", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false );

  ml_sol.AddSolution ( "Um", LAGRANGE, SECOND, 0, false );
  ml_sol.AddSolution ( "Vm", LAGRANGE, SECOND, 0, false );
  if ( !dimension2D ) ml_sol.AddSolution ( "Wm", LAGRANGE, SECOND, 0, false );

  // ******* Initialize solution *******
  ml_sol.Initialize ( "All" );

  if ( simulation == 0 ) {
    ml_sol.AttachSetBoundaryConditionFunction ( SetBoundaryConditionAorticBifurcation );
  }
  else if ( simulation == 1 ) {
    ml_sol.AttachSetBoundaryConditionFunction ( SetBoundaryConditionTubo3D );
  }
  else if ( simulation == 2 ) {
    ml_sol.AttachSetBoundaryConditionFunction ( SetBoundaryConditionCarotidBifurcation );
  }
  else if ( simulation == 3 ) {
    ml_sol.AttachSetBoundaryConditionFunction ( SetBoundaryConditionMagneticStents );
  }

  // ******* Set boundary conditions *******
  ml_sol.GenerateBdc ( "DX", "Steady" );
  ml_sol.GenerateBdc ( "DY", "Steady" );
  if ( !dimension2D ) ml_sol.GenerateBdc ( "DZ", "Steady" );

  if ( simulation == 0 ) {
    ml_sol.GenerateBdc ( "U", "Steady" );
    ml_sol.GenerateBdc ( "V", "Time_dependent" );
  }
  else if ( simulation == 1 ) {
    ml_sol.GenerateBdc ( "U", "Time_dependent" );
    ml_sol.GenerateBdc ( "V", "Steady" );
    ml_sol.GenerateBdc ( "W", "Steady" );
  }
  else if ( simulation == 2 ) {
    ml_sol.GenerateBdc ( "U", "Steady" );
    ml_sol.GenerateBdc ( "V", "Steady" );
    ml_sol.GenerateBdc ( "W", "Time_dependent" );
  }
  else  if ( simulation == 3 ) {
    ml_sol.GenerateBdc ( "U", "Time_dependent" );
    ml_sol.GenerateBdc ( "V", "Steady" );
  }

  ml_sol.GenerateBdc ( "PS", "Steady" );
  ml_sol.GenerateBdc ( "PF", "Steady" );


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
  if ( twoPressure ) system.AddSolutionToSystemPDE ( "PF" );

  // ******* System Fluid-Structure-Interaction Assembly *******
  system.SetAssembleFunction ( FSITimeDependentAssemblySupgNew2 );

  // ******* set MG-Solver *******
  system.SetMgType ( F_CYCLE );

  system.SetNonLinearConvergenceTolerance ( 1.e-7 );
  system.SetMaxNumberOfNonLinearIterations ( 20 );
  if ( dimension2D ) {
//     system.SetMaxNumberOfResidualUpdatesForNonlinearIteration(4);
//     system.SetResidualUpdateConvergenceTolerance(1.e-15);
    system.SetMaxNumberOfLinearIterations ( 1 );
    system.SetAbsoluteLinearConvergenceTolerance ( 1.e-50 );
  }
  else {
//     system.SetMaxNumberOfLinearIterations(4);
//     system.SetAbsoluteLinearConvergenceTolerance(1.e-15);
    system.SetMaxNumberOfLinearIterations ( 1 );
    system.SetAbsoluteLinearConvergenceTolerance ( 1.e-50 );
  }
  system.SetNumberPreSmoothingStep ( 1 );
  system.SetNumberPostSmoothingStep ( 1 );

  // ******* Set Preconditioner *******

  system.SetLinearEquationSolverType ( FEMuS_ASM );

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids ( RICHARDSON );
  //system.SetSolverFineGrids(GMRES);
  system.SetRichardsonScaleFactor ( 0.4 );

  system.SetPreconditionerFineGrids ( MLU_PRECOND );

  if ( dimension2D ) system.SetTolerances ( 1.e-10, 1.e-9, 1.e+50, 40, 40 );
  else system.SetTolerances ( 1.e-10, 1.e-50, 1.e+50, 40, 40 );

  // ******* Add variables to be solved *******
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved ( "All" );

  // ******* Set the last (1) variables in system (i.e. P) to be a schur variable *******
  system.SetNumberOfSchurVariables ( 2 );

  // ******* Set block size for the ASM smoothers *******
  if ( dimension2D ) system.SetElementBlockNumber ( 4 );
  else system.SetElementBlockNumber ( 3 );

  // ******* Print solution *******
  ml_sol.SetWriter ( VTK );


  std::vector<std::string> mov_vars;
  mov_vars.push_back ( "DX" );
  mov_vars.push_back ( "DY" );
  if ( !dimension2D ) mov_vars.push_back ( "DZ" );
  ml_sol.GetWriter()->SetMovingMesh ( mov_vars );

  std::vector<std::string> print_vars;
  print_vars.push_back ( "All" );

  ml_sol.GetWriter()->SetDebugOutput ( true );
  ml_sol.GetWriter()->Write ( DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0 );


  // ******* Solve *******
  std::cout << std::endl;
  std::cout << " *********** Fluid-Structure-Interaction ************  " << std::endl;


  //BEGIN INITIALIZE PARTICLES
  unsigned pSize = 50;
  double PI = acos ( -1. );
  std::vector < std::vector < double > > x ( pSize );
  std::vector < MarkerType > markerType ( pSize );

  if ( simulation == 0 ) { //for aortic_bifurcation

    for ( unsigned j = 0; j < pSize; j++ ) {
      x[j].resize ( 2 );
      x[j][0] = -0.008 + 0.016 * j / ( pSize - 1 );
      x[j][1] = 0.11;
      markerType[j] = VOLUME;
    }
  }

  if ( simulation == 1 ) { //for 3D tube
    unsigned theta_intervals = 100;
    unsigned radius_intervals = 9;
    pSize = radius_intervals * theta_intervals;
    x.resize ( pSize );
    markerType.resize ( pSize );
    srand ( 2 );
    for ( unsigned j = 0; j < pSize; j++ ) {
      double r_rad = static_cast <double> ( rand() ) / RAND_MAX;
      r_rad = 0.0034 * sqrt ( r_rad );
      double r_theta = static_cast <double> ( rand() ) / RAND_MAX * 2 * PI;
      x[j].resize ( 3 );
      x[j][0] = -0.035;
      x[j][1] = 0.0196 + r_rad * sin ( r_theta );
      x[j][2] = r_rad * cos ( r_theta );
    }
  }

  if ( simulation == 2 ) { //for carotidBifurcation
    unsigned theta_intervals = 100;
    unsigned radius_intervals = 9;
    pSize = radius_intervals * theta_intervals;
    x.resize ( pSize );
    markerType.resize ( pSize );
    srand ( 2 );
    for ( unsigned j = 0; j < pSize; j++ ) {
      double r_rad = static_cast <double> ( rand() ) / RAND_MAX;
      r_rad = 0.0034 * sqrt ( r_rad );
      double r_theta = static_cast <double> ( rand() ) / RAND_MAX * 2 * PI;
      x[j].resize ( 3 );
      x[j][0] = r_rad * cos ( r_theta );
      x[j][1] = 0.006 + r_rad * sin ( r_theta );
      x[j][2] = -0.06;
    }
//     unsigned counter = 0;
//     for (unsigned k = 1; k < radius_intervals + 1 ; k++) {
//       for (unsigned j = 0; j < theta_intervals; j++) {
//         x[counter].resize(3);
//         x[counter][0] = 0.00033 * k * sin(2.*PI / theta_intervals * j);
//         x[counter][1] = 0.006 + 0.00033 * k * cos(2.*PI / theta_intervals * j);
//         x[counter][2] = -0.06;
//         counter++;
//       }
//     }
  }

  if ( simulation == 3 ) { //for magnetic_stents
    pSize = 8. * 1.e3;
    x.resize ( pSize );
    double a = -0.0004 + 0.000002;
    double b = 0.0004 - 0.000002;
    double h = ( b - a ) / pSize;
    markerType.resize ( pSize );
    for ( unsigned j = 0; j < pSize; j++ ) {
      x[j].resize ( 2 );
      x[j][0] = 0.000002;
      x[j][1] = a + j * h;
    }
  }

//   std::cout << "SUUUUUUUUUUUUUUUUUUUUUUUUUUUUUCA" <<std::endl;

  //END INITIALIZE PARTICLES

  unsigned itPeriod = 32 /*for steady state with dt =100*/ /* for unsteady with dt =1/32, itPeriod = 32*/;
  unsigned confNumber;
  unsigned partSimMax;
  if ( simulation == 1 ) {
    confNumber = 4;
    partSimMax = 21;
  }
  else if ( simulation == 2 ) {
    confNumber = 6;
    partSimMax = 8;
  }
  else if ( simulation == 3 ) {
    confNumber = 1;
    partSimMax = 1;
  }
  else {
    confNumber = 1;
    partSimMax = 1;
  }


  std::vector < std::vector < std::vector < double > > > streamline ( pSize );
  std::vector < std::vector < std::vector < Line* > > > linea ( confNumber );

  for ( configuration = 0; configuration < confNumber; configuration++ ) {
    linea[configuration].resize ( partSimMax );
    for ( partSim = 0; partSim < partSimMax; partSim++ ) {
      linea[configuration][partSim].resize ( 1 );
      linea[configuration][partSim][0] =  new Line ( x, markerType, ml_sol.GetLevel ( numberOfUniformRefinedMeshes - 1 ), 2 );
      linea[configuration][partSim][0]->GetStreamLine ( streamline, 0 );
      linea[configuration][partSim][0]->GetStreamLine ( streamline, 1 );

      std::ostringstream output_path;

      double diam;
      if ( simulation == 1 ) diam = ( partSim + 1. ) * 0.1 * 1.e-6;
      else if ( simulation == 2 ) diam = ( partSim + 1. ) * 0.5 * 1.e-6;
      else if ( simulation == 3 ) diam = ( partSim + 1. ) * 2 * 0.43 * 1.e-6;
      else diam = 1.;

      output_path << "./output/particles-" << configuration << "-" << diam;

      PrintLine ( output_path.str(), streamline, true, 0 );
    }
  }

//      std::cout << "SUUUUUUUUUUUUUUUUUUUUUUUUUUUUUCA" <<std::endl;

  // time loop parameter
  system.AttachGetTimeIntervalFunction ( SetVariableTimeStep );

  const unsigned int n_timesteps = ( simulation == 1 ) ? 352 : 288 ;

  unsigned count_inside;
  unsigned count_out;
  unsigned count_tot = pSize;
  std::vector < std::vector <  double > > efficiencyVector ( confNumber );

  for ( unsigned time_step = 0; time_step < n_timesteps; time_step++ ) {

    for ( unsigned level = 0; level < numberOfUniformRefinedMeshes; level++ ) {
      SetLambdaNew ( ml_sol, level , SECOND, ELASTICITY );
    }
    if ( time_step > 0 )
      system.SetMgType ( V_CYCLE );
    system.CopySolutionToOldSolution();
    system.MGsolve();
    StoreMeshVelocity ( ml_prob );
    ml_sol.GetWriter()->Write ( DEFAULT_OUTPUTDIR, "biquadratic", print_vars, time_step + 1 );

    //if (time_step >= itPeriod) {
    for ( configuration = 0; configuration < confNumber; configuration++ ) {
      efficiencyVector[configuration].resize ( partSimMax );
      for ( partSim = 0; partSim < partSimMax; partSim++ ) {

        count_out = 0;

        if ( time_step >= 1. / itPeriod /*2.5*/ * itPeriod ) {
          for ( int i = 0; i < linea[configuration][partSim].size(); i++ ) {
            if ( simulation == 1 ) {
              linea[configuration][partSim][i]->AdvectionParallel ( 20, 1. / itPeriod, 4, MagneticForceWire );
            }
            else if ( simulation == 0 || simulation == 2 ) {
              linea[configuration][partSim][i]->AdvectionParallel ( 150, 1. / itPeriod, 4, MagneticForceSC );
              //linea[configuration][partSim][i]->AdvectionParallel(75, 1. / itPeriod, 4, MagneticForceSC);
            }
            else if ( simulation == 3 ) {
              linea[configuration][partSim][i]->AdvectionParallel ( 150, 1. / itPeriod, 4, MagneticForceStents );
            }
            count_out += linea[configuration][partSim][i]->NumberOfParticlesOutsideTheDomain();
          }
        }
//           if (time_step < 4 * itPeriod + itPeriod) {
//             count_tot += pSize;
//             linea.resize(time_step - itPeriod + 2);
//             linea[time_step - itPeriod + 1] =  new Line(x, markerType, ml_sol.GetLevel(numberOfUniformRefinedMeshes - 1), 2);
//           }

        count_inside = count_tot - count_out;

        efficiencyVector[configuration][partSim] =  static_cast< double > ( count_inside ) / count_tot;

        double diam;
        if ( simulation == 1 ) diam = ( partSim + 1. ) * 0.1 * 1.e-6;
        else if ( simulation == 2 ) diam = ( partSim + 1. ) * 0.5 * 1.e-6;
        else if ( simulation == 3 ) diam = ( partSim + 1. ) * 2 * 0.43 * 1.e-6;
        else diam = 0;

        std::cout << "configuration = " << configuration << std::endl;
        std::cout << "diameter = " << std::setw ( 11 ) << std::setprecision ( 12 ) << std::fixed << diam << std::endl;
        std::cout << "time_step = " << time_step + 1 << std::endl;
        std::cout << "particle inside = " << count_inside << std::endl;
        std::cout << "particle outside = " << count_out << std::endl;
        std::cout << "capture efficiency = " << efficiencyVector[configuration][partSim] << std::endl;

        linea[configuration][partSim][0]->GetStreamLine ( streamline, 0 );
        for ( int i = 0; i < linea[configuration][partSim].size(); i++ ) {
          linea[configuration][partSim][i]->GetStreamLine ( streamline, i + 1 );
        }

        std::ostringstream output_path;
        output_path << "./output/particles-" << configuration << "-" << diam;

        PrintLine ( output_path.str(), streamline, true, time_step + 1 );
      }
    }
    //}

  }


  // ******* Clear all systems *******
  for ( configuration = 0; configuration < confNumber; configuration++ ) {
    for ( partSim = 0; partSim < partSimMax; partSim++ ) {
      for ( unsigned i = 0; i < linea[configuration][partSim].size(); i++ ) {
        delete linea[configuration][partSim][i];
      }
    }
  }


  for ( unsigned j = 0; j < confNumber; j++ ) {
    std::cout << " CONFIGURATION " << j << std::endl;
    for ( unsigned i = 0; i < partSimMax; i++ ) {
      std::cout << efficiencyVector[j][i] * 100 << std::endl;
    }
    std::cout << " ---------------------------------------------------------------------------- " << std::endl;
  }

  ml_prob.clear();
  return 0;
}


//-------------------------------------------------------------------------------------------------------------------

double SetVariableTimeStep ( const double time )
{
  //double dt = 1./(64*1.4);

  double dt = 1. / 32; //1./16;
  //double dt = 1. / 4;

  //double dt = 60;
  return dt;
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionAorticBifurcation ( const std::vector < double >& x, const char name[], double & value, const int facename, const double time )
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos ( -1. );

  double ramp = ( time < 1 ) ? sin ( PI / 2 * time ) : 1.;

  if ( !strcmp ( name, "V" ) ) {
    if ( 1 == facename ) {
      double r2 = ( x[0] * 100. ) * ( x[0] * 100. );
      value = -0.15 / .81 * ( .81 - r2 ) * ( 1. + 0.25 * sin ( 2.*PI * time ) ) * ramp; //inflow
    }
    if ( 2 == facename || 3 == facename || 7 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "U" ) ) {
    if ( 2 == facename || 3 == facename || 7 == facename ) {
//       test = 0;
//       value = (10000 + 2500 * sin(2 * PI * time)) * ramp;;
//     }
//     else if (7 == facename) {
      test = 0;
      value = 0;
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
    if ( 7 == facename ) {
      test = 0;
      value = 0;
    }
  }
  else if ( !strcmp ( name, "DY" ) ) {
    if ( 7 == facename ) {
      test = 0;
      value = 0;
    }
  }

  return test;
}

//--------------------------------------------------------------------------------------------------

bool SetBoundaryConditionTubo3D ( const std::vector < double > & x, const char name[], double & value, const int facename, const double time )
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos ( -1. );
  double ramp = ( time < 1 ) ? sin ( PI / 2 * time ) : 1.;

  if ( !strcmp ( name, "U" ) ) {
    if ( 2 == facename ) {
      double r2 = ( ( x[1] - 0.0196 ) * ( x[1] - 0.0196 ) + ( x[2] * x[2] ) ) / ( 0.0035 * 0.0035 );
      value = 2 * 0.1 * ( 1. - r2 ) * ( 1. + 0.25 * sin ( 2.*PI * time ) ) * ramp; //inflow
      //value = 2 * 0.1 * (1. - r2) * ramp; //inflow
    }
    else if ( 1 == facename || 5 == facename ) {
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
    if ( 1 == facename ) {
      value = ( 12500 + 2500 * sin ( 2 * PI * time ) ) * ramp;
      //value = 13335 * ramp;
    }
  }
  else if ( !strcmp ( name, "PF" ) ) {
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

//----------------------------------------------------------------------------------------------

bool SetBoundaryConditionCarotidBifurcation ( const std::vector < double > & x, const char name[], double & value, const int facename, const double time )
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos ( -1. );
  double ramp = ( time < 1 ) ? sin ( PI / 2 * time ) : 1.;

  if ( !strcmp ( name, "W" ) ) {
    if ( 1 == facename ) {
      double r2 = ( ( x[0] * x[0] ) + ( x[1] - 0.006 ) * ( x[1] - 0.006 ) ) / ( 0.0035 * 0.0035 );
      value = 2 * 0.1 * ( 1. - r2 ) * ( 1. + 0.25 * sin ( 2.*PI * time ) ) * ramp; //inflow
      //value = 1.3 * 0.194 * (1. - r2) * (1. + 0.25 * sin(2.*PI * time)) * ramp; //inflow
//       double q;
//       double t = time - floor(time);
//       if (t >= 0. || t < 0.15) {
//         q = 946.67 * t * t * t + 130 * t * t + 0.9333 * t + 3.94;
//       }
//       else if (t >= 0.15 || t < 0.45) {
//         q = - 424889 * t * t * t * t * t * t + 876000 * t * t * t * t * t - 732289 * t * t * t * t + 317210 * t * t * t - 74927 * t * t + 9100.3 * t - 430.54;
//       }
//       else if (t >= 0.45 || t < 1.) {
//         q = - 601.31  * t * t * t * t * t * t + 3582.2 * t * t * t * t * t - 8384.8 * t * t * t * t + 10028 * t * t * t - 6510.4 * t * t + 2177.1 * t - 286.61;
//       }
//       value = 2 * q * 1.e-6 * (1. - r2) / (PI * 0.0035 * 0.0035) * ramp; //inflow
      //std::cout << "velocity we would like to have = " << value1 << " " << "time t is : " << t << std::endl;
    }
    else if ( 2 == facename || 3 == facename || 7 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "U" ) ) {
    if ( 7 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "V" ) ) {
    if ( 7 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "PS" ) ) {
    test = 0;
    value = 0.;
    if ( 2 == facename || 3 == facename ) {
      value = ( 5000 + 2500 * sin ( 2 * PI * time ) ) * ramp;
      //value = 5000 * ramp;//13332
    }
  }
  else if ( !strcmp ( name, "PF" ) ) {
    test = 0;
    value = 0.;
  }
  else if ( !strcmp ( name, "DX" ) || !strcmp ( name, "DY" ) || !strcmp ( name, "DZ" ) ) {
    if ( 7 == facename ) {
      test = 0;
      value = 0;
    }
  }

  return test;
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionMagneticStents ( const std::vector < double >& x, const char name[], double & value, const int facename, const double time )
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos ( -1. );

  double ramp = ( time < 1 ) ? sin ( PI / 2 * time ) : 1.;

  double u0 = 4.68; //0.58, 1.17, 2.34, 4.68, 11.7, 23.4, 52.6 cm/s

  if ( !strcmp ( name, "U" ) ) {
    if ( 1 == facename ) {
      double r2 = ( x[1] * x[1] ) / ( 0.0004 * 0.0004 );
      value = 1.5 * ( u0 * 1.e-2 ) * ( 1. - r2 ) * ramp; //steady-state inflow
      //value = -0.15 / .81 * (.81 - r2) * (1. + 0.25 * sin(2.*PI * time)) * ramp; //time dependent inflow
    }
    if ( 2 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "V" ) ) {
    if ( 2 == facename ) {
      test = 0;
      value = 0;
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

  return test;
}

//-------------------------------------------------------------------------------

void UpdateMeshCoordinates ( MultiLevelMesh & mlMesh, MultiLevelSolution & mlSol )
{

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* solution  = mlSol.GetSolutionLevel ( level );
  Mesh* msh0 = mlSol._mlMesh->GetLevel ( level );

  Mesh* msh = mlMesh.GetLevel ( level );

  const unsigned dim = msh->GetDimension();

  const char varname[3][3] = {"DX", "DY", "DZ"};
  vector <unsigned> indVAR ( dim );

  for ( unsigned k = 0; k < dim; k++ ) {
    indVAR[k] = mlSol.GetIndex ( &varname[k][0] );
  }


  for ( unsigned k = 0; k < dim; k++ ) {

    ( *msh->_topology->_Sol[k] ).zero();
    ( *msh->_topology->_Sol[k] ).close();

    ( *msh->_topology->_Sol[k] ).add ( ( *msh0->_topology->_Sol[k] ) );
    ( *msh->_topology->_Sol[k] ).close();

    ( *msh->_topology->_Sol[k] ).add ( ( *solution->_Sol[indVAR[k]] ) );
    ( *msh->_topology->_Sol[k] ).close();

  }

}

//-----------------------------------------------------------------------------------

void MagneticForceWire ( const std::vector <double> & xMarker, std::vector <double> &Fm, const unsigned & material )
{

  //infinitely long wire with a current I flowing modelled by the line identified by x and v

  //BEGIN magnetic and electric parameters

  double PI = acos ( -1. );
  //double I = 1.e5; // electric current intensity
  //double I = 0.5 * 1.e5;
  double I = 2 * 1.e5;
  double Msat = 1.e6;  //  magnetic saturation
  double  chi = 3.; //magnetic susceptibility
  double mu0 = 4 * PI * 1.e-7;  //magnetic permeability of the vacuum
  double H;
  std::vector <double> gradH ( 3 );
  std::vector <double> gradHSquared ( 3 );
  std::vector<double> vectorH ( 3 );
  double H0 = Msat / chi;

  //END


  //BEGIN fluid viscosity

  double muf = ( material == 2 ) ? 3.5 * 1.0e-3 : 1.0e100; // fluid viscosity

  //END


  //BEGIN geometric parameters

  double D = ( partSim + 1. ) * 0.1 * 1.e-6;      //diameter of the particle //rule con partSim

  std::vector <double> v ( 3 ); //direction vector of the line that identifies the infinite wire
  std::vector <double> x ( 3 ); //point that with v identifies the line of the wire

  if ( configuration == 0 ) {

    x[0] = 0.02093036072;
    x[1] = 0.02093036072;
    x[2] = 0.;

    v[0] = 0.;
    v[1] = 0.;
    v[2] = -1.;

  }

  else if ( configuration == 1 ) {

    x[0] = 0.006788225;
    x[1] = 0.006788225;
    x[2] = 0.;

    v[0] = 0.;
    v[1] = 0.;
    v[2] = -1.;

  }

  else if ( configuration == 2 ) {

    x[0] = 0.0196 * sqrt ( 2 ) / 2;
    x[1] = 0.0196 * sqrt ( 2 ) / 2;
    x[2] = 0.01;

    v[0] = sqrt ( 2 ) / 2.;
    v[1] = sqrt ( 2 ) / 2.;
    v[2] = 0.;

  }


  else if ( configuration == 3 ) {

    x[0] = 0.0196 * sqrt ( 2 ) / 2;
    x[1] = 0.0196 * sqrt ( 2 ) / 2;
    x[2] = 0.01;

    v[0] = - sqrt ( 2 ) / 2.;
    v[1] = sqrt ( 2 ) / 2.;
    v[2] = 0.;

  }


  //END


  //BEGIN extraction of the coordinates of the particle

  std::vector <double> xM ( 3 );
  xM[0] = xMarker[0];
  xM[1] = xMarker[1];
  xM[2] = ( xMarker.size() == 3 ) ? xMarker[2] : 0. ;

  //END


  //BEGIN evaluate H

  double Gamma;
  double Omega;
  std::vector<double> gradOmega ( 3 );

  Gamma = sqrt ( v[0] * v[0] + v[1] * v[1] + v[2] * v[2] );
  Omega = ( v[2] * ( x[1] - xM[1] ) - v[1] * ( x[2] - xM[2] ) ) * ( v[2] * ( x[1] - xM[1] ) - v[1] * ( x[2] - xM[2] ) ) +
          ( v[2] * ( x[0] - xM[0] ) - v[0] * ( x[2] - xM[2] ) ) * ( v[2] * ( x[0] - xM[0] ) - v[0] * ( x[2] - xM[2] ) ) +
          ( v[1] * ( x[0] - xM[0] ) - v[0] * ( x[1] - xM[1] ) ) * ( v[1] * ( x[0] - xM[0] ) - v[0] * ( x[1] - xM[1] ) ) ;

  gradOmega[0] = 2 * v[2] * ( v[2] * ( x[0] - xM[0] ) - v[0] * ( x[2] - xM[2] ) ) + 2 * v[1] * ( v[1] * ( x[0] - xM[0] ) - v[0] * ( x[1] - xM[1] ) );

  gradOmega[1] = 2 * v[2] * ( v[2] * ( x[1] - xM[1] ) - v[1] * ( x[2] - xM[2] ) ) - 2 * v[0] * ( v[1] * ( x[0] - xM[0] ) - v[0] * ( x[1] - xM[1] ) );

  gradOmega[2] = - 2 * v[1] * ( v[2] * ( x[1] - xM[1] ) - v[1] * ( x[2] - xM[2] ) ) - 2 * v[0] * ( v[2] * ( x[0] - xM[0] ) - v[0] * ( x[2] - xM[2] ) );

  H = ( I / ( 2 * PI ) ) * Gamma * ( 1 / sqrt ( Omega ) );

  gradHSquared[0] = ( I / ( 2 * PI ) ) * ( I / ( 2 * PI ) ) * ( ( -1 ) * Gamma * Gamma * ( 1 / ( Omega * Omega ) ) ) * gradOmega[0];
  gradHSquared[1] = ( I / ( 2 * PI ) ) * ( I / ( 2 * PI ) ) * ( ( -1 ) * Gamma * Gamma * ( 1 / ( Omega * Omega ) ) ) * gradOmega[1];
  gradHSquared[2] = ( I / ( 2 * PI ) ) * ( I / ( 2 * PI ) ) * ( ( -1 ) * Gamma * Gamma * ( 1 / ( Omega * Omega ) ) ) * gradOmega[2];

  gradH[0] = ( I / ( 2 * PI ) ) * ( ( -0.5 ) * Gamma * gradOmega[0] ) / sqrt ( Omega * Omega * Omega );
  gradH[1] = ( I / ( 2 * PI ) ) * ( ( -0.5 ) * Gamma * gradOmega[1] ) / sqrt ( Omega * Omega * Omega );
  gradH[2] = ( I / ( 2 * PI ) ) * ( ( -0.5 ) * Gamma * gradOmega[2] ) / sqrt ( Omega * Omega * Omega );


  //BEGIN evaluate Fm

  if ( true || H < H0 ) {


    double C1 = ( PI * D * D * D * mu0 * chi ) / 12;
    // printf("%g\n",C1);
    Fm[0] = C1 * gradHSquared[0];
    Fm[1] = C1 * gradHSquared[1];
    Fm[2] = C1 * gradHSquared[2];
  }

  else {


    double C2 = ( PI * D * D * D * mu0 * Msat ) / 6;

    //printf("%g\n",C2);
    Fm[0] = C2 * gradH[0];
    Fm[1] = C2 * gradH[1];
    Fm[2] = C2 * gradH[2];

  }


  for ( unsigned i = 0 ; i < Fm.size(); i++ ) {
    Fm[i] = Fm[i] / ( 3 * PI * D * muf ) ;
  }

  //END


  //BEGIN cheating to have attractive force

  for ( unsigned i = 0 ; i < Fm.size(); i++ ) {
    Fm[i] = - Fm[i] ;
    //printf("%g ",Fm[i]);
  }


  //END cheating

}

//---------------------------------------------------------------------------------

void MagneticForceSC ( const std::vector <double> & xMarker, std::vector <double> &Fm, const unsigned & material )
{

  //current loop of radius a, center x and for z-axis the line identified by x and v

  //BEGIN magnetic and electric parameters

  double PI = acos ( -1. );
  double I = 1.857e5; // electric current intensity
  double Msat = 1.e6;  //  magnetic saturation
  double  chi = 3.; //magnetic susceptibility
  double mu0 = 4 * PI * 1.e-7;  //magnetic permeability of the vacuum
  double H;
  std::vector <double> gradH ( 3 );
  std::vector <double> gradHSquared ( 3 );
  std::vector<double> vectorH ( 3 );
  double H0 = Msat / chi;

  //END


  //BEGIN fluid viscosity

  double muf = ( material == 2 ) ? 3.5 * 1.0e-3 : 1.0e100; // fluid viscosity

  //END


  //BEGIN geometric parameters

  double D = ( partSim + 1 ) * 0.5 * 1.e-6;     //diameter of the particle

  double a = 0.02; //radius of the circular current loop in m

  std::vector <double> v ( 3 );
  std::vector <double> x ( 3 );

  if ( configuration == 0 ) { //conf1 with z = 1.75 cm
    v[0] = 0.;
    v[1] = -1.;
    v[2] = 0.;

    x[0] = -0.002065;
    x[1] = -0.019932;
    x[2] = 0.000169;
  }
  else if ( configuration == 1 ) { //conf2 with z =1.75 cm
    v[0] = -1;
    v[1] = 0.;
    v[2] = 0.;

    x[0] = 0.0175;
    x[1] = 0.0045;
    x[2] = 0.;
  }
  else if ( configuration == 2 ) { //conf1 with z = 2.75 cm
    v[0] = 0.;
    v[1] = -1.;
    v[2] = 0.;

    x[0] = -0.002065;
    x[1] = -0.029932;
    x[2] = 0.000169;
  }
  else if ( configuration == 3 ) { //conf2 with z = 2.75 cm
    v[0] = -1;
    v[1] = 0.;
    v[2] = 0.;

    x[0] = 0.0275;
    x[1] = 0.0045;
    x[2] = 0.;
  }

  else if ( configuration == 4 ) { //conf 3 with z = 1.75 cm
    v[0] = 1;
    v[1] = 0.;
    v[2] = 0.;

    x[0] = -0.0285;
    x[1] = -0.001;
    x[2] = 0.;
  }

  else if ( configuration == 5 ) { //conf 3 with z = 2.75 cm
    v[0] = 1;
    v[1] = 0.;
    v[2] = 0.;

    x[0] = -0.0385;
    x[1] = -0.001;
    x[2] = 0.;
  }



  //END


  //BEGIN extraction of the coordinates of the particle

  std::vector <double> xM ( 3 );
  xM[0] = xMarker[0];
  xM[1] = xMarker[1];
  xM[2] = ( xMarker.size() == 3 ) ? xMarker[2] : 0. ;

  //END


  //BEGIN evaluate H


  //BEGIN bulid the rotation Matrix;
  double v2 = sqrt ( v[0] * v[0] + v[1] * v[1] + v[2] * v[2] );
  v[0] /= v2;
  v[1] /= v2;
  v[2] /= v2;

  std::vector <double> u ( 3 ); // will be the x'-axis
  std::vector <double> w ( 3 ); // will be the y'-axis

  unsigned imax = 0;
  double vmax = fabs ( v[0] );
  for ( unsigned i = 1; i < 3; i++ ) {
    if ( fabs ( v[i] ) > vmax ) {
      imax = i;
      vmax = fabs ( v[i] );
    }
  }

  u[0] = ( imax == 0 ) ? - ( v[1] + v[2] ) / v[0]  : 1;
  u[1] = ( imax == 1 ) ? - ( v[2] + v[0] ) / v[1]  : 1;
  u[2] = ( imax == 2 ) ? - ( v[0] + v[1] ) / v[2]  : 1;

  double u2 = sqrt ( u[0] * u[0] + u[1] * u[1] + u[2] * u[2] );
  u[0] /= u2;
  u[1] /= u2;
  u[2] /= u2;

  w[0] = v[1] * u[2] - v[2] * u[1];
  w[1] = v[2] * u[0] - v[0] * u[2];
  w[2] = v[0] * u[1] - v[1] * u[0];

  double w2 = sqrt ( w[0] * w[0] + w[1] * w[1] + w[2] * w[2] );
  if ( fabs ( w2 - 1. ) > 1.0e-12 ) {
    std::cout << "ERRORE AAAAAAAAAAAAAAAAAAA" << std::endl;
  }

  std::vector< std::vector <double> > R ( 3 );
  for ( unsigned i = 0; i < 3; i++ ) {
    R[i].resize ( 3 );
    R[i][0] = u[i];
    R[i][1] = w[i];
    R[i][2] = v[i];
  }
  //END bulid the rotation Matrix;

//           std::cout<<std::endl;
//           for(unsigned i=0;i<3;i++){
// 	    for(unsigned j=0;j<3;j++){
// 	      std::cout << R[i][j] <<  " ";
// 	    }
// 	    std::cout<<std::endl;
// 	  }


  //BEGIN find marker local coordinates
  for ( unsigned i = 0; i < 3 ; i++ ) {
    xM[i] -= x[i];
  }




  std::vector<double> xM1 ( 3, 0. );
  for ( unsigned i = 0; i < 3; i++ ) {
    for ( unsigned j = 0; j < 3; j++ ) {
      xM1[i] += R[j][i] * xM[j];
    }
  }
  xM = xM1;


  //END find marker local coordinates

  double rhoSquared = xM[0] * xM[0] + xM[1] * xM[1];
  double rSquared = rhoSquared + xM[2] * xM[2];
  double alphaSquared = a * a + rSquared - 2 * a * sqrt ( rhoSquared );
  double betaSquared = a * a + rSquared + 2 * a * sqrt ( rhoSquared );
  double kSquared = 1 - ( alphaSquared / betaSquared );
  double gamma = xM[0] * xM[0] - xM[1] * xM[1];
  double C = ( I * mu0 ) / PI;



  std::vector < double > vectorHp ( 3 );
  std::vector < std::vector < double > > jacobianVectorHp ( 3 );
  for ( unsigned i = 0; i < 3; i++ ) {
    jacobianVectorHp[i].resize ( 3 );
  }

  bool nearZAxis = ( rhoSquared < 1e-5 * a * a ) ? true : false;

  vectorHp[0] = ( nearZAxis == true ) ?  0 : ( 1 / mu0 ) * ( ( C * xM[0] * xM[2] ) / ( 2 * alphaSquared * sqrt ( betaSquared ) * rhoSquared ) ) * ( ( a * a + rSquared ) * ellint_2 ( kSquared ) - alphaSquared * ellint_1 ( kSquared ) );
  vectorHp[1] = ( nearZAxis == true ) ?  0 : ( 1 / mu0 ) * ( ( C * xM[1] * xM[2] ) / ( 2 * alphaSquared * sqrt ( betaSquared ) * rhoSquared ) ) * ( ( a * a + rSquared ) * ellint_2 ( kSquared ) - alphaSquared * ellint_1 ( kSquared ) );;
  vectorHp[2] = ( 1 / mu0 ) * ( C / ( 2 * alphaSquared * sqrt ( betaSquared ) ) ) * ( ( a * a - rSquared ) * ellint_2 ( kSquared ) + alphaSquared * ellint_1 ( kSquared ) );


  jacobianVectorHp[0][0] = ( nearZAxis == true ) ?  0 : ( ( C * xM[2] ) / ( 2 * alphaSquared * alphaSquared * sqrt ( betaSquared * betaSquared * betaSquared ) * rhoSquared * rhoSquared * mu0 ) ) * (
                             ( a * a * a * a * ( - gamma * ( 3 * xM[2] * xM[2] + a * a ) + rhoSquared * ( 8 * xM[0] * xM[0] - xM[1] * xM[1] ) ) -
                               a * a * ( rhoSquared * rhoSquared * ( 5 * xM[0] * xM[0] + xM[1] * xM[1] ) -
                                         2 * rhoSquared * xM[2] * xM[2] * ( 2 * xM[0] * xM[0] + xM[1] * xM[1] ) + 3 * xM[2] * xM[2] * xM[2] * xM[2] * gamma ) -
                               rSquared * rSquared * ( 2 * xM[0] * xM[0] * xM[0] * xM[0] + gamma * ( xM[1] * xM[1] + xM[2] * xM[2] ) ) )  * ellint_2 ( kSquared ) +
                             ( a * a * ( gamma * ( a * a + 2 * xM[2] * xM[2] ) - rhoSquared * ( 3 * xM[0] * xM[0] - 2 * xM[1] * xM[1] ) ) +
                               rSquared * ( 2 * xM[0] * xM[0] * xM[0] * xM[0] + gamma * ( xM[1] * xM[1] + xM[2] * xM[2] ) ) ) * alphaSquared * ellint_1 ( kSquared ) );


  jacobianVectorHp[0][1] = ( nearZAxis == true ) ?  0 : ( ( C * xM[0] * xM[1] * xM[2] ) / ( 2 * alphaSquared * alphaSquared * sqrt ( betaSquared * betaSquared * betaSquared ) * rhoSquared * rhoSquared * mu0 ) ) * (
                             ( 3 * a * a * a * a * ( 3 * rhoSquared - 2 * xM[2] * xM[2] ) - rSquared * rSquared * ( 2 * rSquared + rhoSquared ) - 2 * a * a * a * a * a * a  -
                               2 * a * a * ( 2 * rhoSquared * rhoSquared - rhoSquared * xM[2] * xM[2] + 3 * xM[2] * xM[2] * xM[2] * xM[2] ) ) * ellint_2 ( kSquared ) +
                             ( rSquared * ( 2 * rSquared + rhoSquared ) - a * a * ( 5 * rhoSquared - 4 * xM[2] * xM[2] ) + 2 * a * a * a * a ) * alphaSquared * ellint_1 ( kSquared ) );


  jacobianVectorHp[0][2] = ( nearZAxis == true ) ?  0 : ( ( C * xM[0] ) / ( 2 * alphaSquared * alphaSquared * sqrt ( betaSquared * betaSquared * betaSquared ) * rhoSquared * mu0 ) ) * (
                             ( ( rhoSquared - a * a ) * ( rhoSquared - a * a ) * ( rhoSquared + a * a ) + 2 * xM[2] * xM[2] * ( a * a * a * a - 6 * a * a * rhoSquared + rhoSquared * rhoSquared ) +
                               xM[2] * xM[2] * xM[2] * xM[2] * ( a * a + rhoSquared ) ) *  ellint_2 ( kSquared ) -
                             ( ( rhoSquared - a * a ) * ( rhoSquared - a * a ) + xM[2] * xM[2] * ( rhoSquared + a * a ) ) * alphaSquared *  ellint_1 ( kSquared ) );


  jacobianVectorHp[1][0] = jacobianVectorHp[0][1];


  jacobianVectorHp[1][1] = ( nearZAxis == true ) ?  0 : ( ( C * xM[2] ) / ( 2 * alphaSquared * alphaSquared * sqrt ( betaSquared * betaSquared * betaSquared ) * rhoSquared * rhoSquared * mu0 ) ) * (
                             ( a * a * a * a * ( gamma * ( 3 * xM[2] * xM[2] + a * a ) + rhoSquared * ( 8 * xM[1] * xM[1] - xM[0] * xM[0] ) ) -
                               a * a * ( rhoSquared * rhoSquared * ( 5 * xM[1] * xM[1] + xM[0] * xM[0] ) - 2 * rhoSquared * xM[2] * xM[2] * ( 2 * xM[1] * xM[1] + xM[0] * xM[0] ) -
                                         3 * xM[2] * xM[2] * xM[2] * xM[2] * gamma )  - rSquared * rSquared * ( 2 * xM[1] * xM[1] * xM[1] * xM[1] - gamma * ( xM[0] * xM[0] + xM[2] * xM[2] ) ) ) * ellint_2 ( kSquared ) +
                             ( a * a * ( - gamma * ( a * a + 2 * xM[2] * xM[2] )  - rhoSquared * ( 3 * xM[1] * xM[1] - 2 * xM[0] * xM[0] ) ) + rSquared * ( 2 * xM[1] * xM[1] * xM[1] * xM[1] -
                                 gamma * ( xM[0] * xM[0] + xM[2] * xM[2] ) ) ) * alphaSquared * ellint_1 ( kSquared ) );


  jacobianVectorHp[1][2] = ( nearZAxis == true ) ?  0 : ( ( C * xM[1] ) / ( 2 * alphaSquared * alphaSquared * sqrt ( betaSquared * betaSquared * betaSquared ) * rhoSquared * mu0 ) ) * (
                             ( ( rhoSquared - a * a ) * ( rhoSquared - a * a ) * ( rhoSquared + a * a ) + 2 * xM[2] * xM[2] * ( a * a * a * a - 6 * a * a * rhoSquared + rhoSquared * rhoSquared ) +
                               xM[2] * xM[2] * xM[2] * xM[2] * ( a * a + rhoSquared ) ) *  ellint_2 ( kSquared ) -
                             ( ( rhoSquared - a * a ) * ( rhoSquared - a * a ) + xM[2] * xM[2] * ( rhoSquared + a * a ) ) * alphaSquared *  ellint_1 ( kSquared ) );


  jacobianVectorHp[2][0] = ( nearZAxis == true ) ?  0 : jacobianVectorHp[0][2];


  jacobianVectorHp[2][1] = ( nearZAxis == true ) ?  0 : jacobianVectorHp[1][2];


  jacobianVectorHp[2][2] = ( ( C * xM[2] ) / ( 2 * alphaSquared * alphaSquared * sqrt ( betaSquared * betaSquared * betaSquared ) * mu0 ) ) * (
                             ( 6 * a * a * ( rhoSquared - xM[2] * xM[2] ) - 7 * a * a * a * a + rSquared * rSquared ) * ellint_2 ( kSquared ) +
                             ( a * a - rSquared ) * alphaSquared * ellint_1 ( kSquared ) );


// // 	  std::cout<<std::endl;
// //           for(unsigned i=0;i<3;i++){
// // 	    for(unsigned j=0;j<3;j++){
// // 	      std::cout << jacobianVectorHp[i][j] <<  " ";
// // 	    }
// // 	    std::cout<<std::endl;
// // 	  }
// // 	  std::cout<<std::endl;
// 	  for(unsigned i=0;i<3;i++){
// 	    std::cout << vectorHp[i] <<  " ";
// 	  }
//           std::cout<<std::endl;  std::cout<<std::endl;

  std::vector < std::vector <double> > jacobianVectorH ( 3 );
  for ( unsigned i = 0; i < 3; i++ ) {
    jacobianVectorH[i].resize ( 3 );
  }

  for ( unsigned i = 0; i < 3; i++ ) {
    vectorH[i] = 0.;
    for ( unsigned j = 0; j < 3; j++ ) {
      vectorH[i] += R[i][j] * vectorHp[j];
      jacobianVectorH[i][j] = 0.;
      // std::cout << R[i][j] <<" ";
      // std::cout << jacobianVectorHp[i][j] <<" ";
      for ( unsigned k = 0; k < 3; k++ ) {
        jacobianVectorH[i][j] += R[i][k] * jacobianVectorHp[k][j];
      }
    }
  }


  for ( unsigned i = 0; i < 3; i++ ) {
    for ( unsigned j = 0; j < 3; j++ ) {
      jacobianVectorHp[i][j] = 0. ;
      for ( unsigned k = 0; k < 3; k++ ) {
        jacobianVectorHp[i][j] += jacobianVectorH[i][k] * R[j][k];
      }
    }
  }

  jacobianVectorH = jacobianVectorHp;


// 	  std::cout<<std::endl;
//           for(unsigned i=0;i<3;i++){
// 	    for(unsigned j=0;j<3;j++){
// 	      std::cout << jacobianVectorH[i][j] <<  " ";
// 	    }
// 	    std::cout<<std::endl;
// 	  }
// 	  std::cout<<std::endl;
// 	  for(unsigned i=0;i<3;i++){
// 	    std::cout << vectorH[i] <<  " ";
// 	  }
//           std::cout<<std::endl;  std::cout<<std::endl;

// 	  std::cout<<std::endl;
//           for(unsigned i=0;i<3;i++){
// 	    for(unsigned j=0;j<3;j++){
// 	      std::cout << jacobianVectorH[i][j] <<  " ";
// 	    }
// 	    std::cout<<std::endl;
// 	  }
// 	  std::cout<<std::endl;
// 	  for(unsigned i=0;i<3;i++){
// 	    std::cout << vectorH[i] <<  " ";
// 	  }
//           std::cout<<std::endl;  std::cout<<std::endl;


  H = 0.;
  for ( unsigned  i = 0; i < 3; i++ ) {
    H += vectorH[i] * vectorH[i];
  }
  H = sqrt ( H );

  for ( unsigned  i = 0; i < 3; i++ ) {
    gradHSquared[i] = 0.;
    for ( unsigned  j = 0; j < 3; j++ ) {
      gradHSquared[i] += 2 * vectorH[j] * jacobianVectorH[j][i];
    }
  }

  for ( unsigned i = 0 ; i < Fm.size(); i++ ) {
    Fm[i] = - Fm[i] ;
    //printf("%g ",Fm[i]);
  }
//           gradHSquared[0] = 2 * vectorH[0] * jacobianVectorH[0][0] + 2 * vectorH[1] * jacobianVectorH[1][0] + 2 * vectorH[2] * jacobianVectorH[2][0];
//           gradHSquared[1] = 2 * vectorH[0] * jacobianVectorH[0][1] + 2 * vectorH[1] * jacobianVectorH[1][1] + 2 * vectorH[2] * jacobianVectorH[2][1];
//           gradHSquared[2] = 2 * vectorH[0] * jacobianVectorH[0][2] + 2 * vectorH[1] * jacobianVectorH[1][2] + 2 * vectorH[2] * jacobianVectorH[2][2];

  gradH[0] = 0.5 * ( 1 / H ) * gradHSquared[0];
  gradH[1] = 0.5 * ( 1 / H ) * gradHSquared[1];
  gradH[2] = 0.5 * ( 1 / H ) * gradHSquared[2];

//           for(unsigned i=0;i<3;i++){
// 	    std::cout <<  gradHSquared[i]<<" "<< gradH[i] <<  " ";
// 	  }
//           std::cout<<std::endl;  std::cout<<std::endl;


  //END valuate H


  //BEGIN evaluate Fm

  if ( H < H0 ) {


    double C1 = ( PI * D * D * D * mu0 * chi ) / 12;
    // printf("%g\n",C1);
    Fm[0] = C1 * gradHSquared[0];
    Fm[1] = C1 * gradHSquared[1];
    Fm[2] = C1 * gradHSquared[2];
  }

  else {

    double C2 = ( PI * D * D * D * mu0 * Msat ) / 6;

    //printf("%g\n",C2);
    Fm[0] = C2 * gradH[0];
    Fm[1] = C2 * gradH[1];
    Fm[2] = C2 * gradH[2];

  }


  for ( unsigned i = 0 ; i < Fm.size(); i++ ) {
    Fm[i] = Fm[i] / ( 3 * PI * D * muf ) ;
  }

  //END


  //BEGIN cheating to have attractive force

//   for (unsigned i = 0 ; i < Fm.size(); i++) {
//     Fm[i] = - Fm[i] ;
//     //printf("%g ",Fm[i]);
//   }


  //END cheating

}

//-------------------------------------------------------------------------------------

void MagneticForceStents ( const std::vector <double> & xMarker, std::vector <double> &Fm, const unsigned & material )
{

  //magnetic force when there are magnezable stents and the magnetic field is constant

  //BEGIN magnetic parameters
  double PI = acos ( -1. );
  double MsatWire = 1261 * 1.e3;  //  magnetic saturation
  double  chiWire = 1000.; //magnetic susceptibility
  double mu0 = 4 * PI * 1.e-7;  //magnetic permeability of the vacuum
  double H0 = 0.15 / mu0; //magnetic field intensity
  double theta = PI * 0.5;
  unsigned numberOfWires = 20;
  double fattore1 = chiWire / ( 2 + chiWire );
  double fattore2 = 0.5 * ( MsatWire / H0 );
  double alphaWire = ( fattore1 < fattore2 ) ? fattore1 : fattore2;
  double RWire = 62.5 * 1.e-6;
  double wireSpace = 0.2 * 1.e-2;
  std::vector <std::vector<double> > wiresCenter ( numberOfWires );
  for ( unsigned i = 0; i < wiresCenter.size(); i++ ) {
    wiresCenter[i].resize ( 2 );
    if ( i < numberOfWires / 2 ) {
      wiresCenter[i][0] = 0.15 * 1.e-2 + i * wireSpace;
      wiresCenter[i][1] = 3 * 1.e-4;
//       std::cout << "wire:" << i << " " << "x=" << wiresCenter[i][0] << " , " << wiresCenter[i][1] << std::endl;
    }
    else {
      wiresCenter[i][0] = 0.05 * 1.e-2 + ( i - numberOfWires / 2 ) * wireSpace;
      wiresCenter[i][1] = - 3 * 1.e-4;
//       std::cout << "wire:" << i << " " << "x=" << wiresCenter[i][0] << " , " << wiresCenter[i][1] << std::endl;
    }
  }
  //END

  std::vector <std::vector<double> > relPosition ( numberOfWires ); //relative position of the wires
  for ( unsigned i = 0; i < relPosition.size(); i++ ) {
    relPosition[i].resize ( 2 );
    for ( unsigned j = 0; j < 2; j++ ) {
      relPosition[i][j] = xMarker[j] - wiresCenter[i][j]; //particle position wrt the ith reference frame, i=1,...,numberOfWires
    }
  }

  //BEGIN GRADIENT PHI and HESSIAN
  // we compute the gradient of each phi and then sum them up to obtain sumGradPhi
  // sumHessPhi is the negative sum of the Hessian matrix of each phi
  std::vector<double> sumGradPhi ( 2, 0. );
  std::vector< std::vector <double> > sumHessPhi ( 2 );
  for ( unsigned i = 0; i < 2; i++ ) {
    sumHessPhi[i].assign ( 2, 0.);
  }

  double value = H0 * RWire * RWire * alphaWire;

  for ( unsigned i = 0; i < numberOfWires; i++ ) {
    
    double x = relPosition[i][0];
    double y = relPosition[i][1];

    double r2 =  x * x + y * y;
    double r4 = r2 * r2;
    
    double phi = value * ( x * cos(theta) + y * sin(theta) ) / r2;
    double phi_x = value * cos ( theta ) / r2 - phi * 2. * x / r2;
    double phi_y = value * sin ( theta ) / r2 - phi * 2. * y / r2;
    
    sumGradPhi[0] += phi_x;
    
    sumGradPhi[1] += phi_y;
    
    sumHessPhi[0][0] += value * cos(theta) * (- 2. * x ) / r4 - ( phi_x * 2. * x / r2 + phi * 2. / r2 + phi * 2. * x * (- 2. * x) / r4 );

    sumHessPhi[1][1] += value * sin(theta) * (- 2. * y ) / r4 - ( phi_y * 2. * y / r2 + phi * 2. / r2 + phi * 2. * y * (- 2. * y) / r4 );
    
    sumHessPhi[0][1] += value * cos(theta) * (- 2. * y ) / r4 - ( phi_y * 2. * x / r2 + phi * 2. * x * (- 2. * y) / r4 );
    
  }


  sumHessPhi[1][0] = sumHessPhi[0][1];
  //END



  //BEGIN B and JACOBIAN of B
  std::vector<double> B ( 2, 0. );
  B[0] = mu0 * ( H0 * cos ( theta )  - sumGradPhi[0] );
  B[1] = mu0 * ( H0 * sin ( theta )  - sumGradPhi[1] );

//   std::cout << "B1 = " << B[0] << " , " << "B2 = " << B[1] << std::endl;

  std::vector< std::vector <double> > JacB ( 2 );
  for ( unsigned i = 0; i < 2; i++ ) {
    JacB[i].resize ( 2 );
    for ( unsigned j = 0; j < 2; j++ ) {
      JacB[i][j] = - mu0 * sumHessPhi[i][j];
    }
  }
//END


  //BEGIN MAGNITUDE of B and its GRADIENT
  double magnitudeB = sqrt ( B[0] * B[0] + B[1] * B[1] );

//   std::vector <double> gradMagnitudeB ( 2, 0. ) ;
//
//   gradMagnitudeB[0] = mu0 * ( B[0] * sumHessPhi[0][0] + B[1] * sumHessPhi[0][1] ) / magnitudeB;
//   gradMagnitudeB[1] = mu0 * ( B[0] * sumHessPhi[0][1] + B[1] * sumHessPhi[1][1] ) / magnitudeB;
//END


  //BEGIN DIVERGENCE of MAGNETIC MOMENT
  double Dp = ( partSim + 1. ) * 2 * 0.43 * 1.e-6;  //particle diameter //rule con partSim

  double m_fm_p = 2.03 * 1.e-19; // magnitude of the magnetic moment of the magnetite in the particle
  double k = 1.38 * 1.e-23; //Boltzman constant
  double T = 300.; //absolute temperature
  double C1 = m_fm_p / ( k * T ) ;

  double beta = C1 * magnitudeB; //this is to highlight that beta depends on B and so on x and y

  double cothBeta =  cosh(beta) / sinh(beta) /*( exp ( beta ) + exp ( - beta ) ) / ( exp ( beta ) - exp ( - beta ) )*/;
  double Langevin = cothBeta - 1. / beta;

  double w_fm_p = 6.4;
  double Vp = ( 4. / 3 ) * PI * ( Dp * 0.5 ) * ( Dp * 0.5 ) * ( Dp * 0.5 ); //particle spherical volume
  double M_fm_p_s = 351.9 * 1.e3;
  double C2 = w_fm_p * Vp * M_fm_p_s;

  std::vector<double> magnMom ( 2, 0. );

  magnMom[0] = C2 * Langevin * B[0] / magnitudeB;
  magnMom[1] = C2 * Langevin * B[1] / magnitudeB;

  /*
    double dMagnMom1_dx = C2 * ( ( 1. - cothBeta * cothBeta ) * C1 * gradMagnitudeB[0] * B[0] / magnitudeB +
                                 B[0] * gradMagnitudeB[0] / pow ( magnitudeB, 3 ) +
                                 Langevin * ( ( JacB[0][0] * magnitudeB - B[0] * gradMagnitudeB[0] ) / ( magnitudeB * magnitudeB ) )
                               );
    double dMagnMom2_dy = C2 * ( ( 1. - cothBeta * cothBeta ) * C1 * gradMagnitudeB[1] * B[1] / magnitudeB +
                                 B[1] * gradMagnitudeB[1] / pow ( magnitudeB, 3 ) +
                                 Langevin * ( ( JacB[1][1] * magnitudeB - B[1] * gradMagnitudeB[1] ) / ( magnitudeB * magnitudeB ) )
                               );

    double divMagnMom = dMagnMom1_dx + dMagnMom2_dy;*/
//END



//BEGIN MAGNETIC FORCE
  double muf = ( material == 2 ) ? 3.5 * 1.0e-3 : 1.0e100; // fluid viscosity

  for ( unsigned i = 0 ; i < 2; i++ ) {
    Fm[i] = 0.;
    for ( unsigned j = 0 ; j < 2; j++ ) {
      Fm[i] += magnMom[j] * JacB[i][j] ;
    }
    Fm[i] /= ( 3. * PI * Dp * muf );
  }
  Fm[2] = 0.;
//END

}





