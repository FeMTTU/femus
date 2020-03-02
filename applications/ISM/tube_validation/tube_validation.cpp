#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"
#include "FemusInit.hpp"
#include "SparseMatrix.hpp"
#include "FElemTypeEnum.hpp"
#include "Files.hpp"
// #include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "TransientSystem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "XDMFWriter.hpp"
#include "Line.hpp"
#include "adept.h"
//<<<<<<< HEAD
#include "include/CFDSteadyStateAssembly.hpp"
//=======
//#include "../../FSI/vascular/include/FSITimeDependentAssemblySupg_OLD.hpp"
//>>>>>>> sara
#include <cmath>


#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

using boost::math::ellint_1;
using boost::math::ellint_2;

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

bool SetBoundaryConditionTubo3D(const std::vector < double >& x, const char name[],
                                double &value, const int facename, const double time);

bool SetBoundaryConditionCarotidBifurcation(const std::vector < double >& x, const char name[],
    double &value, const int facename, const double time);

void GetSolutionNorm(MultiLevelSolution& mlSol, const unsigned & group, std::vector <double> &data);

void UpdateMeshCoordinates(MultiLevelMesh &mlMesh, MultiLevelSolution& mlSol);

void MagneticForceWire(const std::vector <double> & xMarker, std::vector <double> &Fm, const unsigned &material);
void MagneticForceSC(const std::vector <double> & xMarker, std::vector <double> &Fm, const unsigned &material);
//------------------------------------------------------------------------------------------------------------------

unsigned partSim;
unsigned configuration;

int main(int argc, char **args)
{

  // ******* Init Petsc-MPI communicator *******
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  unsigned simulation = 0;

  if (argc >= 2) {
    if (!strcmp("0", args[1])) {   /** FSI Turek2D no stent */
      simulation = 0;
    }
    else if (!strcmp("1", args[1])) {    /** FSI Turek porous */
      simulation = 1;
    }
    else if (!strcmp("2", args[1])) {   /** FSI Turek stents 60 micron */
      simulation = 2;
    }
    else if (!strcmp("3", args[1])) {   /** FSI Turek 11 stents 60 micron */
      simulation = 3;
    }
    else if (!strcmp("4", args[1])) {   /** FSI AAA thrombus 2D */
      simulation = 4;
    }
    else if (!strcmp("5", args[1])) {   /** FSI Aortic Bifurcation 2D*/
      simulation = 5;
    }
    else if (!strcmp("6", args[1])) {   /** FSI Tubo 3D */
      simulation = 6;
    }
    else if (!strcmp("7", args[1])) {   /** FSI Carotid Bifurcation 3D*/
      simulation = 7;
    }
  }

  //Files files;
  //files.CheckIODirectories();
  //files.RedirectCout();

  // ******* Extract the problem dimension and simulation identifier based on the inline input *******


  // ******* Extract the mesh.neu file name based on the simulation identifier *******
//   std::string infile = "./input/aneurysm_Sara_5.neu";
  //std::string infile = "./input/Turek_porous_60micron.neu";
  //std::string infile = "./input/Turek_stents_60micron.neu";
  //std::string infile = "./input/Turek_11stents_60micron.neu";
  std::string infile;
  bool dimension2D = true;

  if (simulation == 0) {
    infile = "./input/Turek.neu";
  }
  else if (simulation == 1) {
    infile = "./input/Turek_porous_60micron.neu";
  }
  else if (simulation == 2) {
    infile = "./input/Turek_stents_60micron.neu";
  }
  else if (simulation == 3) {
    infile = "./input/Turek_11stents_60micron.neu";
  }
  else if (simulation == 4) {
    infile = "./input/AAA_thrombus_2D.neu";
  }
  else if (simulation == 5) {
    infile = "./input/aortic_bifurcation.neu";
  }
  else if (simulation == 6) {
    //infile = "./input/tubo3D.neu";
    infile = "./input/tubo3D_CFD.neu";
    dimension2D = false;
  }
  else if (simulation == 7) {
    infile = "./input/carotid_bifurcation_3D.neu";
    dimension2D = false;
  }



  // ******* Set physics parameters *******
   double Lref, Uref, rhof, muf, rhos, ni, E, E1;
// 
   Lref = 1.;
   Uref = 1.;
// 
//   rhof = 1035.;
//   muf = 3.5 * 1.0e-3; //wrong=3.38*1.0e-4*rhof, note:3.38*1.0e-6*rhof=3.5*1.0e-3
//   rhos = 1120;
//   ni = 0.5;
//   if (simulation == 5) {
//     E = 100 * 1.e6;
//   }
//   else if (simulation == 6) {
//     //E = 1000;
//     E = 1. * 1.e6 ;
//   }
//   else if (simulation == 7) { //carotide
//     E = 1000 * 1.e6;
//   }
//   else {
//     E = 1000000 * 1.e0; //turek: 1000000 * 1.e0;
//   }
//   E1 = 100000;
// 
//   Parameter par(Lref, Uref);

  // Generate Solid Object
//   Solid solid;
//   solid = Solid(par, E, ni, rhos, "Mooney-Rivlin");
// 
//   Solid solid1;
//   solid1 = Solid(par, E1, ni, rhos, "Mooney-Rivlin");
// 
//   cout << "Solid properties: " << endl;
//   cout << solid << endl;
// 
//   // Generate Fluid Object
//   Fluid fluid(par, muf, rhof, "Newtonian");
//   cout << "Fluid properties: " << endl;
//   cout << fluid << endl;

  // ******* Init multilevel mesh from mesh.neu file *******
  unsigned short numberOfUniformRefinedMeshes, numberOfAMRLevels;

  numberOfUniformRefinedMeshes = 2;
  numberOfAMRLevels = 0;

  //std::cout << 0 << std::endl;

  MultiLevelMesh ml_msh(numberOfUniformRefinedMeshes + numberOfAMRLevels, numberOfUniformRefinedMeshes,
                        infile.c_str(), "fifth", Lref, NULL);

//   MultiLevelMesh ml_msh1(numberOfUniformRefinedMeshes + numberOfAMRLevels, numberOfUniformRefinedMeshes,
//                          infile.c_str(), "fifth", Lref, NULL);

//   ml_msh.EraseCoarseLevels(numberOfUniformRefinedMeshes - 1);
//   numberOfUniformRefinedMeshes = 1;

  ml_msh.PrintInfo();

  // mark Solid nodes
// ml_msh.MarkStructureNode();

  // ******* Init multilevel solution ******
  MultiLevelSolution ml_sol(&ml_msh);

  ml_sol.SetIfFSI(false);

  // ******* Add solution variables to multilevel solution and pair them *******
//   ml_sol.AddSolution("DX", LAGRANGE, SECOND, 2);
//   ml_sol.AddSolution("DY", LAGRANGE, SECOND, 2);
//   if (!dimension2D) ml_sol.AddSolution("DZ", LAGRANGE, SECOND, 2);

  ml_sol.AddSolution("U", LAGRANGE, SECOND, 2); 
  ml_sol.AddSolution("V", LAGRANGE, SECOND, 2);
  if (!dimension2D) ml_sol.AddSolution("W", LAGRANGE, SECOND, 2);

  // Pair each velocity variable with the corresponding displacement variable
//   ml_sol.PairSolution("U", "DX"); // Add this line
//   ml_sol.PairSolution("V", "DY"); // Add this line
//   if (!dimension2D) ml_sol.PairSolution("W", "DZ"); // Add this line

  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  ml_sol.AddSolution("P", DISCONTINUOUS_POLYNOMIAL, FIRST, 2);
  ml_sol.AssociatePropertyToSolution("P", "Pressure", false); // Add this line

  //ml_sol.AddSolution("lmbd", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

  // ******* Initialize solution *******
  ml_sol.Initialize("All");

  if (simulation == 0 || simulation == 1 || simulation == 2 || simulation == 3) {
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionTurek2D);
  }
  else if (simulation == 4) {
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionThrombus2D);
  }
  else if (simulation == 5) {
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionAorticBifurcation);
  }
  else if (simulation == 6) {
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionTubo3D);
  }
  else if (simulation == 7) {
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionCarotidBifurcation);
  }

  // ******* Set boundary conditions *******
//   ml_sol.GenerateBdc("DX", "Steady");
//   ml_sol.GenerateBdc("DY", "Steady");
//   if (!dimension2D) ml_sol.GenerateBdc("DZ", "Steady");

  if (simulation == 4 || simulation == 5) {
    ml_sol.GenerateBdc("U", "Steady");
    ml_sol.GenerateBdc("V", "Time_dependent");
  }
  else if (simulation == 7) {
    ml_sol.GenerateBdc("U", "Steady");
    ml_sol.GenerateBdc("V", "Steady");
  }
  else {
    ml_sol.GenerateBdc("U");
    ml_sol.GenerateBdc("V");
  }

  if (!dimension2D && simulation != 7) ml_sol.GenerateBdc("W");
  if (simulation == 7) {
    ml_sol.GenerateBdc("W", "Time_dependent");
  }
  ml_sol.GenerateBdc("P");

//   for(unsigned level = 0; level < numberOfUniformRefinedMeshes; level++ ){
//     SetLambda(ml_sol, level , SECOND, ELASTICITY);
//   }

  // ******* Define the FSI Multilevel Problem *******

  MultiLevelProblem ml_prob(&ml_sol);
//   // Add fluid object
//   ml_prob.parameters.set<Fluid>("Fluid") = fluid;
//   // Add Solid Object
//   ml_prob.parameters.set<Solid>("Solid") = solid;
//   ml_prob.parameters.set<Solid>("Solid1") = solid1;

  // ******* Add FSI system to the MultiLevel problem *******
  //TransientMonolithicFSINonlinearImplicitSystem & system = ml_prob.add_system<TransientMonolithicFSINonlinearImplicitSystem> ("Fluid-Structure-Interaction");
  
  NonLinearImplicitSystem& system = ml_prob.add_system < NonLinearImplicitSystem > ("NS");
  
//   system.AddSolutionToSystemPDE("DX");
//   system.AddSolutionToSystemPDE("DY");
//   if (!dimension2D) system.AddSolutionToSystemPDE("DZ");
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  if (!dimension2D) system.AddSolutionToSystemPDE("W");
  system.AddSolutionToSystemPDE("P");

  // ******* System Assembly *******
  system.SetAssembleFunction(CFDSteadyStateAssembly);

  // ******* set MG-Solver *******
  system.SetMgType(F_CYCLE);

  system.SetNonLinearConvergenceTolerance(1.e-10);
  system.SetMaxNumberOfNonLinearIterations(4);
  if (dimension2D) {
    system.SetMaxNumberOfResidualUpdatesForNonlinearIteration(4);
    system.SetResidualUpdateConvergenceTolerance(1.e-15);
  }
  else {
//     system.SetMaxNumberOfResidualUpdatesForNonlinearIteration(10);
//     system.SetResidualUpdateConvergenceTolerance(1.e-15);
        
    system.SetMaxNumberOfLinearIterations(4);
    system.SetAbsoluteLinearConvergenceTolerance(1.e-15);
  }
  system.SetNumberPreSmoothingStep(0);
  system.SetNumberPostSmoothingStep(2);

  // ******* Set Preconditioner *******

  system.SetLinearEquationSolverType(FEMuS_ASM);

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids(RICHARDSON);
  //system.SetSolverFineGrids(GMRES);

  system.SetPreconditionerFineGrids(ILU_PRECOND);

  system.SetTolerances(1.e-12, 1.e-20, 1.e+50, 20, 10);

  // ******* Add variables to be solved *******
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");

  // ******* Set the last (1) variables in system (i.e. P) to be a schur variable *******
  system.SetNumberOfSchurVariables(1);

  // ******* Set block size for the ASM smoothers *******
  system.SetElementBlockNumber(2);
  system.SetMgType(V_CYCLE);
  system.MGsolve();
  
  // ******* Print solution *******
  ml_sol.SetWriter(VTK);

//   std::vector<std::string> mov_vars;
//   mov_vars.push_back("DX");
//   mov_vars.push_back("DY");
//   mov_vars.push_back("DZ");
//   ml_sol.GetWriter()->SetMovingMesh(mov_vars);
// 
//   std::vector<std::string> print_vars;
//   print_vars.push_back("All");
// 
//   ml_sol.GetWriter()->SetDebugOutput(true);
//   ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&ml_sol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  // print mesh info
  ml_msh.PrintInfo();

  // ******* Solve *******
//   std::cout << std::endl;
//   std::cout << " *********** Fluid-Structure-Interaction ************  " << std::endl;


  //BEGIN INITIALIZE PARTICLES
  unsigned pSize = 50;
  double PI = acos(-1.);
  std::vector < std::vector < double > > x(pSize);
  std::vector < MarkerType > markerType(pSize);


  if (simulation == 0) { //for Turek
    double y0 = 0.007 - 0.0019 / 2;
    double dy = 0.0019 / pSize;

    for (unsigned j = 0; j < pSize; j++) {
      x[j].resize(2);
      x[j][0] = -0.00005;
      x[j][1] = y0 + j * dy;
      markerType[j] = VOLUME;
    }
  }

  if (simulation == 5) { //for aorticBifurcation

    for (unsigned j = 0; j < pSize; j++) {
      x[j].resize(2);
      x[j][0] = -0.008 + 0.016 * j / (pSize - 1);
      x[j][1] = 0.11;
      markerType[j] = VOLUME;
    }
  }

  if (simulation == 6) { //for 3D tube
    unsigned theta_intervals = 100;
    unsigned radius_intervals = 9;
    pSize = radius_intervals * theta_intervals;
    x.resize(pSize);
    markerType.resize(pSize);
    srand(2);
    for (unsigned j = 0; j < pSize; j++) {
      double r_rad = static_cast <double> (rand()) / RAND_MAX;
      r_rad = 0.0034 * sqrt(r_rad);
      double r_theta = static_cast <double> (rand()) / RAND_MAX * 2 * PI;
      x[j].resize(3);
      x[j][0] = -0.035;
      x[j][1] = 0.0196 + r_rad * sin(r_theta);
      x[j][2] = r_rad * cos(r_theta);
    }
  }

  if (simulation == 7) { //for carotidBifurcation
    unsigned theta_intervals = 100;
    unsigned radius_intervals = 9;
    pSize = radius_intervals * theta_intervals;
    x.resize(pSize);
    markerType.resize(pSize);
    srand(2);
    for (unsigned j = 0; j < pSize; j++) {
      double r_rad = static_cast <double> (rand()) / RAND_MAX;
      r_rad = 0.0034 * sqrt(r_rad);
      double r_theta = static_cast <double> (rand()) / RAND_MAX * 2 * PI;
      x[j].resize(3);
      x[j][0] = r_rad * cos(r_theta);
      x[j][1] = 0.006 + r_rad * sin(r_theta);
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

  //END INITIALIZE PARTICLES

  unsigned itPeriod = 32 /*for steady state with dt =100*/ /* for unsteady with dt =1/32, itPeriod = 32*/;
  unsigned confNumber;
  unsigned partSimMax;
  if (simulation == 6) {
    confNumber = 4;
    partSimMax = 21;
  }
  else if (simulation == 7) {
    confNumber = 6;
    partSimMax = 8;
  }
  else {
    confNumber = 1;
    partSimMax = 1;
  }


  std::vector < std::vector < std::vector < double > > > streamline(pSize);
  std::vector < std::vector < std::vector < Line* > > > linea(confNumber);

  for (configuration = 0; configuration < confNumber; configuration++) {
    linea[configuration].resize(partSimMax);
    for (partSim = 0; partSim < partSimMax; partSim++) {
      linea[configuration][partSim].resize(1);
      linea[configuration][partSim][0] =  new Line(x, markerType, ml_sol.GetLevel(numberOfUniformRefinedMeshes - 1), 2);
      linea[configuration][partSim][0]->GetStreamLine(streamline, 0);
      linea[configuration][partSim][0]->GetStreamLine(streamline, 1);

      std::ostringstream output_path;

      double diam;
      if (simulation == 6) diam = (partSim + 1.) * 0.1 * 1.e-6;
      else if (simulation == 7) diam = (partSim + 1.) * 0.5 * 1.e-6;
      else diam = 1.;

      output_path << "./output/particles-" << configuration << "-" << diam;

      PrintLine(output_path.str(), "line", streamline, 0);
    }
  }

  // time loop parameter
  //system.AttachGetTimeIntervalFunction(SetVariableTimeStep);

  const unsigned int n_timesteps = (simulation == 6) ? 352 : 288 ;

  std::vector < std::vector <double> > data(n_timesteps);

  unsigned count_inside;
  unsigned count_out;
  unsigned count_tot = pSize;
  std::vector < std::vector <  double > > efficiencyVector(confNumber);

  for (unsigned time_step = 0; time_step < n_timesteps; time_step++) {
    
    data[time_step].resize(5);

//     for (unsigned level = 0; level < numberOfUniformRefinedMeshes; level++) {
//       SetLambda(ml_sol, level , SECOND, ELASTICITY);
//     }
//     if (time_step > 0)
//       system.SetMgType(V_CYCLE);
//     system.CopySolutionToOldSolution();
//     system.MGsolve();
//     ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, time_step + 1);

    //if (time_step >= itPeriod) {
    for (configuration = 0; configuration < confNumber; configuration++) {
      efficiencyVector[configuration].resize(partSimMax);
      for (partSim = 0; partSim < partSimMax; partSim++) {

        count_out = 0;

//         if (time_step >= 2 * itPeriod) {
          for (int i = 0; i < linea[configuration][partSim].size(); i++) {
            if (simulation == 6) {
              linea[configuration][partSim][i]->AdvectionParallel(20, 1. / itPeriod, 4, MagneticForceWire);
            }
            else if (simulation == 5 || simulation == 7) {
              linea[configuration][partSim][i]->AdvectionParallel(15, 1. / itPeriod, 4, MagneticForceSC);
            }
            count_out += linea[configuration][partSim][i]->NumberOfParticlesOutsideTheDomain();
          }
//         }

//           if (time_step < 4 * itPeriod + itPeriod) {
//             count_tot += pSize;
//             linea.resize(time_step - itPeriod + 2);
//             linea[time_step - itPeriod + 1] =  new Line(x, markerType, ml_sol.GetLevel(numberOfUniformRefinedMeshes - 1), 2);
//           }

        count_inside = count_tot - count_out;

        efficiencyVector[configuration][partSim] =  static_cast< double >(count_inside) / count_tot;

        double diam;
        if (simulation == 6) diam = (partSim + 1.) * 0.1 * 1.e-6;
        else if (simulation == 7) diam = (partSim + 1.) * 0.5 * 1.e-6;
        else diam = 0;

        std::cout << "configuration = " << configuration << std::endl;
        std::cout << "diameter = " << std::setw(11) << std::setprecision(12) << std::fixed << diam << std::endl;
        std::cout << "time_step = " << time_step + 1 << std::endl;
        std::cout << "particle inside = " << count_inside << std::endl;
        std::cout << "particle outside = " << count_out << std::endl;
        std::cout << "capture efficiency = " << efficiencyVector[configuration][partSim] << std::endl;

        linea[configuration][partSim][0]->GetStreamLine(streamline, 0);
        for (int i = 0; i < linea[configuration][partSim].size(); i++) {
          linea[configuration][partSim][i]->GetStreamLine(streamline, i + 1);
        }

        std::ostringstream output_path;
        output_path << "./output/particles-" << configuration << "-" << diam;

        PrintLine(output_path.str(), "streamline", streamline, time_step + 1);

        data[time_step][0] = time_step / 32.;
        //data[time_step][0] = time_step / (64*1.4);
        if (simulation == 0 || simulation == 1 || simulation == 2 || simulation == 3) {
          GetSolutionNorm(ml_sol, 9, data[time_step]);
        }
        else if (simulation == 4) {   //AAA_thrombus, 15=thrombus
          GetSolutionNorm(ml_sol, 7, data[time_step]);
        }
      }
    }
    //}

  }


  int  iproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  if (iproc == 0) {
    std::ofstream outf;
    if (simulation == 0) {
      outf.open("DataPrint_Turek.txt");
    }
    else if (simulation == 1) {
      outf.open("DataPrint_TurekPorous.txt");
    }
    else if (simulation == 2) {
      outf.open("DataPrint_TurekStents.txt");
    }
    else if (simulation == 3) {
      outf.open("DataPrint_Turek11Stents.txt");
    }
    else if (simulation == 4) {
      outf.open("DataPrint_AAA_thrombus_2D.txt");
    }


    if (!outf) {
      std::cout << "Error in opening file DataPrint.txt";
      return 1;
    }
    for (unsigned k = 0; k < n_timesteps; k++) {
      outf << data[k][0] << "\t" << data[k][1] << "\t" << data[k][2] << "\t" << data[k][3] << "\t" << data[k][4] << std::endl;
    }
    outf.close();
  }



  // ******* Clear all systems *******
  for (configuration = 0; configuration < confNumber; configuration++) {
    for (partSim = 0; partSim < partSimMax; partSim++) {
      for (unsigned i = 0; i < linea[configuration][partSim].size(); i++) {
        delete linea[configuration][partSim][i];
      }
    }
  }


  for (unsigned j = 0; j < confNumber; j++) {
    std::cout << " CONFIGURATION " << j << std::endl;
    for (unsigned i = 0; i < partSimMax; i++) {
      std::cout << efficiencyVector[j][i] * 100 << std::endl;
    }
    std::cout << " ---------------------------------------------------------------------------- " << std::endl;
  }

  ml_prob.clear();
  return 0;
}

double SetVariableTimeStep(const double time)
{
  //double dt = 1./(64*1.4);

  double dt = 1. / 32;
  //double dt = 1. / 4;

  //double dt = 60;

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

bool SetBoundaryConditionTurek2D(const std::vector < double >& x, const char name[], double & value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;


//   std::ifstream inf;
//   inf.open("./input/womersleyProfile64_R0p001_f84.txt");
//   if(!inf) {
//     std::cout << "velocity file ./input/womersleyProfile64_R0p001_f84.txt can not be opened\n";
//     exit(0);
//   }

//   std::vector<double> vel(64);
//
//   for(unsigned i = 0; i < 64; i++) {
//     inf >> vel[i];
//   }
//   inf.close();

//   double period = 1. / 1.4;
//   double dt = period / 64;
//
//   double time1 = time - floor(time / period) * period;
//
//   unsigned j = static_cast < unsigned >(floor(time1 / dt));

  //git pstd::cout<< name << " " << time <<" "<< j <<" "<<  vel[j] << std::endl;

  double PI = acos(-1.);
  if (!strcmp(name, "U")) {

    if (1 == facename) {
      double ramp = (time < 1) ? sin(PI / 2 * time) : 1.;
      //double ramp = (time < period) ? sin(PI/2 * time / period) : 1.;
      value = 0.05 * (x[1] * 1000 - 6) * (x[1] * 1000 - 8) * (1. + 0.75 * sin(2.*PI * time)) * ramp; //inflow
      //value = (x[1] * 1000 - 6) * ( x[1] * 1000 - 8) * vel[j] * ramp; //inflow
    }
    else if (2 == facename || 5 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if (!strcmp(name, "V")) {
    if (2 == facename || 5 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if (!strcmp(name, "P")) {
    test = 0;
    value = 0.;
  }
  else if (!strcmp(name, "DX")) {
    //if(2 == facename || 4 == facename || 5 == facename || 6 == facename) {
    if (5 == facename || 6 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if (!strcmp(name, "DY")) {
    //if(1 == facename || 3 == facename || 5 == facename || 6 == facename) {
    if (5 == facename || 6 == facename) {
      test = 0;
      value = 0;
    }
  }

  return test;

}


bool SetBoundaryConditionThrombus2D(const std::vector < double >& x, const char name[], double & value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos(-1.);

  double ramp = (time < 1) ? sin(PI / 2 * time) : 1.;

  if (!strcmp(name, "V")) {
    if (1 == facename) {
      double r2 = (x[0] * 100.) * (x[0] * 100.);
      //value = -0.01/.9 * (.9 - r2); //inflow
      value = -0.04 / .81 * (.81 - r2) * (1. + 0.75 * sin(2.*PI * time)) * ramp; //inflow
    }
    if (2 == facename || 5 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if (!strcmp(name, "U")) {
    if (2 == facename) {
      test = 0;
      value = (10000 + 2500 * sin(2 * PI * time)) * ramp;;
    }
    else if (5 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if (!strcmp(name, "P")) {
    test = 0;
    value = 0.;
  }
  else if (!strcmp(name, "DX")) {
    if (5 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if (!strcmp(name, "DY")) {
    if (5 == facename) {
      test = 0;
      value = 0;
    }
  }

  return test;
}


bool SetBoundaryConditionAorticBifurcation(const std::vector < double >& x, const char name[], double & value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos(-1.);

  double ramp = (time < 1) ? sin(PI / 2 * time) : 1.;

  if (!strcmp(name, "V")) {
    if (1 == facename) {
      double r2 = (x[0] * 100.) * (x[0] * 100.);
      value = -0.2 / .81 * (.81 - r2) * (1. + 0.25 * sin(2.*PI * time)) * ramp; //inflow
    }
    if (2 == facename || 3 == facename || 7 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if (!strcmp(name, "U")) {
    if (2 == facename || 3 == facename) {
      test = 0;
      value = (10000 + 2500 * sin(2 * PI * time)) * ramp;;
    }
    else if (7 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if (!strcmp(name, "P")) {
    test = 0;
    value = 0.;
  }
  else if (!strcmp(name, "DX")) {
    if (7 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if (!strcmp(name, "DY")) {
    if (7 == facename) {
      test = 0;
      value = 0;
    }
  }

  return test;
}


bool SetBoundaryConditionTubo3D(const std::vector < double > & x, const char name[], double & value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos(-1.);
  //double ramp = (time < 1) ? sin(PI / 2 * time) : 1.;

  if (!strcmp(name, "U")) {
    if (2 == facename) {
      double r2 = ((x[1] - 0.0196) * (x[1] - 0.0196) + (x[2] * x[2])) / (0.0035 * 0.0035);
      //value = 2 * 0.1 * (1. - r2) * (1. + 0.25 * sin(2.*PI * time)) * ramp; //inflow
      value = 2 * 0.1 * (1. - r2); //inflow
    }
    else if (1 == facename || 5 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if (!strcmp(name, "V") ){
    if (1 == facename || 5 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if (!strcmp(name, "W") ){
    if (5 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if (!strcmp(name, "P")) {
    test = 0;
    value = 0.;
    if (1 == facename) {
      //value = (12500 + 2500 * sin(2 * PI * time)) * ramp;
      value = 0;
    }
  }
//   else if (!strcmp(name, "DX") || !strcmp(name, "DY") || !strcmp(name, "DZ")) {
//     if (5 == facename) {
//       test = 0;
//       value = 0;
//     }
//   }

  return test;
}

bool SetBoundaryConditionCarotidBifurcation(const std::vector < double > & x, const char name[], double & value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos(-1.);
  double ramp = (time < 1) ? sin(PI / 2 * time) : 1.;

  if (!strcmp(name, "W")) {
    if (1 == facename) {
      double r2 = ((x[0] * x[0]) + (x[1] - 0.006) * (x[1] - 0.006)) / (0.0035 * 0.0035);
      value = 2 * 0.1 * (1. - r2) * (1. + 0.25 * sin(2.*PI * time)) * ramp; //inflow
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
    else if (2 == facename || 3 == facename || 7 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if (!strcmp(name, "U")) {
    if (7 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if (!strcmp(name, "V")) {
    if (7 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if (!strcmp(name, "P")) {
    test = 0;
    value = 0.;
    if (2 == facename || 3 == facename) {
      //value = (12500 + 2500 * sin(2 * PI * time)) * ramp;
      value = 5000 * ramp;//13332
    }
  }
  else if (!strcmp(name, "DX") || !strcmp(name, "DY") || !strcmp(name, "DZ")) {
    if (7 == facename) {
      test = 0;
      value = 0;
    }
  }

  return test;
}

void GetSolutionNorm(MultiLevelSolution & mlSol, const unsigned & group, std::vector <double> &data)
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

  if (nprocs == 1) {
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
  for (unsigned d = 0; d < dim; d++) {
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
  if (dim == 3) solVIndex[2] = mlSol.GetIndex("W");      // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol.GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  vector < unsigned > solDIndex(dim);
  solDIndex[0] = mlSol.GetIndex("DX");    // get the position of "U" in the ml_sol object
  solDIndex[1] = mlSol.GetIndex("DY");    // get the position of "V" in the ml_sol object
  if (dim == 3) solDIndex[2] = mlSol.GetIndex("DZ");      // get the position of "V" in the ml_sol object

  unsigned solDType = mlSol.GetSolutionType(solDIndex[0]);

  unsigned solPIndex;
  solPIndex = mlSol.GetIndex("P");
  unsigned solPType = mlSol.GetSolutionType(solPIndex);

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    if (msh->GetElementGroup(iel) == group) {
      short unsigned ielt = msh->GetElementType(iel);
      unsigned ndofV = msh->GetElementDofNumber(iel, solVType);
      unsigned ndofP = msh->GetElementDofNumber(iel, solPType);
      unsigned ndofD = msh->GetElementDofNumber(iel, solDType);
      // resize

      phiV.resize(ndofV);
      gradphiV.resize(ndofV * dim);
      nablaphiV.resize(ndofV * (3 * (dim - 1) + !(dim - 1)));

      solP.resize(ndofP);
      for (int d = 0; d < dim; d++) {
        solV[d].resize(ndofV);
        x0[d].resize(ndofD);
        x[d].resize(ndofD);
      }
      // get local to global mappings
      for (unsigned i = 0; i < ndofD; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solDType);
        for (unsigned d = 0; d < dim; d++) {
          x0[d][i] = (*msh->_topology->_Sol[d])(idof);

          x[d][i] = (*msh->_topology->_Sol[d])(idof) +
                    (*solution->_Sol[solDIndex[d]])(idof);
        }
      }

      for (unsigned i = 0; i < ndofV; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof
        for (unsigned  d = 0; d < dim; d++) {
          solV[d][i] = (*solution->_Sol[solVIndex[d]])(idof);      // global extraction and local storage for the solution
        }
      }



      for (unsigned i = 0; i < ndofP; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solPType);
        solP[i] = (*solution->_Sol[solPIndex])(idof);
      }


      for (unsigned ig = 0; ig < mlSol._mlMesh->_finiteElement[ielt][solVType]->GetGaussPointNumber(); ig++) {
        // *** get Jacobian and test function and test function derivatives ***
        msh->_finiteElement[ielt][solVType]->Jacobian(x0, ig, weight0, phiV, gradphiV, nablaphiV);
        msh->_finiteElement[ielt][solVType]->Jacobian(x, ig, weight, phiV, gradphiV, nablaphiV);
        phiP = msh->_finiteElement[ielt][solPType]->GetPhi(ig);

        vol0->add(iproc, weight0);
        vol->add(iproc, weight);

        std::vector < double> SolV2(dim, 0.);
        for (unsigned i = 0; i < ndofV; i++) {
          for (unsigned d = 0; d < dim; d++) {
            SolV2[d] += solV[d][i] * phiV[i];
          }
        }

        double V2 = 0.;
        for (unsigned d = 0; d < dim; d++) {
          V2 += SolV2[d] * SolV2[d];
        }
        v2->add(iproc, V2 * weight);

        double P2 = 0;
        for (unsigned i = 0; i < ndofP; i++) {
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
  p2_l2 = sqrt(p2_l2);
  double v2_l2 = v2->l1_norm();
  v2_l2 = sqrt(v2_l2);
  double VOL0 = vol0->l1_norm();
  double VOL = vol->l1_norm();

  std::cout.precision(14);
  std::scientific;
  std::cout << " vol0 = " << VOL0 << std::endl;
  std::cout << " vol = " << VOL << std::endl;
  std::cout << " (vol-vol0)/vol0 = " << (VOL - VOL0) / VOL0 << std::endl;
  std::cout << " p_l2 norm / vol = " << p2_l2 / VOL  << std::endl;
  std::cout << " v_l2 norm / vol = " << v2_l2 / VOL  << std::endl;

  data[1] = VOL0;
  data[2] = VOL;
  data[3] = p2_l2;
  data[4] = v2_l2;

  delete p2;
  delete v2;
  delete vol;

}

void UpdateMeshCoordinates(MultiLevelMesh & mlMesh, MultiLevelSolution & mlSol)
{

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* solution  = mlSol.GetSolutionLevel(level);
  Mesh* msh0 = mlSol._mlMesh->GetLevel(level);

  Mesh* msh = mlMesh.GetLevel(level);

  const unsigned dim = msh->GetDimension();

  const char varname[3][3] = {"DX", "DY", "DZ"};
  vector <unsigned> indVAR(dim);

  for (unsigned k = 0; k < dim; k++) {
    indVAR[k] = mlSol.GetIndex(&varname[k][0]);
  }


  for (unsigned k = 0; k < dim; k++) {

    (*msh->_topology->_Sol[k]).zero();
    (*msh->_topology->_Sol[k]).close();

    (*msh->_topology->_Sol[k]).add((*msh0->_topology->_Sol[k]));
    (*msh->_topology->_Sol[k]).close();

    (*msh->_topology->_Sol[k]).add((*solution->_Sol[indVAR[k]]));
    (*msh->_topology->_Sol[k]).close();

  }

}

void MagneticForceWire(const std::vector <double> & xMarker, std::vector <double> &Fm, const unsigned & material)
{

  //infinitely long wire with a current I flowing modelled by the line identified by x and v

  //BEGIN magnetic and electric parameters

  double PI = acos(-1.);
  double I = 1.e5; // electric current intensity
  double Msat = 1.e6;  //  magnetic saturation
  double  chi = 3.; //magnetic susceptibility
  double mu0 = 4 * PI * 1.e-7;  //magnetic permeability of the vacuum
  double H;
  std::vector <double> gradH(3);
  std::vector <double> gradHSquared(3);
  std::vector<double> vectorH(3);
  double H0 = Msat / chi;

  //END


  //BEGIN fluid viscosity

  double muf = (material == 2) ? 3.5 * 1.0e-3 : 1.0e100; // fluid viscosity

  //END


  //BEGIN geometric parameters

  double D =  (partSim + 1.) * 0.1 * 1.e-6;       //diameter of the particle //rule con partSim

  std::vector <double> v(3);    //direction vector of the line that identifies the infinite wire
  std::vector <double> x(3);   //point that with v identifies the line of the wire

  if (configuration == 0 ) {

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

    x[0] = 0.0196 * sqrt(2) / 2;
    x[1] = 0.0196 * sqrt(2) / 2;
    x[2] = 0.01;

    v[0] = sqrt(2) / 2.;
    v[1] = sqrt(2) / 2.;
    v[2] = 0.;

  }


  else if ( configuration == 3 ) {

    x[0] = 0.0196 * sqrt(2) / 2;
    x[1] = 0.0196 * sqrt(2) / 2;
    x[2] = 0.01;

    v[0] = - sqrt(2) / 2.;
    v[1] = sqrt(2) / 2.;
    v[2] = 0.;

  }


  //END


  //BEGIN extraction of the coordinates of the particle

  std::vector <double> xM(3);
  xM[0] = xMarker[0];
  xM[1] = xMarker[1];
  xM[2] = (xMarker.size() == 3) ? xMarker[2] : 0. ;

  //END


  //BEGIN evaluate H

  double Gamma;
  double Omega;
  std::vector<double> gradOmega(3);

  Gamma = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  Omega = (v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2])) * (v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2])) +
          (v[2] * (x[0] - xM[0]) - v[0] * (x[2] - xM[2])) * (v[2] * (x[0] - xM[0]) - v[0] * (x[2] - xM[2])) +
          (v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1])) * (v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1])) ;

  gradOmega[0] = 2 * v[2] * (v[2] * (x[0] - xM[0]) - v[0] * (x[2] - xM[2])) + 2 * v[1] * (v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1]));

  gradOmega[1] = 2 * v[2] * (v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2])) - 2 * v[0] * (v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1]));

  gradOmega[2] = - 2 * v[1] * (v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2])) - 2 * v[0] * (v[2] * (x[0] - xM[0]) - v[0] * (x[2] - xM[2]));

  H = (I / (2 * PI)) * Gamma * (1 / sqrt(Omega));

  gradHSquared[0] = (I / (2 * PI)) * (I / (2 * PI)) * ((-1) * Gamma * Gamma * (1 / (Omega * Omega))) * gradOmega[0];
  gradHSquared[1] = (I / (2 * PI)) * (I / (2 * PI)) * ((-1) * Gamma * Gamma * (1 / (Omega * Omega))) * gradOmega[1];
  gradHSquared[2] = (I / (2 * PI)) * (I / (2 * PI)) * ((-1) * Gamma * Gamma * (1 / (Omega * Omega))) * gradOmega[2];

  gradH[0] = (I / (2 * PI)) * ((-0.5) * Gamma * gradOmega[0]) / sqrt(Omega * Omega * Omega);
  gradH[1] = (I / (2 * PI)) * ((-0.5) * Gamma * gradOmega[1]) / sqrt(Omega * Omega * Omega);
  gradH[2] = (I / (2 * PI)) * ((-0.5) * Gamma * gradOmega[2]) / sqrt(Omega * Omega * Omega);


  //BEGIN evaluate Fm

  if (true || H < H0 ) {


    double C1 = (PI * D * D * D * mu0 * chi) / 12;
    // printf("%g\n",C1);
    Fm[0] = C1 * gradHSquared[0];
    Fm[1] = C1 * gradHSquared[1];
    Fm[2] = C1 * gradHSquared[2];
  }

  else {


    double C2 = (PI * D * D * D * mu0 * Msat) / 6;

    //printf("%g\n",C2);
    Fm[0] = C2 * gradH[0];
    Fm[1] = C2 * gradH[1];
    Fm[2] = C2 * gradH[2];

  }


  for (unsigned i = 0 ; i < Fm.size(); i++) {
    Fm[i] = Fm[i] / (3 * PI * D * muf) ;
  }

  //END


  //BEGIN cheating to have attractive force

  for (unsigned i = 0 ; i < Fm.size(); i++) {
    Fm[i] = - Fm[i] ;
    //printf("%g ",Fm[i]);
  }


  //END cheating

}


void MagneticForceSC(const std::vector <double> & xMarker, std::vector <double> &Fm, const unsigned & material)
{

  //current loop of radius a, center x and for z-axis the line identified by x and v

  //BEGIN magnetic and electric parameters

  double PI = acos(-1.);
  double I = 1.857e5; // electric current intensity
  double Msat = 1.e6;  //  magnetic saturation
  double  chi = 3.; //magnetic susceptibility
  double mu0 = 4 * PI * 1.e-7;  //magnetic permeability of the vacuum
  double H;
  std::vector <double> gradH(3);
  std::vector <double> gradHSquared(3);
  std::vector<double> vectorH(3);
  double H0 = Msat / chi;

  //END


  //BEGIN fluid viscosity

  double muf = (material == 2) ? 3.5 * 1.0e-3 : 1.0e100; // fluid viscosity

  //END


  //BEGIN geometric parameters

  double D = (partSim + 1) * 0.5 * 1.e-6;       //diameter of the particle

  double a = 0.04; //radius of the circular current loop in m

  std::vector <double> v(3);
  std::vector <double> x(3);

  if (configuration == 0) { //conf1 with z = 1.75 cm
    v[0] = 0.;
    v[1] = -1.;
    v[2] = 0.;

    x[0] = -0.002065;
    x[1] = -0.019932;
    x[2] = 0.000169;
  }
  else if (configuration == 1) { //conf2 with z =1.75 cm
    v[0] = -1;
    v[1] = 0.;
    v[2] = 0.;

    x[0] = 0.0175;
    x[1] = 0.0045;
    x[2] = 0.;
  }
  else if (configuration == 2) { //conf1 with z = 2.75 cm
    v[0] = 0.;
    v[1] = -1.;
    v[2] = 0.;

    x[0] = -0.002065;
    x[1] = -0.029932;
    x[2] = 0.000169;
  }
  else if (configuration == 3) { //conf2 with z = 2.75 cm
    v[0] = -1;
    v[1] = 0.;
    v[2] = 0.;

    x[0] = 0.0275;
    x[1] = 0.0045;
    x[2] = 0.;
  }
  
    else if (configuration == 4) { //conf 3 with z = 1.75 cm
    v[0] = 1;
    v[1] = 0.;
    v[2] = 0.;

    x[0] = -0.0285;
    x[1] = -0.001;
    x[2] = 0.;
  }
  
    else if (configuration == 5) { //conf 3 with z = 2.75 cm
    v[0] = 1;
    v[1] = 0.;
    v[2] = 0.;

    x[0] = -0.0385;
    x[1] = -0.001;
    x[2] = 0.;
  }
  
  

  //END


  //BEGIN extraction of the coordinates of the particle

  std::vector <double> xM(3);
  xM[0] = xMarker[0];
  xM[1] = xMarker[1];
  xM[2] = (xMarker.size() == 3) ? xMarker[2] : 0. ;

  //END


  //BEGIN evaluate H


  //BEGIN bulid the rotation Matrix;
  double v2 = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] /= v2;
  v[1] /= v2;
  v[2] /= v2;

  std::vector <double> u(3); // will be the x'-axis
  std::vector <double> w(3); // will be the y'-axis

  unsigned imax = 0;
  double vmax = fabs(v[0]);
  for (unsigned i = 1; i < 3; i++) {
    if (fabs(v[i]) > vmax) {
      imax = i;
      vmax = fabs(v[i]);
    }
  }

  u[0] = (imax == 0) ? - (v[1] + v[2]) / v[0]  : 1;
  u[1] = (imax == 1) ? - (v[2] + v[0]) / v[1]  : 1;
  u[2] = (imax == 2) ? - (v[0] + v[1]) / v[2]  : 1;

  double u2 = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  u[0] /= u2;
  u[1] /= u2;
  u[2] /= u2;

  w[0] = v[1] * u[2] - v[2] * u[1];
  w[1] = v[2] * u[0] - v[0] * u[2];
  w[2] = v[0] * u[1] - v[1] * u[0];

  double w2 = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
  if (fabs(w2 - 1.) > 1.0e-12) {
    std::cout << "ERRORE AAAAAAAAAAAAAAAAAAA" << std::endl;
  }

  std::vector< std::vector <double> > R(3);
  for (unsigned i = 0; i < 3; i++) {
    R[i].resize(3);
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
  for (unsigned i = 0; i < 3 ; i++) {
    xM[i] -= x[i];
  }




  std::vector<double> xM1(3, 0.);
  for (unsigned i = 0; i < 3; i++) {
    for (unsigned j = 0; j < 3; j++) {
      xM1[i] += R[j][i] * xM[j];
    }
  }
  xM = xM1;


  //END find marker local coordinates

  double rhoSquared = xM[0] * xM[0] + xM[1] * xM[1];
  double rSquared = rhoSquared + xM[2] * xM[2];
  double alphaSquared = a * a + rSquared - 2 * a * sqrt(rhoSquared);
  double betaSquared = a * a + rSquared + 2 * a * sqrt(rhoSquared);
  double kSquared = 1 - (alphaSquared / betaSquared);
  double gamma = xM[0] * xM[0] - xM[1] * xM[1];
  double C = (I * mu0) / PI;



  std::vector < double > vectorHp(3);
  std::vector < std::vector < double > > jacobianVectorHp(3);
  for (unsigned i = 0; i < 3; i++) {
    jacobianVectorHp[i].resize(3);
  }

  bool nearZAxis = (rhoSquared < 1e-5 * a * a) ? true : false;

  vectorHp[0] = (nearZAxis == true) ?  0 : (1 / mu0) * ((C * xM[0] * xM[2]) / (2 * alphaSquared * sqrt(betaSquared) * rhoSquared)) * ((a * a + rSquared) * ellint_2(kSquared) - alphaSquared * ellint_1(kSquared));
  vectorHp[1] = (nearZAxis == true) ?  0 : (1 / mu0) * ((C * xM[1] * xM[2]) / (2 * alphaSquared * sqrt(betaSquared) * rhoSquared)) * ((a * a + rSquared) * ellint_2(kSquared) - alphaSquared * ellint_1(kSquared));;
  vectorHp[2] = (1 / mu0) * (C / (2 * alphaSquared * sqrt(betaSquared))) * ((a * a - rSquared) * ellint_2(kSquared) + alphaSquared * ellint_1(kSquared));


  jacobianVectorHp[0][0] = (nearZAxis == true) ?  0 : ((C * xM[2]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * rhoSquared * rhoSquared * mu0)) * (
                             (a * a * a * a * (- gamma * (3 * xM[2] * xM[2] + a * a) + rhoSquared * (8 * xM[0] * xM[0] - xM[1] * xM[1])) -
                              a * a * (rhoSquared * rhoSquared * (5 * xM[0] * xM[0] + xM[1] * xM[1]) -
                                       2 * rhoSquared * xM[2] * xM[2] * (2 * xM[0] * xM[0] + xM[1] * xM[1]) + 3 * xM[2] * xM[2] * xM[2] * xM[2] * gamma) -
                              rSquared * rSquared * (2 * xM[0] * xM[0] * xM[0] * xM[0] + gamma * (xM[1] * xM[1] + xM[2] * xM[2])))  * ellint_2(kSquared) +
                             (a * a * (gamma * (a * a + 2 * xM[2] * xM[2]) - rhoSquared * (3 * xM[0] * xM[0] - 2 * xM[1] * xM[1])) +
                              rSquared * (2 * xM[0] * xM[0] * xM[0] * xM[0] + gamma * (xM[1] * xM[1] + xM[2] * xM[2]))) * alphaSquared * ellint_1(kSquared));


  jacobianVectorHp[0][1] = (nearZAxis == true) ?  0 : ((C * xM[0] * xM[1] * xM[2]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * rhoSquared * rhoSquared * mu0)) * (
                             (3 * a * a * a * a * (3 * rhoSquared - 2 * xM[2] * xM[2]) - rSquared * rSquared * (2 * rSquared + rhoSquared) - 2 * a * a * a * a * a * a  -
                              2 * a * a * (2 * rhoSquared * rhoSquared - rhoSquared * xM[2] * xM[2] + 3 * xM[2] * xM[2] * xM[2] * xM[2])) * ellint_2(kSquared) +
                             (rSquared * (2 * rSquared + rhoSquared) - a * a * (5 * rhoSquared - 4 * xM[2] * xM[2]) + 2 * a * a * a * a) * alphaSquared * ellint_1(kSquared));


  jacobianVectorHp[0][2] = (nearZAxis == true) ?  0 : ((C * xM[0]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * rhoSquared * mu0)) * (
                             ((rhoSquared - a * a) * (rhoSquared - a * a) * (rhoSquared + a * a) + 2 * xM[2] * xM[2] * (a * a * a * a - 6 * a * a * rhoSquared + rhoSquared * rhoSquared) +
                              xM[2] * xM[2] * xM[2] * xM[2] * (a * a + rhoSquared)) *  ellint_2(kSquared) -
                             ((rhoSquared - a * a) * (rhoSquared - a * a) + xM[2] * xM[2] * (rhoSquared + a * a)) * alphaSquared *  ellint_1(kSquared));


  jacobianVectorHp[1][0] = jacobianVectorHp[0][1];


  jacobianVectorHp[1][1] = (nearZAxis == true) ?  0 : ((C * xM[2]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * rhoSquared * rhoSquared * mu0)) * (
                             (a * a * a * a * (gamma * (3 * xM[2] * xM[2] + a * a) + rhoSquared * (8 * xM[1] * xM[1] - xM[0] * xM[0])) -
                              a * a * (rhoSquared * rhoSquared * (5 * xM[1] * xM[1] + xM[0] * xM[0]) - 2 * rhoSquared * xM[2] * xM[2] * (2 * xM[1] * xM[1] + xM[0] * xM[0]) -
                                       3 * xM[2] * xM[2] * xM[2] * xM[2] * gamma)  - rSquared * rSquared * (2 * xM[1] * xM[1] * xM[1] * xM[1] - gamma * (xM[0] * xM[0] + xM[2] * xM[2]))) * ellint_2(kSquared) +
                             (a * a * (- gamma * (a * a + 2 * xM[2] * xM[2])  - rhoSquared * (3 * xM[1] * xM[1] - 2 * xM[0] * xM[0])) + rSquared * (2 * xM[1] * xM[1] * xM[1] * xM[1] -
                                 gamma * (xM[0] * xM[0] + xM[2] * xM[2]))) * alphaSquared * ellint_1(kSquared));


  jacobianVectorHp[1][2] = (nearZAxis == true) ?  0 : ((C * xM[1]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * rhoSquared * mu0)) * (
                             ((rhoSquared - a * a) * (rhoSquared - a * a) * (rhoSquared + a * a) + 2 * xM[2] * xM[2] * (a * a * a * a - 6 * a * a * rhoSquared + rhoSquared * rhoSquared) +
                              xM[2] * xM[2] * xM[2] * xM[2] * (a * a + rhoSquared)) *  ellint_2(kSquared) -
                             ((rhoSquared - a * a) * (rhoSquared - a * a) + xM[2] * xM[2] * (rhoSquared + a * a)) * alphaSquared *  ellint_1(kSquared));


  jacobianVectorHp[2][0] = (nearZAxis == true) ?  0 : jacobianVectorHp[0][2];


  jacobianVectorHp[2][1] = (nearZAxis == true) ?  0 : jacobianVectorHp[1][2];


  jacobianVectorHp[2][2] = ((C * xM[2]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * mu0)) * (
                             (6 * a * a * (rhoSquared - xM[2] * xM[2]) - 7 * a * a * a * a + rSquared * rSquared) * ellint_2(kSquared) +
                             (a * a - rSquared) * alphaSquared * ellint_1(kSquared));


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

  std::vector < std::vector <double> > jacobianVectorH(3);
  for (unsigned i = 0; i < 3; i++) {
    jacobianVectorH[i].resize(3);
  }

  for (unsigned i = 0; i < 3; i++) {
    vectorH[i] = 0.;
    for (unsigned j = 0; j < 3; j++) {
      vectorH[i] += R[i][j] * vectorHp[j];
      jacobianVectorH[i][j] = 0.;
      // std::cout << R[i][j] <<" ";
      // std::cout << jacobianVectorHp[i][j] <<" ";
      for (unsigned k = 0; k < 3; k++) {
        jacobianVectorH[i][j] += R[i][k] * jacobianVectorHp[k][j];
      }
    }
  }


  for (unsigned i = 0; i < 3; i++) {
    for (unsigned j = 0; j < 3; j++) {
      jacobianVectorHp[i][j] = 0. ;
      for (unsigned k = 0; k < 3; k++) {
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
  for (unsigned  i = 0; i < 3; i++) {
    H += vectorH[i] * vectorH[i];
  }
  H = sqrt(H);

  for (unsigned  i = 0; i < 3; i++) {
    gradHSquared[i] = 0.;
    for (unsigned  j = 0; j < 3; j++) {
      gradHSquared[i] += 2 * vectorH[j] * jacobianVectorH[j][i];
    }
  }

  for (unsigned i = 0 ; i < Fm.size(); i++) {
    Fm[i] = - Fm[i] ;
    //printf("%g ",Fm[i]);
  }
//           gradHSquared[0] = 2 * vectorH[0] * jacobianVectorH[0][0] + 2 * vectorH[1] * jacobianVectorH[1][0] + 2 * vectorH[2] * jacobianVectorH[2][0];
//           gradHSquared[1] = 2 * vectorH[0] * jacobianVectorH[0][1] + 2 * vectorH[1] * jacobianVectorH[1][1] + 2 * vectorH[2] * jacobianVectorH[2][1];
//           gradHSquared[2] = 2 * vectorH[0] * jacobianVectorH[0][2] + 2 * vectorH[1] * jacobianVectorH[1][2] + 2 * vectorH[2] * jacobianVectorH[2][2];

  gradH[0] = 0.5 * (1 / H) * gradHSquared[0];
  gradH[1] = 0.5 * (1 / H) * gradHSquared[1];
  gradH[2] = 0.5 * (1 / H) * gradHSquared[2];

//           for(unsigned i=0;i<3;i++){
// 	    std::cout <<  gradHSquared[i]<<" "<< gradH[i] <<  " ";
// 	  }
//           std::cout<<std::endl;  std::cout<<std::endl;


  //END valuate H


  //BEGIN evaluate Fm

  if (H < H0) {


    double C1 = (PI * D * D * D * mu0 * chi) / 12;
    // printf("%g\n",C1);
    Fm[0] = C1 * gradHSquared[0];
    Fm[1] = C1 * gradHSquared[1];
    Fm[2] = C1 * gradHSquared[2];
  }

  else {

    double C2 = (PI * D * D * D * mu0 * Msat) / 6;

    //printf("%g\n",C2);
    Fm[0] = C2 * gradH[0];
    Fm[1] = C2 * gradH[1];
    Fm[2] = C2 * gradH[2];

  }


  for (unsigned i = 0 ; i < Fm.size(); i++) {
    Fm[i] = Fm[i] / (3 * PI * D * muf) ;
  }

  //END


  //BEGIN cheating to have attractive force

//   for (unsigned i = 0 ; i < Fm.size(); i++) {
//     Fm[i] = - Fm[i] ;
//     //printf("%g ",Fm[i]);
//   }


  //END cheating

}







