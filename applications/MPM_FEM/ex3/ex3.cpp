
#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "Marker.hpp"
#include "Line.hpp"

#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"

#include "NumericVector.hpp"
#include "adept.h"

#include "../include/mpmFem.hpp"

using namespace femus;

double SetVariableTimeStep(const double time) {
  double dt =  0.02;
  return dt;
}

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int facename, const double time) {
  bool test = 1; //dirichlet
  value = 0.;

  if(!strcmp(name, "DY")) {
    if(3 == facename || 4 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "DX")) {
    if(2 == facename) {
      test = 0;
      value = 0;
    }
  }

  return test;

}

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfUniformLevels = 4; //for refinement in 3D
  //unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;

  double Lref = 1.;
  double Uref = 1.;
  
  //initialize parameters for rolling ball (MPM)
  double rho_MPM = 1000.;
  double nu_MPM = 0.4;
  double E_MPM = 4.2 * 1.e6;
  
    //initialize parameters for plate (FEM)
  double rho_FEM = 10000.;
  double nu_FEM = 0.4;
  double E_FEM = 4.2 * 1.e7;

  Parameter par(Lref, Uref);

  // Generate Solid Object
  Solid solidMPM;
  Solid solidFEM;
  
  solidMPM = Solid(par, E_MPM, nu_MPM, rho_MPM, "Neo-Hookean");
  solidFEM = Solid(par, E_FEM, nu_FEM, rho_FEM, "Neo-Hookean");

  mlMsh.ReadCoarseMesh("../input/inclined_plane_coupled_2D.neu", "fifth", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);

  unsigned dim = mlMsh.GetDimension();

  MultiLevelSolution mlSol(&mlMsh);
  // add variables to mlSol
  mlSol.AddSolution("DX", LAGRANGE, SECOND, 2);
  if(dim > 1) mlSol.AddSolution("DY", LAGRANGE, SECOND, 2);
  if(dim > 2) mlSol.AddSolution("DZ", LAGRANGE, SECOND, 2);

  mlSol.AddSolution("VX", LAGRANGE, SECOND, 2);
  if(dim > 1) mlSol.AddSolution("VY", LAGRANGE, SECOND, 2);
  if(dim > 2) mlSol.AddSolution("VZ", LAGRANGE, SECOND, 2);

  mlSol.AddSolution("AX", LAGRANGE, SECOND, 2);
  if(dim > 1) mlSol.AddSolution("AY", LAGRANGE, SECOND, 2);
  if(dim > 2) mlSol.AddSolution("AZ", LAGRANGE, SECOND, 2);

  mlSol.AddSolution("M", LAGRANGE, FIRST, 2);
  mlSol.AddSolution("Mat", DISCONTINOUS_POLYNOMIAL, ZERO, 0, false);

  mlSol.Initialize("All");

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc("DX", "Steady");
  if(dim > 1) mlSol.GenerateBdc("DY", "Steady");
  if(dim > 2) mlSol.GenerateBdc("DZ", "Steady");
  mlSol.GenerateBdc("M", "Steady");

  MultiLevelProblem ml_prob(&mlSol);

  ml_prob.parameters.set<Solid> ("SolidMPM") = solidMPM;
  ml_prob.parameters.set<Solid> ("SolidFEM") = solidFEM;

  // ******* Add MPM system to the MultiLevel problem *******
  TransientNonlinearImplicitSystem& system = ml_prob.add_system < TransientNonlinearImplicitSystem > ("MPM_FEM");
  system.AddSolutionToSystemPDE("DX");
  if(dim > 1)system.AddSolutionToSystemPDE("DY");
  if(dim > 2) system.AddSolutionToSystemPDE("DZ");

  // ******* System MPM Assembly *******
  system.SetAssembleFunction(AssembleMPMSys);
  //system.SetAssembleFunction(AssembleFEM);
  // ******* set MG-Solver *******
  system.SetMgType(V_CYCLE);


  system.SetAbsoluteLinearConvergenceTolerance(1.0e-10);
  system.SetMaxNumberOfLinearIterations(1);
  system.SetNonLinearConvergenceTolerance(1.e-9);
  system.SetMaxNumberOfNonLinearIterations(20);

  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);

  // ******* Set Preconditioner *******
  system.SetMgSmoother(GMRES_SMOOTHER);

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids(GMRES);

  system.SetPreconditionerFineGrids(ILU_PRECOND);

  system.SetTolerances(1.e-10, 1.e-15, 1.e+50, 40, 40);


  unsigned size = 1;
  std::vector < std::vector < double > > x; // marker
  x.resize(size);
  x[0].resize(dim, 0.);
  x[0][1] = 0.05;

  double R = 1.6;
  double PI = acos(-1.);
  unsigned NR = 600;
  unsigned NL = NR / (2 * PI);
  double DL = R / NL;

  for(unsigned i = 0; i < NL; i++) {
    double  r = R - i * DL;
    unsigned Nr = static_cast <unsigned>(ceil(NR * r / R));
//     std::cout << r << " " << Nr << " ";
    double dtheta = 2 * PI / Nr;
    unsigned sizeOld = x.size();
    x.resize(sizeOld + Nr);
    for(unsigned s = sizeOld; s < x.size(); s++) {
      x[s].resize(dim);
    }
//     std::cout << x.size() << " ";
    for(unsigned j = 0; j < Nr; j++) {
      x[sizeOld + j][0] = r * cos(j * dtheta);
      x[sizeOld + j][1] = 0.05 + r * sin(j * dtheta);
    }
  }

  size = x.size();
  std::vector < MarkerType > markerType;
  markerType.resize(size);

  for(unsigned j = 0; j < size; j++) {
    markerType[j] = VOLUME;
  }

  std::vector < std::vector < std::vector < double > > > line(1);
  std::vector < std::vector < std::vector < double > > > line0(1);

  unsigned solType = 2;
  linea = new Line(x, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);

  double diskArea = PI * R * R;
  linea->SetParticlesMass(diskArea, rho_MPM);

  linea->GetLine(line0[0]);
  PrintLine(DEFAULT_OUTPUTDIR, line0, false, 0);


  linea->GetParticlesToGridMaterial();

  // ******* Print solution *******
  mlSol.SetWriter(VTK);

  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  mov_vars.push_back("DZ");
  mlSol.GetWriter()->SetMovingMesh(mov_vars);

  std::vector<std::string> print_vars;
  print_vars.push_back("All");


  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);

  gravity[0] = 9.81 * sqrt(2.) / 2.;
  gravity[1] = -9.81 * sqrt(2.) / 2.;

  system.AttachGetTimeIntervalFunction(SetVariableTimeStep);
  unsigned n_timesteps = 300;
  for(unsigned time_step = 1; time_step <= n_timesteps; time_step++) {

    system.CopySolutionToOldSolution();

    system.MGsolve();

    // ******* Print solution *******
    mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, time_step);

    GridToParticlesProjection(ml_prob, *linea);

    linea->GetLine(line[0]);
    PrintLine(DEFAULT_OUTPUTDIR, line, false, time_step);

  }

  delete linea;
  return 0;

} //end main


