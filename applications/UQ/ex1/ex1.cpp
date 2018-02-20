
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

#include "../include/sfem_assembly.hpp"

using namespace femus;


bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int facename, const double time) {
  bool test = 1; //dirichlet
  value = 0.;

  if(!strcmp(name, "UX")) {
    if(3 == facename || 4 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "UY")) {
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

  mlMsh.ReadCoarseMesh("../input/square.neu", "fifth", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);

  unsigned dim = mlMsh.GetDimension();

  MultiLevelSolution mlSol(&mlMsh);
  // add variables to mlSol
  mlSol.AddSolution("UX", LAGRANGE, SECOND, 2);
  if(dim > 1) mlSol.AddSolution("UY", LAGRANGE, SECOND, 2);
  if(dim > 2) mlSol.AddSolution("UZ", LAGRANGE, SECOND, 2);

  mlSol.Initialize("All");

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc("UX", "Steady");
  if(dim > 1) mlSol.GenerateBdc("UY", "Steady");
  if(dim > 2) mlSol.GenerateBdc("UZ", "Steady");

  MultiLevelProblem ml_prob(&mlSol);

  // ******* Add MPM system to the MultiLevel problem *******
  TransientNonlinearImplicitSystem& system = ml_prob.add_system < TransientNonlinearImplicitSystem > ("UQ");
  system.AddSolutionToSystemPDE("UX");
  if(dim > 1)system.AddSolutionToSystemPDE("UY");
  if(dim > 2) system.AddSolutionToSystemPDE("UZ");

  // ******* System MPM Assembly *******
  system.SetAssembleFunction(AssembleUQSys);
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


  system.MGsolve();

  // ******* Print solution *******
  mlSol.SetWriter(VTK);
  std::vector<std::string> print_vars;
  print_vars.push_back("All");
  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);

  return 0;

} //end main


