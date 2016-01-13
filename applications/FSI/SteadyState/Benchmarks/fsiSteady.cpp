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
#include "../include/FSISteadyStateAssembly.hpp"
#include <dlfcn.h>

using namespace std;
using namespace femus;


int main(int argc,char **args) {

  // ******* Init Petsc-MPI communicator *******
  FemusInit mpinit(argc,args,MPI_COMM_WORLD);

  // process options
  int dimension=2;
  char infile[256] = "";
  char outer_ksp_solver[256] = "gmres";
  size_t len_infile_name = 256;
  double Lref=1., Uref=1., rhof=1., muf=1., rhos=1., ni=0., E=1.;
  int numofmeshlevels = 1;
  int numofrefinements = 1;
  std::string gauss_integration_order = "fifth";
  char bdcfilename[256] = "";
  int numlineariter = 1;
  int numnonlineariter = 15;
  double lin_tol = 1.e-04;
  double alin_tol = 1.e-20;
  double div_tol = 1.e+10;
  double nonlin_tol = 1.e-08;
  int asm_block = 2;
  int npre = 0;
  int npost = 2;
  int max_outer_solver_iter = 40;
  PetscBool equation_pivoting = PETSC_TRUE;

  // ******* reading input parameters *******
  PetscOptionsBegin(PETSC_COMM_WORLD, "", "FSI steady problem options", "Unstructured mesh");

  cout << " Reading flags:" << endl;

  PetscOptionsInt("-nlevel", "The number of mesh levels", "fsiSteady.cpp", numofmeshlevels , &numofmeshlevels, NULL);
  printf(" nlevel: %i\n", numofmeshlevels);

  PetscOptionsInt("-nrefinement", "The number of refinements", "fsiSteady.cpp", numofrefinements , &numofrefinements, NULL);
  printf(" nrefinement: %i\n", numofrefinements);

  PetscOptionsString("-input", "The name of the input file", "fsiSteady.cpp", "./mesh.neu", infile, len_infile_name, NULL);
  printf(" input: %s\n", infile);

  PetscOptionsString("-ic_bdc", "The name of the file with bdc and ic functions", "fsiSteady.cpp", "", bdcfilename, len_infile_name, NULL);
  printf(" ic_bdc: %s\n", bdcfilename);

  PetscOptionsReal("-rhof", "The density of the fluid", "fsiSteady.cpp", rhof, &rhof, NULL);
  printf(" rhof: %f\n", rhof);

  PetscOptionsReal("-rhof", "The density of the solid", "fsiSteady.cpp", rhos, &rhos, NULL);
  printf(" rhos: %f\n", rhos);

  PetscOptionsReal("-E", "The young module of the solid", "fsiSteady.cpp", E, &E, NULL);
  printf(" E: %f\n", E);

  PetscOptionsReal("-muf", "The viscosity of the fluid", "fsiSteady.cpp", muf, &muf, NULL);
  printf(" muf: %f\n", muf);

  PetscOptionsReal("-ni", "The Poisson coefficient of the Solid", "fsiSteady.cpp", ni, &ni, NULL);
  printf(" ni: %f\n", ni);

//   PetscOptionsInt("-nlin_iter", "The number of linear iteration", "fsiSteady.cpp", numlineariter , &numlineariter, NULL);
//   printf(" nlin_iter: %i\n", numlineariter);

  PetscOptionsBool("-equation_pivoting", "Set equation pivoting during assembly", "fsiSteady.cpp", equation_pivoting , &equation_pivoting, NULL);
  printf(" equation_pivoting: %i\n", equation_pivoting);

  PetscOptionsInt("-nnonlin_iter", "The number of non-linear iteration", "fsiSteady.cpp", numnonlineariter, &numnonlineariter, NULL);
  printf(" nnonlin_iter: %i\n", numnonlineariter);

  PetscOptionsReal("-lin_tol", "The linear solver tolerance", "fsiSteady.cpp", lin_tol, &lin_tol, NULL);
  printf(" lin_tol: %g\n", lin_tol);

  PetscOptionsReal("-alin_tol", "The abs linear solver tolerance", "fsiSteady.cpp", alin_tol, &alin_tol, NULL);
  printf(" alin_tol: %g\n", alin_tol);

  PetscOptionsReal("-nonlin_tol", "The nonlinear solver tolerance", "fsiSteady.cpp", nonlin_tol, &nonlin_tol, NULL);
  printf(" nonlin_tol: %g\n", nonlin_tol);

  PetscOptionsInt("-asm_block", "The asm block dimension", "fsiSteady.cpp", asm_block, &asm_block, NULL);
  printf(" asm_block: %i\n", asm_block);

  PetscOptionsString("-outer_ksp_solver", "The outer ksp solver", "fsiSteady.cpp", "gmres", outer_ksp_solver, len_infile_name, NULL);
  printf(" outer_ksp_solver: %s\n", outer_ksp_solver);

  PetscOptionsInt("-max_outer_solver_iter", "The maximum outer solver iterations", "fsiSteady.cpp", max_outer_solver_iter, &max_outer_solver_iter, NULL);
  printf(" max_outer_solver_iter: %i\n", max_outer_solver_iter);

  printf("\n");

  PetscOptionsEnd();

  // loading external functions
  cout << " Loading lib bdcfilename...";
  void *handle = dlopen (bdcfilename, RTLD_LAZY);
  if (!handle) {
    cerr << "Cannot open library: " << dlerror() << '\n';
    return 1;
  }
  else {
    cout << " done" << endl;
  }

  // load the symbols
  cout << " Loading symbol InitalValueU...";
  typedef double (*InitalValueU_t)(const std::vector < double >& x);

  // reset errors
  dlerror();
  InitalValueU_t InitalValueU = (InitalValueU_t) dlsym(handle, "InitalValueU");
  const char *dlsym_error = dlerror();
  if (dlsym_error) {
      cerr << "Cannot load symbol 'InitalValueU': " << dlsym_error << '\n';
      dlclose(handle);
      return 1;
  }
  else {
    cout << " done" << endl;
  }

  cout << " Loading symbol BdcFunction...";
  typedef bool (*BdcFunction_t)(const std::vector < double >& x,const char name[], double &value, const int facename, const double time);

  // reset errors
  dlerror();
  BdcFunction_t BdcFunction = (BdcFunction_t) dlsym(handle, "BdcFunction");
  dlsym_error = dlerror();
  if (dlsym_error) {
      cerr << "Cannot load symbol 'BdcFunction': " << dlsym_error << '\n';
      dlclose(handle);
      return 1;
  }
  else {
    cout << " done" << endl;
  }

  // ******* Init multilevel mesh from mesh.neu file *******

  MultiLevelMesh ml_msh(numofrefinements, numofrefinements, infile, gauss_integration_order.c_str(), Lref, NULL);

  ml_msh.EraseCoarseLevels(numofrefinements - numofmeshlevels);

  ml_msh.PrintInfo();
  
  dimension = ml_msh.GetLevel(0)->GetDimension();

  // mark Solid nodes
  ml_msh.MarkStructureNode();


  // ******* Fluid and Solid Parameters *******
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


  // ******* Init multilevel solution ******
  MultiLevelSolution ml_sol(&ml_msh);

  // ******* Add solution variables to multilevel solution and pair them *******
  ml_sol.AddSolution("DX",LAGRANGE,SECOND,1);
  ml_sol.AddSolution("DY",LAGRANGE,SECOND,1);
  if (dimension==3) ml_sol.AddSolution("DZ",LAGRANGE,SECOND,1);

  ml_sol.AddSolution("U",LAGRANGE,SECOND,1);
  ml_sol.AddSolution("V",LAGRANGE,SECOND,1);
  if (dimension==3) ml_sol.AddSolution("W",LAGRANGE,SECOND,1);

  // Pair each velocity variable with the corresponding displacement variable
  ml_sol.PairSolution("U","DX"); // Add this line
  ml_sol.PairSolution("V","DY"); // Add this line
  if (dimension==3) ml_sol.PairSolution("W","DZ"); // Add this line

  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  ml_sol.AddSolution("P",DISCONTINOUS_POLYNOMIAL,FIRST,1);
  ml_sol.AssociatePropertyToSolution("P","Pressure",false); // Add this line

  // ******* Initialize solution *******
  ml_sol.Initialize("All");
  ml_sol.Initialize("U",InitalValueU); //only for turek problem

  // ******* Set boundary functions *******
  ml_sol.AttachSetBoundaryConditionFunction(BdcFunction);


  // ******* Set boundary conditions *******
  ml_sol.GenerateBdc("DX","Steady");
  ml_sol.GenerateBdc("DY","Steady");
  if (dimension==3) ml_sol.GenerateBdc("DZ","Steady");
  ml_sol.GenerateBdc("U","Steady");
  ml_sol.GenerateBdc("V","Steady");
  if (dimension==3) ml_sol.GenerateBdc("W","Steady");
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
  if (dimension==3) system.AddSolutionToSystemPDE("DZ");
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  if (dimension==3) system.AddSolutionToSystemPDE("W");
  system.AddSolutionToSystemPDE("P");

  // ******* System Fluid-Structure-Interaction Assembly *******
  if ( equation_pivoting == PETSC_TRUE){
    system.SetAssembleFunction(FSISteadyStateAssembly);
  }
  else {
    system.SetAssembleFunction(FSISteadyStateAssemblyWithNoPivoting);
  }


  // Solver settings
  // ******* set Non-linear MG-Solver *******
  system.SetMaxNumberOfLinearIterations(numlineariter);
  system.SetAbsoluteLinearConvergenceTolerance(lin_tol);

  system.SetMaxNumberOfNonLinearIterations(numnonlineariter);
  system.SetNonLinearConvergenceTolerance(nonlin_tol);

  // Set type of multigrid preconditioner (V-cycle, F-cycle)
  system.SetMgType(F_CYCLE);
  system.SetNumberPreSmoothingStep(npre);
  system.SetNumberPostSmoothingStep(npost);

  // ******* Set Smoother *******
  // Set Preconditioner of the smoother (name to be changed)
  system.SetMgSmoother(ASM_SMOOTHER);

  // System init
  system.init();

  // Set the preconditioner for each ASM block
  system.SetPreconditionerFineGrids(ILU_PRECOND);
  // Set block size for the ASM smoother
  system.SetElementBlockNumber(asm_block);
  // Set Solver of the smoother (name to be changed)
  system.SetSolverFineGrids(RICHARDSON);

  // set the tolerances for the GMRES outer solver
  system.SetTolerances(lin_tol,alin_tol,div_tol,max_outer_solver_iter);
  system.SetOuterKSPSolver(outer_ksp_solver);

  // ******* Add variables to be solved *******
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");

  // ******* Set the last (1) variables in system (i.e. P) to be a schur variable *******
  system.SetNumberOfSchurVariables(1);


  // ******* Solve *******
  std::cout << std::endl;

  std::cout << " *********** Solving... ************  " << std::endl;

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

  ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR,"biquadratic",print_vars);

  // ******* Clear all systems *******
  ml_prob.clear();

  // close the library
  dlclose(handle);

  return 0;
}

