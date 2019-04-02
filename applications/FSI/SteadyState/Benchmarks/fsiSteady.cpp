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

void PrintMumpsInfo(char *stdOutfile, char* infile, const unsigned &numofrefinements);
void PrintConvergenceInfo(char *stdOutfile, char* infile, const unsigned &numofrefinements);
void PrintMultigridTime(char *stdOutfile, char* infile, const unsigned &numofrefinements);

int main(int argc,char **args) {

  // ******* Init Petsc-MPI communicator *******
  FemusInit mpinit(argc,args,MPI_COMM_WORLD);

  // process options
  int dimension=2;
  char infile[256] = "";
  char stdOutfile[256] = "";
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
  int npre = 4;
  int npost = 4;
  int max_outer_solver_iter = 40;
  int ksp_restart = 10;
  PetscBool equation_pivoting = PETSC_TRUE;
  PetscBool mem_infos = PETSC_FALSE;
  PetscLogDouble memory_current_usage, memory_maximum_usage;

  PetscMemorySetGetMaximumUsage();
  if(mem_infos) {
    PetscMemoryGetCurrentUsage(&memory_current_usage);
    PetscPrintf(PETSC_COMM_WORLD, "0: Memory current usage at beginning: %g M\n", (double)(memory_current_usage)/(1024.*1024.));
  }

  // ******* reading input parameters *******
  PetscOptionsBegin(PETSC_COMM_WORLD, "", "FSI steady problem options", "Unstructured mesh");

  cout << " Reading flags:" << endl;

  PetscOptionsBool("-mem_infos", "Print memory infos", "fsiSteady.cpp", mem_infos, &mem_infos, NULL);
  printf(" mem_infos: %d\n", mem_infos);

  PetscOptionsInt("-nlevel", "The number of mesh levels", "fsiSteady.cpp", numofmeshlevels , &numofmeshlevels, NULL);
  printf(" nlevel: %i\n", numofmeshlevels);

  PetscOptionsInt("-nrefinement", "The number of refinements", "fsiSteady.cpp", numofrefinements , &numofrefinements, NULL);
  printf(" nrefinement: %i\n", numofrefinements);

  PetscOptionsString("-input", "The name of the input file", "fsiSteady.cpp", "./mesh.neu", infile, len_infile_name, NULL);
  printf(" input: %s\n", infile);

  PetscOptionsString("-std_output", "The name of the redirected standard output file", "fsiSteady.cpp", "", stdOutfile, len_infile_name, NULL);
  printf(" redirected standard output: %s\n", stdOutfile);

  PetscOptionsString("-ic_bdc", "The name of the file with bdc and ic functions", "fsiSteady.cpp", "", bdcfilename, len_infile_name, NULL);
  printf(" ic_bdc: %s\n", bdcfilename);

  PetscOptionsReal("-rhof", "The density of the fluid", "fsiSteady.cpp", rhof, &rhof, NULL);
  printf(" rhof: %f\n", rhof);

  PetscOptionsReal("-rhos", "The density of the solid", "fsiSteady.cpp", rhos, &rhos, NULL);
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

  PetscOptionsInt("-npre", "The number of presmoothing step", "fsiSteady.cpp", npre, &npre, NULL);
  printf(" npre: %i\n", npre);

  PetscOptionsInt("-npost", "The number of postmoothing step", "fsiSteady.cpp", npost, &npost, NULL);
  printf(" npost: %i\n", npost);

  PetscOptionsInt("-ksp_restart", "The number of ksp linear step before restarting", "fsiSteady.cpp", ksp_restart, &ksp_restart, NULL);
  printf(" ksp_restart: %i\n", ksp_restart);

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
  ml_sol.AddSolution("P",DISCONTINUOUS_POLYNOMIAL,FIRST,1);
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
  system.SetLinearEquationSolverType(FEMuS_ASM);

  if(mem_infos) {
    PetscMemoryGetCurrentUsage(&memory_current_usage);
    PetscPrintf(PETSC_COMM_WORLD, "1: Memory current usage before system init: %g M\n", (double)(memory_current_usage)/(1024.*1024.));
  }

  // System init
  system.init();

  if(mem_infos) {
    PetscMemoryGetCurrentUsage(&memory_current_usage);
    PetscPrintf(PETSC_COMM_WORLD, "2: Memory current usage after system init: %g M\n", (double)(memory_current_usage)/(1024.*1024.));
  }

  // Set the preconditioner for each ASM block
  system.SetPreconditionerFineGrids(ILU_PRECOND);
  // Set block size for the ASM smoother
  system.SetElementBlockNumber(asm_block);
  // Set Solver of the smoother (name to be changed)
  system.SetSolverFineGrids(RICHARDSON);

  // set the tolerances for the GMRES outer solver
  system.SetTolerances(lin_tol,alin_tol,div_tol,max_outer_solver_iter,ksp_restart);
  system.SetOuterSolver(GMRES);//system.SetOuterKSPSolver(outer_ksp_solver);


  // ******* Add variables to be solved *******
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");

  // ******* Set the last (1) variables in system (i.e. P) to be a schur variable *******
  system.SetNumberOfSchurVariables(1);


  // ******* Solve *******
  std::cout << std::endl;

  std::cout << " *********** Solving... ************  " << std::endl;

  system.PrintSolverInfo(true);

  if(mem_infos) {
    PetscMemoryGetCurrentUsage(&memory_current_usage);
    PetscPrintf(PETSC_COMM_WORLD, "3: Memory current usage before solve: %g M\n", (double)(memory_current_usage)/(1024.*1024.));
  }

  system.MGsolve();

  if(mem_infos) {
    PetscMemoryGetCurrentUsage(&memory_current_usage);
    PetscPrintf(PETSC_COMM_WORLD, "4: Memory current usage after solve: %g M\n", (double)(memory_current_usage)/(1024.*1024.));
    PetscMemoryGetMaximumUsage(&memory_maximum_usage);
    PetscPrintf(PETSC_COMM_WORLD, "4: Memory maximum usage after solve: %g M\n", (double)(memory_maximum_usage)/(1024.*1024.));
  }

  // ******* Print solution *******
  ml_sol.SetWriter(VTK);


  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  mov_vars.push_back("DZ");
  ml_sol.GetWriter()->SetMovingMesh(mov_vars);

  std::vector<std::string> print_vars;
  print_vars.push_back("All");

  ml_sol.GetWriter()->SetDebugOutput(true);
  ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR,"biquadratic",print_vars);

  if(mem_infos) {
    PetscMemoryGetCurrentUsage(&memory_current_usage);
    PetscPrintf(PETSC_COMM_WORLD, "5: Memory current usage before clear: %g M\n", (double)(memory_current_usage)/(1024.*1024.));
    PetscMemoryGetMaximumUsage(&memory_maximum_usage);
    PetscPrintf(PETSC_COMM_WORLD, "5: Memory maximum usage before clear: %g M\n", (double)(memory_maximum_usage)/(1024.*1024.));
  }

  // ******* Clear all systems *******
  ml_prob.clear();

  // close the library
  dlclose(handle);
  if(strcmp (stdOutfile,"") != 0){
    PrintMumpsInfo(stdOutfile, infile, numofrefinements);
    PrintConvergenceInfo(stdOutfile, infile, numofrefinements);
    PrintMultigridTime(stdOutfile, infile, numofrefinements);
  }

  return 0;
}

void PrintMumpsInfo(char *stdOutfile, char* infile, const unsigned &numofrefinements){

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
  if(strcmp (infile,"./input/turek_FSI1.neu") == 0){
    sprintf(outFileName, "turek_hron_FSI1_mumps_info.txt");
  }
  else if(strcmp (infile,"./input/richter3d.neu") == 0){
    sprintf(outFileName, "richter3d_mumps_info.txt");
  }
  else{
    sprintf(outFileName, "generic_mumps_info.txt");
  }

  outf.open(outFileName, std::ofstream::app);
  outf << std::endl << std::endl;
  outf << "Number_of_refinements="<<numofrefinements<<std::endl;
  outf << "Nonlinear_Iteration,RINFOG(7),RINFOG(8),RINFOG(9),RINFOG(10),RINFOG(11),INFOG(19),E5*B5*E5/D5+F5*C5*F5/D5";

  std::string str1;
  inf >> str1;
  while (str1.compare("END_COMPUTATION") != 0) {
    if (str1.compare("Nonlinear") == 0) {
      inf >> str1;
      if (str1.compare("iteration") == 0) {
        inf >> str1;
        outf << std::endl << str1;
      }
    }
    else if (str1.compare("RINFOG(7),RINFOG(8)") == 0) {
      inf >> str1 >> str1 >> str1 >> str1;
      outf <<","<< str1;
      inf >> str1;
      outf << str1;
    }
    else if (str1.compare("RINFOG(9)") == 0) {
      inf >> str1 >> str1 >> str1;
      outf <<","<< str1;
    }
    else if (str1.compare("RINFOG(10),RINFOG(11)(condition") == 0) {
      inf >> str1>> str1;
      outf <<","<< str1;
      inf >> str1;
      outf << str1;
    }
    else if (str1.compare("INFOG(19)") == 0){
      inf >> str1 >> str1 >> str1>> str1 >> str1 >> str1;
      inf >> str1 >> str1 >> str1>> str1 >> str1 >> str1 >> str1;
      inf >> str1;
      outf <<","<< str1;
    }
    inf >> str1;
  }

  outf.close();
  inf.close();

};



void PrintConvergenceInfo(char *stdOutfile, char* infile, const unsigned &numofrefinements){

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
  if(strcmp (infile,"./input/turek_FSI1.neu") == 0){
    sprintf(outFileName, "turek_hron_FSI1_convergence_info.txt");
  }
  else if(strcmp (infile,"./input/richter3d.neu") == 0){
    sprintf(outFileName, "richter3d_convergence_info.txt");
  }
  else{
    sprintf(outFileName, "generic_convergence_info.txt");
  }

  outf.open(outFileName, std::ofstream::app);
  outf << std::endl << std::endl;
  outf << "Number_of_refinements="<<numofrefinements<<std::endl;
  outf << "Nonlinear_Iteration,resid_norm0,resid_normN,N,convergence";

  std::string str1;
  inf >> str1;
  while (str1.compare("END_COMPUTATION") != 0) {

    if (str1.compare("Nonlinear") == 0) {
      inf >> str1;
      if (str1.compare("iteration") == 0) {
        inf >> str1;
        outf << std::endl << str1;
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

void PrintMultigridTime(char *stdOutfile, char* infile, const unsigned &numofrefinements){

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

  if(strcmp (infile,"./input/turek_FSI1.neu") == 0){
    sprintf(outFileName, "turek_hron_FSI1_multigrid_time.txt");
  }
  else if(strcmp (infile,"./input/richter3d.neu") == 0){
    sprintf(outFileName, "richter3d_multigrid_time.txt");
  }
  else{
    sprintf(outFileName, "generic_multigrid_time.txt");
  }

  outf.open(outFileName, std::ofstream::app);
  outf << std::endl;
  outf << "\nLevel_Max,Average_Time";

  int counter = 0;
  double ave_lin_solver_time = 0.;

  std::string str1;
  inf >> str1;
  while (str1.compare("END_COMPUTATION") != 0) {

    if (str1.compare("Start") == 0){
      inf >> str1;
      if (str1.compare("Level") == 0){
        inf >> str1;
        if (str1.compare("Max") == 0){
          inf >> str1;
          outf <<"\n"<< str1 <<",";
          counter = 0;
          ave_lin_solver_time = 0.;
        }
      }
    }
    else if (str1.compare("MG") == 0) {
      inf >> str1;
      if (str1.compare("linear") == 0) {
        inf >> str1;
        if (str1.compare("solver") == 0) {
	  inf >> str1;
          if (str1.compare("time:") == 0) {
            double value;
            inf >> value;
	    ++counter;
	    ave_lin_solver_time += value;
	  }
        }
      }
    }
    if (str1.compare("End") == 0){
      inf >> str1;
      if (str1.compare("Level") == 0){
        inf >> str1;
        if (str1.compare("Max") == 0){
          outf << ave_lin_solver_time / counter;
        }
      }
    }
    inf >> str1;
  }

//   outf << ave_lin_solver_time / counter;

  outf.close();
  inf.close();

}

