#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"
#include "FemusInit.hpp"
#include "SparseMatrix.hpp"
#include "FElemTypeEnum.hpp"
#include "Files.hpp"
#include "TransientSystem.hpp"
#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "../include/FSITimeDependentAssembly.hpp"
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
  size_t len_infile_name = 256;
  char infile[256] = "";
  char stdOutfile[256] = "";
  char outer_ksp_solver[256] = "gmres";
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
  int ksp_restart = 30;
  int n_timesteps = 1;
  double time_step = 0.01;
  char restart_file_name[256] = "";
  int autosave_time_interval = 1;
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

  PetscOptionsBool("-mem_infos", "Print memory infos", "fsiTimeDependent.cpp", mem_infos, &mem_infos, NULL);
  printf(" mem_infos: %d\n", mem_infos);

  PetscOptionsInt("-n_timesteps", "The number of time steps", "fsiTimeDependent.cpp", n_timesteps, &n_timesteps, NULL);
  printf(" n_timesteps: %d\n", n_timesteps);

  PetscOptionsReal("-time_step", "The time step", "fsiTimeDependent.cpp", time_step, &time_step, NULL);
  printf(" time_step: %f\n", time_step);

  PetscOptionsString("-restart_file_name", "The name of the file for restart", "fsiTimeDependent.cpp", "", restart_file_name, len_infile_name, NULL);
  printf(" restart_file_name: %s\n", restart_file_name);

  PetscOptionsInt("-autosave_time_interval", "The autosave interval time", "fsiTimeDependent.cpp", autosave_time_interval, &autosave_time_interval, NULL);
  printf(" autosave_time_interval: %d\n", autosave_time_interval);

  PetscOptionsInt("-nlevel", "The number of mesh levels", "fsiTimeDependent.cpp", numofmeshlevels , &numofmeshlevels, NULL);
  printf(" nlevel: %i\n", numofmeshlevels);

  PetscOptionsInt("-nrefinement", "The number of refinements", "fsiTimeDependent.cpp", numofrefinements , &numofrefinements, NULL);
  printf(" nrefinement: %i\n", numofrefinements);

  PetscOptionsString("-input", "The name of the input file", "fsiTimeDependent.cpp", "./mesh.neu", infile, len_infile_name, NULL);
  printf(" input: %s\n", infile);

  PetscOptionsString("-std_output", "The name of the redirected standard output file", "fsiTimeDependent.cpp", "", stdOutfile, len_infile_name, NULL);
  printf(" redirected standard output: %s\n", stdOutfile);

  PetscOptionsString("-ic_bdc", "The name of the file with bdc and ic functions", "fsiTimeDependent.cpp", "", bdcfilename, len_infile_name, NULL);
  printf(" ic_bdc: %s\n", bdcfilename);

  PetscOptionsReal("-rhof", "The density of the fluid", "fsiTimeDependent.cpp", rhof, &rhof, NULL);
  printf(" rhof: %f\n", rhof);

  PetscOptionsReal("-rhos", "The density of the solid", "fsiTimeDependent.cpp", rhos, &rhos, NULL);
  printf(" rhos: %f\n", rhos);

  PetscOptionsReal("-E", "The young module of the solid", "fsiTimeDependent.cpp", E, &E, NULL);
  printf(" E: %f\n", E);

  PetscOptionsReal("-muf", "The viscosity of the fluid", "fsiTimeDependent.cpp", muf, &muf, NULL);
  printf(" muf: %f\n", muf);

  PetscOptionsReal("-ni", "The Poisson coefficient of the Solid", "fsiTimeDependent.cpp", ni, &ni, NULL);
  printf(" ni: %f\n", ni);

//   PetscOptionsInt("-nlin_iter", "The number of linear iteration", "fsiTimeDependent.cpp", numlineariter , &numlineariter, NULL);
//   printf(" nlin_iter: %i\n", numlineariter);

//   PetscOptionsBool("-equation_pivoting", "Set equation pivoting during assembly", "fsiTimeDependent.cpp", equation_pivoting , &equation_pivoting, NULL);
//   printf(" equation_pivoting: %i\n", equation_pivoting);

  PetscOptionsInt("-nnonlin_iter", "The number of non-linear iteration", "fsiTimeDependent.cpp", numnonlineariter, &numnonlineariter, NULL);
  printf(" nnonlin_iter: %i\n", numnonlineariter);

  PetscOptionsReal("-lin_tol", "The linear solver tolerance", "fsiTimeDependent.cpp", lin_tol, &lin_tol, NULL);
  printf(" lin_tol: %g\n", lin_tol);

  PetscOptionsReal("-alin_tol", "The abs linear solver tolerance", "fsiTimeDependent.cpp", alin_tol, &alin_tol, NULL);
  printf(" alin_tol: %g\n", alin_tol);

  PetscOptionsReal("-nonlin_tol", "The nonlinear solver tolerance", "fsiTimeDependent.cpp", nonlin_tol, &nonlin_tol, NULL);
  printf(" nonlin_tol: %g\n", nonlin_tol);

  PetscOptionsInt("-asm_block", "The asm block dimension", "fsiTimeDependent.cpp", asm_block, &asm_block, NULL);
  printf(" asm_block: %i\n", asm_block);

  PetscOptionsString("-outer_ksp_solver", "The outer ksp solver", "fsiTimeDependent.cpp", "gmres", outer_ksp_solver, len_infile_name, NULL);
  printf(" outer_ksp_solver: %s\n", outer_ksp_solver);

  PetscOptionsInt("-npre", "The number of presmoothing step", "fsiTimeDependent.cpp", npre, &npre, NULL);
  printf(" npre: %i\n", npre);

  PetscOptionsInt("-npost", "The number of postmoothing step", "fsiTimeDependent.cpp", npost, &npost, NULL);
  printf(" npost: %i\n", npost);

  PetscOptionsInt("-ksp_restart", "The number of ksp linear step before restarting", "fsiTimeDependent.cpp", ksp_restart, &ksp_restart, NULL);
  printf(" ksp_restart: %i\n", ksp_restart);

  PetscOptionsInt("-max_outer_solver_iter", "The maximum outer solver iterations", "fsiTimeDependent.cpp", max_outer_solver_iter, &max_outer_solver_iter, NULL);
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
  const char *dlsym_error = dlerror();

//   cout << " Loading symbol InitalValueU...";
//   typedef double (*InitalValueU_t)(const std::vector < double >& x);

// reset errors
//   dlerror();
//   InitalValueU_t InitalValueU = (InitalValueU_t) dlsym(handle, "InitalValueU");
//   if (dlsym_error) {
//       cerr << "Cannot load symbol 'InitalValueU': " << dlsym_error << '\n';
//       dlclose(handle);
//       return 1;
//   }
//   else {
//     cout << " done" << endl;
//   }

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

  cout << " Loading symbol TimeStepFunction...";
  typedef double (*TimeStepFunction_t)(const double time);

  // reset errors
  dlerror();
  TimeStepFunction_t TimeStepFunction = (TimeStepFunction_t) dlsym(handle, "TimeStepFunction");
  dlsym_error = dlerror();
  if (dlsym_error) {
      cerr << "Cannot load symbol 'TimeStepFunction': " << dlsym_error << '\n';
      cout << " skip" << endl;
      TimeStepFunction = NULL;
//       return 1;
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
  ml_sol.AddSolution("DX",LAGRANGE,SECOND,2);
  ml_sol.AddSolution("DY",LAGRANGE,SECOND,2);
  if (dimension==3) ml_sol.AddSolution("DZ",LAGRANGE,SECOND,2);

  ml_sol.AddSolution("U",LAGRANGE,SECOND,2);
  ml_sol.AddSolution("V",LAGRANGE,SECOND,2);
  if (dimension==3) ml_sol.AddSolution("W",LAGRANGE,SECOND,2);

  // Pair each velocity variable with the corresponding displacement variable
  ml_sol.PairSolution("U","DX"); // Add this line
  ml_sol.PairSolution("V","DY"); // Add this line
  if (dimension==3) ml_sol.PairSolution("W","DZ"); // Add this line

  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  ml_sol.AddSolution("P",DISCONTINUOUS_POLYNOMIAL,FIRST,2);
  ml_sol.AssociatePropertyToSolution("P","Pressure",false); // Add this line

  // ******* Initialize solution *******
  ml_sol.Initialize("All");

  // ******* Set boundary functions *******
  ml_sol.AttachSetBoundaryConditionFunction(BdcFunction);


  // ******* Set boundary conditions *******
  ml_sol.GenerateBdc("DX","Steady");
  ml_sol.GenerateBdc("DY","Steady");
  if (dimension==3) ml_sol.GenerateBdc("DZ","Steady");
  // TODO: time dep or steady cannot be hardcoded
  ml_sol.GenerateBdc("U","Time_dependent");
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
  TransientMonolithicFSINonlinearImplicitSystem & system = ml_prob.add_system<TransientMonolithicFSINonlinearImplicitSystem> ("Fluid-Structure-Interaction");
  system.AddSolutionToSystemPDE("DX");
  system.AddSolutionToSystemPDE("DY");
  if (dimension==3) system.AddSolutionToSystemPDE("DZ");
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  if (dimension==3) system.AddSolutionToSystemPDE("W");
  system.AddSolutionToSystemPDE("P");

  // ******* System Fluid-Structure-Interaction Assembly *******
  system.SetAssembleFunction(FSITimeDependentAssembly);


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
  system.SetPreconditionerFineGrids(MLU_PRECOND);
  // Set block size for the ASM smoother
  system.SetElementBlockNumber(asm_block);
  // Set Solver of the smoother (name to be changed)
  system.SetSolverFineGrids(RICHARDSON);

  // set the tolerances for the GMRES outer solver
  system.SetTolerances(lin_tol,alin_tol,div_tol,max_outer_solver_iter,ksp_restart);

  system.SetOuterSolver(GMRES);


  // ******* Add variables to be solved *******
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");

  // ******* Set the last (1) variables in system (i.e. P) to be a schur variable *******
  system.SetNumberOfSchurVariables(1);

  // time loop parameter
  if(TimeStepFunction) {
    system.AttachGetTimeIntervalFunction(TimeStepFunction);
  } else
  {
    system.SetIntervalTime(time_step);
  }

  // TODO cannot be hardcoded
  if(strcmp (restart_file_name,"") != 0) {
    ml_sol.LoadSolution(restart_file_name);
    std::string s = restart_file_name;
    std::string delimiter = "_time";
    size_t pos = 0;
    pos = s.find(delimiter);
    std::cout << "substring: " << s.substr(pos + delimiter.length(), s.length()) << std::endl;
    double time = atof(s.substr(pos + delimiter.length(), s.length()).c_str());
    system.SetTime(time);
  }

  // ******* Solve *******
  std::cout << std::endl;

  std::cout << " *********** Solving... ************  " << std::endl;

  system.PrintSolverInfo(true);

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



  for (unsigned i_time_step = 0; i_time_step < n_timesteps; i_time_step++) {

    if( i_time_step > 0 || strcmp (restart_file_name,"") != 0 )
      system.SetMgType(V_CYCLE);

    system.CopySolutionToOldSolution();

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

    //system.UpdateSolution();

    if( (i_time_step+1)%autosave_time_interval == 0) {
      ml_sol.SaveSolution("run", system.GetTime());
      std::cout << " it: " << i_time_step + 1 << " store save solution for restart at time " << system.GetTime() << std::endl;
    }

    ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR,"biquadratic",print_vars, i_time_step+1);
  }


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
  if(strcmp (infile,"./input/turek_FSI2.neu") == 0){
    sprintf(outFileName, "turek_hron_FSI2_mumps_info.txt");
  }
  else if(strcmp (infile,"./input/turek_FSI3.neu") == 0){
    sprintf(outFileName, "turek_hron_FSI3_mumps_info.txt");
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
  outf << "Simulation_Time,Nonlinear_Iteration,RINFOG(7),RINFOG(8),RINFOG(9),RINFOG(10),RINFOG(11),INFOG(19),F5*C5*F5/E5+G5*D5*G5/E5";

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
        outf << std::endl <<simulationTime<<","<<str1;
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
  if(strcmp (infile,"./input/turek_FSI2.neu") == 0){
    sprintf(outFileName, "turek_hron_FSI2_convergence_info.txt");
  }
  else if(strcmp (infile,"./input/turek_FSI3.neu") == 0){
    sprintf(outFileName, "turek_hron_FSI3_convergence_info.txt");
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

  if(strcmp (infile,"./input/turek_FSI2.neu") == 0){
    sprintf(outFileName, "turek_hron_FSI2_multigrid_time.txt");
  }
  else if(strcmp (infile,"./input/turek_FSI3.neu") == 0){
    sprintf(outFileName, "turek_hron_FSI3_multigrid_info.txt");
  }
  else if(strcmp (infile,"./input/richter3d.neu") == 0){
    sprintf(outFileName, "richter3d_multigrid_time.txt");
  }
  else{
    sprintf(outFileName, "generic_multigrid_time.txt");
  }

  outf.open(outFileName, std::ofstream::app);
  outf << std::endl;
  outf << "\nLevel_Max,Simulation_Time,Average_Time";

  int counter = 0;
  double ave_lin_solver_time = 0.;

  std::string str1;
  inf >> str1;
  double simulationTime = 0.;
  while (str1.compare("END_COMPUTATION") != 0) {

    if (str1.compare("Simulation") == 0){
      inf >> str1;
      if (str1.compare("Time:") == 0){
        inf >> simulationTime;
      }
    }
    else if (str1.compare("Start") == 0){
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
          outf << simulationTime <<","<<ave_lin_solver_time / counter;
        }
      }
    }
    inf >> str1;
  }

//   outf << ave_lin_solver_time / counter;

  outf.close();
  inf.close();

}

