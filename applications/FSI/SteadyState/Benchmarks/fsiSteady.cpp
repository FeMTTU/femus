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
#include "CurrentElem.hpp"
#include "../include/FSISteadyStateAssembly.hpp"
#include <dlfcn.h>





using namespace std;
using namespace femus;

void PrintMumpsInfo      (const std::string output_path, const char *stdOutfile, const char* mesh_file, const unsigned &numofrefinements);
void PrintConvergenceInfo(const std::string output_path, const char *stdOutfile, const char* mesh_file, const unsigned &numofrefinements);
void PrintMultigridTime  (const std::string output_path, const char *stdOutfile, const char* mesh_file, const unsigned &numofrefinements);


///@todo I believe I have to review how the computation of the normal is performed
///@todo let the files be read from each run folder




int main(int argc, char **args) {

  // ******* Init Petsc-MPI communicator *******
  FemusInit mpinit(argc,args,MPI_COMM_WORLD);

  // ===============================
  // process options
#define MAX_CHAR_LENGTH 256
  char output_time[MAX_CHAR_LENGTH] = "";
  int dimension=2;
  char mesh_file[MAX_CHAR_LENGTH] = "";
  char stdOutfile[MAX_CHAR_LENGTH] = "";
  char outer_ksp_solver[MAX_CHAR_LENGTH] = "gmres";
  size_t len_mesh_file_name = MAX_CHAR_LENGTH;
  double Lref=1., Uref=1., rhof=1., muf=1., rhos=1., ni=0., E=1.;
  int numofmeshlevels = 1;
  int numofrefinements = 1;
  char bdcfilename[MAX_CHAR_LENGTH] = "";
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
  std::string third_arg_options = "fsiSteady.cpp";
  
  PetscOptionsBegin(PETSC_COMM_WORLD, "", "FSI steady problem options", "Unstructured mesh");

  cout << " Reading flags:" << endl;

  PetscOptionsBool("-mem_infos", "Print memory infos", third_arg_options.c_str(), mem_infos, &mem_infos, NULL);
  printf(" mem_infos: %d\n", mem_infos);

  PetscOptionsInt("-nlevel", "The number of mesh levels", third_arg_options.c_str(), numofmeshlevels , &numofmeshlevels, NULL);
  printf(" nlevel: %i\n", numofmeshlevels);

  PetscOptionsInt("-nrefinement", "The number of refinements", third_arg_options.c_str(), numofrefinements , &numofrefinements, NULL);
  printf(" nrefinement: %i\n", numofrefinements);

  PetscOptionsString("-input", "The name of the input file", third_arg_options.c_str(), "./mesh.neu", mesh_file, len_mesh_file_name, NULL);
  printf(" input: %s\n", mesh_file);

  PetscOptionsString("-std_output", "The name of the redirected standard output file", third_arg_options.c_str(), "", stdOutfile, len_mesh_file_name, NULL);
  printf(" redirected standard output: %s\n", stdOutfile);

  PetscOptionsString("-ic_bdc", "The name of the file with bdc and ic functions", third_arg_options.c_str(), "", bdcfilename, len_mesh_file_name, NULL);
  printf(" ic_bdc: %s\n", bdcfilename);

  PetscOptionsReal("-rhof", "The density of the fluid", third_arg_options.c_str(), rhof, &rhof, NULL);
  printf(" rhof: %f\n", rhof);

  PetscOptionsReal("-rhos", "The density of the solid", third_arg_options.c_str(), rhos, &rhos, NULL);
  printf(" rhos: %f\n", rhos);

  PetscOptionsReal("-E", "The young module of the solid", third_arg_options.c_str(), E, &E, NULL);
  printf(" E: %f\n", E);

  PetscOptionsReal("-muf", "The viscosity of the fluid", third_arg_options.c_str(), muf, &muf, NULL);
  printf(" muf: %f\n", muf);

  PetscOptionsReal("-ni", "The Poisson coefficient of the Solid", third_arg_options.c_str(), ni, &ni, NULL);
  printf(" ni: %f\n", ni);

//   PetscOptionsInt("-nlin_iter", "The number of linear iteration", third_arg_options.c_str(), numlineariter , &numlineariter, NULL);
//   printf(" nlin_iter: %i\n", numlineariter);

  PetscOptionsBool("-equation_pivoting", "Set equation pivoting during assembly", third_arg_options.c_str(), equation_pivoting , &equation_pivoting, NULL);
  printf(" equation_pivoting: %i\n", equation_pivoting);

  PetscOptionsInt("-nnonlin_iter", "The number of non-linear iteration", third_arg_options.c_str(), numnonlineariter, &numnonlineariter, NULL);
  printf(" nnonlin_iter: %i\n", numnonlineariter);

  PetscOptionsReal("-lin_tol", "The linear solver tolerance", third_arg_options.c_str(), lin_tol, &lin_tol, NULL);
  printf(" lin_tol: %g\n", lin_tol);

  PetscOptionsReal("-alin_tol", "The abs linear solver tolerance", third_arg_options.c_str(), alin_tol, &alin_tol, NULL);
  printf(" alin_tol: %g\n", alin_tol);

  PetscOptionsReal("-nonlin_tol", "The nonlinear solver tolerance", third_arg_options.c_str(), nonlin_tol, &nonlin_tol, NULL);
  printf(" nonlin_tol: %g\n", nonlin_tol);

  PetscOptionsInt("-asm_block", "The asm block dimension", third_arg_options.c_str(), asm_block, &asm_block, NULL);
  printf(" asm_block: %i\n", asm_block);

  PetscOptionsString("-outer_ksp_solver", "The outer ksp solver", third_arg_options.c_str(), "gmres", outer_ksp_solver, len_mesh_file_name, NULL);
  printf(" outer_ksp_solver: %s\n", outer_ksp_solver);

  PetscOptionsInt("-npre", "The number of presmoothing step", third_arg_options.c_str(), npre, &npre, NULL);
  printf(" npre: %i\n", npre);

  PetscOptionsInt("-npost", "The number of postmoothing step", third_arg_options.c_str(), npost, &npost, NULL);
  printf(" npost: %i\n", npost);

  PetscOptionsInt("-ksp_restart", "The number of ksp linear step before restarting", third_arg_options.c_str(), ksp_restart, &ksp_restart, NULL);
  printf(" ksp_restart: %i\n", ksp_restart);

  PetscOptionsInt("-max_outer_solver_iter", "The maximum outer solver iterations", third_arg_options.c_str(), max_outer_solver_iter, &max_outer_solver_iter, NULL);
  printf(" max_outer_solver_iter: %i\n", max_outer_solver_iter);

   PetscOptionsString("-output_time", "The name of the redirected standard output file", third_arg_options.c_str(), "", output_time, len_mesh_file_name, NULL);
  printf(" Output time folder: %s\n", output_time);
  
  printf("\n");

 PetscOptionsEnd();
  

  
  // ======= Files ========================
  std::cout << "This is the command line" << std::endl;
  
  for (int i = 0; i < argc; ++i)         std::cout << args[i] << " "; 
  std::cout << std::endl << std::endl;
  
  Files files; 
//         files.CheckIODirectories();
//         files.RedirectCout();
// std::string output_path = files.GetOutputPath();
//   const std::string output_file_to_parse = files.GetOutputTime() + stdOutfile;
// std::string output_path = DEFAULT_OUTPUTDIR;
std::string output_path = output_time;
output_path.append("/");
  const std::string output_file_to_parse = output_path + stdOutfile;
  std::cout << output_file_to_parse << std::endl;

  // ======= Quad Rule ========================
  std::string fe_quad_rule = "fifth";

  

  // *********** loading external functions *******************
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
  // *********** loading external functions - end *******************
  
  

  // ******* Init multilevel mesh from mesh.neu file *******
  const std::string mesh_file_folder = "./input/";
  const std::string mesh_file_path = mesh_file_folder + mesh_file;
  MultiLevelMesh ml_msh(numofrefinements, numofrefinements, mesh_file_path.c_str(), fe_quad_rule.c_str(), Lref, NULL);

  
  ml_msh.EraseCoarseLevels(numofrefinements - numofmeshlevels);

  ml_msh.PrintInfo();

  dimension = ml_msh.GetLevel(0)->GetDimension();

 // ==================================
//   This is done at the beginning, so the mesh is in its original position and I can do a check with if''s on coordinates correctly
 // ==================================
 //How do I handle the computation of the NORMAL for inner edges?
 // I need to make sure that the normal is in the right direction
 // What I'll do is the following: I will restrict the loop over Volume Elements 
 // only to those that belong to Group 6 !
 // ==================================
  //If you want this field to be defined at all levels, this operation has to be performed BEFORE we remove all levels
  //If you want this field to be defined only at the finest level, then do it here
  //here I have to create another structure that only aims at isolating the interface faces inside the mesh
  
  // - either I do it through a check on the coordinates (it has to be on the fixed mesh, the initial one)
  // - or I do a structure on the mesh skeleton that sets a -1 to all faces and then sets a different number to the others

  //it will be similar to elementnearface
  //si potrebbe fare anche svincolata dagli elementi di volume ma supponiamo sia insieme
  //faccio una funzione che non contiene informazione sugli elementi di volume, solo facce.
  
  //indice POSITIVO per le facce NON FLAGGATE (+1) nel mesh file 
  //indice NEGATIVO per le facce FLAGGATE nel mesh file (o di Boundary, o di qualsiasi Interfaccia) (parte da -2 mi sa, per conformita')
  
  //faccio questa struttura con la stessa allocazione iniziale di elementnearface

  std::vector< MyMatrix <int> > _element_faces(numofmeshlevels);  //@todo is this about the faces of each element, only
  //I need to allocate this at all levels! or better, at least at the finest level, which requires the coarser
  //it seems like it is first allocated at the coarse level, and then with the refinement it goes to all levels
  
   allocate_element_faces(_element_faces, ml_msh);

   const std::vector < int >  all_face_flags =  {WET_RIGID,   WET_DEFORMABLE /*DRY_RIGID_DEFORMABLE*/};
   
   const int group_outside_solid_to_which_all_faces_belong = GROUP_VOL_ELEMS; 
   

   fill_element_faces(_element_faces, ml_msh, all_face_flags, group_outside_solid_to_which_all_faces_belong);

   
   
  
  // ******* Fluid and Solid Parameters *******
//    const std::string solid_model = "Neo-Hookean";
   const std::string solid_model = "Mooney-Rivlin";
   
   const std::string fluid_model = "Newtonian";
   
   
  Parameter par(Lref, Uref);

  // Generate Solid Object
  Solid solid;
  solid = Solid(par, E, ni, rhos, solid_model.c_str() );

  cout << "Solid properties: " << endl;
  cout << solid << endl;

  // Generate Fluid Object

  Fluid fluid(par, muf, rhof, fluid_model.c_str() );
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;


  // ******* Init multilevel solution ******
  MultiLevelSolution ml_sol(&ml_msh);

  // ******* Add solution variables to multilevel solution and pair them *******
  ml_sol.AddSolution("DX", LAGRANGE, SECOND, 1);
  ml_sol.AddSolution("DY", LAGRANGE, SECOND, 1);
  if (dimension == 3) ml_sol.AddSolution("DZ", LAGRANGE, SECOND, 1);

  ml_sol.AddSolution("U", LAGRANGE, SECOND, 1);
  ml_sol.AddSolution("V", LAGRANGE, SECOND, 1);
  if (dimension == 3) ml_sol.AddSolution("W", LAGRANGE, SECOND, 1);

  // Pair each velocity variable with the corresponding displacement variable
  ml_sol.PairSolution("U", "DX"); // Add this line
  ml_sol.PairSolution("V", "DY"); // Add this line
  if (dimension == 3) ml_sol.PairSolution("W","DZ"); // Add this line

  
//   const FEFamily pressure_fe_fam   = LAGRANGE;
//   const FEOrder  pressure_fe_order = FIRST;
  const FEFamily pressure_fe_fam   = DISCONTINUOUS_POLYNOMIAL;
  const FEOrder  pressure_fe_order = FIRST;
//   const FEFamily pressure_fe_fam   = DISCONTINUOUS_POLYNOMIAL;
//   const FEOrder  pressure_fe_order = ZERO;
  
  
  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  ml_sol.AddSolution("P", pressure_fe_fam, pressure_fe_order, 1);
  ml_sol.AssociatePropertyToSolution("P", "Pressure", false); // Add this line

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
  
  ml_prob.SetFilesHandler(&files);
  ml_prob.SetQuadratureRuleAllGeomElems(fe_quad_rule);
  ml_prob.set_all_abstract_fe();
  
  
  // Add fluid object
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;
  // Add Solid Object
  ml_prob.parameters.set<Solid>("Solid") = solid;

//     std::cout << ml_prob.parameters.get<Fluid>("Fluid") << std::endl;
//     std::cout << fluid << std::endl;

  
  const bool solve_system = true;
  
  if (solve_system) {
      
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

  
  } //end solve system
  
  
  // ******* Print solution *******
  ml_sol.SetWriter(VTK);

  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  if (dimension==3) mov_vars.push_back("DZ");
  ml_sol.GetWriter()->SetMovingMesh(mov_vars);

  std::vector<std::string> print_vars;
  print_vars.push_back("All");

  ml_sol.GetWriter()->SetDebugOutput(true);
  ml_sol.GetWriter()->Write(output_path, "biquadratic", print_vars);

  if(mem_infos) {
    PetscMemoryGetCurrentUsage(&memory_current_usage);
    PetscPrintf(PETSC_COMM_WORLD, "5: Memory current usage before clear: %g M\n", (double)(memory_current_usage)/(1024.*1024.));
    PetscMemoryGetMaximumUsage(&memory_maximum_usage);
    PetscPrintf(PETSC_COMM_WORLD, "5: Memory maximum usage before clear: %g M\n", (double)(memory_maximum_usage)/(1024.*1024.));
  }

  // ******* Postprocessing *******
  
  Compute_normal_stress_interface(ml_prob, numofmeshlevels - 1, all_face_flags, 0, _element_faces, NULL);

  if(strcmp (output_file_to_parse.c_str(), "") != 0) {
    PrintMumpsInfo      (output_path, output_file_to_parse.c_str(), mesh_file, numofrefinements);
    PrintConvergenceInfo(output_path, output_file_to_parse.c_str(), mesh_file, numofrefinements);
    PrintMultigridTime  (output_path, output_file_to_parse.c_str(), mesh_file, numofrefinements);
  }

  // ******* Clear all systems *******
  ml_prob.clear();

  // ******* close the library *******
  dlclose(handle);
  
  return 0;
}

//these routines are printing another file by processing the run output file
void PrintMumpsInfo(const std::string output_path, const char *stdOutfile, const char* mesh_file, const unsigned &numofrefinements){

  const std::string routine_file_suffix = "_mumps_info.txt";
  
  char output_path_char[MAX_CHAR_LENGTH]; strcpy(output_path_char, output_path.c_str());
  
  std::cout << "END_COMPUTATION\n" << std::flush;

  std::ifstream inf;
  inf.open(stdOutfile);
  if (!inf) {
    std::cout << "Redirected standard output file not found\n";
    std::cout << "add option -std_output std_out_filename > std_out_filename\n";
    return;
  }

  std::ofstream outf;
  char outFileName[MAX_CHAR_LENGTH];
  
  std::ostringstream info_file_stream; info_file_stream << mesh_file << routine_file_suffix;
  sprintf(outFileName, strcat(output_path_char, info_file_stream.str().c_str()));
  
  outf.open(outFileName, std::ofstream::app);
  outf << std::endl << std::endl;
  outf << "Number_of_refinements = " << numofrefinements << std::endl;
  outf << "Nonlinear_Iteration, RINFOG(7), RINFOG(8), RINFOG(9), RINFOG(10), RINFOG(11), INFOG(19), E5*B5*E5/D5+F5*C5*F5/D5";

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



void PrintConvergenceInfo(const std::string output_path, const char *stdOutfile, const char* mesh_file, const unsigned &numofrefinements){

    const std::string routine_file_suffix = "_convergence_info.txt";
    
   char output_path_char[MAX_CHAR_LENGTH]; strcpy(output_path_char, output_path.c_str());

   std::cout << "END_COMPUTATION\n" << std::flush;

  std::ifstream inf;
  inf.open(stdOutfile);
  if (!inf) {
    std::cout<<"Redirected standard output file not found\n";
    std::cout<<"add option -std_output std_out_filename > std_out_filename\n";
    return;
  }

  std::ofstream outf;
  char outFileName[MAX_CHAR_LENGTH];
  
  std::ostringstream info_file_stream; info_file_stream << mesh_file << routine_file_suffix;
  sprintf(outFileName, strcat(output_path_char, info_file_stream.str().c_str()));


  outf.open(outFileName, std::ofstream::app);
  outf << std::endl << std::endl;
  outf << "Number_of_refinements = " << numofrefinements << std::endl;
  outf << "Nonlinear_Iteration, resid_norm0, resid_normN, N, convergence";

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
            outf << "," << normN;
            if(counter != 0){
              outf << "," << counter << "," << pow(normN/norm0,1./counter);
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

void PrintMultigridTime(const std::string output_path, const char *stdOutfile, const char* mesh_file, const unsigned &numofrefinements){

  const std::string routine_file_suffix = "_multigrid_time.txt";
    
  char output_path_char[MAX_CHAR_LENGTH]; strcpy(output_path_char, output_path.c_str());

  std::cout << "END_COMPUTATION\n" << std::flush;

  std::ifstream inf;
  inf.open(stdOutfile);
  if (!inf) {
    std::cout << "Redirected standard output file not found\n";
    std::cout << "add option -std_output std_out_filename > std_out_filename\n";
    return;
  }

  std::ofstream outf;
  char outFileName[MAX_CHAR_LENGTH];

  std::ostringstream info_file_stream; info_file_stream << mesh_file << routine_file_suffix;
  sprintf(outFileName, strcat(output_path_char, info_file_stream.str().c_str()));

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





