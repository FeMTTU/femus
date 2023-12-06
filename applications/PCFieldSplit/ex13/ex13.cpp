
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"



#include "03_navier_stokes.hpp"



using namespace femus;

bool SetBoundaryCondition(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
    value = 0.;
  if (!strcmp(SolName, "V")) {
    if (facename == 2) {
      if(x[1]>-0.5 && x[1]<0.5){
	value = 1.;
  } 
    }
  } else if (!strcmp(SolName, "P")) {
    dirichlet = false;
    value = 0.;
  }

  return dirichlet;
}



int main(int argc, char** args) {



  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  const std::string relative_path_to_build_directory =  "../../../";
  const std::string mesh_file = relative_path_to_build_directory + Files::mesh_folder_path() + "01_gambit/2d/square/minus0p5-plus0p5_minus0p5-plus0p5/square_2x2_quad_Three_boundary_groups.neu";

  mlMsh.ReadCoarseMesh (mesh_file.c_str(), "seventh", scalingFactor );
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 6;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  //mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);
  
  
  
  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);
  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND);
  
  //  mlSol.AddSolution("P", LAGRANGE, FIRST);
  mlSol.AddSolution("P",  DISCONTINUOUS_POLYNOMIAL, FIRST);
 
  mlSol.AssociatePropertyToSolution("P", "Pressure");
  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.FixSolutionAtOnePoint("P");
  
  mlSol.GenerateBdc("All", "Steady", & mlProb);


  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("NS");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  system.AddSolutionToSystemPDE("P");

  if (dim == 3) system.AddSolutionToSystemPDE("W");

// BEGIN 
  std::vector < unsigned > solutionTypeUV(2);
  solutionTypeUV[0] = mlSol.GetSolutionType("U");
  solutionTypeUV[1] = mlSol.GetSolutionType("V");

  std::vector < unsigned > fieldUV(2);
  fieldUV[0] = system.GetSolPdeIndex("U");
  fieldUV[1] = system.GetSolPdeIndex("V");
  FieldSplitTree FS_UV( PREONLY, ASM_PRECOND, fieldUV, solutionTypeUV,  "Velocity");
  FS_UV.SetAsmBlockSize(4);
  FS_UV.SetAsmNumeberOfSchurVariables(0);
// END

// BEGIN 
  std::vector < unsigned > solutionTypeP(1);
  solutionTypeP[0] = mlSol.GetSolutionType("P");
  std::vector < unsigned > fieldP(1);
  fieldP[0] = system.GetSolPdeIndex("P");
  //FieldSplitTree FS_P(GMRES, LSC_PRECOND, fieldP, "Pressure");// It works, but it is slower than ILU_PRECOND
  FieldSplitTree FS_P(PREONLY, ASM_PRECOND, fieldP, solutionTypeP, "Pressure");// It works, but it is slower that Vanka-ASM
  FS_P.SetAsmBlockSize(4);
  FS_P.SetAsmNumeberOfSchurVariables(0);
// END

  
// BEGIN 
  std::vector < FieldSplitTree *> FS1;
  FS1.reserve(2);
  FS1.push_back(&FS_UV);
  FS1.push_back(&FS_P);
// END


// BEGIN 
  FieldSplitTree FS_NS(RICHARDSON, FIELDSPLIT_SCHUR_PRECOND, FS1, "Navier-Stokes");
// END


  // system.SetLinearEquationSolverType(FEMuS_DEFAULT);
  //system.SetLinearEquationSolverType(FEMuS_ASM); // Additive Swartz Method
  system.SetLinearEquationSolverType(FEMuS_FIELDSPLIT);

  // attach the assembling function to system
  system.SetAssembleFunction( femus::AssembleNavierStokes_AD );

  system.SetMaxNumberOfNonLinearIterations(10);
  system.SetNonLinearConvergenceTolerance(1.e-8);
  system.SetMaxNumberOfResidualUpdatesForNonlinearIteration(10);
  system.SetResidualUpdateConvergenceTolerance(1.e-12);

  system.SetMgType(F_CYCLE);

  system.SetNumberPreSmoothingStep(0);
  system.SetNumberPostSmoothingStep(2);
  // initilaize and solve the system
  system.init();

  system.SetSolverFineGrids(RICHARDSON);
  system.SetPreconditionerFineGrids(ILU_PRECOND);

  system.SetTolerances(1.e-10, 1.e-20, 1.e+50, 20, 5);

  system.SetFieldSplitTree(&FS_NS);


  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");
  system.SetNumberOfSchurVariables(1);
  system.SetElementBlockNumber(4);

  system.MGsolve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(Files::_application_output_directory, fe_fams_for_files[ FILES_CONTINUOUS_BIQUADRATIC ], variablesToBePrinted);

  return 0;
}




