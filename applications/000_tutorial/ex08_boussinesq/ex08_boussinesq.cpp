/** \file Ex7.cpp
 *  \brief This example shows how to set and solve the weak form
 *   of the Boussinesq approximation of the Navier-Stokes Equation
 *
 *  \f{eqnarray*}
 *  && \mathbf{V} \cdot \nabla T - \nabla \cdot\alpha \nabla T = 0 \\
 *  && \mathbf{V} \cdot \nabla \mathbf{V} - \nabla \cdot \nu (\nabla \mathbf{V} +(\nabla \mathbf{V})^T)
 *  +\nabla P = \beta T \mathbf{j} \\
 *  && \nabla \cdot \mathbf{V} = 0
 *  \f}
 *  in a unit box domain (in 2D and 3D) with given temperature 0 and 1 on
 *  the left and right walls, respectively, and insulated walls elsewhere.
 *  \author Eugenio Aulisa
 */

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"


#include "04_boussinesq.hpp"


using namespace femus;

bool SetBoundaryCondition(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;

  if (!strcmp(SolName, "T")) {
    if (facename == 2) {
      value = 1.;
    } else if (facename == 3) {
      dirichlet = false; //Neumann
    }
  } else if (!strcmp(SolName, "P")) {
    dirichlet = false;
  }

  return dirichlet;
}




int main(int argc, char** args) {

  
  boussinesq_equation::Pr = 1.;
  boussinesq_equation::Ra = 1000.;
  
  
  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;

  const std::string relative_path_to_build_directory =  "../../../";
  const std::string mesh_file = relative_path_to_build_directory + Files::mesh_folder_path() + "01_gambit/02_2d/square/minus0p5-plus0p5_minus0p5-plus0p5/square_2x2_quad_Three_boundary_groups.neu";
  mlMsh.ReadCoarseMesh(mesh_file.c_str(), "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 7;
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
  mlSol.AddSolution("T", LAGRANGE, FIRST);
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);

  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND);

  //mlSol.AddSolution("P", LAGRANGE, FIRST);
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
  system.AddSolutionToSystemPDE("T");
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");

  if (dim == 3) system.AddSolutionToSystemPDE("W");

  system.AddSolutionToSystemPDE("P");

  // system.SetLinearEquationSolverType(FEMuS_ASM); // GMRES with ADDITIVE SWRARTZ METHOD (domain decomposition)
  system.SetLinearEquationSolverType(FEMuS_DEFAULT); // GMRES
  // attach the assembling function to system
  system.SetAssembleFunction( femus::boussinesq_equation::AssembleBoussinesqApproximation_AD );

  system.SetMaxNumberOfNonLinearIterations(20);
  system.SetNonLinearConvergenceTolerance(1.e-8);

  system.SetMaxNumberOfLinearIterations(3);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-12);

  system.SetMgType(F_CYCLE);
  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);


  // initilaize and solve the system
  system.init();

  system.SetSolverFineGrids(GMRES);

  system.SetPreconditionerFineGrids(ILU_PRECOND);

  system.SetTolerances(1.e-20, 1.e-20, 1.e+50, 40); //PETSC GMRES tolerances

  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");

  system.SetNumberOfSchurVariables(1); // option for FEMuS_ASM fot "P"
  system.SetElementBlockNumber(4); // option for FEMuS_ASM

  system.MGsolve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.Write(Files::_application_output_directory, fe_fams_for_files[ FILES_CONTINUOUS_BIQUADRATIC ], variablesToBePrinted);

  return 0;
}

