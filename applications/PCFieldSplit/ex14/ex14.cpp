/** \file Ex7.cpp
 *  \brief This example shows how to set and solve the weak form
 *   of the Boussinesq appoximation of the Navier-Stokes Equation
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
#include "NumericVector.hpp"
#include "SparseMatrix.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "LinearEquationSolver.hpp"
#include "FieldSplitTree.hpp"

#include "adept.h"


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
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  system.AddSolutionToSystemPDE("P");

  if (dim == 3) system.AddSolutionToSystemPDE("W");
 
  std::vector < unsigned > fieldUV(2);
  fieldUV[0] = system.GetSolPdeIndex("U");
  fieldUV[1] = system.GetSolPdeIndex("V");
  FieldSplitTree FS_UV( PREONLY, ILU_PRECOND, fieldUV , "Velocity");
  FS_UV.SetTolerances(1.e-3,1.e-20,1.e+50, 1); // changed by Guoyi Ke 

  std::vector < unsigned > fieldP(1);
  fieldP[0] = system.GetSolPdeIndex("P");

  //FS_P.SetFieldSplitSchurFactType{PC_FIELDSPLIT_SCHUR_FACT_LOWER}; //changed by Guoyi Ke
  FieldSplitTree FS_P(PREONLY, ILU_PRECOND, fieldP, "Pressure");
  FS_P.SetTolerances(1.e-3,1.e-20,1.e+50, 1); //changed by Guoyi Ke
  
  std::vector < FieldSplitTree *> FS1;
  FS1.reserve(2);
  FS1.push_back(&FS_UV);
  FS1.push_back(&FS_P);
  
  FieldSplitTree FS_NS(GMRES, FIELDSPLIT_SCHUR_PRECOND, FS1, "Navier-Stokes");
  FS_NS.SetSchurFactorizationType(SCHUR_FACT_UPPER); // SCHUR_FACT_UPPER, SCHUR_FACT_LOWER,SCHUR_FACT_FULL; how to use if FS_SCHUR_PRECOND? Guoyike
  FS_NS.SetSchurPreType(SCHUR_PRE_SELFP);// SCHUR_PRE_SELF, SCHUR_PRE_SELFP, SCHUR_PRE_USER, SCHUR_PRE_A11,SCHUR_PRE_FULL;

  //system.SetLinearEquationSolverType(FEMuS_DEFAULT);
  //system.SetLinearEquationSolverType(FEMuS_ASM); // Additive Swartz Method
  system.SetLinearEquationSolverType(FEMuS_FIELDSPLIT); // Additive Swartz Method

  // attach the assembling function to system
  system.SetAssembleFunction( femus::AssembleNavierStokes_AD );

  system.SetMaxNumberOfNonLinearIterations(20);
  system.SetMaxNumberOfLinearIterations(3);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-12);
  system.SetNonLinearConvergenceTolerance(1.e-8);
  system.SetMgType(F_CYCLE);

  system.SetNumberPreSmoothingStep(2);
  system.SetNumberPostSmoothingStep(2);
  // initilaize and solve the system
  system.init();

  system.SetSolverFineGrids(GMRES);
  system.SetPreconditionerFineGrids(ILU_PRECOND);
  system.SetFieldSplitTree(&FS_NS);

  system.SetTolerances(1.e-3, 1.e-20, 1.e+50, 5);


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
  vtkIO.Write(Files::_application_output_directory, "biquadratic", variablesToBePrinted);

  return 0;
}


