
#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "MultiLevelSolution.hpp"

#include "NumericVector.hpp"
#include "adept.h"

#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"
#include "FieldSplitTree.hpp"

#include "slepceps.h"

using namespace femus;

//2D NONLOCAL DOMAIN DECOMPOSITION WITH FETI: nonlocal diffusion using a nonlocal version of FETI

unsigned numberOfUniformLevels = 3;

bool linearMu = true; //decides the sol type of the Lagrange multipliers mu (if false then it is piece-wise constant)

// solver specifics (default is direct solver (MUMPS))
bool Schur = true;

#include "../include/nonlocal_assembly_2D_FETI.hpp"

double InitalValueU (const std::vector < double >& x) {
  double value;

  value = 0.;
//     value =  x[0];
//     value =  x[0] * x[0];
//     value =  x[0] * x[0] * x[0] ;
  return value;
}

void GetL2Norm (MultiLevelSolution & mlSol, MultiLevelSolution & mlSolGlobal);

bool SetBoundaryCondition (const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {

  bool dirichlet = true;
  value = 0.;
//     value = x[0];
//     value = x[0] * x[0];
//     value = x[0] * x[0] * x[0] ;

  return dirichlet;
}

int main (int argc, char** argv) {

  clock_t total_time = clock();

  // init Petsc-MPI communicator
  FemusInit mpinit (argc, argv, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;

  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.ReadCoarseMesh ("../input/FETI_domain.neu", "second", scalingFactor);
//   mlMsh.ReadCoarseMesh ("../input/FETI_domain_small_delta.neu", "second", scalingFactor);
  mlMsh.RefineMesh (numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);

  mlMsh.EraseCoarseLevels (numberOfUniformLevels - 1);

  unsigned dim = mlMsh.GetDimension();

  MultiLevelSolution mlSol (&mlMsh);
  MultiLevelSolution mlSolGlobal (&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution ("u1", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("u2", LAGRANGE, FIRST, 2);
  if (linearMu) mlSol.AddSolution ("mu", LAGRANGE, FIRST, 2);
  else mlSol.AddSolution ("mu", DISCONTINUOUS_POLYNOMIAL, ZERO, 2);

  mlSolGlobal.AddSolution ("u", LAGRANGE, FIRST, 2);

//   mlSol.AddSolution ("u_exact", LAGRANGE, FIRST, 2);

  mlSol.AddSolution ("u1Flag", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("u2Flag", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("muFlag", LAGRANGE, FIRST, 2);

  mlSol.Initialize ("All");
  mlSolGlobal.Initialize ("All");

//   mlSol.Initialize ("u_exact", InitalValueU);

  // ******* Set boundary conditions *******
  mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
  mlSolGlobal.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
  mlSol.GenerateBdc ("All");
  mlSolGlobal.GenerateBdc ("All");

  // ******* Set volume constraints for the nonlocal *******
  std::vector < unsigned > volumeConstraintFlags (3);

  volumeConstraintFlags[0] = 5;
  volumeConstraintFlags[1] = 6;
  volumeConstraintFlags[2] = 7;

  unsigned solu1Index = mlSol.GetIndex ("u1");
  mlSol.GenerateBdcOnVolumeConstraint (volumeConstraintFlags, solu1Index, 0);

  unsigned solu2Index = mlSol.GetIndex ("u2");
  mlSol.GenerateBdcOnVolumeConstraint (volumeConstraintFlags, solu2Index, 0);

  unsigned solmuIndex = mlSol.GetIndex ("mu");
  mlSol.GenerateBdcOnVolumeConstraint (volumeConstraintFlags, solmuIndex, 0);

  unsigned soluIndex = mlSolGlobal.GetIndex ("u");
  mlSolGlobal.GenerateBdcOnVolumeConstraint (volumeConstraintFlags, soluIndex, 0);

  //BEGIN assemble and solve nonlocal FETI problem
  MultiLevelProblem ml_prob (&mlSol);

  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("NonLocal_FETI");
  system.AddSolutionToSystemPDE ("u1");
  system.AddSolutionToSystemPDE ("u2");
  system.AddSolutionToSystemPDE ("mu");

  //BEGIN FIELD SPLIT
  std::vector < unsigned > solutionTypeU1U2 (2);
  solutionTypeU1U2[0] = mlSol.GetSolutionType ("u1");
  solutionTypeU1U2[1] = mlSol.GetSolutionType ("u2");

  std::vector < unsigned > fieldU1U2 (2);
  fieldU1U2[0] = system.GetSolPdeIndex ("u1");
  fieldU1U2[1] = system.GetSolPdeIndex ("u2");
  FieldSplitTree FS_U1U2 (PREONLY, MLU_PRECOND, fieldU1U2, solutionTypeU1U2,  "u1u2");
  //FS_U1U2.SetTolerances (1.e-3, 1.e-20, 1.e+50, 1); // by Guoyi Ke
  //FS_U1U2.PrintFieldSplitTree();


  std::vector < unsigned > solutionTypeMu (1);
  solutionTypeMu[0] = mlSol.GetSolutionType ("mu");
  std::vector < unsigned > fieldMu (1);
  fieldMu[0] = system.GetSolPdeIndex ("mu");
  FieldSplitTree FS_MU (PREONLY, MLU_PRECOND, fieldMu, solutionTypeMu, "mu");
  //FS_MU.SetTolerances (1.e-3, 1.e-20, 1.e+50, 1); // by Guoyi Ke
  //FS_MU.PrintFieldSplitTree();

  std::vector < FieldSplitTree *> FS1;
  FS1.reserve (2);
  FS1.push_back (&FS_U1U2);
  FS1.push_back (&FS_MU);
  FieldSplitTree FS_Nonlocal (PREONLY, FIELDSPLIT_SCHUR_PRECOND, FS1, "Nonlocal_FETI");

  FS_Nonlocal.PrintFieldSplitTree();

  //FS_Nonlocal.SetSchurFactorizationType (SCHUR_FACT_UPPER); // SCHUR_FACT_UPPER, SCHUR_FACT_LOWER,SCHUR_FACT_FULL;
  FS_Nonlocal.SetSchurPreType (SCHUR_PRE_FULL); // SCHUR_PRE_SELF, SCHUR_PRE_SELFP, SCHUR_PRE_USER, SCHUR_PRE_A11, SCHUR_PRE_FULL;
  //END FIELD SPLIT

  // ******* System FEM Assembly *******
  system.SetAssembleFunction (AssembleNonLocalSysFETI);
  system.SetMaxNumberOfLinearIterations (1);
  // ******* set MG-Solver *******
  system.SetMgType (V_CYCLE);

  system.SetAbsoluteLinearConvergenceTolerance (1.e-50);
  //   system.SetNonLinearConvergenceTolerance(1.e-9);
//   system.SetMaxNumberOfNonLinearIterations(20);

  system.SetNumberPreSmoothingStep (1);
  system.SetNumberPostSmoothingStep (1);

  // ******* Set Preconditioner *******
  system.SetLinearEquationSolverType (FEMuS_DEFAULT);
  system.SetOuterSolver (PREONLY);
  if (Schur) {
    system.SetOuterSolver (RICHARDSON);
    system.SetLinearEquationSolverType (FEMuS_FIELDSPLIT, INCLUDE_COARSE_LEVEL_TRUE);
  }
  system.SetSparsityPatternMinimumSize (5000u);   //TODO tune

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids (RICHARDSON);
//   system.SetRichardsonScaleFactor(1.);

  system.SetPreconditionerFineGrids (ILU_PRECOND);

  if (Schur) system.SetFieldSplitTree (&FS_Nonlocal);

  system.SetTolerances (1.e-20, 1.e-20, 1.e+50, 100);

// ******* Solution *******

  system.MGsolve();

  //END assemble and solve nonlocal problem

  //BEGIN assemble and solve global nonlocal problem
  MultiLevelProblem ml_probGlobal (&mlSolGlobal);

  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& systemGlobal = ml_probGlobal.add_system < LinearImplicitSystem > ("NonLocal");
  systemGlobal.AddSolutionToSystemPDE ("u");

  // ******* System FEM Assembly *******
  systemGlobal.SetAssembleFunction (AssembleNonLocalSys);
  systemGlobal.SetMaxNumberOfLinearIterations (1);
  // ******* set MG-Solver *******
  systemGlobal.SetMgType (V_CYCLE);

  systemGlobal.SetAbsoluteLinearConvergenceTolerance (1.e-50);

  systemGlobal.SetNumberPreSmoothingStep (1);
  systemGlobal.SetNumberPostSmoothingStep (1);

  // ******* Set Preconditioner *******
  systemGlobal.SetLinearEquationSolverType (FEMuS_DEFAULT);

  systemGlobal.SetSparsityPatternMinimumSize (5000u);   //TODO tune

  systemGlobal.init();

  // ******* Set Smoother *******
  systemGlobal.SetSolverFineGrids (RICHARDSON);
//   systemGlobal.SetRichardsonScaleFactor(0.7);

  systemGlobal.SetPreconditionerFineGrids (ILU_PRECOND);

  systemGlobal.SetTolerances (1.e-20, 1.e-20, 1.e+50, 100);

// ******* Solution *******

  systemGlobal.MGsolve();

  //END assemble and solve global nonlocal problem


  //BEGIN compute errors
  GetL2Norm (mlSol, mlSolGlobal);
  //END compute errors

  // ******* Print solution *******
  mlSol.SetWriter (VTK);
  std::vector<std::string> print_vars;
  print_vars.push_back ("All");
  mlSol.GetWriter()->SetDebugOutput (true);
  mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "FETI", print_vars, 0);

  mlSolGlobal.SetWriter (VTK);
  std::vector<std::string> print_vars2;
  print_vars2.push_back ("All");
  mlSolGlobal.GetWriter()->SetDebugOutput (true);
  mlSolGlobal.GetWriter()->Write (DEFAULT_OUTPUTDIR, "global", print_vars2, 0);

//   std::cout << std::endl << " total CPU time : " << std::setw (11) << std::setprecision (6) << std::fixed
//             << static_cast<double> ( (clock() - total_time)) / CLOCKS_PER_SEC << " s" << std::endl;

  return 0;

} //end main


void GetL2Norm (MultiLevelSolution & mlSol, MultiLevelSolution & mlSolGlobal) {

  //using the same mesh for FETI and global solutions
  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel (level);

  Solution* sol  = mlSol.GetSolutionLevel (level);
  Solution* solGlobal = mlSolGlobal.GetSolutionLevel (level);

  const unsigned  dim = msh->GetDimension();

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  double L2_error_FETI_global = 0.;

  unsigned solu1Index;
  solu1Index = mlSol.GetIndex ("u1");
  unsigned solu1Type = mlSol.GetSolutionType (solu1Index);

  unsigned solu2Index;
  solu2Index = mlSol.GetIndex ("u2");
  unsigned solu2Type = mlSol.GetSolutionType (solu2Index);

  unsigned soluIndex;
  soluIndex = mlSolGlobal.GetIndex ("u");
  unsigned soluType = mlSolGlobal.GetSolutionType (soluIndex);

  unsigned    iproc = msh->processor_id();
  unsigned    nprocs = msh->n_processors();

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    short unsigned ielGroup = msh->GetElementGroup (iel);
    unsigned nDofu  = msh->GetElementDofNumber (iel, soluType);
    unsigned nDofx = msh->GetElementDofNumber (iel, xType);

    vector < vector < double > > x1 (dim);

    for (int i = 0; i < dim; i++) {
      x1[i].resize (nDofx);
    }

    vector < double >  solu (nDofu);
    vector < double >  solu1 (nDofu);
    vector < double >  solu2 (nDofu);

    for (unsigned i = 0; i < nDofu; i++) {
      unsigned solDofu = msh->GetSolutionDof (i, iel, soluType);
      unsigned solDofu1 = msh->GetSolutionDof (i, iel, solu1Type);
      unsigned solDofu2 = msh->GetSolutionDof (i, iel, solu2Type);
      solu[i] = (*solGlobal->_Sol[soluIndex]) (solDofu);
      solu1[i] = (*sol->_Sol[solu1Index]) (solDofu1);
      solu2[i] = (*sol->_Sol[solu2Index]) (solDofu2);
    }

    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x1[jdim][i] = (*msh->_topology->_Sol[jdim]) (xDof);
      }
    }

    vector <double> phi;
    vector <double> phi_x;
    double weight;

    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielGeom][soluType]->Jacobian (x1, ig, weight, phi, phi_x);

      double solu_gss = 0.;
      double solu1_gss = 0.;
      double solu2_gss = 0.;
      double x_gss = 0.;
      double y_gss = 0.;

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i];
        solu1_gss += phi[i] * solu1[i];
        solu2_gss += phi[i] * solu2[i];
        x_gss += phi[i] * x1[0][i]; // this is x at the Gauss point
        y_gss += phi[i] * x1[1][i]; // this is y at the Gauss point
      }

      double uFETI_gss;
      if (ielGroup == 5 || ielGroup == 6 || ielGroup == 8 || ielGroup == 9) uFETI_gss =  solu1_gss;
      else  uFETI_gss = solu2_gss;

      L2_error_FETI_global += (solu_gss - uFETI_gss) * (solu_gss - uFETI_gss) * weight;

    }
  }

  double norm2 = 0.;
  MPI_Allreduce (&L2_error_FETI_global, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double norm = sqrt (norm2);
  std::cout.precision (16);
  std::cout << "L2 norm of ERROR: Global - FETI = " << norm << std::endl;


}








