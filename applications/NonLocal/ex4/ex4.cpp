
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

#include "../include/nonlocal_assembly_2D_FETI_4domains.hpp"

//2D NONLOCAL DOMAIN DECOMPOSITION WITH FETI AND 4 SUBDOMAINS: nonlocal diffusion using a nonlocal version of FETI

// solver specifics (default is direct solver (MUMPS))
bool Schur = false;

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

unsigned numberOfUniformLevels = 3;

int main (int argc, char** argv) {

  clock_t total_time = clock();

  // init Petsc-MPI communicator
  FemusInit mpinit (argc, argv, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;

  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.ReadCoarseMesh ("../input/FETI_4domains.neu", "second", scalingFactor);
  mlMsh.RefineMesh (numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);

  mlMsh.EraseCoarseLevels (numberOfUniformLevels - 1);

  MultiLevelSolution mlSol (&mlMsh);
  MultiLevelSolution mlSolGlobal (&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution ("u1", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("u2", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("u3", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("u4", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("mu", LAGRANGE, FIRST, 2); //use for all the constraints not on group 13
  mlSol.AddSolution ("muExtra", LAGRANGE, FIRST, 2); // use to impose u2=u3 on group 13
  mlSol.AddSolution ("muExtra2", LAGRANGE, FIRST, 2); // use to impose u3=u4 on group 13
  mlSol.AddSolution ("muExtra3", LAGRANGE, FIRST, 2); // use to impose u1=u2 on group 13

  mlSolGlobal.AddSolution ("u", LAGRANGE, FIRST, 2);

//   mlSol.AddSolution ("u_exact", LAGRANGE, FIRST, 2);

  mlSol.AddSolution ("u1Flag", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("u2Flag", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("u3Flag", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("u4Flag", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("muFlag", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("muExtraFlag", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("muExtra2Flag", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("muExtra3Flag", LAGRANGE, FIRST, 2);

  mlSol.Initialize ("All");
  mlSolGlobal.Initialize ("All");

//   mlSol.Initialize ("u_exact", InitalValueU);

  // ******* Set boundary conditions *******
  mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
  mlSolGlobal.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
  mlSol.GenerateBdc ("All");
  mlSolGlobal.GenerateBdc ("All");

  // ******* Set volume constraints for the nonlocal *******
  std::vector<unsigned> volumeConstraintFlags (8);
  volumeConstraintFlags[0] = 5;
  volumeConstraintFlags[1] = 6;
  volumeConstraintFlags[2] = 7;
  volumeConstraintFlags[3] = 11;
  volumeConstraintFlags[4] = 15;
  volumeConstraintFlags[5] = 16;
  volumeConstraintFlags[6] = 17;
  volumeConstraintFlags[7] = 18;

  unsigned solu1Index = mlSol.GetIndex ("u1");
  mlSol.GenerateBdcOnVolumeConstraint (volumeConstraintFlags, solu1Index, 0);

  unsigned solu2Index = mlSol.GetIndex ("u2");
  mlSol.GenerateBdcOnVolumeConstraint (volumeConstraintFlags, solu2Index, 0);

  unsigned solu3Index = mlSol.GetIndex ("u3");
  mlSol.GenerateBdcOnVolumeConstraint (volumeConstraintFlags, solu3Index, 0);

  unsigned solu4Index = mlSol.GetIndex ("u4");
  mlSol.GenerateBdcOnVolumeConstraint (volumeConstraintFlags, solu4Index, 0);

  unsigned solmuIndex = mlSol.GetIndex ("mu");
  mlSol.GenerateBdcOnVolumeConstraint (volumeConstraintFlags, solmuIndex, 0);

  unsigned solmuExtraIndex = mlSol.GetIndex ("muExtra");
  mlSol.GenerateBdcOnVolumeConstraint (volumeConstraintFlags, solmuExtraIndex, 0);

  unsigned solmuExtra2Index = mlSol.GetIndex ("muExtra2");
  mlSol.GenerateBdcOnVolumeConstraint (volumeConstraintFlags, solmuExtra2Index, 0);

  unsigned solmuExtra3Index = mlSol.GetIndex ("muExtra3");
  mlSol.GenerateBdcOnVolumeConstraint (volumeConstraintFlags, solmuExtra3Index, 0);

  unsigned soluIndex = mlSolGlobal.GetIndex ("u");
  mlSolGlobal.GenerateBdcOnVolumeConstraint (volumeConstraintFlags, soluIndex, 0);

  //BEGIN assemble and solve nonlocal FETI problem
  MultiLevelProblem ml_prob (&mlSol);

  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("NonLocal_FETI");
  system.AddSolutionToSystemPDE ("u1");
  system.AddSolutionToSystemPDE ("u2");
  system.AddSolutionToSystemPDE ("u3");
  system.AddSolutionToSystemPDE ("u4");
  system.AddSolutionToSystemPDE ("mu");
  system.AddSolutionToSystemPDE ("muExtra");
  system.AddSolutionToSystemPDE ("muExtra2");
  system.AddSolutionToSystemPDE ("muExtra3");

  //BEGIN FIELD SPLIT
  std::vector < unsigned > solutionTypeUjs (4);
  solutionTypeUjs[0] = mlSol.GetSolutionType ("u1");
  solutionTypeUjs[1] = mlSol.GetSolutionType ("u2");
  solutionTypeUjs[2] = mlSol.GetSolutionType ("u3");
  solutionTypeUjs[3] = mlSol.GetSolutionType ("u4");

  std::vector < unsigned > fieldUjs (4);
  fieldUjs[0] = system.GetSolPdeIndex ("u1");
  fieldUjs[1] = system.GetSolPdeIndex ("u2");
  fieldUjs[2] = system.GetSolPdeIndex ("u3");
  fieldUjs[3] = system.GetSolPdeIndex ("u4");
  FieldSplitTree FS_Ujs (PREONLY, MLU_PRECOND, fieldUjs, solutionTypeUjs,  "u1u2u3u4");
  //FS_Ujs.SetTolerances (1.e-3, 1.e-20, 1.e+50, 1); // by Guoyi Ke
  //FS_Ujs.PrintFieldSplitTree();


  std::vector < unsigned > solutionTypeMus (4);
  solutionTypeMus[0] = mlSol.GetSolutionType ("mu");
  solutionTypeMus[1] = mlSol.GetSolutionType ("muExtra");
  solutionTypeMus[2] = mlSol.GetSolutionType ("muExtra2");
  solutionTypeMus[3] = mlSol.GetSolutionType ("muExtra3");
  std::vector < unsigned > fieldMus (4);
  fieldMus[0] = system.GetSolPdeIndex ("mu");
  fieldMus[1] = system.GetSolPdeIndex ("muExtra");
  fieldMus[2] = system.GetSolPdeIndex ("muExtra2");
  fieldMus[3] = system.GetSolPdeIndex ("muExtra3");
  FieldSplitTree FS_Mus (PREONLY, MLU_PRECOND, fieldMus, solutionTypeMus, "mus");
  //FS_Mus.SetTolerances (1.e-3, 1.e-20, 1.e+50, 1); // by Guoyi Ke
  //FS_Mus.PrintFieldSplitTree();

  std::vector < FieldSplitTree *> FS1;
  FS1.reserve (2);
  FS1.push_back (&FS_Ujs);
  FS1.push_back (&FS_Mus);
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
//   system.SetRichardsonScaleFactor(0.7);

  system.SetPreconditionerFineGrids (ILU_PRECOND);

  if (Schur) system.SetFieldSplitTree (&FS_Nonlocal);

  system.SetTolerances (1.e-20, 1.e-20, 1.e+50, 100);

// ******* Solution *******

  system.MGsolve();

  //END assemble and solve nonlocal FETI problem


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
  //   systemGlobal.SetNonLinearConvergenceTolerance(1.e-9);
//   systemGlobal.SetMaxNumberOfNonLinearIterations(20);

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

  unsigned solu3Index;
  solu3Index = mlSol.GetIndex ("u3");
  unsigned solu3Type = mlSol.GetSolutionType (solu3Index);

  unsigned solu4Index;
  solu4Index = mlSol.GetIndex ("u4");
  unsigned solu4Type = mlSol.GetSolutionType (solu4Index);

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
    vector < double >  solu3 (nDofu);
    vector < double >  solu4 (nDofu);

    for (unsigned i = 0; i < nDofu; i++) {
      unsigned solDofu = msh->GetSolutionDof (i, iel, soluType);
      unsigned solDofu1 = msh->GetSolutionDof (i, iel, solu1Type);
      unsigned solDofu2 = msh->GetSolutionDof (i, iel, solu2Type);
      unsigned solDofu3 = msh->GetSolutionDof (i, iel, solu3Type);
      unsigned solDofu4 = msh->GetSolutionDof (i, iel, solu4Type);
      solu[i] = (*solGlobal->_Sol[soluIndex]) (solDofu);
      solu1[i] = (*sol->_Sol[solu1Index]) (solDofu1);
      solu2[i] = (*sol->_Sol[solu2Index]) (solDofu2);
      solu3[i] = (*sol->_Sol[solu3Index]) (solDofu3);
      solu4[i] = (*sol->_Sol[solu4Index]) (solDofu4);
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
      double solu3_gss = 0.;
      double solu4_gss = 0.;
      double x_gss = 0.;
      double y_gss = 0.;

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i];
        solu1_gss += phi[i] * solu1[i];
        solu2_gss += phi[i] * solu2[i];
        solu3_gss += phi[i] * solu3[i];
        solu4_gss += phi[i] * solu4[i];
        x_gss += phi[i] * x1[0][i]; // this is x at the Gauss point
        y_gss += phi[i] * x1[1][i]; // this is y at the Gauss point
      }

      double uFETI_gss;
      if (ielGroup == 5 || ielGroup == 6 || ielGroup == 8 || ielGroup == 9 || ielGroup == 11 || ielGroup == 12 || ielGroup == 13) {
        uFETI_gss =  solu1_gss;
      }
      else if (ielGroup == 7 || ielGroup == 10 || ielGroup == 14 || ielGroup == 15) {
        uFETI_gss = solu2_gss;
      }
      else if (ielGroup == 16 || ielGroup == 17 || ielGroup == 19 || ielGroup == 20) {
        uFETI_gss = solu3_gss;
      }
      else uFETI_gss = solu4_gss;

      L2_error_FETI_global += (solu_gss - uFETI_gss) * (solu_gss - uFETI_gss) * weight;

    }
  }

  double norm2 = 0.;
  MPI_Allreduce (&L2_error_FETI_global, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double norm = sqrt (norm2);
  std::cout.precision (16);
  std::cout << "L2 norm of ERROR: Global - FETI = " << norm << std::endl;


}







