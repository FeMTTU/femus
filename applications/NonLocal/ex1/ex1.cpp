
#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"

#include "NumericVector.hpp"
#include "adept.h"

#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

#include "slepceps.h"

#include "../include/nonlocal_assembly.hpp"


//1D NONLOCAL EX : nonlocal diffusion for a body with different material properties

using namespace femus;

double InitalValueU (const std::vector < double >& x) {

  double u1 = a1 + b1 * x[0] - 1. / (2. * kappa1) * x[0] * x[0];
  double u2 = a2 + b2 * x[0] - 1. / (2. * kappa2) * x[0] * x[0];

  double value = (x[0] < 0.) ? u1 : u2;

//     return x[0] + 0. * ( 0.51 * 0.51 - x[0] * x[0] ) * ( 0.51 * 0.51 - x[1] * x[1] );
//     return x[0];
//     return x[0] * x[0];
//     return x[0] * x[0] * x[0];
//     return x[0] * x[0] * x[0] * x[0] + 0.1 * x[0] * x[0]; //this is x^4 + delta x^2
//     return x[0] * x[0] * x[0] * x[0]; //this is x^4
//     return 2 * x[0] + x[0] * x[0] * x[0] * x[0] * x[0]; //this is 2x + x^5

  return value;
}

void GetL2Norm (MultiLevelSolution &mlSol, MultiLevelSolution &mlSolFine);

void GetL2NormLocalConnectedNonlocalFine (MultiLevelSolution & mlSol, MultiLevelSolution & mlSolFine);

void PutADoubleNodeAtTheInterface (MultiLevelMesh &mlMsh, const double &meshSize, double &leftBound, double &rightBound);

void ShiftTheExtrema (MultiLevelMesh &mlMsh, const double &meshSize, const double &delta1Shift, const double &delta2Shift, double &leftBound, double &rightBound);

void BuildElementSkipFlags (MultiLevelMesh &mlMsh, std::vector<unsigned> &elementSkipFlags);

bool SetBoundaryCondition (const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {

  bool dirichlet = true;

  double u1 = a1 + b1 * x[0] - 1. / (2. * kappa1) * x[0] * x[0];
  double u2 = a2 + b2 * x[0] - 1. / (2. * kappa2) * x[0] * x[0];

  value = (x[0] < 0.) ? u1 : u2;

//     value = 0.;
//     value = x[0];
//     value = x[0] * x[0];
//     value = x[0] * x[0] * x[0];
//     value = x[0] * x[0] * x[0] * x[0] + 0.1 * x[0] * x[0]; //this is x^4 + delta x^2
//     value = x[0] * x[0] * x[0] * x[0];
//     value =  2 * x[0] + x[0] * x[0] * x[0] * x[0] * x[0]; //this is 2x + x^5

  return dirichlet;
}

unsigned numberOfUniformLevels = 1;

int main (int argc, char** argv) {

  // init Petsc-MPI communicator
  FemusInit mpinit (argc, argv, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;

  //coarse soln mesh
  double xMinCoarseBox = leftBound - delta1Mesh;
  double xMaxCoarseBox = rightBound + delta2Mesh;

  mlMsh.GenerateCoarseBoxMesh (numberOfElements, 0, 0, xMinCoarseBox, xMaxCoarseBox, 0., 0., 0., 0., EDGE3, "fifth");

  if (doubleIntefaceNode) PutADoubleNodeAtTheInterface (mlMsh, desiredMeshSize, leftBound, rightBound);

  if (shiftExternalNodes) ShiftTheExtrema (mlMsh, desiredMeshSize, delta1Shift, delta2Shift, leftBound, rightBound);

//   mlMsh.RefineMesh (numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);

  BuildElementSkipFlags (mlMsh, elementSkipFlags);

  mlMsh.EraseCoarseLevels (numberOfUniformLevels - 1);
//     numberOfUniformLevels = 1;

  //fine soln mesh
  MultiLevelMesh mlMshFine;

  double xMinCoarseBoxFine = leftBoundFine - delta1Mesh;
  double xMaxCoarseBoxFine = rightBoundFine + delta2Mesh;

  mlMshFine.GenerateCoarseBoxMesh (numberOfElementsFine, 0, 0, xMinCoarseBoxFine, xMaxCoarseBoxFine, 0., 0., 0., 0., EDGE3, "fifth");

  if (doubleIntefaceNodeFine) PutADoubleNodeAtTheInterface (mlMshFine, desiredMeshSizeFine, leftBoundFine, rightBoundFine);

  if (shiftExternalNodes) ShiftTheExtrema (mlMshFine, desiredMeshSizeFine, delta1ShiftFine, delta2ShiftFine, leftBoundFine, rightBoundFine);

//   mlMshFine.RefineMesh (numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);

  BuildElementSkipFlags (mlMshFine, elementSkipFlagsFine);

  mlMshFine.EraseCoarseLevels (numberOfUniformLevels - 1);
//     numberOfUniformLevels = 1;

  unsigned dim = mlMsh.GetDimension();

  MultiLevelSolution mlSol (&mlMsh);
  MultiLevelSolution mlSolFine (&mlMshFine);

  // add variables to mlSol
  mlSol.AddSolution ("u", LAGRANGE, FIRST, 2);
  mlSolFine.AddSolution ("u_fine", LAGRANGE, FIRST, 2);

  mlSol.AddSolution ("u_local", LAGRANGE, FIRST, 2);

  mlSol.AddSolution ("u_exact", LAGRANGE, FIRST, 2);

  mlSol.Initialize ("All");
  mlSolFine.Initialize ("All");

  mlSol.Initialize ("u_exact", InitalValueU);

  mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
  mlSolFine.AttachSetBoundaryConditionFunction (SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc ("All");
  mlSolFine.GenerateBdc ("All");


  //BEGIN assemble and solve nonlocal problem
  MultiLevelProblem ml_prob (&mlSol);

  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("NonLocal");
  system.AddSolutionToSystemPDE ("u");

  // ******* System FEM Assembly *******
  system.SetAssembleFunction (AssembleNonLocalSys);
  system.SetMaxNumberOfLinearIterations (1);
  // ******* set MG-Solver *******
  system.SetMgType (V_CYCLE);

  system.SetAbsoluteLinearConvergenceTolerance (1.e-50);
  //   system.SetNonLinearConvergenceTolerance(1.e-9);
//   system.SetMaxNumberOfNonLinearIterations(20);

  system.SetNumberPreSmoothingStep (1);
  system.SetNumberPostSmoothingStep (1);

  system.SetLinearEquationSolverType ( FEMuS_DEFAULT );

  system.SetSparsityPatternMinimumSize (10000u);   //TODO tune was 500

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids (GMRES);

  system.SetPreconditionerFineGrids (ILU_PRECOND);

  system.SetTolerances (1.e-20, 1.e-20, 1.e+50, 100);

// ******* Solution *******

  system.MGsolve(); //TODO

  //END assemble and solve nonlocal problem

  //BEGIN assemble and solve local problem
  MultiLevelProblem ml_prob2 (&mlSol);

  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& system2 = ml_prob2.add_system < LinearImplicitSystem > ("Local");
  system2.AddSolutionToSystemPDE ("u_local");

  // ******* System FEM Assembly *******
  system2.SetAssembleFunction (AssembleLocalSys);
  system2.SetMaxNumberOfLinearIterations (1);
  // ******* set MG-Solver *******
  system2.SetMgType (V_CYCLE);

  system2.SetAbsoluteLinearConvergenceTolerance (1.e-50);

  system2.SetNumberPreSmoothingStep (1);
  system2.SetNumberPostSmoothingStep (1);

  system2.SetLinearEquationSolverType ( FEMuS_DEFAULT );

  system2.init();

  // ******* Set Smoother *******
  system2.SetSolverFineGrids (GMRES);

  system2.SetPreconditionerFineGrids (ILU_PRECOND);

  system2.SetTolerances (1.e-20, 1.e-20, 1.e+50, 100);

// ******* Solution *******

  system2.MGsolve();

  //END assemble and solve local problem

  //BEGIN assemble and solve fine nonlocal problem
  MultiLevelProblem ml_probFine (&mlSolFine);

  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& systemFine = ml_probFine.add_system < LinearImplicitSystem > ("NonLocalFine");
  systemFine.AddSolutionToSystemPDE ("u_fine");

  // ******* System FEM Assembly *******
  systemFine.SetAssembleFunction (AssembleNonLocalSysFine);
  systemFine.SetMaxNumberOfLinearIterations (1);
  // ******* set MG-Solver *******
  systemFine.SetMgType (V_CYCLE);

  systemFine.SetAbsoluteLinearConvergenceTolerance (1.e-50);

  systemFine.SetNumberPreSmoothingStep (1);
  systemFine.SetNumberPostSmoothingStep (1);

  // ******* Set Preconditioner *******
  systemFine.SetLinearEquationSolverType ( FEMuS_DEFAULT );

  systemFine.SetSparsityPatternMinimumSize (10000u);   //TODO tune

  systemFine.init();

  // ******* Set Smoother *******
  systemFine.SetSolverFineGrids (GMRES);

  systemFine.SetPreconditionerFineGrids (ILU_PRECOND);

  systemFine.SetTolerances (1.e-20, 1.e-20, 1.e+50, 100);

// ******* Solution *******

//   systemFine.MGsolve();  //TODO uncomment

  //END assemble and solve fine nonlocal problem


  clock_t L2norm_time = clock();

  //BEGIN compute errors
  GetL2Norm (mlSol, mlSolFine);
  GetL2NormLocalConnectedNonlocalFine (mlSol, mlSolFine);
  //END compute errors

  std::cout << std::endl << " L2 norm CPU time : " << std::setw (11) << std::setprecision (6) << std::fixed
            << static_cast<double> ( (clock() - L2norm_time)) / CLOCKS_PER_SEC << " s" << std::endl;

  // ******* Print solution *******
  mlSol.SetWriter (VTK);
  std::vector<std::string> print_vars;
  print_vars.push_back ("All");
  mlSol.GetWriter()->SetDebugOutput (true);
  mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "nonlocal_local_exact", print_vars, 0);

  mlSolFine.SetWriter (VTK);
  std::vector<std::string> print_vars2;
  print_vars2.push_back ("All");
  mlSolFine.GetWriter()->SetDebugOutput (true);
  mlSolFine.GetWriter()->Write (DEFAULT_OUTPUTDIR, "fine", print_vars2, 0);

  std::cout.precision (16);
  std::cout << "Mesh size h = " << (xMaxCoarseBox - xMinCoarseBox) / (numberOfElements * pow (2, numberOfUniformLevels - 1)) << std::endl;
  std::cout << "Mesh size fine h = " << (xMaxCoarseBoxFine - xMinCoarseBoxFine) / (numberOfElementsFine * pow (2, numberOfUniformLevels - 1)) << std::endl;

  return 0;

} //end main


void GetL2Norm (MultiLevelSolution & mlSol, MultiLevelSolution & mlSolFine) {

  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel (level);
  Solution* sol  = mlSol.GetSolutionLevel (level);

  const unsigned levelFine = mlSolFine._mlMesh->GetNumberOfLevels() - 1;
  Mesh* mshFine = mlSolFine._mlMesh->GetLevel (levelFine);
  Solution* solFine  = mlSolFine.GetSolutionLevel (levelFine);

  const unsigned  dim = msh->GetDimension();

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  double error_solExact_norm2 = 0.;

  double error_solExact_local_norm2 = 0.;

  double error_solLocal_norm2 = 0.;

  double error_NonLocCoarse_NonLocFine_norm2 = 0.;

  double solNonlocal_norm2 = 0.;

  double solNonlocalFine_norm2 = 0.;

  double solLocal_norm2 = 0.;

  double sol_exact_norm2 = 0.;

  unsigned soluIndex;
  soluIndex = mlSol.GetIndex ("u");
  unsigned soluType = mlSol.GetSolutionType (soluIndex);

  unsigned soluIndexLocal;
  soluIndexLocal = mlSol.GetIndex ("u_local");

  unsigned soluIndexFine;
  soluIndexFine = mlSolFine.GetIndex ("u_fine");

  unsigned    iproc = msh->processor_id();
  unsigned    nprocs = msh->n_processors();

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    if (elementSkipFlags[iel] == 0) {

      short unsigned ielGeom = msh->GetElementType (iel);
      unsigned nDofu  = msh->GetElementDofNumber (iel, soluType);

      vector < vector < double > > x1 (dim);

      for (int i = 0; i < dim; i++) {
        x1[i].resize (nDofu);
      }

      vector < double >  soluNonLoc (nDofu);
      vector < double >  soluLoc (nDofu);

      for (unsigned i = 0; i < nDofu; i++) {
        unsigned solDof = msh->GetSolutionDof (i, iel, soluType);
        soluNonLoc[i] = (*sol->_Sol[soluIndex]) (solDof);
        soluLoc[i] = (*sol->_Sol[soluIndexLocal]) (solDof);
      }

      for (unsigned i = 0; i < nDofu; i++) {
        unsigned xDof  = msh->GetSolutionDof (i, iel, xType);

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          x1[jdim][i] = (*msh->_topology->_Sol[jdim]) (xDof);
        }
      }

      vector <double> phi;  // local test function
      vector <double> phi_x; // local test function first order partial derivatives
      double weight; // gauss point weight

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[ielGeom][soluType]->Jacobian (x1, ig, weight, phi, phi_x);
        double soluNonLoc_gss = 0.;
        double soluLoc_gss = 0.;
        double soluExact_gss = 0.;
        double x_gss = 0.;


        for (unsigned i = 0; i < nDofu; i++) {
          soluNonLoc_gss += phi[i] * soluNonLoc[i];
          soluLoc_gss += phi[i] * soluLoc[i];
          x_gss += phi[i] * x1[0][i]; // this is x at the Gauss point

        }

        double u1 = a1 + b1 * x_gss - 1. / (2. * kappa1) * x_gss * x_gss;
        double u2 = a2 + b2 * x_gss - 1. / (2. * kappa2) * x_gss * x_gss;

        soluExact_gss = (x_gss < 0.) ? u1 : u2;


//       soluExact_gss = soluExact_gss * soluExact_gss * soluExact_gss * soluExact_gss + 0.1 * soluExact_gss * soluExact_gss; // this is x^4 + delta * x^2

//       soluExact_gss = soluExact_gss * soluExact_gss; // this is x^2

//       soluExact_gss = soluExact_gss * soluExact_gss * soluExact_gss; // this is x^3

//       soluExact_gss = soluExact_gss * soluExact_gss * soluExact_gss * soluExact_gss; // this is x^4

//       soluExact_gss = 2 * soluExact_gss  + soluExact_gss * soluExact_gss * soluExact_gss * soluExact_gss * soluExact_gss ; // this is 2x + x^5

        error_solExact_norm2 += (soluNonLoc_gss - soluExact_gss) * (soluNonLoc_gss - soluExact_gss) * weight;

        error_solExact_local_norm2 += (soluLoc_gss - soluExact_gss) * (soluLoc_gss - soluExact_gss) * weight;

        error_solLocal_norm2 += (soluNonLoc_gss - soluLoc_gss) * (soluNonLoc_gss - soluLoc_gss) * weight;

        solNonlocal_norm2 += soluNonLoc_gss * soluNonLoc_gss * weight;

        solLocal_norm2 += soluLoc_gss * soluLoc_gss * weight;

        sol_exact_norm2 += soluExact_gss * soluExact_gss * weight;
      }
    }
  }

  double norm2 = 0.;
  MPI_Allreduce (&error_solExact_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double norm = sqrt (norm2);
  std::cout.precision (16);
  std::cout << "L2 norm of ERROR: Nonlocal - exact = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce (&error_solExact_local_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt (norm2);
  std::cout.precision (16);
  std::cout << "L2 norm of ERROR: Local - exact = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce (&error_solLocal_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt (norm2);
  std::cout.precision (16);
  std::cout << "L2 norm of ERROR: Nonlocal - local = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce (&solNonlocal_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt (norm2);
  std::cout.precision (16);
  std::cout << "L2 norm of NONLOCAL soln = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce (&solLocal_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt (norm2);
  std::cout.precision (16);
  std::cout << "L2 norm of LOCAL soln = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce (&sol_exact_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt (norm2);
  std::cout.precision (16);
  std::cout << "L2 norm of EXACT soln = " << norm << std::endl;


  //BEGIN computation of the l2 and linfinity norms

//   double littleL2norm = 0.;
//   std::vector<double> littleLInfinitynorm (nprocs, 0.);
//
//   for (unsigned i =  msh->_dofOffset[soluType][iproc]; i <  msh->_dofOffset[soluType][iproc + 1]; i++) {
//
//     double nonLocalNodalValue = (*sol->_Sol[soluIndex]) (i);
//     double LocalNodalValue = (*sol->_Sol[soluIndexLocal]) (i);
//
//     double difference = fabs (nonLocalNodalValue - LocalNodalValue);
//
//     if (difference > littleLInfinitynorm[iproc]) littleLInfinitynorm[iproc] = difference;
//
//     littleL2norm += difference * difference;
//
//   }
//
//   norm2 = 0.;
//   MPI_Allreduce (&littleL2norm, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//   norm = sqrt (norm2);
//   std::cout.precision (16);
//   std::cout << "l2 norm of ERROR: Nonlocal - local = " << norm << std::endl;
//
//   for (int kproc = 0; kproc < nprocs; kproc++) {
//     MPI_Bcast (&littleLInfinitynorm[iproc], 1, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
//   }
//
//   double littleLInfinityNorm = littleLInfinitynorm[0];
//
//   for (unsigned kproc = 0; kproc < nprocs; kproc++) {
//     if (littleLInfinitynorm[kproc] > littleLInfinityNorm) littleLInfinityNorm = littleLInfinitynorm[kproc];
//   }
//
//   std::cout.precision (16);
//   std::cout << "linfinity norm of ERROR: Nonlocal - local = " << littleLInfinityNorm << std::endl;

  //END

  //BEGIN nonlocal fine - coarse L2 norm on fine grid

  std::cout << "------------------------------------- " << std::endl;

  unsigned    iprocFine = mshFine->processor_id(); // get the process_id (for parallel computation)

  for (int ielFine = solFine->GetMesh()->_elementOffset[iprocFine]; ielFine < solFine->GetMesh()->_elementOffset[iprocFine + 1]; ielFine ++) {

    if (elementSkipFlagsFine[ielFine] == 0) {

      short unsigned ielFineGeom = mshFine->GetElementType (ielFine);
      unsigned nDofFine  = mshFine->GetElementDofNumber (ielFine, soluType);

      vector < double >  soluNonLocFine (nDofFine);

      std::vector < std::vector <double> > xFine (dim);

      for (int k = 0; k < dim; k++) {
        xFine[k].assign (nDofFine, 0.);
      }

      for (unsigned i = 0; i < nDofFine; i++) {
        unsigned solDof = mshFine->GetSolutionDof (i, ielFine, soluType);
        soluNonLocFine[i] = (*solFine->_Sol[soluIndexFine]) (solDof);
      }

      for (unsigned i = 0; i < nDofFine; i++) {
        unsigned xDof  = mshFine->GetSolutionDof (i, ielFine, xType);

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          xFine[jdim][i] = (*mshFine->_topology->_Sol[jdim]) (xDof);
        }
      }

      vector <double> phi;  // local test function
      vector <double> phi_x; // local test function first order partial derivatives
      double weight; // gauss point weight

      // *** Gauss point loop ***
      unsigned igNumber = mshFine->_finiteElement[ielFineGeom][soluType]->GetGaussPointNumber();
//     unsigned igNumber = femQuadrature->GetGaussPointNumber();

      for (unsigned ig = 0; ig < igNumber; ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        mshFine->_finiteElement[ielFineGeom][soluType]->Jacobian (xFine, ig, weight, phi, phi_x);
//       femQuadrature->Jacobian (xFine, ig, weight, phi, phi_x);

        double soluNonLocCoarse_gss = 0.;
        double soluNonLocFine_gss = 0.;
        double x_gss_fine = 0.;


        for (unsigned i = 0; i < nDofFine; i++) {
          soluNonLocFine_gss += phi[i] * soluNonLocFine[i];
          x_gss_fine += phi[i] * xFine[0][i]; // this is x at the Gauss point

        }

        //BEGIN computation of the coarse solution at the fine Gauss point

        for (int ielCoarse = sol->GetMesh()->_elementOffset[iproc]; ielCoarse < sol->GetMesh()->_elementOffset[iproc + 1]; ielCoarse++) {

          if (elementSkipFlags[ielCoarse] == 0) {

            short unsigned ielGeomCoarse = msh->GetElementType (ielCoarse);
            unsigned nDofCoarse  = msh->GetElementDofNumber (ielCoarse, soluType);

            vector < double >  soluNonLocCoarse (nDofCoarse);

            std::vector < std::vector <double> > xCoarse (dim);

            for (int k = 0; k < dim; k++) {
              xCoarse[k].assign (nDofCoarse, 0.);
            }

            unsigned  xLeftDof = sol->GetMesh()->GetSolutionDof (0, ielCoarse, 2);
            unsigned  xRightDof = sol->GetMesh()->GetSolutionDof (1, ielCoarse, 2);

            xCoarse[0][0] = (*sol->GetMesh()->_topology->_Sol[0]) (xLeftDof);
            xCoarse[0][1] = (*sol->GetMesh()->_topology->_Sol[0]) (xRightDof);

            if ( (x_gss_fine > xCoarse[0][0] && x_gss_fine < xCoarse[0][1]) || fabs (x_gss_fine - xCoarse[0][0]) < 1.e-10 || fabs (x_gss_fine - xCoarse[0][1]) < 1.e-10) {

//           std::cout.precision (16);
//           std::cout << "FOUND ----------> " << xCoarse[0][0] << " , " << xCoarse[0][1] << " , " << x_gss_fine << std::endl;

              for (unsigned jdof = 0; jdof < nDofCoarse; jdof++) {
                unsigned solDof = msh->GetSolutionDof (jdof, ielCoarse, soluType);
                soluNonLocCoarse[jdof] = (*sol->_Sol[soluIndex]) (solDof);
              }

              std::vector<double> x_gss_fine_Local (1, 0.);
              x_gss_fine_Local[0] = - 1. + 2. * (x_gss_fine - xCoarse[0][0]) / (xCoarse[0][1] - xCoarse[0][0]);

//           std::cout << "=====-------------------==== x_gss_fine_Local[0] = " << x_gss_fine_Local[0] << std::endl;

              vector <double> phi2;  // local test function
              vector <double> phi_x2; // local test function first order partial derivatives
              double weight2; // gauss point weight
              msh->_finiteElement[ielGeomCoarse][soluType]->Jacobian (xCoarse, x_gss_fine_Local, weight2, phi2, phi_x2);
//           femQuadrature->Jacobian (xCoarse, x_gss_fine_Local, weight2, phi2, phi_x2);

              for (unsigned jdof = 0; jdof < nDofCoarse; jdof++) {
                soluNonLocCoarse_gss += soluNonLocCoarse[jdof] * phi2[jdof];
              }

              break;

            }

          }

        }

        //END computation of the fine solution at the coarse Gauss point


        error_NonLocCoarse_NonLocFine_norm2 += (soluNonLocCoarse_gss - soluNonLocFine_gss) * (soluNonLocCoarse_gss - soluNonLocFine_gss) * weight;

        solNonlocalFine_norm2 += soluNonLocFine_gss * soluNonLocFine_gss * weight;

      }

    }

  }

  norm2 = 0.;
  MPI_Allreduce (&error_NonLocCoarse_NonLocFine_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt (norm2);
  std::cout.precision (16);
  std::cout << "L2 norm of ERROR: Nonlocal - Nonlocal Fine = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce (&solNonlocalFine_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt (norm2);
  std::cout.precision (16);
  std::cout << "L2 norm of NONLOCAL FINE soln = " << norm << std::endl;

  //END

//   for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
//
//     unsigned xDof = msh->GetSolutionDof (0, iel, soluType);
//
//     double x = (*msh->_topology->_Sol[0]) (xDof);
//
//     double u = (*sol->_Sol[soluIndex]) (xDof);
//
//     double u_local = (*sol->_Sol[soluIndexLocal]) (xDof);
//
//     double u_exact = x * x ;
//
//     std::cout << x << " " << u << " " << u_local << " " << u_exact << std::endl;
//
//   }

//   for (unsigned idof = msh->_dofOffset[0][iproc]; idof < msh->_dofOffset[0][iproc + 1]; idof++) {
//
//     double x = (*msh->_topology->_Sol[0]) (idof);
//
//     double u = (*sol->_Sol[soluIndex]) (idof);
//
//     double u_local = (*sol->_Sol[soluIndexLocal]) (idof);
//
//     double u1 = a1 + b1 * x - 1. / (2. * kappa1) * x * x;
//     double u2 = a2 + b2 * x - 1. / (2. * kappa2) * x * x;
//
//     double u_exact = (x < 0.) ? u1 : u2;
//
//
//     std::cout << x << " " << u << " " << u_local << " " << u_exact << std::endl;
//
//   }

}



void PutADoubleNodeAtTheInterface (MultiLevelMesh & mlMsh, const double & meshSize, double & leftBound, double & rightBound) {

  unsigned level = 0;

  unsigned xType = 2;

  Mesh* msh = mlMsh.GetLevel (level);

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the noumber of processes (for parallel computation)

  //BEGIN TO REMOVE
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned xMinDof  = msh->GetSolutionDof (0, iel, xType);
    unsigned xMaxDof  = msh->GetSolutionDof (1, iel, xType);
    unsigned xMidDof  = msh->GetSolutionDof (2, iel, xType);


    double xMin = (*msh->_topology->_Sol[0]) (xMinDof);
    double xMax = (*msh->_topology->_Sol[0]) (xMaxDof);
    double xMid = (*msh->_topology->_Sol[0]) (xMidDof);

    std::cout.precision (16);
    std::cout << "xMin prima del move= " << xMin << " , " << "xMid = " << xMid << " , " << "xMax = " << xMax << std::endl;

  }

  //END

  unsigned numberOfNodes = msh->GetNumberOfNodes();

  std::vector<unsigned> nodeShiftFlags (numberOfNodes, 0);

  unsigned leftDofsIproc = msh->_dofOffset[xType][iproc];

  unsigned rightDofsIproc = msh->_dofOffset[xType][iproc + 1];

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned xMinDof  = msh->GetSolutionDof (0, iel, xType);
    unsigned xMaxDof  = msh->GetSolutionDof (1, iel, xType);
    unsigned xMidDof  = msh->GetSolutionDof (2, iel, xType);

    double xMin = (*msh->_topology->_Sol[0]) (xMinDof);
    double xMax = (*msh->_topology->_Sol[0]) (xMaxDof);
    double xMid = (*msh->_topology->_Sol[0]) (xMidDof);

//     if (fabs (xMid) <= 1.e-14) {
//       elementToSkip = iel;
//       elementToSkipFound = true;
//       procWhoFoundIt = iproc;
//     }

    bool iprocOwnsXmin = (leftDofsIproc <= xMinDof < rightDofsIproc) ? true : false;
    bool iprocOwnsXmax = (leftDofsIproc <= xMaxDof < rightDofsIproc) ? true : false;

    if (nodeShiftFlags[xMinDof] == 0 && iprocOwnsXmin) {

      if (xMin < 0.) msh->_topology->_Sol[0]->set (xMinDof, xMin + 0.5 * meshSize);

      else msh->_topology->_Sol[0]->set (xMinDof, xMin - 0.5 * meshSize);

      nodeShiftFlags[xMinDof] = 1;

    }

    if (nodeShiftFlags[xMaxDof] == 0 && iprocOwnsXmax) {

      if (xMax < 0.) msh->_topology->_Sol[0]->set (xMaxDof, xMax + 0.5 * meshSize);

      else msh->_topology->_Sol[0]->set (xMaxDof, xMax - 0.5 * meshSize);

      nodeShiftFlags[xMaxDof] = 1;

    }

  }

  msh->_topology->_Sol[0]->close();

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned xMinDof  = msh->GetSolutionDof (0, iel, xType);
    unsigned xMaxDof  = msh->GetSolutionDof (1, iel, xType);
    unsigned xMidDof  = msh->GetSolutionDof (2, iel, xType);


    double xMin = (*msh->_topology->_Sol[0]) (xMinDof);
    double xMax = (*msh->_topology->_Sol[0]) (xMaxDof);

    msh->_topology->_Sol[0]->set (xMidDof, 0.5 * (xMin + xMax));

  }

  msh->_topology->_Sol[0]->close();

  leftBound += 0.5 * meshSize;
  rightBound -= 0.5 * meshSize;

//         //BEGIN TO REMOVE
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned xMinDof  = msh->GetSolutionDof (0, iel, xType);
    unsigned xMaxDof  = msh->GetSolutionDof (1, iel, xType);
    unsigned xMidDof  = msh->GetSolutionDof (2, iel, xType);


    double xMin = (*msh->_topology->_Sol[0]) (xMinDof);
    double xMax = (*msh->_topology->_Sol[0]) (xMaxDof);
    double xMid = (*msh->_topology->_Sol[0]) (xMidDof);

    std::cout.precision (16);
    std::cout << "xMin mesh spostato= " << xMin << " , " << "xMid = " << xMid << " , " << "xMax = " << xMax << std::endl;

  }

  //END

//         END

}


void ShiftTheExtrema (MultiLevelMesh & mlMsh, const double & meshSize, const double & delta1Shift, const double & delta2Shift, double & leftBound, double & rightBound) {

  unsigned level = 0;

  unsigned xType = 2;

  Mesh* msh = mlMsh.GetLevel (level);

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the noumber of processes (for parallel computation)

//   //         //BEGIN TO REMOVE
//   for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
//
//     unsigned xMinDof  = msh->GetSolutionDof (0, iel, xType);
//     unsigned xMaxDof  = msh->GetSolutionDof (1, iel, xType);
//     unsigned xMidDof  = msh->GetSolutionDof (2, iel, xType);
//
//
//     double xMin = (*msh->_topology->_Sol[0]) (xMinDof);
//     double xMax = (*msh->_topology->_Sol[0]) (xMaxDof);
//     double xMid = (*msh->_topology->_Sol[0]) (xMidDof);
//
//     std::cout << "xMin prima dello shift= " << xMin << " , " << "xMid = " << xMid << " , " << "xMax = " << xMax << std::endl;
//
//   }
//   //END

  if (!doubleIntefaceNode) {
    leftBound += 0.5 * meshSize;
    rightBound -= 0.5 * meshSize;
  }

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned xMinDof  = msh->GetSolutionDof (0, iel, xType);
    unsigned xMaxDof  = msh->GetSolutionDof (1, iel, xType);
    unsigned xMidDof  = msh->GetSolutionDof (2, iel, xType);

    double xMin = (*msh->_topology->_Sol[0]) (xMinDof);
    double xMax = (*msh->_topology->_Sol[0]) (xMaxDof);

    if (fabs (xMin - (leftBound - delta1Mesh)) <= 1.e-10) {
      msh->_topology->_Sol[0]->set (xMinDof, xMin + delta1Shift);
      msh->_topology->_Sol[0]->set (xMidDof, 0.5 * (xMin + delta1Shift + xMax));
    }

    if (fabs (xMax - (rightBound + delta2Mesh)) <= 1.e-10) {
      msh->_topology->_Sol[0]->set (xMaxDof, xMax - delta2Shift);
      msh->_topology->_Sol[0]->set (xMidDof, 0.5 * (xMin + xMax - delta2Shift));
    }

  }

  msh->_topology->_Sol[0]->close();

//   //         //BEGIN TO REMOVE
//   for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
//
//     unsigned xMinDof  = msh->GetSolutionDof (0, iel, xType);
//     unsigned xMaxDof  = msh->GetSolutionDof (1, iel, xType);
//     unsigned xMidDof  = msh->GetSolutionDof (2, iel, xType);
//
//
//     double xMin = (*msh->_topology->_Sol[0]) (xMinDof);
//     double xMax = (*msh->_topology->_Sol[0]) (xMaxDof);
//     double xMid = (*msh->_topology->_Sol[0]) (xMidDof);
//
//     std::cout << "xMin dopo lo shift= " << xMin << " , " << "xMid = " << xMid << " , " << "xMax = " << xMax << std::endl;
//
//   }

}

void BuildElementSkipFlags (MultiLevelMesh & mlMsh, std::vector<unsigned> &elementSkipFlags) {

  const unsigned level = mlMsh.GetNumberOfLevels() - 1;

  unsigned xType = 2;

  Mesh* msh = mlMsh.GetLevel (level);

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the noumber of processes (for parallel computation)

  unsigned numberOfElements = msh->GetNumberOfElements();

  elementSkipFlags.assign (numberOfElements, 0);

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned xMinDof  = msh->GetSolutionDof (0, iel, xType);
    unsigned xMaxDof  = msh->GetSolutionDof (1, iel, xType);

    double xMin = (*msh->_topology->_Sol[0]) (xMinDof);
    double xMax = (*msh->_topology->_Sol[0]) (xMaxDof);

    if (fabs (xMax - xMin) <= 1.e-14) elementSkipFlags[iel] = 1;

  }

  for (unsigned iel = 0; iel < numberOfElements; iel++) {
    unsigned flag = 0;
    MPI_Allreduce (&elementSkipFlags[iel], &flag, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    elementSkipFlags[iel] = flag;
  }

  //BEGIN TO REMOVE (use with serial)
//   for (unsigned iel = 0; iel < numberOfElements; iel++) {
//
//     unsigned xMinDof  = msh->GetSolutionDof (0, iel, xType);
//     unsigned xMaxDof  = msh->GetSolutionDof (1, iel, xType);
//
//     double xMin = (*msh->_topology->_Sol[0]) (xMinDof);
//     double xMax = (*msh->_topology->_Sol[0]) (xMaxDof);
//
//     std::cout << xMin << " , " << xMax << " , " << "flag = " << elementSkipFlags[iel] << std::endl ;
//
//
//   }
  //END TO REMOVE

}

void GetL2NormLocalConnectedNonlocalFine (MultiLevelSolution & mlSol, MultiLevelSolution & mlSolFine) {

  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel (level);
  Solution* sol  = mlSol.GetSolutionLevel (level);

  const unsigned levelFine = mlSolFine._mlMesh->GetNumberOfLevels() - 1;
  Mesh* mshFine = mlSolFine._mlMesh->GetLevel (levelFine);
  Solution* solFine  = mlSolFine.GetSolutionLevel (levelFine);

  const unsigned  dim = msh->GetDimension();

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  unsigned soluIndexLocal;
  soluIndexLocal = mlSol.GetIndex ("u_local");
  unsigned soluType = mlSol.GetSolutionType (soluIndexLocal);

  unsigned soluIndexFine;
  soluIndexFine = mlSolFine.GetIndex ("u_fine");

  unsigned    iproc = msh->processor_id();
  unsigned    nprocs = msh->n_processors();

  double error_Loc_NonLocFine_norm2 = 0.;

  //BEGIN ||local on connected grid - nonlocal on nonconnected grid||_L2 on nonconnected grid

  std::cout << "------------------------------------- " << std::endl;

  unsigned    iprocFine = mshFine->processor_id(); // get the process_id (for parallel computation)

  for (int kproc = 0; kproc < nprocs; kproc++) {
    for (int ielLocal = sol->GetMesh()->_elementOffset[kproc]; ielLocal < sol->GetMesh()->_elementOffset[kproc + 1]; ielLocal++) {

      short unsigned ielGeomLocal;
      unsigned nDofLocal;

      if (iproc == kproc) {
        ielGeomLocal = msh->GetElementType (ielLocal);
        nDofLocal  = msh->GetElementDofNumber (ielLocal, soluType);
      }

      MPI_Bcast (&ielGeomLocal, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD);
      MPI_Bcast (&nDofLocal, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);


      std::vector < std::vector <double> > xLocal (dim);
      vector < double >  soluLocLocal (nDofLocal);

      for (int k = 0; k < dim; k++) {
        xLocal[k].resize (nDofLocal);
      }

      if (iproc == kproc) {
        for (unsigned jdof = 0; jdof < nDofLocal; jdof++) {
          unsigned xDof  = msh->GetSolutionDof (jdof, ielLocal, xType);
          unsigned solDof = msh->GetSolutionDof (jdof, ielLocal, soluType);
          soluLocLocal[jdof] = (*sol->_Sol[soluIndexLocal]) (solDof);
          for (unsigned k = 0; k < dim; k++) {
            xLocal[k][jdof] = (*sol->GetMesh()->_topology->_Sol[k]) (xDof);
          }
        }

        std::vector<int> dofsTemp (nDofLocal);
//         ReorderElement (dofsTemp, soluLocLocal, xLocal);
      }

      MPI_Bcast (&soluLocLocal[0], nDofLocal, MPI_DOUBLE, kproc, MPI_COMM_WORLD);

      for (unsigned k = 0; k < dim; k++) {
        MPI_Bcast (& xLocal[k][0], nDofLocal, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }

      std::vector<double> xMinAndMax (nDofLocal);
      xMinAndMax.assign (2, 0.);

      for (unsigned i = 0; i < nDofLocal; i++) {
        xMinAndMax[i] = xLocal[0][i];
      }

      for (unsigned i = 0; i < nDofLocal; i++) {
        if (xLocal[0][i] < xMinAndMax[0]) xMinAndMax[0] = xLocal[0][i];

        if (xLocal[0][i] > xMinAndMax[1]) xMinAndMax[1] = xLocal[0][i];
      }

      for (int ielFine = solFine->GetMesh()->_elementOffset[iprocFine]; ielFine < solFine->GetMesh()->_elementOffset[iprocFine + 1]; ielFine ++) {

        short unsigned ielFineGeom = mshFine->GetElementType (ielFine);
        unsigned nDofFine  = mshFine->GetElementDofNumber (ielFine, soluType);

        vector < double >  soluNonLocFine (nDofFine);

        std::vector < std::vector <double> > xFine (dim);

        for (int k = 0; k < dim; k++) {
          xFine[k].assign (nDofFine, 0.);
        }

        for (unsigned i = 0; i < nDofFine; i++) {
          unsigned solDof = mshFine->GetSolutionDof (i, ielFine, soluType);
          soluNonLocFine[i] = (*solFine->_Sol[soluIndexFine]) (solDof);
          unsigned xDof  = mshFine->GetSolutionDof (i, ielFine, xType);
          for (unsigned jdim = 0; jdim < dim; jdim++) {
            xFine[jdim][i] = (*mshFine->_topology->_Sol[jdim]) (xDof);
          }
        }

        vector <double> phi;  // local test function
        vector <double> phi_x; // local test function first order partial derivatives
        double weight; // gauss point weight

        unsigned igNumberFine = mshFine->_finiteElement[ielFineGeom][soluType]->GetGaussPointNumber();
//     unsigned igNumberFine = femQuadrature->GetGaussPointNumber();

        // *** Gauss point loop ***
        for (unsigned ig = 0; ig < igNumberFine; ig++) {
          // *** get gauss point weight, test function and test function partial derivatives ***
          mshFine->_finiteElement[ielFineGeom][soluType]->Jacobian (xFine, ig, weight, phi, phi_x);
//       femQuadrature->Jacobian (xFine, ig, weight, phi, phi_x);

          double soluLoc_gss = 0.;
          double soluNonLocFine_gss = 0.;
          std::vector <double> x_gss_fine (dim, 0.);


          for (unsigned i = 0; i < nDofFine; i++) {
            soluNonLocFine_gss += phi[i] * soluNonLocFine[i];
            for (unsigned jdim = 0; jdim < dim; jdim++) {
              x_gss_fine[jdim] += phi[i] * xFine[jdim][i];
            }
          }

          unsigned fineGaussPointInLocalMeshElem = 0;

          if ( (x_gss_fine[0] > xMinAndMax[0] && x_gss_fine[0] < xMinAndMax[1]) || fabs (x_gss_fine[0] - xMinAndMax[0]) < 1.e-10 || fabs (x_gss_fine[0] - xMinAndMax[1]) < 1.e-10) {
            fineGaussPointInLocalMeshElem++;
          }

          if (fineGaussPointInLocalMeshElem == dim) {


            std::vector<double> x_gss_fine_Local (dim, 0.);
            for (unsigned ii = 0; ii < dim; ii++) {
              x_gss_fine_Local[ii] = - 1. + 2. * (x_gss_fine[ii] - xLocal[ii][ii]) / (xLocal[ii][ii + 1] - xLocal[ii][ii]);
            }

            vector <double> phi2;  // local test function
            vector <double> phi_x2; // local test function first order partial derivatives
            double weight2; // gauss point weight
            msh->_finiteElement[ielGeomLocal][soluType]->Jacobian (xLocal, x_gss_fine_Local, weight2, phi2, phi_x2);
//           femQuadrature->Jacobian (xLocal, x_gss_fine_Local, weight2, phi2, phi_x2);

            for (unsigned jdof = 0; jdof < nDofLocal; jdof++) {
              soluLoc_gss += soluLocLocal[jdof] * phi2[jdof];
            }

            error_Loc_NonLocFine_norm2 += (soluLoc_gss - soluNonLocFine_gss) * (soluLoc_gss - soluNonLocFine_gss) * weight;

          }
        }
      }

    }

  }

  double norm2 = 0.;
  MPI_Allreduce (&error_Loc_NonLocFine_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double norm = sqrt (norm2);
  std::cout.precision (14);
  std::cout << "L2 norm of ERROR: Local (connected grid) - Nonlocal (nonconnected grid) = " << norm << std::endl;

  //END

}






