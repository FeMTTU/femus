
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
//     return x[0] + 0. * ( 0.51 * 0.51 - x[0] * x[0] ) * ( 0.51 * 0.51 - x[1] * x[1] );
//     return x[0];
  return x[0] * x[0];
//     return x[0] * x[0] * x[0];
//     return x[0] * x[0] * x[0] * x[0] + 0.1 * x[0] * x[0]; //this is x^4 + delta x^2
//     return x[0] * x[0] * x[0] * x[0]; //this is x^4
//        return 2 * x[0] + x[0] * x[0] * x[0] * x[0] * x[0]; //this is 2x + x^5
}

void GetL2Norm (MultiLevelSolution &mlSol);

void PutADoubleNodeAtTheInterface (MultiLevelMesh &mlMsh);

void ShiftTheExtrema (MultiLevelMesh &mlMsh);

bool SetBoundaryCondition (const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {

  bool dirichlet = true;
//     value = 0.;
//     value = x[0];
  value = x[0] * x[0];
//     value =  x[0] * x[0] * x[0];
//     value = x[0] * x[0] * x[0] * x[0] + 0.1 * x[0] * x[0]; //this is x^4 + delta x^2
//     value = x[0] * x[0] * x[0] * x[0];
//        value =  2 * x[0] + x[0] * x[0] * x[0] * x[0] * x[0]; //this is 2x + x^5

  if (fabs (x[0]) <= 1.e-10) {
    dirichlet = false;
    value = 0.;
  }

  return dirichlet;
}

unsigned numberOfUniformLevels = 1; //keep = 1

int main (int argc, char** argv) {

  // init Petsc-MPI communicator
  FemusInit mpinit (argc, argv, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;

  double xMinCoarseBox = leftBound - delta1Mesh;
  double xMaxCoarseBox = rightBound + delta2Mesh;

  mlMsh.GenerateCoarseBoxMesh (numberOfElements, 0, 0, xMinCoarseBox, xMaxCoarseBox, 0., 0., 0., 0., EDGE3, "fifth");

  if (doubleIntefaceNode) PutADoubleNodeAtTheInterface (mlMsh);

  if (shiftExternalNodes) ShiftTheExtrema (mlMsh);

  mlMsh.RefineMesh (numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);

  mlMsh.EraseCoarseLevels (numberOfUniformLevels - 1);
//     numberOfUniformLevels = 1;

  unsigned dim = mlMsh.GetDimension();

  MultiLevelSolution mlSol (&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution ("u", LAGRANGE, FIRST, 2);

  mlSol.AddSolution ("u_local", LAGRANGE, FIRST, 2);

  mlSol.AddSolution ("u_exact", LAGRANGE, FIRST, 2);

  mlSol.Initialize ("All");

  mlSol.Initialize ("u_exact", InitalValueU);

  mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc ("All");


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

  // ******* Set Preconditioner *******
  system.SetMgSmoother (GMRES_SMOOTHER);

  system.SetSparsityPatternMinimumSize (100u);   //TODO tune

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids (GMRES);

  system.SetPreconditionerFineGrids (ILU_PRECOND);

  system.SetTolerances (1.e-20, 1.e-20, 1.e+50, 100);

// ******* Solution *******

  system.MGsolve();

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

  // ******* Set Preconditioner *******
  system2.SetMgSmoother (GMRES_SMOOTHER);

  system2.init();

  // ******* Set Smoother *******
  system2.SetSolverFineGrids (GMRES);

  system2.SetPreconditionerFineGrids (ILU_PRECOND);

  system2.SetTolerances (1.e-20, 1.e-20, 1.e+50, 100);

// ******* Solution *******

  system2.MGsolve();

  //END assemble and solve local problem


  //BEGIN compute errors
  GetL2Norm (mlSol);
  //END compute errors

  // ******* Print solution *******
  mlSol.SetWriter (VTK);
  std::vector<std::string> print_vars;
  print_vars.push_back ("All");
  mlSol.GetWriter()->SetDebugOutput (true);
  mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);

  std::cout << "Mesh size h = " << (xMaxCoarseBox - xMinCoarseBox) / (numberOfElements * pow (2, numberOfUniformLevels - 1)) << std::endl;

  return 0;

} //end main


void GetL2Norm (MultiLevelSolution &mlSol) {

  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel (level);
  elem*                     el = msh->el;
  Solution* sol  = mlSol.GetSolutionLevel (level);

  const unsigned  dim = msh->GetDimension();

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  double error_solExact_norm2 = 0.;

  double error_solExact_local_norm2 = 0.;

  double error_solLocal_norm2 = 0.;

  double solNonlocal_norm2 = 0.;

  double solLocal_norm2 = 0.;

  double sol_exact_norm2 = 0.;

  unsigned soluIndex;
  soluIndex = mlSol.GetIndex ("u");
  unsigned soluType = mlSol.GetSolutionType (soluIndex);

  unsigned soluIndexLocal;
  soluIndexLocal = mlSol.GetIndex ("u_local");

  unsigned    iproc = msh->processor_id();
  unsigned    nprocs = msh->n_processors();

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

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


      for (unsigned i = 0; i < nDofu; i++) {
        soluNonLoc_gss += phi[i] * soluNonLoc[i];
        soluLoc_gss += phi[i] * soluLoc[i];
        soluExact_gss += phi[i] * x1[0][i]; // this is x at the Gauss point

      }

//             soluExact_gss = soluExact_gss * soluExact_gss * soluExact_gss * soluExact_gss + 0.1 * soluExact_gss * soluExact_gss; // this is x^4 + delta * x^2

      soluExact_gss = soluExact_gss * soluExact_gss; // this is x^2

//             soluExact_gss = soluExact_gss * soluExact_gss * soluExact_gss; // this is x^3

//             soluExact_gss = soluExact_gss * soluExact_gss * soluExact_gss * soluExact_gss; // this is x^4

//             soluExact_gss = 2 * soluExact_gss  + soluExact_gss * soluExact_gss * soluExact_gss * soluExact_gss * soluExact_gss ; // this is 2x + x^5

      error_solExact_norm2 += (soluNonLoc_gss - soluExact_gss) * (soluNonLoc_gss - soluExact_gss) * weight;

      error_solExact_local_norm2 += (soluLoc_gss - soluExact_gss) * (soluLoc_gss - soluExact_gss) * weight;

      error_solLocal_norm2 += (soluNonLoc_gss - soluLoc_gss) * (soluNonLoc_gss - soluLoc_gss) * weight;

      solNonlocal_norm2 += soluNonLoc_gss * soluNonLoc_gss * weight;

      solLocal_norm2 += soluLoc_gss * soluLoc_gss * weight;

      sol_exact_norm2 += soluExact_gss * soluExact_gss * weight;
    }
  }

  double norm2 = 0.;
  MPI_Allreduce (&error_solExact_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double norm = sqrt (norm2);
  std::cout.precision (14);
  std::cout << "L2 norm of ERROR: Nonlocal - exact = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce (&error_solExact_local_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt (norm2);
  std::cout.precision (14);
  std::cout << "L2 norm of ERROR: Local - exact = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce (&error_solLocal_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt (norm2);
  std::cout.precision (14);
  std::cout << "L2 norm of ERROR: Nonlocal - local = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce (&solNonlocal_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt (norm2);
  std::cout.precision (14);
  std::cout << "L2 norm of NONLOCAL soln = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce (&solLocal_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt (norm2);
  std::cout.precision (14);
  std::cout << "L2 norm of LOCAL soln = " << norm << std::endl;

  norm2 = 0.;
  MPI_Allreduce (&sol_exact_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt (norm2);
  std::cout.precision (14);
  std::cout << "L2 norm of EXACT soln = " << norm << std::endl;


  double littleL2norm = 0.;
  std::vector<double> littleLInfinitynorm (nprocs, 0.);

  for (unsigned i =  msh->_dofOffset[soluType][iproc]; i <  msh->_dofOffset[soluType][iproc + 1]; i++) {

    double nonLocalNodalValue = (*sol->_Sol[soluIndex]) (i);
    double LocalNodalValue = (*sol->_Sol[soluIndexLocal]) (i);

    double difference = fabs (nonLocalNodalValue - LocalNodalValue);

    if (difference > littleLInfinitynorm[iproc]) littleLInfinitynorm[iproc] = difference;

    littleL2norm += difference * difference;

  }

  norm2 = 0.;
  MPI_Allreduce (&littleL2norm, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  norm = sqrt (norm2);
  std::cout.precision (14);
  std::cout << "l2 norm of ERROR: Nonlocal - local = " << norm << std::endl;

  for (int kproc = 0; kproc < nprocs; kproc++) {
    MPI_Bcast (&littleLInfinitynorm[iproc], 1, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
  }

  double littleLInfinityNorm = littleLInfinitynorm[0];

  for (unsigned kproc = 0; kproc < nprocs; kproc++) {
    if (littleLInfinitynorm[kproc] > littleLInfinityNorm) littleLInfinityNorm = littleLInfinitynorm[kproc];
  }

  std::cout.precision (14);
  std::cout << "linfinity norm of ERROR: Nonlocal - local = " << littleLInfinityNorm << std::endl;


  for (unsigned idof = msh->_dofOffset[0][iproc]; idof < msh->_dofOffset[0][iproc + 1]; idof++) {

    double x = (*msh->_topology->_Sol[0]) (idof);

    double u = (*sol->_Sol[soluIndex]) (idof);

    double u_local = (*sol->_Sol[soluIndexLocal]) (idof);

    double u_exact = x * x ;

    std::cout << x << " " << u << " " << u_local << " " << x * x << std::endl;

  }


}


void PutADoubleNodeAtTheInterface (MultiLevelMesh &mlMsh) {

  unsigned level = 0;

  unsigned xType = 2;

  Mesh* msh = mlMsh.GetLevel (level);

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the noumber of processes (for parallel computation)

//   //BEGIN TO REMOVE
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
//     std::cout << "xMin prima del move= " << xMin << " , " << "xMid = " << xMid << " , " << "xMax = " << xMax << std::endl;
//
//   }
//
//   //END


  unsigned x0Dof  = msh->GetSolutionDof (0, msh->_elementOffset[iproc], xType);
  unsigned x1Dof  = msh->GetSolutionDof (1, msh->_elementOffset[iproc], xType);

  double x0 = (*msh->_topology->_Sol[0]) (x0Dof);
  double x1 = (*msh->_topology->_Sol[0]) (x1Dof);

  meshSize = fabs (x1 - x0);

//   std::cout<<"MEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEESH " << meshSize << std::endl;

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

    if (fabs (xMid) <= 1.e-14) {
      elementToSkip = iel;
      elementToSkipFound = true;
      procWhoFoundIt = iproc;
    }

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

// //         //BEGIN TO REMOVE
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
//     std::cout << "xMin mesh spostato= " << xMin << " , " << "xMid = " << xMid << " , " << "xMax = " << xMax << std::endl;
//
//   }
//
//   //END

  if (elementToSkipFound) {
    for (unsigned kproc = 0; kproc < nprocs; kproc++) {
      if (kproc != iproc) MPI_Send (&elementToSkip, 1, MPI_UNSIGNED, kproc, 1 , PETSC_COMM_WORLD);
    }
  }

  else MPI_Recv (&elementToSkip, 1, MPI_UNSIGNED, procWhoFoundIt, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);

//         END

}


void ShiftTheExtrema (MultiLevelMesh &mlMsh) {

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

    if (fabs (xMin - (leftBound - delta1Mesh)) <= 1.e-10){
      msh->_topology->_Sol[0]->set (xMinDof, xMin + delta1Shift);
      msh->_topology->_Sol[0]->set (xMidDof, 0.5* (xMin + xMax));
    }

    if (fabs (xMax - (rightBound + delta2Mesh)) <= 1.e-10){
      msh->_topology->_Sol[0]->set (xMaxDof, xMax - delta2Shift);
      msh->_topology->_Sol[0]->set (xMidDof, 0.5* (xMin + xMax));
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
