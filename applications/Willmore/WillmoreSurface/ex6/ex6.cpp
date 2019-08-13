/* This example details the full implementation of the p-Willmore flow
 *   algorithm, which involves three nonlinear systems.
 *
 *   System0 AssembleInit computes the initial curvatures given mesh positions.
 *   System AssemblePWillmore solves the flow equations.
 *   System2 AssembleConformalMinimization "reparametrizes" the surface to
 *   correct the mesh. */

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "TransientSystem.hpp"
#include "adept.h"
#include <cstdlib>
#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

/* Vector option for P (to handle polynomials).
 * ap is the coefficient in front of the power of H. */
const unsigned P[3] = {2, 3, 4};
const double ap[3] = {1, 0., 0.};

using namespace femus;

// Toggle for setting volume and area constraints, as well as sign of N.
bool volumeConstraint = true;
bool areaConstraint = false;
const double normalSign = -1.;

// Penalty parameter for conformal minimization (eps).
// Trick for system0 (delta).
// Trick for system2 (timederiv).
const double eps = 1e-4;
const double delta = 0.005;
const double timederiv = 0.;

// Declaration of systems.
void CopyDisplacement (MultiLevelSolution &mlSol,  const bool &forward);
void AssembleInit (MultiLevelProblem&);
void AssemblePWillmore (MultiLevelProblem&);
void AssembleConformalMinimization (MultiLevelProblem&);  //stable and not bad
void AssembleShearMinimization (MultiLevelProblem&);  //vastly inferior
void AssembleO2ConformalMinimization (MultiLevelProblem&);  //vastly superior.. when convergent
void AssembleSphereConformalMinimization (MultiLevelProblem&);  //inferior

// Function to control the time stepping.
double GetTimeStep (const double t) {
  //if(time==0) return 1.0e-10;
  //return 0.0001;
  double dt0 = 0.1;
  //double dt0 = 0.00000032; // spot
  double s = 1.;
  double n = 2;
  return dt0 * pow (1. + t / pow (dt0, s), n);
}

// IBVs.  No boundary, and IVs set to sphere (just need something).
bool SetBoundaryCondition (const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = false;
  value = 0.;
  return dirichlet;
}
double InitalValueY1 (const std::vector < double >& x) {
  return -2. * x[0];
}
double InitalValueY2 (const std::vector < double >& x) {
  return -2. * x[1];
}
double InitalValueY3 (const std::vector < double >& x) {
  return -2. * x[2];
}
double InitalValueW1 (const std::vector < double >& x) {
  return -2. * P[2] * pow (2., P[2] - 2) * x[0];
}
double InitalValueW2 (const std::vector < double >& x) {
  return -2. * P[2] * pow (2., P[2] - 2) * x[1];
}
double InitalValueW3 (const std::vector < double >& x) {
  return -2. * P[2] * pow (2., P[2] - 2) * x[2];
}


// Main program starts here.
int main (int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit (argc, args, MPI_COMM_WORLD);


  // define multilevel mesh
  unsigned maxNumberOfMeshes;
  MultiLevelMesh mlMsh;

  // Read coarse level mesh and generate finer level meshes.
  double scalingFactor = 1.;

  //mlMsh.ReadCoarseMesh("../input/torus.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/sphere.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/ellipsoidRef3.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/ellipsoidV1.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/genusOne.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/knot.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/cube.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/horseShoe.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/tiltedTorus.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/dog.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/virus3.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("../input/ellipsoidSphere.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("../input/CliffordTorus.neu", "seventh", scalingFactor);

  mlMsh.ReadCoarseMesh ("../input/duckpoint85.med", "seventh", scalingFactor);

  // Set number of mesh levels.
  unsigned numberOfUniformLevels = 2;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh (numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // Erase all the coarse mesh levels.
  mlMsh.EraseCoarseLevels (numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  // Define the multilevel solution and attach the mlMsh object to it.
  MultiLevelSolution mlSol (&mlMsh);

  // Add variables X,Y,W to mlSol.
  mlSol.AddSolution ("Dx1", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("Dx2", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("Dx3", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("W1", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("W2", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("W3", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("Y1", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("Y2", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("Y3", LAGRANGE, FIRST, 2);

  // Add variable "Lambda" based on constraint choice.
  if (volumeConstraint || areaConstraint) {
    mlSol.AddSolution ("Lambda", DISCONTINUOUS_POLYNOMIAL, ZERO, 0);
  }

  // Add variables "nDx" and "Lambda1" for the conformal system.
  mlSol.AddSolution ("nDx1", LAGRANGE, FIRST, 0);
  mlSol.AddSolution ("nDx2", LAGRANGE, FIRST, 0);
  mlSol.AddSolution ("nDx3", LAGRANGE, FIRST, 0);
  mlSol.AddSolution ("Lambda1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0);

  // Initialize the variables and attach boundary conditions.
  mlSol.Initialize ("All");
  mlSol.Initialize ("W1", InitalValueW1);
  mlSol.Initialize ("W2", InitalValueW2);
  mlSol.Initialize ("W3", InitalValueW3);
  mlSol.Initialize ("Y1", InitalValueY1);
  mlSol.Initialize ("Y2", InitalValueY2);
  mlSol.Initialize ("Y3", InitalValueY3);

  mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
  mlSol.GenerateBdc ("All");

  MultiLevelProblem mlProb (&mlSol);

  //Add system0 in mlProb to compute the initial curvature data.
  NonLinearImplicitSystem& system0 = mlProb.add_system < NonLinearImplicitSystem > ("Init");

  // Add solutions Y,W to system0.
  system0.AddSolutionToSystemPDE ("Y1");
  system0.AddSolutionToSystemPDE ("Y2");
  system0.AddSolutionToSystemPDE ("Y3");
  system0.AddSolutionToSystemPDE ("W1");
  system0.AddSolutionToSystemPDE ("W2");
  system0.AddSolutionToSystemPDE ("W3");

  // Parameters for convergence and # of iterations.
  system0.SetMaxNumberOfNonLinearIterations (5);
  system0.SetNonLinearConvergenceTolerance (1.e-12);

  // Add the assembling function to system0 and initialize.
  system0.SetAssembleFunction (AssembleInit);
  system0.init();

  // Add system P-Willmore in mlProb as a time-dependent system.
  TransientNonlinearImplicitSystem& system = mlProb.add_system < TransientNonlinearImplicitSystem > ("PWillmore");

  // Add solutions X, Y, W to P-Willmore system.
  system.AddSolutionToSystemPDE ("Dx1");
  system.AddSolutionToSystemPDE ("Dx2");
  system.AddSolutionToSystemPDE ("Dx3");
  system.AddSolutionToSystemPDE ("Y1");
  system.AddSolutionToSystemPDE ("Y2");
  system.AddSolutionToSystemPDE ("Y3");
  system.AddSolutionToSystemPDE ("W1");
  system.AddSolutionToSystemPDE ("W2");
  system.AddSolutionToSystemPDE ("W3");

  // Add solution Lambda to system based on constraint choice.
  if (volumeConstraint || areaConstraint) {
    system.AddSolutionToSystemPDE ("Lambda");
    system.SetNumberOfGlobalVariables (volumeConstraint + areaConstraint);
  }

  // Parameters for convergence and # of iterations for Willmore.
  system.SetMaxNumberOfNonLinearIterations (15);
  system.SetNonLinearConvergenceTolerance (1.e-12);

  // Attach the assembling function to P-Willmore system.
  system.SetAssembleFunction (AssemblePWillmore);

  // Attach time step function to P-Willmore sysyem.
  system.AttachGetTimeIntervalFunction (GetTimeStep);

  // Initialize the P-Willmore system.
  system.init();
  system.SetMgType (V_CYCLE);

  // Add system2 Conformal Minimization in mlProb.
  NonLinearImplicitSystem& system2 = mlProb.add_system < NonLinearImplicitSystem > ("nProj");

  // Add solutions newDX, Lambda1 to system2.
  system2.AddSolutionToSystemPDE ("nDx1");
  system2.AddSolutionToSystemPDE ("nDx2");
  system2.AddSolutionToSystemPDE ("nDx3");
  system2.AddSolutionToSystemPDE ("Lambda1");

  // Parameters for convergence and # of iterations.
  system2.SetMaxNumberOfNonLinearIterations (40);
  system2.SetNonLinearConvergenceTolerance (1.e-10);

  // Attach the assembling function to system2 and initialize.
  system2.SetAssembleFunction (AssembleConformalMinimization);
  system2.init();

  mlSol.SetWriter (VTK);
  std::vector<std::string> mov_vars;
  mov_vars.push_back ("Dx1");
  mov_vars.push_back ("Dx2");
  mov_vars.push_back ("Dx3");
  mlSol.GetWriter()->SetMovingMesh (mov_vars);

  // and this?
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back ("All");

  // and this?
  mlSol.GetWriter()->SetDebugOutput (false);
  mlSol.GetWriter()->Write ("./output1",
                            "linear", variablesToBePrinted, 0);

  // First, solve system2 to "conformalize" the initial mesh.
  CopyDisplacement (mlSol, true);
  system2.MGsolve();

  // Then, solve system0 to compute initial curvatures.
  CopyDisplacement (mlSol, false);
  system0.MGsolve();

  // mlSol.SetWriter (VTK);
  // std::vector<std::string> mov_vars;
  // mov_vars.push_back ("Dx1");
  // mov_vars.push_back ("Dx2");
  // mov_vars.push_back ("Dx3");
  // mlSol.GetWriter()->SetMovingMesh (mov_vars);
  //
  // std::vector < std::string > variablesToBePrinted;
  // variablesToBePrinted.push_back ("All");
  //
  // mlSol.GetWriter()->SetDebugOutput (true);
  mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "linear", variablesToBePrinted, 0);

  // Parameters for the main algorithm loop.
  unsigned numberOfTimeSteps = 10000u;
  unsigned printInterval = 1u;

  // Main algorithm loop.
  for (unsigned time_step = 0; time_step < numberOfTimeSteps; time_step++) {
    system.CopySolutionToOldSolution();
    system.MGsolve();

    if (time_step % 1 == 0) {
      // mlSol.GetWriter()->Write ("./output1", "linear",
      //                           variablesToBePrinted, (time_step + 1) / printInterval);

      CopyDisplacement (mlSol, true);
      system2.MGsolve();

      CopyDisplacement (mlSol, false);
      system0.MGsolve();
    }

    if ( (time_step + 1) % printInterval == 0)
      mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "linear", variablesToBePrinted, (time_step + 1) / printInterval);
  }
  return 0;
}


/* @@@@@@@@@@@@@@@@@@@@@ ASSEMBLY FUNCTIONS BELOW @@@@@@@@@@@@@@@@@@@@@@ */

// Eugenio's standard FEMuS function.
void CopyDisplacement (MultiLevelSolution &mlSol,  const bool &forward) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* solution  = mlSol.GetSolutionLevel (level);
  Mesh* msh = mlSol._mlMesh->GetLevel (level);

  unsigned DIM = 3u;
  vector < unsigned > solDxIndex (DIM);
  solDxIndex[0] = mlSol.GetIndex ("Dx1");
  solDxIndex[1] = mlSol.GetIndex ("Dx2");
  solDxIndex[2] = mlSol.GetIndex ("Dx3");

  vector < unsigned > solNDxIndex (DIM);
  solNDxIndex[0] = mlSol.GetIndex ("nDx1");
  solNDxIndex[1] = mlSol.GetIndex ("nDx2");
  solNDxIndex[2] = mlSol.GetIndex ("nDx3");

  if (forward) {
    for (unsigned i = 0; i < DIM; i++) {
      * (solution->_Sol[solNDxIndex[i]]) = * (solution->_Sol[solDxIndex[i]]);
    }
  }
  else {
    for (unsigned i = 0; i < DIM; i++) {
      * (solution->_Sol[solDxIndex[i]]) = * (solution->_Sol[solNDxIndex[i]]);
    }
  }
}


// Build system0 to compute initial curvature data.
void AssembleInit (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  // Call the adept stack object.
  adept::Stack& s = FemusInit::_adeptStack;

  //  Extract pointers to the several objects that we are going to use.
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem> ("Init");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Point to the mesh and element objects.
  Mesh *msh = ml_prob._ml_msh->GetLevel (level);
  elem *el = msh->el;

  // Point to mlSol, solution (level), and equation (level) objects.
  MultiLevelSolution *mlSol = ml_prob._ml_sol;
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel (level);
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level];

  // Point to the global stiffness mtx and residual vectors in pdeSys (level).
  SparseMatrix *KK = pdeSys->_KK;
  NumericVector *RES = pdeSys->_RES;

  // Convenience variables to encode the dimension.
  const unsigned  dim = 2;
  const unsigned  DIM = 3;

  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // Extract the solution vector; get solDx positions in the ml_sol object.
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol->GetIndex ("Dx1");
  solDxIndex[1] = mlSol->GetIndex ("Dx2");
  solDxIndex[2] = mlSol->GetIndex ("Dx3");

  // Extract the finite element type for solx.
  unsigned solxType;
  solxType = mlSol->GetSolutionType (solDxIndex[0]);  // get the finite element type for "U"

  // Define solx (only a double).
  std::vector < double > solx[DIM];

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get positions of Y in the ml_sol object.
  unsigned solYIndex[DIM];
  solYIndex[0] = mlSol->GetIndex ("Y1");
  solYIndex[1] = mlSol->GetIndex ("Y2");
  solYIndex[2] = mlSol->GetIndex ("Y3");

  // Extract the finite element type for Y.
  unsigned solYType;
  solYType = mlSol->GetSolutionType (solYIndex[0]);

  // Get positions of Y in the pdeSys object.
  unsigned solYPdeIndex[DIM];
  solYPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("Y1");
  solYPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("Y2");
  solYPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("Y3");

  // Define local Y solution.
  std::vector < adept::adouble > solY[DIM];

  // Get positions of W in the ml_sol object.
  unsigned solWIndex[DIM];
  solWIndex[0] = mlSol->GetIndex ("W1");
  solWIndex[1] = mlSol->GetIndex ("W2");
  solWIndex[2] = mlSol->GetIndex ("W3");

  // Extract the finite element type for W.
  unsigned solWType;
  solWType = mlSol->GetSolutionType (solWIndex[0]);

  // Get positions of W in the pdeSys object.
  unsigned solWPdeIndex[DIM];
  solWPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("W1");
  solWPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("W2");
  solWPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("W3");

  // Define local W solution.
  std::vector < adept::adouble > solW[DIM];

  // Local-to-global pdeSys dofs.
  std::vector< unsigned > SYSDOF;

  // Define local residual vectors.
  vector< double > Res;
  std::vector< adept::adouble > aResY[3];
  std::vector< adept::adouble > aResW[3];

  // Local (column-ordered) Jacobian matrix.
  vector < double > Jac;

  // Zero all entries of global matrix and residual vectors.
  KK->zero();
  RES->zero();

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Number of solution element dofs.
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solxType);
    unsigned nYDofs  = msh->GetElementDofNumber (iel, solYType);
    unsigned nWDofs  = msh->GetElementDofNumber (iel, solWType);

    // Resize solution vectors.
    for (unsigned K = 0; K < DIM; K++) {
      solx[K].resize (nxDofs);
      solY[K].resize (nYDofs);
      solW[K].resize (nWDofs);
    }

    // Resize local arrays.
    SYSDOF.resize (DIM * (nYDofs + nWDofs));
    Res.resize (DIM * (nYDofs + nWDofs));

    for (unsigned K = 0; K < DIM; K++) {
      aResY[K].assign (nYDofs, 0.);  //resize and zet to zero
      aResW[K].assign (nWDofs, 0.);  //resize and zet to zero
    }

    // Local storage of global X mapping and solution.
    for (unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof (i, iel, solxType);
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);

      for (unsigned K = 0; K < DIM; K++) {
        solx[K][i] = (*msh->_topology->_Sol[K]) (iXDof) + (*sol->_Sol[solDxIndex[K]]) (iDDof);
      }
    }

    // Local storage of global Y mapping and solution.
    for (unsigned i = 0; i < nYDofs; i++) {
      unsigned iYDof = msh->GetSolutionDof (i, iel, solYType);

      // Global-to-local mapping between solution node and solution dof.
      for (unsigned K = 0; K < DIM; K++) {
        solY[K][i] = (*sol->_Sol[solYIndex[K]]) (iYDof);

        // Global-to-global mapping between solution node and pdeSys dof.
        SYSDOF[ K * nYDofs + i] = pdeSys->GetSystemDof (solYIndex[K], solYPdeIndex[K], i, iel);
      }
    }


    // Local storage of global W mapping and solution
    for (unsigned i = 0; i < nWDofs; i++) {
      unsigned iWDof = msh->GetSolutionDof (i, iel, solWType);

      // Global-to-local mapping between solution node and solution dof.
      for (unsigned K = 0; K < DIM; K++) {
        solW[K][i] = (*sol->_Sol[solWIndex[K]]) (iWDof);

        // Global-to-global mapping between solution node and pdeSys dof
        SYSDOF[ DIM * nYDofs + K * nWDofs + i] = pdeSys->GetSystemDof (solWIndex[K], solWPdeIndex[K], i, iel);
      }
    }

    // Start a new recording of all the operations involving adept::adouble variables.
    s.new_recording();

    // begin GAUSS POINT LOOP
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solxType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phix_uv[dim]; // first order derivatives in (u,v)

      const double *phiY;  // local test function
      const double *phiY_uv[dim]; // first order derivatives in (u,v)

      const double *phiW;  // local test function

      double weight; // gauss point weight

      //Extract Gauss point weight, test functions, and their partial derivatives.
      // "0" is derivative in u, "1" is derivative in v.
      weight = msh->_finiteElement[ielGeom][solxType]->GetGaussWeight (ig);

      phix = msh->_finiteElement[ielGeom][solxType]->GetPhi (ig);
      phix_uv[0] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDXi (ig);
      phix_uv[1] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDEta (ig);

      phiY = msh->_finiteElement[ielGeom][solYType]->GetPhi (ig);
      phiY_uv[0] = msh->_finiteElement[ielGeom][solYType]->GetDPhiDXi (ig);
      phiY_uv[1] = msh->_finiteElement[ielGeom][solYType]->GetDPhiDEta (ig);

      phiW = msh->_finiteElement[ielGeom][solWType]->GetPhi (ig);

      // Initialize fist order derivatives of X,Y.
      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solY_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      // Initialize solutions X,Y,W at Gauss points.
      double solxg[3] = {0., 0., 0.};
      adept::adouble solYg[3] = {0., 0., 0.};
      adept::adouble solWg[3] = {0., 0., 0.};

      // Compute values of the above at the Gauss points.
      for (unsigned K = 0; K < DIM; K++) {

        for (unsigned i = 0; i < nxDofs; i++) {
          solxg[K] += phix[i] * solx[K][i];
        }

        for (unsigned i = 0; i < nYDofs; i++) {
          solYg[K] += phiY[i] * solY[K][i];
        }

        for (unsigned i = 0; i < nWDofs; i++) {
          solWg[K] += phiW[i] * solW[K][i];
        }

        for (int j = 0; j < dim; j++) {

          for (unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j] += phix_uv[j][i] * solx[K][i];
          }
        }

        for (int j = 0; j < dim; j++) {

          for (unsigned i = 0; i < nYDofs; i++) {
            solY_uv[K][j] += phiY_uv[j][i] * solY[K][i];
          }
        }
      }

      // Compute the metric, metric determinant, and area element.
      double g[dim][dim] = {{0., 0.}, {0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      double Area = weight * sqrt (detg);

      // Compute the components of the unit normal N.
      adept::adouble normal[DIM];
      normal[0] = normalSign * (solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) / sqrt (detg);
      normal[1] = normalSign * (solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);
      normal[2] = normalSign * (solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);

      // Initialize and compute Y.N and |Y|^2 -- essentially 2H and 4H^2.
      adept::adouble YdotN = 0.;
      adept::adouble YdotY = 0.;

      for (unsigned K = 0; K < DIM; K++) {
        YdotN += solYg[K] * normal[K];
        YdotY += solYg[K] * solYg[K];
      }
      double signYdotN = (YdotN.value() >= 0.) ? 1. : -1.;

      adept::adouble sumP1 = 0.;
      for (unsigned p = 0; p < 3; p++) {
        double signP = (P[p] % 2u == 0) ? 1. : signYdotN;
        sumP1 += signP * ap[p] * P[p] * pow (YdotY, (P[p] - 2.) / 2.);
      }

      // Compute the metric inverse.
      double gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;

      // Extra code for debugging.
      // for (unsigned i = 0; i < dim; i++) {
      //   for (unsigned j = 0; j < dim; j++) {
      //     double id = 0.;
      //     for (unsigned j2 = 0; j2 < dim; j2++) {
      //       id +=  g[i][j2] * gi[j2][j];
      //     }
      //     if (i == j && fabs (id - 1.) > 1.0e-10) std::cout << id << " error ";
      //     else if (i != j && fabs (id) > 1.0e-10) std::cout << id << " error ";
      //   }
      // }

      // Compute the "reduced Jacobian" g^{ij}X_j .
      double Jir[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            Jir[i][J] += gi[i][k] * solx_uv[J][k];
          }
        }
      }

      // Initialize and compute tangential gradients of X,Y.
      adept::adouble solx_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solY_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            solx_Xtan[I][J] += solx_uv[I][k] * Jir[k][J];
            solY_Xtan[I][J] += solY_uv[I][k] * Jir[k][J];
          }
        }
      }

      // Initiailize and compute tangential gradients of test functions.
      std::vector < double > phix_Xtan[DIM];
      std::vector < double > phiY_Xtan[DIM];

      for (unsigned J = 0; J < DIM; J++) {
        phix_Xtan[J].assign (nxDofs, 0.);
        phiY_Xtan[J].assign (nYDofs, 0.);

        for (unsigned inode  = 0; inode < nxDofs; inode++) {
          for (unsigned k = 0; k < dim; k++) {
            phix_Xtan[J][inode] += phix_uv[k][inode] * Jir[k][J];
          }
        }

        for (unsigned inode  = 0; inode < nYDofs; inode++) {
          for (unsigned k = 0; k < dim; k++) {
            phiY_Xtan[J][inode] += phiY_uv[k][inode] * Jir[k][J];
          }
        }
      }

      // Implement the equations for Y and W.
      for (unsigned K = 0; K < DIM; K++) {

        for (unsigned i = 0; i < nYDofs; i++) {
          adept::adouble term1 = 0.;
          adept::adouble term2 = 0.;

          for (unsigned J = 0; J < DIM; J++) {
            term1 +=  solx_Xtan[K][J] * phiY_Xtan[J][i];
            term2 +=  solY_Xtan[K][J] * phiY_Xtan[J][i];
          }
          // This is a trick to smooth the initial data.
          aResY[K][i] += (solYg[K] * phiY[i] + delta * term2 + term1) * Area;
        }

        //  Solve the equation W = |Y|^{p-2}Y .
        for (unsigned i = 0; i < nWDofs; i++) {
          aResW[K][i] += (solWg[K] - sumP1 * solYg[K]) * phiW[i] * weight;
        }
      }
    } // end GAUSS POINT LOOP

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    //copy the value of the adept::adouble aRes in double Res and store

    for (int K = 0; K < DIM; K++) {
      for (int i = 0; i < nYDofs; i++) {
        Res[ K * nYDofs + i] = -aResY[K][i].value();
      }
    }
    for (int K = 0; K < DIM; K++) {
      for (int i = 0; i < nWDofs; i++) {
        Res[ DIM * nYDofs + K * nWDofs + i] = -aResW[K][i].value();
      }
    }

    RES->add_vector_blocked (Res, SYSDOF);

    // Resize Jacobian.
    Jac.resize (DIM * (nYDofs + nWDofs) * DIM * (nYDofs + nWDofs));

    // Define the dependent variables.
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResY[K][0], nYDofs);
    }
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResW[K][0], nWDofs);
    }

    // Define the independent variables.
    for (int K = 0; K < DIM; K++) {
      s.independent (&solY[K][0], nYDofs);
    }
    for (int K = 0; K < DIM; K++) {
      s.independent (&solW[K][0], nWDofs);
    }

    // Extract the Jacobian matrix (ordered by row).
    s.jacobian (&Jac[0], true);

    KK->add_matrix_blocked (Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } // end ELEMENT LOOP for each process.

  RES->close();
  KK->close();

} // end AssembleInit.


// Building the P-Willmore assembly function.
void AssemblePWillmore (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  // Call the adept stack object.
  adept::Stack& s = FemusInit::_adeptStack;

  // Extract pointers to the several objects that we are going to use.
  TransientNonlinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientNonlinearImplicitSystem> ("PWillmore");   // pointer to the linear implicit system named "Poisson"

  // Define level and time variable.
  double dt = mlPdeSys->GetIntervalTime();
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Point to the mesh and element objects.
  Mesh *msh = ml_prob._ml_msh->GetLevel (level);
  elem *el = msh->el;

  // Point to mlSol, solution (level), and equation (level) objects.
  MultiLevelSolution *mlSol = ml_prob._ml_sol;
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel (level);
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level];

  // Point to the global stiffness mtx and residual vectors in pdeSys (level).
  SparseMatrix *KK = pdeSys->_KK;
  NumericVector *RES = pdeSys->_RES;

  // Convenience variables to encode the dimension.
  const unsigned dim = 2;
  const unsigned DIM = 3;

  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // Extract the solution vector; get solDx positions in the ml_sol object.
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol->GetIndex ("Dx1");
  solDxIndex[1] = mlSol->GetIndex ("Dx2");
  solDxIndex[2] = mlSol->GetIndex ("Dx3");

  // Extract the finite element type for solx.
  unsigned solxType;
  solxType = mlSol->GetSolutionType (solDxIndex[0]);

  // Get positions of solDx in the pdeSys object.
  unsigned solDxPdeIndex[DIM];
  solDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("Dx1");
  solDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("Dx2");
  solDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("Dx3");

  // Define solx and solxOld.
  std::vector < adept::adouble > solx[DIM];
  std::vector < double > solxOld[DIM];

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get positions of Y in the ml_sol object.
  unsigned solYIndex[DIM];
  solYIndex[0] = mlSol->GetIndex ("Y1");
  solYIndex[1] = mlSol->GetIndex ("Y2");
  solYIndex[2] = mlSol->GetIndex ("Y3");

  // Extract the finite element type for Y.
  unsigned solYType;
  solYType = mlSol->GetSolutionType (solYIndex[0]);

  // Get positions of Y in the pdeSys object.
  unsigned solYPdeIndex[DIM];
  solYPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("Y1");
  solYPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("Y2");
  solYPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("Y3");

  // Define solY and solYOld.
  std::vector < adept::adouble > solY[DIM];
  std::vector < double > solYOld[DIM];

  // Get positions of W in the ml_sol object.
  unsigned solWIndex[DIM];
  solWIndex[0] = mlSol->GetIndex ("W1");
  solWIndex[1] = mlSol->GetIndex ("W2");
  solWIndex[2] = mlSol->GetIndex ("W3");

  // Extract the finite element type for W.
  unsigned solWType;
  solWType = mlSol->GetSolutionType (solWIndex[0]);

  // Get positions of W in the pdeSys object.
  unsigned solWPdeIndex[DIM];
  solWPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("W1");
  solWPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("W2");
  solWPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("W3");

  // Define local W, WOld solutions.
  std::vector < adept::adouble > solW[DIM];
  std::vector < double > solWOld[DIM];

  // Local-to-global pdeSys dofs.
  std::vector< unsigned > SYSDOF;

  // Define local residual vectors.
  vector < double > Res;
  std::vector< adept::adouble > aResx[3];
  std::vector< adept::adouble > aResY[3];
  std::vector< adept::adouble > aResW[3];

  // Local (column-ordered) Jacobian matrix
  vector < double > Jac;

  //MatSetOption ( ( static_cast<PetscMatrix*> ( KK ) )->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE );

  KK->zero();  // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual


  // Setting up solLambda1 (vol) and solLambda2 (area).
  unsigned solLambaPdeIndex;

  adept::adouble solLambda1 = 0.;
  adept::adouble aResLambda1;
  unsigned lambda1PdeDof;

  adept::adouble solLambda2 = 0.;
  adept::adouble aResLambda2;
  unsigned lambda2PdeDof;

  if (volumeConstraint || areaConstraint) {
    unsigned solLambdaIndex;
    solLambdaIndex = mlSol->GetIndex ("Lambda");
    solLambaPdeIndex = mlPdeSys->GetSolPdeIndex ("Lambda");

    if (volumeConstraint) {
      double lambda1;
      if (iproc == 0) {
        lambda1 = (*sol->_Sol[solLambdaIndex]) (0); // global to local solution
        lambda1PdeDof = pdeSys->GetSystemDof (solLambdaIndex, solLambaPdeIndex, 0, 0);
      }
      MPI_Bcast (&lambda1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast (&lambda1PdeDof, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      solLambda1 = lambda1;
    }

    if (areaConstraint) {
      double lambda2;
      if (iproc == 0) {
        lambda2 = (*sol->_Sol[solLambdaIndex]) (volumeConstraint); // global to local solution
        lambda2PdeDof = pdeSys->GetSystemDof (solLambdaIndex, solLambaPdeIndex, 0, volumeConstraint);
      }
      MPI_Bcast (&lambda2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast (&lambda2PdeDof, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      solLambda2 = lambda2;
    }

    std::vector < double > value (2);
    std::vector < int > row (1);
    std::vector < int > columns (2);
    value[0] = 1;
    value[1] = -1;
    columns[1] = (volumeConstraint) ? lambda1PdeDof : lambda2PdeDof;

    // For equations other than Lagrange multiplier:
    for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
      if (iel > volumeConstraint * areaConstraint) {
        row[0] = pdeSys->GetSystemDof (solLambdaIndex, solLambaPdeIndex, 0, iel);
        columns[0] = row[0];
        KK->add_matrix_blocked (value, row, columns);
      }
    }
  }

  // Initialize area, volume, P-Willmore energy.
  double surface = 0.;
  double volume = 0.;
  double energy = 0.;



  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Number of solution element dofs.
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solxType);
    unsigned nYDofs  = msh->GetElementDofNumber (iel, solYType);
    unsigned nWDofs  = msh->GetElementDofNumber (iel, solWType);

    // Resize solution vectors.
    for (unsigned K = 0; K < DIM; K++) {
      solx[K].resize (nxDofs);
      solxOld[K].resize (nxDofs);
      solY[K].resize (nYDofs);
      solYOld[K].resize (nYDofs);
      solW[K].resize (nWDofs);
      solWOld[K].resize (nWDofs);
    }

    // Convenience variable for keeping track of problem size.
    unsigned sizeAll = DIM * (nxDofs + nYDofs +  nWDofs) + volumeConstraint + areaConstraint;

    // Resize local arrays.
    SYSDOF.resize (sizeAll);
    Res.resize (sizeAll);

    for (unsigned K = 0; K < DIM; K++) {
      aResx[K].assign (nxDofs, 0.);  //resize and set to zero
      aResY[K].assign (nYDofs, 0.);  //resize and set to zero
      aResW[K].assign (nWDofs, 0.);  //resize and zet to zero
    }
    aResLambda1 = 0.;
    aResLambda2 = 0.;

    // Loop which handles local storage of global mapping and solution X.
    for (unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof (i, iel, solxType);
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);

      for (unsigned K = 0; K < DIM; K++) {
        solxOld[K][i] = (*msh->_topology->_Sol[K]) (iXDof) + (*sol->_SolOld[solDxIndex[K]]) (iDDof);
        solx[K][i] = (*msh->_topology->_Sol[K]) (iXDof) + (*sol->_Sol[solDxIndex[K]]) (iDDof);

        // Global-to-global mapping between solution node and pdeSys dof.
        SYSDOF[K * nxDofs + i] = pdeSys->GetSystemDof (solDxIndex[K], solDxPdeIndex[K], i, iel);
      }
    }

    // Loop which handles local storage of global mapping and solution Y.
    for (unsigned i = 0; i < nYDofs; i++) {

      // Global-to-local mapping between solution node and solution dof.
      unsigned iYDof = msh->GetSolutionDof (i, iel, solYType);
      for (unsigned K = 0; K < DIM; K++) {

        // Global-to-local solutions.
        solYOld[K][i] = (*sol->_SolOld[solYIndex[K]]) (iYDof);
        solY[K][i] = (*sol->_Sol[solYIndex[K]]) (iYDof);

        // Global-to-global mapping between solution node and pdeSys dof.
        SYSDOF[DIM * nxDofs + K * nYDofs + i] =
          pdeSys->GetSystemDof (solYIndex[K], solYPdeIndex[K], i, iel);
      }
    }

    // Loop which handles local storage of global mapping and solution W.
    for (unsigned i = 0; i < nWDofs; i++) {

      // Global-to-local mapping between solution node and solution dof.
      unsigned iWDof = msh->GetSolutionDof (i, iel, solWType);
      for (unsigned K = 0; K < DIM; K++) {

        // Global-to-local solutions.
        solWOld[K][i] = (*sol->_SolOld[solWIndex[K]]) (iWDof);
        solW[K][i] = (*sol->_Sol[solWIndex[K]]) (iWDof);

        // Global-to-global mapping between solution node and pdeSys dof.
        SYSDOF[DIM * (nxDofs + nYDofs) + K * nWDofs + i] =
          pdeSys->GetSystemDof (solWIndex[K], solWPdeIndex[K], i, iel);
      }
    }

    // Conditions for local storage of global Lagrange multipliers.
    if (volumeConstraint) {
      SYSDOF[sizeAll - 1u - areaConstraint ] = lambda1PdeDof;
    }

    if (areaConstraint) {
      SYSDOF[sizeAll - 1u ] = lambda2PdeDof;
    }

    // Start a new recording of all the operations involving adept variables.
    s.new_recording();

    // begin GAUSS POINT LOOP
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solxType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phix_uv[dim]; // first order derivatives in (u,v)

      const double *phiY;  // local test function

      const double *phiW;  // local test function
      const double *phiW_uv[dim]; // first order derivatives in (u,v)

      double weight; // gauss point weight

      //Extract Gauss point weight, test functions, and their partial derivatives.
      // "0" is derivative in u, "1" is derivative in v.
      weight = msh->_finiteElement[ielGeom][solxType]->GetGaussWeight (ig);

      phix = msh->_finiteElement[ielGeom][solxType]->GetPhi (ig);
      phix_uv[0] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDXi (ig);
      phix_uv[1] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDEta (ig);

      phiY = msh->_finiteElement[ielGeom][solYType]->GetPhi (ig);

      phiW = msh->_finiteElement[ielGeom][solWType]->GetPhi (ig);
      phiW_uv[0] = msh->_finiteElement[ielGeom][solWType]->GetDPhiDXi (ig);
      phiW_uv[1] = msh->_finiteElement[ielGeom][solWType]->GetDPhiDEta (ig);

      // Initialize quantities xNew, xOld, Y, W at the Gauss points.
      adept::adouble solxNewg[3] = {0., 0., 0.};
      double solxOldg[3] = {0., 0., 0.};
      adept::adouble solYg[3] = {0., 0., 0.};
      adept::adouble solWg[3] = {0., 0., 0.};

      // Initialize derivatives of x and W (new, middle, old) at the Gauss points.
      adept::adouble solxNew_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solWNew_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solW_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solxOld_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solWOld_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};


      for (unsigned K = 0; K < DIM; K++) {

        for (unsigned i = 0; i < nxDofs; i++) {
          solxNewg[K] += phix[i] * solx[K][i];
          solxOldg[K] += phix[i] * solxOld[K][i];
        }

        for (unsigned i = 0; i < nYDofs; i++) {
          //solYNewg[K] += phiY[i] * solY[K][i];
          solYg[K] += phiY[i] * 0.5 * (solYOld[K][i] + solY[K][i]);
        }

        for (unsigned i = 0; i < nWDofs; i++) {
          solWg[K] += phiW[i] * 0.5 * (solWOld[K][i] + solW[K][i]);
          //solWOldg[K] += phiW[i] * solWOld[K][i];
        }

        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nxDofs; i++) {
            solxNew_uv[K][j]    += phix_uv[j][i] * solx[K][i];
            solx_uv[K][j]    += phix_uv[j][i] * 0.5 * (solx[K][i] + solxOld[K][i]);
            solxOld_uv[K][j] += phix_uv[j][i] * solxOld[K][i];
          }
        }

        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nWDofs; i++) {
            solWNew_uv[K][j] += phiW_uv[j][i] * solW[K][i];
            solW_uv[K][j] += phiW_uv[j][i] * 0.5 * (solW[K][i] + solWOld[K][i]);
            solWOld_uv[K][j] += phiW_uv[j][i] * solWOld[K][i];
          }
        }
      }

      // Computing the metric, metric determinant, and area element.
      adept::adouble g[dim][dim] = {{0., 0.}, {0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      adept::adouble detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      adept::adouble Area = weight * sqrt (detg);

      // Computing the unit normal vector N.
      adept::adouble normal[DIM];
      normal[0] = normalSign * (solx_uv[1][0] * solx_uv[2][1]
                                - solx_uv[2][0] * solx_uv[1][1]) / sqrt (detg);
      normal[1] = normalSign * (solx_uv[2][0] * solx_uv[0][1]
                                - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);
      normal[2] = normalSign * (solx_uv[0][0] * solx_uv[1][1]
                                - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);

      // Computing Y.N and |Y|^2, which are essentially 2H and 4H^2.
      adept::adouble YdotN = 0.;
      adept::adouble YdotY = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        YdotN += solYg[K] * normal[K];
        YdotY += solYg[K] * solYg[K];
      }
      double signYdotN = (YdotN.value() >= 0.) ? 1. : -1.;

      // Some necessary quantities when working with polynomials.
      adept::adouble sumP1 = 0.;
      adept::adouble sumP2 = 0.;
      adept::adouble sumP3 = 0.;
      for (unsigned p = 0; p < 3; p++) {
        double signP = (P[p] % 2u == 0) ? 1. : signYdotN;
        sumP1 += signP * ap[p] * P[p] * pow (YdotY, (P[p] - 2.) / 2.);
        sumP2 += signP * ap[p] * (1. - P[p]) * pow (YdotY , P[p] / 2.);
        //sumP2 += signP * (ap[p] - ap[p] * P[p]) * pow (YdotY , P[p]/2.);
        sumP3 += signP * ap[p] * pow (YdotY, P[p] / 2.);
      }

      // Computing the metric inverse
      adept::adouble gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;

      // Computing the "reduced Jacobian" g^{ij}X_j .
      adept::adouble Jir[dim][DIM] = {{0., 0., 0.}, {0., 0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            Jir[i][J] += gi[i][k] * solx_uv[J][k];
          }
        }
      }

      // Initializing tangential gradients of X and W (new, middle, old).
      adept::adouble solxNew_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solx_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solxOld_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

      adept::adouble solWNew_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solW_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solWOld_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

      // Computing tangential gradients defined above.
      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            solxNew_Xtan[I][J] += solxNew_uv[I][k] * Jir[k][J];
            solx_Xtan[I][J] += solx_uv[I][k] * Jir[k][J];
            solxOld_Xtan[I][J] += solxOld_uv[I][k] * Jir[k][J];

            solWNew_Xtan[I][J] += solWNew_uv[I][k] * Jir[k][J];
            solW_Xtan[I][J] += solW_uv[I][k] * Jir[k][J];
            solWOld_Xtan[I][J] += solWOld_uv[I][k] * Jir[k][J];
          }
        }
      }

      // Define and compute gradients of test functions for X and W.
      std::vector < adept::adouble > phiW_Xtan[DIM];
      std::vector < adept::adouble > phix_Xtan[DIM];

      for (unsigned J = 0; J < DIM; J++) {
        phix_Xtan[J].assign (nxDofs, 0.);
        phiW_Xtan[J].assign (nWDofs, 0.);

        for (unsigned inode  = 0; inode < nxDofs; inode++) {
          for (unsigned k = 0; k < dim; k++) {
            phix_Xtan[J][inode] += phix_uv[k][inode] * Jir[k][J];
          }
        }

        for (unsigned inode  = 0; inode < nWDofs; inode++) {
          for (unsigned k = 0; k < dim; k++) {
            phiW_Xtan[J][inode] += phiW_uv[k][inode] * Jir[k][J];
          }
        }
      }

      // Implement the curvature equation Y = \Delta X .
      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;

          for (unsigned J = 0; J < DIM; J++) {
            term1 +=  solxNew_Xtan[K][J] * phix_Xtan[J][i];
          }
          aResx[K][i] += (solYg[K] * phix[i] + term1) * Area;
        }

        // Implement the equation relating Y and W.
        for (unsigned i = 0; i < nWDofs; i++) {
          aResY[K][i] += (solWg[K] - sumP1 * solYg[K]) * phiY[i] * Area;
        }

        // Implement the main P-Willmore equation.
        for (unsigned i = 0; i < nWDofs; i++) {
          adept::adouble term0 = 0.;
          adept::adouble term1 = 0.;
          adept::adouble term2 = 0.;
          adept::adouble term3 = 0.;

          adept::adouble term5 = 0.;

          for (unsigned J = 0; J < DIM; J++) {
            term0 +=  solWNew_Xtan[K][J] * phiW_Xtan[J][i];
            term1 +=  solx_Xtan[K][J] * phiW_Xtan[J][i];
            term2 +=  solW_Xtan[J][J];

            adept::adouble term4 = 0.;
            for (unsigned L = 0; L < DIM; L++) {
              term4 += solxOld_Xtan[L][J] * solWOld_Xtan[L][K] + solxOld_Xtan[L][K] * solWOld_Xtan[L][J];
            }
            term3 += phiW_Xtan[J][i] * term4;
          }

          // P-Willmore equation
          aResW[K][i] += ( ( (solLambda1 /*- YdotN * solLambda2*/) * normal[K] + (solxNewg[K] - solxOldg[K])  / dt) * phiW[i]
                           + solLambda2 * term1
                           - term0
                           + sumP2 * term1
                           - term2 * phiW_Xtan[K][i]
                           + term3
                         ) * Area;
        }

        // Lagrange multiplier (volume) equation Dx.N = 0.
        if (volumeConstraint) {
          aResLambda1 += ( (solxNewg[K] - solxOldg[K]) * normal[K]) * Area;
        }

        // Lagrange multiplier (area) equation.
        if (areaConstraint) {
          //aResLambda2 += ( -YdotN * (solxNewg[K] - solxOldg[K]) * normal[K]) * Area;
          adept::adouble term1t = 0.;
          for (unsigned J = 0; J < DIM; J++) {
            term1t +=  solx_Xtan[K][J] * (solxNew_Xtan[K][J] - solxOld_Xtan[K][J]) ;
          }
          aResLambda2 += term1t * Area;

        }
      }

      // Compute new surface area, volume, and P-Willmore energy.
      for (unsigned K = 0; K < DIM; K++) {
        surface += Area.value();
      }
      for (unsigned K = 0; K < DIM; K++) {
        volume += normalSign * (solxNewg[K].value()  * normal[K].value()) * Area.value();
      }
      energy += sumP3.value() * Area.value();

    } // end GAUSS POINT LOOP.

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    //copy the value of the adept::adoube aRes in double Res and store

    for (int K = 0; K < DIM; K++) {
      for (int i = 0; i < nxDofs; i++) {
        Res[ K * nxDofs + i] = -aResx[K][i].value();
      }
    }

    for (int K = 0; K < DIM; K++) {
      for (int i = 0; i < nYDofs; i++) {
        Res[DIM * nxDofs + K * nYDofs + i] = -aResY[K][i].value();
      }
    }

    for (int K = 0; K < DIM; K++) {
      for (int i = 0; i < nWDofs; i++) {
        Res[DIM * (nxDofs + nYDofs) + K * nWDofs + i] = -aResW[K][i].value();
      }
    }

    if (volumeConstraint) {
      Res[sizeAll - 1u - areaConstraint] = - aResLambda1.value();
    }

    if (areaConstraint) {
      Res[sizeAll - 1u] = - aResLambda2.value();
    }

    RES->add_vector_blocked (Res, SYSDOF);

    // Resize Jacobian.
    Jac.resize (sizeAll * sizeAll);

    // Define the dependent variables.
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResx[K][0], nxDofs);
    }
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResY[K][0], nYDofs);
    }
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResW[K][0], nWDofs);
    }
    if (volumeConstraint) {
      s.dependent (&aResLambda1, 1);
    }
    if (areaConstraint) {
      s.dependent (&aResLambda2, 1);
    }

    // Define the independent variables.
    for (int K = 0; K < DIM; K++) {
      s.independent (&solx[K][0], nxDofs);
    }
    for (int K = 0; K < DIM; K++) {
      s.independent (&solY[K][0], nYDofs);
    }
    for (int K = 0; K < DIM; K++) {
      s.independent (&solW[K][0], nWDofs);
    }
    if (volumeConstraint) {
      s.independent (&solLambda1, 1);
    }
    if (areaConstraint) {
      s.independent (&solLambda2, 1);
    }

    // Get the Jacobian matrix (ordered by row).
    s.jacobian (&Jac[0], true);

    KK->add_matrix_blocked (Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } // End ELEMENT LOOP for each process.

  RES->close();
  KK->close();

  // Get data from each process running in parallel.
  double surfaceAll;
  MPI_Reduce (&surface, &surfaceAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "SURFACE = " << surfaceAll << std::endl;

  double volumeAll;
  MPI_Reduce (&volume, &volumeAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "VOLUME = " << volumeAll << std::endl;

  double energyAll;
  MPI_Reduce (&energy, &energyAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "ENERGY = " << energyAll << std::endl;

//   VecView ( (static_cast<PetscVector*> (RES))->vec(),  PETSC_VIEWER_STDOUT_SELF);
//   MatView ( (static_cast<PetscMatrix*> (KK))->mat(), PETSC_VIEWER_STDOUT_SELF);

//     PetscViewer    viewer;
//     PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//     PetscObjectSetName ( (PetscObject) viewer, "PWilmore matrix");
//     PetscViewerPushFormat (viewer, PETSC_VIEWER_DRAW_LG);
//     MatView ( (static_cast<PetscMatrix*> (KK))->mat(), viewer);
//     double a;
//     std::cin >> a;

} // end AssemblePWillmore.


// Building the Conformal Minimization system.
void AssembleConformalMinimization (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  Extract pointers to the several objects that we are going to use.
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem> ("nProj");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Pointers to the mesh (level) object and elem object in mesh (level).
  Mesh *msh = ml_prob._ml_msh->GetLevel (level);
  elem *el = msh->el;

  // Pointers to the multilevel solution, solution (level) and equation (level).
  MultiLevelSolution *mlSol = ml_prob._ml_sol;
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel (level);
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level];

  // Pointers to global stiffness matrix and residual vector in pdeSys (level).
  SparseMatrix *KK = pdeSys->_KK;
  NumericVector *RES = pdeSys->_RES;

  // Convenience variables to keep track of the dimension.
  const unsigned  dim = 2;
  const unsigned  DIM = 3;

  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // Setting the reference elements to be equilateral triangles.
  std::vector < std::vector < double > > xT (2);
  xT[0].resize (3);
  xT[0][0] = -0.5;
  xT[0][1] = 0.5;
  xT[0][2] = 0.;

  xT[1].resize (3);
  xT[1][0] = 0.;
  xT[1][1] = 0.;
  xT[1][2] = sqrt (3.) / 2.;

  std::vector<double> phi_uv0;
  std::vector<double> phi_uv1;

  std::vector< double > stdVectorPhi;
  std::vector< double > stdVectorPhi_uv;

  // Extract positions of Dx in ml_sol object.
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol->GetIndex ("Dx1");
  solDxIndex[1] = mlSol->GetIndex ("Dx2");
  solDxIndex[2] = mlSol->GetIndex ("Dx3");

  // Extract finite element type for the solution.
  unsigned solType;
  solType = mlSol->GetSolutionType (solDxIndex[0]);

  // Local solution vectors for X, Dx, Xhat, XC.
  std::vector < double > solx[DIM];
  std::vector < double > solDx[DIM];
  std::vector < double > xhat[DIM];
  std::vector < double > xc[DIM];

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get the poitions of Y in the ml_sol object.
  unsigned solNDxIndex[DIM];
  solNDxIndex[0] = mlSol->GetIndex ("nDx1");
  solNDxIndex[1] = mlSol->GetIndex ("nDx2");
  solNDxIndex[2] = mlSol->GetIndex ("nDx3");

  // Get the positions of Y in the pdeSys object.
  unsigned solNDxPdeIndex[DIM];
  solNDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("nDx1");
  solNDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("nDx2");
  solNDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("nDx3");

  // Local solution vectors for Nx and NDx.
  std::vector < adept::adouble > solNDx[DIM];
  std::vector < adept::adouble > solNx[DIM];

  // Get the position of "Lambda1" in the ml_sol object.
  unsigned solLIndex;
  solLIndex = mlSol->GetIndex ("Lambda1");

  // Get the finite element type for "Lambda1".
  unsigned solLType;
  solLType = mlSol->GetSolutionType (solLIndex);

  // Get the position of "Lambda1" in the pdeSys object.
  unsigned solLPdeIndex;
  solLPdeIndex = mlPdeSys->GetSolPdeIndex ("Lambda1");

  // Local Lambda1 solution.
  std::vector < adept::adouble > solL;

  // Local-to-global pdeSys dofs.
  std::vector < int > SYSDOF;

  // Local residual vectors.
  vector< double > Res;
  std::vector< adept::adouble > aResNDx[3];
  std::vector< adept::adouble > aResL;

  // Local Jacobian matrix (ordered by column).
  vector < double > Jac;

  KK->zero();  // Zero all the entries of the Global Matrix
  RES->zero(); // Zero all the entries of the Global Residual

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Numer of solution element dofs.
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solType);
    unsigned nLDofs  = msh->GetElementDofNumber (iel, solLType);

    // Resize local arrays.
    for (unsigned K = 0; K < DIM; K++) {

      xhat[K].resize (nxDofs);

      solDx[K].resize (nxDofs);
      solx[K].resize (nxDofs);

      solNDx[K].resize (nxDofs);
      solNx[K].resize (nxDofs);

      solL.resize (nLDofs);
    }

    // Resize local arrays
    SYSDOF.resize (DIM * nxDofs + nLDofs);
    Res.resize (DIM * nxDofs + nLDofs);

    for (unsigned K = 0; K < DIM; K++) {
      aResNDx[K].assign (nxDofs, 0.);
    }
    aResL.assign (nLDofs, 0.);

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof (i, iel, solType);
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);

      for (unsigned K = 0; K < DIM; K++) {
        xhat[K][i] = (*msh->_topology->_Sol[K]) (iXDof);
        solDx[K][i] = (*sol->_Sol[solDxIndex[K]]) (iDDof);
        solx[K][i] = xhat[K][i] + solDx[K][i];
        solNDx[K][i] = (*sol->_Sol[solNDxIndex[K]]) (iDDof);

        // Global-to-global mapping between NDx solution node and pdeSys dof.
        SYSDOF[ K * nxDofs + i] =
          pdeSys->GetSystemDof (solNDxIndex[K], solNDxPdeIndex[K], i, iel);
      }
    }

    // Local storage of global mapping and solution.
    for (unsigned i = 0; i < nLDofs; i++) {

      // Global-to-local mapping between Lambda solution node and solution dof.
      unsigned iLDof = msh->GetSolutionDof (i, iel, solLType);
      solL[i] = (*sol->_Sol[solLIndex]) (iLDof);

      // Global-to-global mapping between Lambda solution node and pdeSys dof.
      SYSDOF[DIM * nxDofs + i] =
        pdeSys->GetSystemDof (solLIndex, solLPdeIndex, i, iel);
    }

    // start a new recording of all the operations involving adept variables.
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phiL;  // local test function
      const double *phix_uv[dim]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // Get Gauss point weight, test function, and first order derivatives.
      if (ielGeom == QUAD) {
        phix = msh->_finiteElement[ielGeom][solType]->GetPhi (ig);
        phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi (ig);

        phix_uv[0] = msh->_finiteElement[ielGeom][solType]->GetDPhiDXi (ig);
        phix_uv[1] = msh->_finiteElement[ielGeom][solType]->GetDPhiDEta (ig);

        weight = msh->_finiteElement[ielGeom][solType]->GetGaussWeight (ig);
      }

      // Special adjustments for triangles.
      else {
        msh->_finiteElement[ielGeom][solType]->Jacobian (xT, ig, weight, stdVectorPhi, stdVectorPhi_uv);

        phix = &stdVectorPhi[0];
///////// WHAT ABOUT PHIL? ////////////////////

        phi_uv0.resize(nxDofs);
        phi_uv1.resize(nxDofs);


        for (unsigned i = 0; i < nxDofs; i++) {
          phi_uv0[i] = stdVectorPhi_uv[i * dim];
          phi_uv1[i] = stdVectorPhi_uv[i * dim + 1];
        }

        phix_uv[0] = &phi_uv0[0];
        phix_uv[1] = &phi_uv1[0];

        phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi (ig);

      }

      // Initialize and compute values of x, Dx, NDx, x_uv at the Gauss points.
      double solDxg[3] = {0., 0., 0.};
      adept::adouble solNDxg[3] = {0., 0., 0.};

      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solNx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          solDxg[K] += phix[i] * solDx[K][i];
          solNDxg[K] += phix[i] * solNDx[K][i];
        }
        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * solx[K][i];
            solNx_uv[K][j]   += phix_uv[j][i] * (xhat[K][i] + solNDx[K][i]);
          }
        }
      }

      ///////// ADDED THIS /////////
      adept::adouble solLg = 0.;
      for (unsigned i = 0; i < nLDofs; i++) {
        solLg += phiL[i] * solL[i];
      }

      // Compute the metric, metric determinant, and area element.
      double g[dim][dim] = {{0., 0.}, {0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      double Area = weight * sqrt (detg);
      double Area2 = weight; // Trick to give equal weight to each element.

      // Compute components of the unit normal N.
      double normal[DIM];
      normal[0] = (solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) / sqrt (detg);
      normal[1] = (solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);
      normal[2] = (solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);

      // Discretize the equation \delta CD = 0 on the basis d/du, d/dv.
      adept::adouble V[DIM];
      V[0] = solNx_uv[0][1] - normal[1] * solNx_uv[2][0] + normal[2] * solNx_uv[1][0];
      V[1] = solNx_uv[1][1] - normal[2] * solNx_uv[0][0] + normal[0] * solNx_uv[2][0];
      V[2] = solNx_uv[2][1] - normal[0] * solNx_uv[1][0] + normal[1] * solNx_uv[0][0];

      adept::adouble W[DIM];
      W[0] = solNx_uv[0][0] + normal[1] * solNx_uv[2][1] - normal[2] * solNx_uv[1][1];
      W[1] = solNx_uv[1][0] + normal[2] * solNx_uv[0][1] - normal[0] * solNx_uv[2][1];
      W[2] = solNx_uv[2][0] + normal[0] * solNx_uv[1][1] - normal[1] * solNx_uv[0][1];

      adept::adouble M[DIM][dim];
      M[0][0] = W[0] - normal[2] * V[1] + normal[1] * V[2];
      M[1][0] = W[1] - normal[0] * V[2] + normal[2] * V[0];
      M[2][0] = W[2] - normal[1] * V[0] + normal[0] * V[1];

      M[0][1] = V[0] + normal[2] * W[1] - normal[1] * W[2];
      M[1][1] = V[1] + normal[0] * W[2] - normal[2] * W[0];
      M[2][1] = V[2] + normal[1] * W[0] - normal[0] * W[1];

      // Compute new X minus old X dot N, for "reparametrization".
      adept::adouble DnXmDxdotN = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        DnXmDxdotN += (solDxg[K] - solNDxg[K]) * normal[K];
      }

      // Compute the metric inverse.
      double gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;

      // Compute the "reduced Jacobian" g^{ij}X_j .
      double Jir[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            Jir[i][J] += gi[i][k] * solx_uv[J][k];
          }
        }
      }

      // Implement the Conformal Minimization equations.
      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;

          for (unsigned j = 0; j < dim; j++) {
            term1 +=  M[K][j] * phix_uv[j][i];
            //term1 +=  X_uv[K][j] * phix_uv[j][i];
          }

          // Conformal energy equation (with trick).
          aResNDx[K][i] += term1 * Area2
                           + timederiv * (solNDxg[K] - solDxg[K]) * phix[i] * Area2
                           + solL[0] * phix[i] * normal[K] * Area;
        }
      }

      // Lagrange multiplier equation (with trick).
      for (unsigned i = 0; i < nLDofs; i++) {
        aResL[i] += phiL[i] * (DnXmDxdotN + eps * solL[i]) * Area; // no2
      }
      //aResL[0] += (DnXmDxdotN + eps * solL[0]) * Area;

    } // end GAUSS POINT LOOP

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    //copy the value of the adept::adoube aRes in double Res and store

    for (int K = 0; K < DIM; K++) {
      for (int i = 0; i < nxDofs; i++) {
        Res[ K * nxDofs + i] = -aResNDx[K][i].value();
      }
    }

    for (int i = 0; i < nLDofs; i++) {
      Res[DIM * nxDofs + i] = - aResL[i].value();
    }
    RES->add_vector_blocked (Res, SYSDOF);

    // Resize Jacobian.
    Jac.resize ( (DIM * nxDofs + nLDofs) * (DIM * nxDofs + nLDofs));

    // Define the dependent variables.
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResNDx[K][0], nxDofs);
    }
    s.dependent (&aResL[0], nLDofs);


    // Define the independent variables.
    for (int K = 0; K < DIM; K++) {
      s.independent (&solNDx[K][0], nxDofs);
    }
    s.independent (&solL[0], nLDofs);

    // Get the jacobian matrix (ordered by row).
    s.jacobian (&Jac[0], true);

    KK->add_matrix_blocked (Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } //end ELEMENT LOOP for each process.

  RES->close();
  KK->close();

} // end AssembleConformalMinimization.


void AssembleShearMinimization (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem> ("nProj");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh *msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem *el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution *mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object

  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix *KK = pdeSys->_KK;  // pointer to the global stiffness matrix object in pdeSys (level)
  NumericVector *RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = 2;
  const unsigned  DIM = 3;
  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol->GetIndex ("Dx1"); // get the position of "DX" in the ml_sol object
  solDxIndex[1] = mlSol->GetIndex ("Dx2"); // get the position of "DY" in the ml_sol object
  solDxIndex[2] = mlSol->GetIndex ("Dx3"); // get the position of "DZ" in the ml_sol object
  unsigned solType;
  solType = mlSol->GetSolutionType (solDxIndex[0]);  // get the finite element type for "U"
  std::vector < double > solx[DIM];  // surface coordinates
  std::vector < double > solDx[DIM];  // surface coordinates
  std::vector < double > solOldDx[DIM];  // surface coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  unsigned solNDxIndex[DIM];
  solNDxIndex[0] = mlSol->GetIndex ("nDx1");   // get the position of "Y1" in the ml_sol object
  solNDxIndex[1] = mlSol->GetIndex ("nDx2");   // get the position of "Y2" in the ml_sol object
  solNDxIndex[2] = mlSol->GetIndex ("nDx3");   // get the position of "Y3" in the ml_sol object
  unsigned solNDxPdeIndex[DIM];
  solNDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("nDx1");   // get the position of "Y1" in the pdeSys object
  solNDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("nDx2");   // get the position of "Y2" in the pdeSys object
  solNDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("nDx3");   // get the position of "Y3" in the pdeSys object
  std::vector < adept::adouble > solNDx[DIM]; // local Y solution


  unsigned solLIndex;
  solLIndex = mlSol->GetIndex ("Lambda1");   // get the position of "lambda" in the ml_sol object
  unsigned solLType;
  solLType = mlSol->GetSolutionType (solLIndex);  // get the finite element type for "lambda"
  unsigned solLPdeIndex;
  solLPdeIndex = mlPdeSys->GetSolPdeIndex ("Lambda1");   // get the position of "lambda" in the pdeSys object
  std::vector < adept::adouble > solL; // local lambda solution

  std::vector< int > SYSDOF; // local to global pdeSys dofs

  vector< double > Res; // local redidual vector
  std::vector< adept::adouble > aResNDx[3]; // local redidual vector
  std::vector< adept::adouble > aResL; // local redidual vector


  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)

  KK->zero();  // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solType);   // number of solution element dofs
    unsigned nLDofs  = msh->GetElementDofNumber (iel, solLType);   // number of solution element dofs


    for (unsigned K = 0; K < DIM; K++) {
      solDx[K].resize (nxDofs);
      solOldDx[K].resize (nxDofs);
      solx[K].resize (nxDofs);
      solNDx[K].resize (nxDofs);
      solL.resize (nLDofs);
    }

    // resize local arrays
    SYSDOF.resize (DIM * nxDofs + nLDofs);
    Res.resize (DIM * nxDofs + nLDofs);       //resize

    for (unsigned K = 0; K < DIM; K++) {
      aResNDx[K].assign (nxDofs, 0.);  //resize and zet to zero
    }
    aResL.assign (nLDofs, 0.);


    // local storage of global mapping and solution
    for (unsigned i = 0; i < nxDofs; i++) {
      unsigned iDDof = msh->GetSolutionDof (i, iel, solType); // global to local mapping between solution node and solution dof
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);
      for (unsigned K = 0; K < DIM; K++) {
        solx[K][i] = (*msh->_topology->_Sol[K]) (iXDof) + (*sol->_Sol[solDxIndex[K]]) (iDDof);
        solDx[K][i] = (*sol->_Sol[solDxIndex[K]]) (iDDof);
        solOldDx[K][i] = (*sol->_SolOld[solDxIndex[K]]) (iDDof);
        solNDx[K][i] = (*sol->_Sol[solNDxIndex[K]]) (iDDof);
        SYSDOF[ K * nxDofs + i] = pdeSys->GetSystemDof (solNDxIndex[K], solNDxPdeIndex[K], i, iel); // global to global mapping between solution node and pdeSys dof
      }
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nLDofs; i++) {
      unsigned iLDof = msh->GetSolutionDof (i, iel, solLType); // global to local mapping between solution node and solution dof
      solL[i] = (*sol->_Sol[solLIndex]) (iLDof); // global to local solution
      SYSDOF[DIM * nxDofs + i] = pdeSys->GetSystemDof (solLIndex, solLPdeIndex, i, iel);  // global to global mapping between solution node and pdeSys dof
    }


    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phix_uv[dim]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // *** get gauss point weight, test function and test function partial derivatives ***
      phix = msh->_finiteElement[ielGeom][solType]->GetPhi (ig);
      phix_uv[0] = msh->_finiteElement[ielGeom][solType]->GetDPhiDXi (ig); //derivative in u
      phix_uv[1] = msh->_finiteElement[ielGeom][solType]->GetDPhiDEta (ig); //derivative in v

      weight = msh->_finiteElement[ielGeom][solType]->GetGaussWeight (ig);

      double solOldDxg[3] = {0., 0., 0.};
      double solDxg[3] = {0., 0., 0.};
      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      adept::adouble solNDx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solDxOld_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      adept::adouble solNDxg[3] = {0., 0., 0.};


      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          solOldDxg[K] += phix[i] * solOldDx[K][i];
          solDxg[K] += phix[i] * solDx[K][i];
          solNDxg[K] += phix[i] * solNDx[K][i];
        }
        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]     += phix_uv[j][i] * solx[K][i];
            solDxOld_uv[K][j] += phix_uv[j][i] * solOldDx[K][i];
            solNDx_uv[K][j]   += phix_uv[j][i] * solNDx[K][i];
          }
        }
      }


      double g[dim][dim] = {{0., 0.}, {0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];

      double normal[DIM];
      normal[0] = (solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) / sqrt (detg);
      normal[1] = (solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);
      normal[2] = (solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);

      adept::adouble DnXmDxdotN = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        DnXmDxdotN += (solDxg[K] - solNDxg[K]) * normal[K];
      }

      double Area = weight * sqrt (detg);

      double gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;


      double Jir[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            Jir[i][J] += gi[i][k] * solx_uv[J][k];
          }
        }
      }

      adept::adouble solNDx_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solDxOld_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            solNDx_Xtan[I][J] += solNDx_uv[I][k] * Jir[k][J];
            solDxOld_Xtan[I][J] += solDxOld_uv[I][k] * Jir[k][J];
          }
        }
      }

      std::vector < adept::adouble > phix_Xtan[DIM];

      for (unsigned J = 0; J < DIM; J++) {
        phix_Xtan[J].assign (nxDofs, 0.);
        for (unsigned inode  = 0; inode < nxDofs; inode++) {
          for (unsigned k = 0; k < dim; k++) {
            phix_Xtan[J][inode] += phix_uv[k][inode] * Jir[k][J];
          }
        }
      }

      for (unsigned K = 0; K < DIM; K++) {
        //Energy equation
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;
          for (unsigned J = 0; J < DIM; J++) {
            if (J != K) {
              term1 += (solNDx_Xtan[K][J] + solNDx_Xtan[J][K]) * phix_Xtan[J][i];
            }
            else {
              term1 += 0.5 * (solNDx_Xtan[K][J] + solNDx_Xtan[J][K]) * phix_Xtan[J][i];
            }
          }
          aResNDx[K][i] += (term1 + solL[0] * phix[i] * normal[K]) * Area;
        }
      }

      aResL[0] += DnXmDxdotN * Area;




    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store


    for (int K = 0; K < DIM; K++) {
      for (int i = 0; i < nxDofs; i++) {
        Res[ K * nxDofs + i] = -aResNDx[K][i].value();
      }
    }

    for (int i = 0; i < nLDofs; i++) {
      Res[DIM * nxDofs + i] = - aResL[i].value();
    }

    RES->add_vector_blocked (Res, SYSDOF);



    Jac.resize ( (DIM * nxDofs + nLDofs) * (DIM * nxDofs + nLDofs));

    // define the dependent variables

    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResNDx[K][0], nxDofs);
    }
    s.dependent (&aResL[0], nLDofs);


    // define the dependent variables

    for (int K = 0; K < DIM; K++) {
      s.independent (&solNDx[K][0], nxDofs);
    }
    s.independent (&solL[0], nLDofs);

    // get the jacobian matrix (ordered by row)
    s.jacobian (&Jac[0], true);

    KK->add_matrix_blocked (Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();

  // ***************** END ASSEMBLY *******************
}

// Building the Conformal Minimization system.
void AssembleO2ConformalMinimization (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  Extract pointers to the several objects that we are going to use.
  NonLinearImplicitSystem* mlPdeSys =
    &ml_prob.get_system< NonLinearImplicitSystem> ("nProj");

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Pointers to the mesh (level) object and elem object in mesh (level).
  Mesh *msh = ml_prob._ml_msh->GetLevel (level);
  elem *el = msh->el;

  // Pointers to the multilevel solution, solution (level) and equation (level).
  MultiLevelSolution *mlSol = ml_prob._ml_sol;
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel (level);
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level];

  // Pointers to global stiffness matrix and residual vector in pdeSys (level).
  SparseMatrix *KK = pdeSys->_KK;
  NumericVector *RES = pdeSys->_RES;

  // Convenience variables to keep track of the dimension.
  const unsigned  dim = 2;
  const unsigned  DIM = 3;

  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // Setting the reference elements to be equilateral triangles.
  std::vector < std::vector < double > > xT (2);

  xT[0].resize (3);
  xT[0][0] = -0.5;
  xT[0][1] = 0.5;
  xT[0][2] = 0.;

  xT[1].resize (3);
  xT[1][0] = 0.;
  xT[1][1] = 0.;
  xT[1][2] = sqrt (3.) / 2.;

  std::vector<double> phi_uv0;
  std::vector<double> phi_uv1;

  std::vector< double > stdVectorPhi;
  std::vector< double > stdVectorPhi_uv;

  // Extract positions of Dx in ml_sol object.
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol->GetIndex ("Dx1");
  solDxIndex[1] = mlSol->GetIndex ("Dx2");
  solDxIndex[2] = mlSol->GetIndex ("Dx3");

  // Extract finite element type for the solution.
  unsigned solType;
  solType = mlSol->GetSolutionType (solDxIndex[0]);

  // Local solution vectors for X, Dx, Xhat, XC.
  std::vector < double > solx[DIM];
  std::vector < double > solDx[DIM];
  std::vector < double > xhat[DIM];
  std::vector < double > xc[DIM];

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get the poitions of Y in the ml_sol object.
  unsigned solNDxIndex[DIM];
  solNDxIndex[0] = mlSol->GetIndex ("nDx1");
  solNDxIndex[1] = mlSol->GetIndex ("nDx2");
  solNDxIndex[2] = mlSol->GetIndex ("nDx3");

  // Get the positions of Y in the pdeSys object.
  unsigned solNDxPdeIndex[DIM];
  solNDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("nDx1");
  solNDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("nDx2");
  solNDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("nDx3");

  // Local solution vectors for Nx and NDx.
  std::vector < adept::adouble > solNDx[DIM];
  std::vector < adept::adouble > solNx[DIM];

  // Get the position of "Lambda1" in the ml_sol object.
  unsigned solLIndex;
  solLIndex = mlSol->GetIndex ("Lambda1");

  // Get the finite element type for "Lambda1".
  unsigned solLType;
  solLType = mlSol->GetSolutionType (solLIndex);

  // Get the position of "Lambda1" in the pdeSys object.
  unsigned solLPdeIndex;
  solLPdeIndex = mlPdeSys->GetSolPdeIndex ("Lambda1");

  // Local Lambda1 solution.
  std::vector < adept::adouble > solL;

  // Local-to-global pdeSys dofs.
  std::vector < int > SYSDOF;

  // Local residual vectors.
  vector< double > Res;
  std::vector< adept::adouble > aResNDx[3];
  std::vector< adept::adouble > aResL;

  // Local Jacobian matrix (ordered by column).
  vector < double > Jac;

  KK->zero();  // Zero all the entries of the Global Matrix
  RES->zero(); // Zero all the entries of the Global Residual

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for (int iel = msh->_elementOffset[iproc];
       iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Numer of solution element dofs.
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solType);
    unsigned nLDofs  = msh->GetElementDofNumber (iel, solLType);

    // Resize local arrays.
    for (unsigned K = 0; K < DIM; K++) {

      xhat[K].resize (nxDofs);
      xc[K].assign (nxDofs, 0.);

      solDx[K].resize (nxDofs);
      solx[K].resize (nxDofs);

      solNDx[K].resize (nxDofs);
      solNx[K].resize (nxDofs);

      solL.resize (nLDofs);
    }

    SYSDOF.resize (DIM * nxDofs + nLDofs);
    Res.resize (DIM * nxDofs + nLDofs);

    for (unsigned K = 0; K < DIM; K++) {
      aResNDx[K].assign (nxDofs, 0.);
    }
    aResL.assign (nLDofs, 0.);


    // local storage of global mapping and solution
    for (unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof (i, iel, solType);
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);

      for (unsigned K = 0; K < DIM; K++) {
        xhat[K][i] = (*msh->_topology->_Sol[K]) (iXDof);
        solDx[K][i] = (*sol->_Sol[solDxIndex[K]]) (iDDof);
        solx[K][i] = xhat[K][i] + solDx[K][i];
        solNDx[K][i] = (*sol->_Sol[solNDxIndex[K]]) (iDDof);

        // Global-to-global mapping between NDx solution node and pdeSys dof.
        SYSDOF[ K * nxDofs + i] =
          pdeSys->GetSystemDof (solNDxIndex[K], solNDxPdeIndex[K], i, iel);
      }
    }

    // Local storage of global mapping and solution.
    for (unsigned i = 0; i < nLDofs; i++) {

      // Global-to-local mapping between Lambda solution node and solution dof.
      unsigned iLDof = msh->GetSolutionDof (i, iel, solLType);
      solL[i] = (*sol->_Sol[solLIndex]) (iLDof);

      // Global-to-global mapping between Lambda solution node and pdeSys dof.
      SYSDOF[DIM * nxDofs + i] =
        pdeSys->GetSystemDof (solLIndex, solLPdeIndex, i, iel);
    }

    // start a new recording of all the operations involving adept variables.
    s.new_recording();

    // begin GAUSS POINT LOOP
    for (unsigned ig = 0; ig < msh->
         _finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phiL;  // local test function
      const double *phix_uv[dim]; // first order derivatives in (u,v)

      double weight; // gauss point weight

      // Get Gauss point weight, test function, and first order derivatives.
      if (ielGeom == QUAD) {
        phix = msh->_finiteElement[ielGeom][solType]->GetPhi (ig);
        phix_uv[0] = msh->_finiteElement[ielGeom][solType]->GetDPhiDXi (ig);
        phix_uv[1] = msh->_finiteElement[ielGeom][solType]->GetDPhiDEta (ig);
        weight = msh->_finiteElement[ielGeom][solType]->GetGaussWeight (ig);
        phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi (ig);
      }

      // Special adjustments for triangles.
      else {
        msh->_finiteElement[ielGeom][solType]->Jacobian (xT, ig,
                                                         weight, stdVectorPhi, stdVectorPhi_uv);
        phix = &stdVectorPhi[0];

        phi_uv0.resize (nxDofs);
        phi_uv1.resize (nxDofs);

        for (unsigned i = 0; i < nxDofs; i++) {
          phi_uv0[i] = stdVectorPhi_uv[i * dim];
          phi_uv1[i] = stdVectorPhi_uv[i * dim + 1];
        }
        phix_uv[0] = &phi_uv0[0];
        phix_uv[1] = &phi_uv1[0];

        phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi (ig);

      }

      // Initialize and compute values of x, Dx, NDx, x_uv at the Gauss points.
      double solDxg[3] = {0., 0., 0.};
      adept::adouble solNDxg[3] = {0., 0., 0.};

      adept::adouble solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solNx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};



      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          solDxg[K] += phix[i] * solDx[K][i];
          solNDxg[K] += phix[i] * solNDx[K][i];
        }

        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * (xhat[K][i] + 0.5 * (solDx[K][i] + solNDx[K][i]));
            solNx_uv[K][j]   += phix_uv[j][i] * (xhat[K][i] + ( solNDx[K][i] ) );
          }
        }
      }

      adept::adouble solLg = 0.;
      for (unsigned i = 0; i < nLDofs; i++) {
        solLg += phiL[i] * solL[i];
      }

      // Compute the metric, metric determinant, and area element.
      adept::adouble g[dim][dim] = {{0., 0.}, {0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      adept::adouble detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      adept::adouble Area = weight * sqrt (detg);
      double Area2 = weight; // Trick to give equal weight to each element.

      // Compute components of the unit normal N.
      adept::adouble normal[DIM];
      normal[0] = (solx_uv[1][0] * solx_uv[2][1]
                   - solx_uv[2][0] * solx_uv[1][1]) / sqrt (detg);
      normal[1] = (solx_uv[2][0] * solx_uv[0][1]
                   - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);
      normal[2] = (solx_uv[0][0] * solx_uv[1][1]
                   - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);

      // Discretize the equation \delta CD = 0 on the basis d/du, d/dv.
      adept::adouble V[DIM];
      V[0] = solNx_uv[0][1] - normal[1] * solNx_uv[2][0] + normal[2] * solNx_uv[1][0];
      V[1] = solNx_uv[1][1] - normal[2] * solNx_uv[0][0] + normal[0] * solNx_uv[2][0];
      V[2] = solNx_uv[2][1] - normal[0] * solNx_uv[1][0] + normal[1] * solNx_uv[0][0];

      adept::adouble W[DIM];
      W[0] = solNx_uv[0][0] + normal[1] * solNx_uv[2][1] - normal[2] * solNx_uv[1][1];
      W[1] = solNx_uv[1][0] + normal[2] * solNx_uv[0][1] - normal[0] * solNx_uv[2][1];
      W[2] = solNx_uv[2][0] + normal[0] * solNx_uv[1][1] - normal[1] * solNx_uv[0][1];

      adept::adouble M[DIM][dim];
      M[0][0] = W[0] - normal[2] * V[1] + normal[1] * V[2];
      M[1][0] = W[1] - normal[0] * V[2] + normal[2] * V[0];
      M[2][0] = W[2] - normal[1] * V[0] + normal[0] * V[1];

      M[0][1] = V[0] + normal[2] * W[1] - normal[1] * W[2];
      M[1][1] = V[1] + normal[0] * W[2] - normal[2] * W[0];
      M[2][1] = V[2] + normal[1] * W[0] - normal[0] * W[1];

      // Compute new X minus old X dot N, for "reparametrization".
      adept::adouble DnXmDxdotN = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        DnXmDxdotN += (solDxg[K] - solNDxg[K]) * normal[K];
      }

      // Implement the Conformal Minimization equations.
      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;

          for (unsigned j = 0; j < dim; j++) {
            term1 +=  M[K][j] * phix_uv[j][i];
          }

          // Conformal energy equation (with trick).
          aResNDx[K][i] += term1 * Area2
                           + timederiv * (solNDxg[K] - solDxg[K]) * phix[i] * Area2
                           + solLg * phix[i] * normal[K] * Area;  //no2
        }
      }

      // Lagrange multiplier equation (with trick).
      for (unsigned i = 0; i < nLDofs; i++) {
        aResL[i] += phiL[i] * (DnXmDxdotN + eps * solL[i]) * Area; // no2
      }

    } // end GAUSS POINT LOOP

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    //copy the value of the adept::adoube aRes in double Res and store

    for (int K = 0; K < DIM; K++) {
      for (int i = 0; i < nxDofs; i++) {
        Res[ K * nxDofs + i] = -aResNDx[K][i].value(); // for X
      }
    }

    for (int i = 0; i < nLDofs; i++) {
      Res[DIM * nxDofs + i] = - aResL[i].value(); // for Lambda
    }

    RES->add_vector_blocked (Res, SYSDOF);

    // Resize Jacobian.
    Jac.resize ( (DIM * nxDofs + nLDofs) * (DIM * nxDofs + nLDofs));

    // Define the dependent variables.
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResNDx[K][0], nxDofs);
    }
    s.dependent (&aResL[0], nLDofs);

    // Define the independent variables.
    for (int K = 0; K < DIM; K++) {
      s.independent (&solNDx[K][0], nxDofs);
    }
    s.independent (&solL[0], nLDofs);

    // Get the jacobian matrix (ordered by row).
    s.jacobian (&Jac[0], true);

    KK->add_matrix_blocked (Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } //end ELEMENT LOOP for each process.

  RES->close();
  KK->close();

} // end AssembleO2ConformalMinimization.


void AssembleSphereConformalMinimization (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  Extract pointers to the several objects that we are going to use.
  NonLinearImplicitSystem* mlPdeSys =
    &ml_prob.get_system< NonLinearImplicitSystem> ("nProj");

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Pointers to the mesh (level) object and elem object in mesh (level).
  Mesh *msh = ml_prob._ml_msh->GetLevel (level);
  elem *el = msh->el;

  // Pointers to the multilevel solution, solution (level) and equation (level).
  MultiLevelSolution *mlSol = ml_prob._ml_sol;
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel (level);
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level];

  // Pointers to global stiffness matrix and residual vector in pdeSys (level).
  SparseMatrix *KK = pdeSys->_KK;
  NumericVector *RES = pdeSys->_RES;

  // Convenience variables to keep track of the dimension.
  const unsigned  dim = 2;
  const unsigned  DIM = 3;

  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // Setting the reference elements to be equilateral triangles.
  std::vector < std::vector < double > > xT (2);

  xT[0].resize (3);
  xT[0][0] = -0.5;
  xT[0][1] = 0.5;
  xT[0][2] = 0.;

  xT[1].resize (3);
  xT[1][0] = 0.;
  xT[1][1] = 0.;
  xT[1][2] = sqrt (3.) / 2.;

  std::vector<double> phi_uv0;
  std::vector<double> phi_uv1;

  std::vector< double > stdVectorPhi;
  std::vector< double > stdVectorPhi_uv;

  // Extract positions of Dx in ml_sol object.
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol->GetIndex ("Dx1");
  solDxIndex[1] = mlSol->GetIndex ("Dx2");
  solDxIndex[2] = mlSol->GetIndex ("Dx3");


  unsigned solYIndex[DIM];
  solYIndex[0] = mlSol->GetIndex ("Y1");
  solYIndex[1] = mlSol->GetIndex ("Y2");
  solYIndex[2] = mlSol->GetIndex ("Y3");

  // Extract finite element type for the solution.
  unsigned solType;
  solType = mlSol->GetSolutionType (solDxIndex[0]);

  // Local solution vectors for X, Dx, Xhat, XC.
  std::vector < double > solY[DIM];
  std::vector < double > solx[DIM];
  std::vector < double > solDx[DIM];
  std::vector < double > xhat[DIM];


  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get the poitions of Y in the ml_sol object.
  unsigned solNDxIndex[DIM];
  solNDxIndex[0] = mlSol->GetIndex ("nDx1");
  solNDxIndex[1] = mlSol->GetIndex ("nDx2");
  solNDxIndex[2] = mlSol->GetIndex ("nDx3");

  // Get the positions of Y in the pdeSys object.
  unsigned solNDxPdeIndex[DIM];
  solNDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("nDx1");
  solNDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("nDx2");
  solNDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("nDx3");

  // Local solution vectors for Nx and NDx.
  std::vector < adept::adouble > solNDx[DIM];
  std::vector < adept::adouble > solNx[DIM];

  // Get the position of "Lambda1" in the ml_sol object.
  unsigned solLIndex;
  solLIndex = mlSol->GetIndex ("Lambda1");

  // Get the finite element type for "Lambda1".
  unsigned solLType;
  solLType = mlSol->GetSolutionType (solLIndex);

  // Get the position of "Lambda1" in the pdeSys object.
  unsigned solLPdeIndex;
  solLPdeIndex = mlPdeSys->GetSolPdeIndex ("Lambda1");

  // Local Lambda1 solution.
  std::vector < adept::adouble > solL;

  // Local-to-global pdeSys dofs.
  std::vector < int > SYSDOF;

  // Local residual vectors.
  vector< double > Res;
  std::vector< adept::adouble > aResNDx[3];
  std::vector< adept::adouble > aResL;

  // Local Jacobian matrix (ordered by column).
  vector < double > Jac;

  KK->zero();  // Zero all the entries of the Global Matrix
  RES->zero(); // Zero all the entries of the Global Residual

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for (int iel = msh->_elementOffset[iproc];
       iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Numer of solution element dofs.
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solType);
    unsigned nLDofs  = msh->GetElementDofNumber (iel, solLType);

    // Resize local arrays.
    for (unsigned K = 0; K < DIM; K++) {

      xhat[K].resize (nxDofs);

      solDx[K].resize (nxDofs);
      solx[K].resize (nxDofs);

      solY[K].resize (nxDofs);

      solNDx[K].resize (nxDofs);
      solNx[K].resize (nxDofs);

      solL.resize (nLDofs);
    }

    SYSDOF.resize (DIM * nxDofs + nLDofs);
    Res.resize (DIM * nxDofs + nLDofs);

    for (unsigned K = 0; K < DIM; K++) {
      aResNDx[K].assign (nxDofs, 0.);
    }
    aResL.assign (nLDofs, 0.);


    // local storage of global mapping and solution
    for (unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof (i, iel, solType);
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);

      for (unsigned K = 0; K < DIM; K++) {
        xhat[K][i] = (*msh->_topology->_Sol[K]) (iXDof);
        solDx[K][i] = (*sol->_Sol[solDxIndex[K]]) (iDDof);

        solY[K][i] = (*sol->_Sol[solYIndex[K]]) (iDDof);

        solx[K][i] = xhat[K][i] + solDx[K][i];
        solNDx[K][i] = (*sol->_Sol[solNDxIndex[K]]) (iDDof);

        // Global-to-global mapping between NDx solution node and pdeSys dof.
        SYSDOF[ K * nxDofs + i] =
          pdeSys->GetSystemDof (solNDxIndex[K], solNDxPdeIndex[K], i, iel);
      }
    }

    // Local storage of global mapping and solution.
    for (unsigned i = 0; i < nLDofs; i++) {

      // Global-to-local mapping between Lambda solution node and solution dof.
      unsigned iLDof = msh->GetSolutionDof (i, iel, solLType);
      solL[i] = (*sol->_Sol[solLIndex]) (iLDof);

      // Global-to-global mapping between Lambda solution node and pdeSys dof.
      SYSDOF[DIM * nxDofs + i] =
        pdeSys->GetSystemDof (solLIndex, solLPdeIndex, i, iel);
    }

    // start a new recording of all the operations involving adept variables.
    s.new_recording();

    // begin GAUSS POINT LOOP
    for (unsigned ig = 0; ig < msh->
         _finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phiL;  // local test function
      const double *phix_uv[dim]; // first order derivatives in (u,v)

      double weight; // gauss point weight

      // Get Gauss point weight, test function, and first order derivatives.
      if (ielGeom == QUAD) {
        phix = msh->_finiteElement[ielGeom][solType]->GetPhi (ig);
        phix_uv[0] = msh->_finiteElement[ielGeom][solType]->GetDPhiDXi (ig);
        phix_uv[1] = msh->_finiteElement[ielGeom][solType]->GetDPhiDEta (ig);
        weight = msh->_finiteElement[ielGeom][solType]->GetGaussWeight (ig);
        phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi (ig);
      }

      // Special adjustments for triangles.
      else {
        msh->_finiteElement[ielGeom][solType]->Jacobian (xT, ig,
                                                         weight, stdVectorPhi, stdVectorPhi_uv);
        phix = &stdVectorPhi[0];

        phi_uv0.resize (nxDofs);
        phi_uv1.resize (nxDofs);

        for (unsigned i = 0; i < nxDofs; i++) {
          phi_uv0[i] = stdVectorPhi_uv[i * dim];
          phi_uv1[i] = stdVectorPhi_uv[i * dim + 1];
        }
        phix_uv[0] = &phi_uv0[0];
        phix_uv[1] = &phi_uv1[0];

        phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi (ig);

      }

      // Initialize and compute values of x, Dx, NDx, x_uv at the Gauss points.
      double solYg[3] = {0., 0., 0.};
      double solxg[3] = {0., 0., 0.};

      double solDxg[3] = {0., 0., 0.};
      adept::adouble solNDxg[3] = {0., 0., 0.};

      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble X_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};



      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          solxg[K] += phix[i] * solx[K][i];
          solYg[K] += phix[i] * solY[K][i];
          solDxg[K] += phix[i] * solDx[K][i];
          solNDxg[K] += phix[i] * solNDx[K][i];
        }

        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * solx[K][i];
            X_uv[K][j]   += phix_uv[j][i] * (xhat[K][i] + solNDx[K][i]);
          }
        }
      }

      adept::adouble solLg = 0.;
      for (unsigned i = 0; i < nLDofs; i++) {
        solLg += phiL[i] * solL[i];
      }

      // Compute the metric, metric determinant, and area element.
      double g[dim][dim] = {{0., 0.}, {0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      double Area = weight * sqrt (detg);
      double Area2 = weight; // Trick to give equal weight to each element.

      // Compute components of the unit normal N.
      double normal[DIM];
      normal[0] = (solx_uv[1][0] * solx_uv[2][1]
                   - solx_uv[2][0] * solx_uv[1][1]) / sqrt (detg);
      normal[1] = (solx_uv[2][0] * solx_uv[0][1]
                   - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);
      normal[2] = (solx_uv[0][0] * solx_uv[1][1]
                   - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);

      double a = 0;
      for (unsigned I = 0; I < DIM; I++) {
        a -= solYg[I] * normal[I];
      }
      double signa = (a > 0) ? 1 : -1;
      a = 2. / (a + 0 * signa * eps);

      //       std::cout << solxg[0] << " " << solxg[1] <<" "<< solxg[2]<<std::endl;
      //       std::cout << normal[0] << " " << normal[1] <<" "<< normal[2]<<std::endl;
      //       std::cout << a<<std::endl;

      // Discretize the equation \delta CD = 0 on the basis d/du, d/dv.
      adept::adouble V[DIM];
      V[0] = X_uv[0][1] - normal[1] * X_uv[2][0] + normal[2] * X_uv[1][0];
      V[1] = X_uv[1][1] - normal[2] * X_uv[0][0] + normal[0] * X_uv[2][0];
      V[2] = X_uv[2][1] - normal[0] * X_uv[1][0] + normal[1] * X_uv[0][0];

      adept::adouble W[DIM];
      W[0] = X_uv[0][0] + normal[1] * X_uv[2][1] - normal[2] * X_uv[1][1];
      W[1] = X_uv[1][0] + normal[2] * X_uv[0][1] - normal[0] * X_uv[2][1];
      W[2] = X_uv[2][0] + normal[0] * X_uv[1][1] - normal[1] * X_uv[0][1];

      adept::adouble M[DIM][dim];
      M[0][0] = W[0] - normal[2] * V[1] + normal[1] * V[2];
      M[1][0] = W[1] - normal[0] * V[2] + normal[2] * V[0];
      M[2][0] = W[2] - normal[1] * V[0] + normal[0] * V[1];

      M[0][1] = V[0] + normal[2] * W[1] - normal[1] * W[2];
      M[1][1] = V[1] + normal[0] * W[2] - normal[2] * W[0];
      M[2][1] = V[2] + normal[1] * W[0] - normal[0] * W[1];

      // Compute new X minus old X dot N, for "reparametrization".
      adept::adouble DnXmDxdotN = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        DnXmDxdotN += (solDxg[K] - solNDxg[K]) * normal[K];
      }

      // Compute the metric inverse.
      double gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;

      // Compute the "reduced Jacobian" g^{ij}X_j .
      double Jir[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            Jir[i][J] += gi[i][k] * solx_uv[J][k];
          }
        }
      }

      // Initializing tangential gradients of X and W (new, middle, old).
      double solx_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble NX_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble XminusXold_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      // adept::adouble solW_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

      // Computing tangential gradients defined above.
      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            NX_Xtan[I][J] += X_uv[I][k] * Jir[k][J];
            solx_Xtan[I][J] += solx_uv[I][k] * Jir[k][J];
            // solW_Xtan[I][J] += solW_uv[I][k] * Jir[k][J];
          }
        }
      }

      // Define and compute gradients of test functions for X and W.
      std::vector < adept::adouble > phiY_Xtan[DIM];
      std::vector < adept::adouble > phix_Xtan[DIM];

      for (unsigned J = 0; J < DIM; J++) {
        phix_Xtan[J].assign (nxDofs, 0.);
        phiY_Xtan[J].assign (nxDofs, 0.);

        for (unsigned inode  = 0; inode < nxDofs; inode++) {
          for (unsigned k = 0; k < dim; k++) {
            phix_Xtan[J][inode] += phix_uv[k][inode] * Jir[k][J];
          }
        }

        for (unsigned inode  = 0; inode < nxDofs; inode++) {
          for (unsigned k = 0; k < dim; k++) {
            phiY_Xtan[J][inode] += phix_uv[k][inode] * Jir[k][J];
          }
        }
      }

      // Implement the Conformal Minimization equations.
      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;

          for (unsigned j = 0; j < dim; j++) {
            term1 +=  M[K][j] * phix_uv[j][i];
            // term3 += solx_uv[K][j] * phix_uv[j][i];
          }


          adept::adouble term2 = 0.;
          //           for (unsigned J = 0; J < DIM; J++) {
          //             term2 += (a * normal[J] - solDxg[J] + solNDxg [J]) * phix[i];
          //           }

          term2 += (a * normal[K] - solDxg[K] + solNDxg [K]) * phix[i];

          // Conformal energy equation (with trick).
          aResNDx[K][i] += term1 * Area2
                           + solLg * term2 * Area;  //no2
        }

        // Lagrange multiplier equation (with trick).
        for (unsigned i = 0; i < nLDofs; i++) {
          adept::adouble term1 = 0.;
          adept::adouble term2 = 0.;
          for (unsigned J = 0; J < DIM; J++) {
            term1 += (a * normal[J] - solDxg[J] + solNDxg [J]) * (a * normal[J] - solDxg[J] + solNDxg [J]) ;
            term2 += (a * normal[J]) * (a * normal[J]) ;
          }
          aResL[i] += (phiL[i] * (term1 - term2) + eps * 0 * solLg * phiL[i]) * Area;    // no2
        }

        //         for(unsigned i = 0; i< nLDofs; i++){
        //           aResL[i] += phiL[i] * (DnXmDxdotN + eps * solL[i]) * Area2; // no2
        //         }

      }

    } // end GAUSS POINT LOOP

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    //copy the value of the adept::adoube aRes in double Res and store

    for (int K = 0; K < DIM; K++) {
      for (int i = 0; i < nxDofs; i++) {
        Res[ K * nxDofs + i] = -aResNDx[K][i].value(); // for X
      }
    }

    for (int i = 0; i < nLDofs; i++) {
      Res[DIM * nxDofs + i] = - aResL[i].value(); // for Lambda
    }

    RES->add_vector_blocked (Res, SYSDOF);

    // Resize Jacobian.
    Jac.resize ( (DIM * nxDofs + nLDofs) * (DIM * nxDofs + nLDofs));

    // Define the dependent variables.
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResNDx[K][0], nxDofs);
    }
    s.dependent (&aResL[0], nLDofs);

    // Define the independent variables.
    for (int K = 0; K < DIM; K++) {
      s.independent (&solNDx[K][0], nxDofs);
    }
    s.independent (&solL[0], nLDofs);

    // Get the jacobian matrix (ordered by row).
    s.jacobian (&Jac[0], true);

    KK->add_matrix_blocked (Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } //end ELEMENT LOOP for each process.

  RES->close();
  KK->close();

} // end AssembleSphereConformalMinimization.
