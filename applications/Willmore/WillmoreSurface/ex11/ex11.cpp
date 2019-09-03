/** tutorial/Ex3
 * This example shows how to set and solve the weak form of the nonlinear problem
 *                     -\Delta^2 u = f(x) \text{ on }\Omega,
 *            u=0 \text{ on } \Gamma,
 *      \Delat u=0 \text{ on } \Gamma,
 * on a box domain $\Omega$ with boundary $\Gamma$,
 * by using a system of second order partial differential equation.
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/

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
#include "PetscMatrix.hpp"

#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

const unsigned P[3] = {2, 3, 4};
const double ap[3] = {1, 0., 0.};

const double normalSign = -1.;
using namespace femus;


const unsigned volumeConstraint = true;
const unsigned areaConstraint = true;


void CopyDisplacement(MultiLevelSolution &mlSol,  const bool &forward);

void AssemblePWillmore (MultiLevelProblem&);

void AssembleInit (MultiLevelProblem&);

void AssembleNProj (MultiLevelProblem&);

double GetTimeStep (const double t) {
  //if(time==0) return 1.0e-10;
  //return 0.0001;

  double dt0 = 0.001;
  double s = 1.;
  double n = 0.3;
  return dt0 * pow (1. + t / pow (dt0, s), n);


}

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

int main (int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit (argc, args, MPI_COMM_WORLD);


  // define multilevel mesh

  unsigned maxNumberOfMeshes;

  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;

  //mlMsh.ReadCoarseMesh("./input/torus.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("./input/sphere.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("./input/ellipsoidRef3.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("./input/ellipsoidV1.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("./input/genusOne.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh ("./input/knot.neu", "seventh", scalingFactor);
  // mlMsh.ReadCoarseMesh ("./input/cube.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("./input/horseShoe.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("./input/tiltedTorus.neu", "seventh", scalingFactor);
  
  //mlMsh.ReadCoarseMesh ("./input/ellipsoidSphere.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/CliffordTorus.neu", "seventh", scalingFactor);

  unsigned numberOfUniformLevels = 2;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh (numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels (numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol (&mlMsh);

  // add variables to mlSol

  mlSol.AddSolution ("Dx1", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("Dx2", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("Dx3", LAGRANGE, FIRST, 2);

  mlSol.AddSolution ("W1", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("W2", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("W3", LAGRANGE, FIRST, 2);

  mlSol.AddSolution ("Y1", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("Y2", LAGRANGE, FIRST, 2);
  mlSol.AddSolution ("Y3", LAGRANGE, FIRST, 2);

  if (volumeConstraint || areaConstraint) {
    mlSol.AddSolution ("Lambda", DISCONTINUOUS_POLYNOMIAL, ZERO, 0);
  }
  
  mlSol.AddSolution ("nDx1", LAGRANGE, FIRST, 0);
  mlSol.AddSolution ("nDx2", LAGRANGE, FIRST, 0);
  mlSol.AddSolution ("nDx3", LAGRANGE, FIRST, 0);
  
  mlSol.AddSolution ("Lambda1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0);

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

  //BEGIN Get initial curbvature data
  NonLinearImplicitSystem& system0 = mlProb.add_system < NonLinearImplicitSystem > ("Init");

  // add solution "X", "Y", "Z" and "H" to the system

  system0.AddSolutionToSystemPDE ("Y1");
  system0.AddSolutionToSystemPDE ("Y2");
  system0.AddSolutionToSystemPDE ("Y3");

  system0.AddSolutionToSystemPDE ("W1");
  system0.AddSolutionToSystemPDE ("W2");
  system0.AddSolutionToSystemPDE ("W3");

  system0.SetMaxNumberOfNonLinearIterations (40);
  system0.SetNonLinearConvergenceTolerance (1.e-9);

  // attach the assembling function to system
  system0.SetAssembleFunction (AssembleInit);

  // initilaize and solve the system
  system0.init();

  system0.MGsolve();

  mlSol.SetWriter (VTK);
  std::vector<std::string> mov_vars;
  mov_vars.push_back ("Dx1");
  mov_vars.push_back ("Dx2");
  mov_vars.push_back ("Dx3");
  mlSol.GetWriter()->SetMovingMesh (mov_vars);

  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back ("All");

  mlSol.GetWriter()->SetDebugOutput (true);
  mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "linear", variablesToBePrinted, 0);
  //END Get initial curbvature data

  // add system Wilmore in mlProb as a Linear Implicit System
  TransientNonlinearImplicitSystem& system = mlProb.add_system < TransientNonlinearImplicitSystem > ("PWillmore");

  // add solution "X", "Y", "Z" and "H" to the system

  system.AddSolutionToSystemPDE ("Dx1");
  system.AddSolutionToSystemPDE ("Dx2");
  system.AddSolutionToSystemPDE ("Dx3");
  system.AddSolutionToSystemPDE ("Y1");
  system.AddSolutionToSystemPDE ("Y2");
  system.AddSolutionToSystemPDE ("Y3");
  system.AddSolutionToSystemPDE ("W1");
  system.AddSolutionToSystemPDE ("W2");
  system.AddSolutionToSystemPDE ("W3");

  if (volumeConstraint || areaConstraint) {
    system.AddSolutionToSystemPDE ("Lambda");
    system.SetNumberOfGlobalVariables (volumeConstraint + areaConstraint);
  }

  system.SetMaxNumberOfNonLinearIterations (20);
  system.SetNonLinearConvergenceTolerance (1.e-9);

  // attach the assembling function to system
  system.SetAssembleFunction (AssemblePWillmore);

  system.AttachGetTimeIntervalFunction (GetTimeStep);

  // initilaize and solve the system
  
  
  system.init();
  system.SetMgType (V_CYCLE);
  
  
  NonLinearImplicitSystem& system2 = mlProb.add_system < NonLinearImplicitSystem > ("nProj");
  
  // add solution "X", "Y", "Z" and "H" to the system
  
  system2.AddSolutionToSystemPDE ("nDx1");
  system2.AddSolutionToSystemPDE ("nDx2");
  system2.AddSolutionToSystemPDE ("nDx3");
  system2.AddSolutionToSystemPDE ("Lambda1");
  
  
  system2.SetMaxNumberOfNonLinearIterations (40);
  system2.SetNonLinearConvergenceTolerance (1.e-9);
  
  // attach the assembling function to system
  system2.SetAssembleFunction (AssembleNProj);
  
  system2.init();
  
  
  
  
  unsigned numberOfTimeSteps = 1000u;
  unsigned printInterval = 1u;
  for (unsigned time_step = 0; time_step < numberOfTimeSteps; time_step++) {
    system.CopySolutionToOldSolution();
    system.MGsolve();
    
    if(time_step%1 == 0){
      CopyDisplacement(mlSol, true);
        
      system2.MGsolve();
    
      CopyDisplacement(mlSol, false);
      
      system0.MGsolve();
    }
    
    if ( (time_step + 1) % printInterval == 0)
      mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "linear", variablesToBePrinted, (time_step + 1) / printInterval);
  }
  return 0;
}

void AssemblePWillmore (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  TransientNonlinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientNonlinearImplicitSystem> ("PWillmore");   // pointer to the linear implicit system named "Poisson"

  double dt = mlPdeSys->GetIntervalTime();

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

  unsigned solxType;
  solxType = mlSol->GetSolutionType (solDxIndex[0]);  // get the finite element type for "U"

  unsigned solDxPdeIndex[DIM];
  solDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("Dx1");   // get the position of "DX" in the pdeSys object
  solDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("Dx2");   // get the position of "DY" in the pdeSys object
  solDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("Dx3");   // get the position of "DZ" in the pdeSys object

  std::vector < adept::adouble > solx[DIM];  // surface coordinates
  std::vector < double > solxOld[DIM];  // surface coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  unsigned solYIndex[DIM];
  solYIndex[0] = mlSol->GetIndex ("Y1");   // get the position of "Y1" in the ml_sol object
  solYIndex[1] = mlSol->GetIndex ("Y2");   // get the position of "Y2" in the ml_sol object
  solYIndex[2] = mlSol->GetIndex ("Y3");   // get the position of "Y3" in the ml_sol object

  unsigned solYType;
  solYType = mlSol->GetSolutionType (solYIndex[0]);  // get the finite element type for "Y"

  unsigned solYPdeIndex[DIM];
  solYPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("Y1");   // get the position of "Y1" in the pdeSys object
  solYPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("Y2");   // get the position of "Y2" in the pdeSys object
  solYPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("Y3");   // get the position of "Y3" in the pdeSys object

  std::vector < adept::adouble > solY[DIM];  // surface coordinates
  std::vector < double > solYOld[DIM];  // surface coordinates

  unsigned solWIndex[DIM];
  solWIndex[0] = mlSol->GetIndex ("W1");   // get the position of "W1" in the ml_sol object
  solWIndex[1] = mlSol->GetIndex ("W2");   // get the position of "W2" in the ml_sol object
  solWIndex[2] = mlSol->GetIndex ("W3");   // get the position of "W3" in the ml_sol object

  unsigned solWType;
  solWType = mlSol->GetSolutionType (solWIndex[0]);  // get the finite element type for "W"

  unsigned solWPdeIndex[DIM];
  solWPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("W1");   // get the position of "W1" in the pdeSys object
  solWPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("W2");   // get the position of "W2" in the pdeSys object
  solWPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("W3");   // get the position of "W3" in the pdeSys object

  std::vector < adept::adouble > solW[DIM]; // local W solution
  std::vector < double > solWOld[DIM];  // surface coordinates
  std::vector< unsigned > SYSDOF; // local to global pdeSys dofs

  vector< double > Res; // local redidual vector
  std::vector< adept::adouble > aResx[3]; // local redidual vector
  std::vector< adept::adouble > aResY[3]; // local redidual vector
  std::vector< adept::adouble > aResW[3]; // local redidual vector

  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)

  //MatSetOption ( ( static_cast<PetscMatrix*> ( KK ) )->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE );

  KK->zero();  // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual



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
    
    if(volumeConstraint){
      double lambda1;
      if (iproc == 0) {
        lambda1 = (*sol->_Sol[solLambdaIndex]) (0); // global to local solution
        lambda1PdeDof = pdeSys->GetSystemDof (solLambdaIndex, solLambaPdeIndex, 0, 0);
      }
      MPI_Bcast (&lambda1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast (&lambda1PdeDof, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      solLambda1 = lambda1;
    }

    if(areaConstraint){
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
    columns[1] = (volumeConstraint)? lambda1PdeDof : lambda2PdeDof;
    for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) { //for equations other than Lagrange multiplier
      if (iel > volumeConstraint * areaConstraint) {
        row[0] = pdeSys->GetSystemDof (solLambdaIndex, solLambaPdeIndex, 0, iel);
        columns[0] = row[0];
        KK->add_matrix_blocked (value, row, columns);
      }
    }
  }
  double surface = 0.;
  double volume = 0.;
  double energy = 0.;
  
  

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solxType);   // number of solution element dofs
    unsigned nYDofs  = msh->GetElementDofNumber (iel, solYType);   // number of solution element dofs
    unsigned nWDofs  = msh->GetElementDofNumber (iel, solWType);   // number of solution element dofs

    for (unsigned K = 0; K < DIM; K++) {
      solx[K].resize (nxDofs);
      solxOld[K].resize (nxDofs);
      solY[K].resize (nYDofs);
      solYOld[K].resize (nYDofs);
      solW[K].resize (nWDofs);
      solWOld[K].resize (nWDofs);
    }

    unsigned sizeAll = DIM * (nxDofs + nYDofs +  nWDofs) + volumeConstraint + areaConstraint;

    // resize local arrays
    SYSDOF.resize (sizeAll);

    Res.resize (sizeAll);       //resize

    for (unsigned K = 0; K < DIM; K++) {
      aResx[K].assign (nxDofs, 0.);  //resize and set to zero
      aResY[K].assign (nYDofs, 0.);  //resize and set to zero
      aResW[K].assign (nWDofs, 0.);  //resize and zet to zero
    }
    aResLambda1 = 0.;
    aResLambda2 = 0.;

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nxDofs; i++) {
      unsigned iDDof = msh->GetSolutionDof (i, iel, solxType); // global to local mapping between solution node and solution dof
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);

      for (unsigned K = 0; K < DIM; K++) {
        solxOld[K][i] = (*msh->_topology->_Sol[K]) (iXDof) + (*sol->_SolOld[solDxIndex[K]]) (iDDof);
        solx[K][i] = (*msh->_topology->_Sol[K]) (iXDof) + (*sol->_Sol[solDxIndex[K]]) (iDDof);
        SYSDOF[K * nxDofs + i] = pdeSys->GetSystemDof (solDxIndex[K], solDxPdeIndex[K], i, iel); // global to global mapping between solution node and pdeSys dof
      }
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nYDofs; i++) {
      unsigned iYDof = msh->GetSolutionDof (i, iel, solYType); // global to local mapping between solution node and solution dof
      for (unsigned K = 0; K < DIM; K++) {
        solYOld[K][i] = (*sol->_SolOld[solYIndex[K]]) (iYDof); // global to local solution
        solY[K][i] = (*sol->_Sol[solYIndex[K]]) (iYDof); // global to local solution
        SYSDOF[DIM * nxDofs + K * nYDofs + i] = pdeSys->GetSystemDof (solYIndex[K], solYPdeIndex[K], i, iel); // global to global mapping between solution node and pdeSys dof
      }
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nWDofs; i++) {
      unsigned iWDof = msh->GetSolutionDof (i, iel, solWType); // global to local mapping between solution node and solution dof
      for (unsigned K = 0; K < DIM; K++) {
        solWOld[K][i] = (*sol->_SolOld[solWIndex[K]]) (iWDof); // global to local solution
        solW[K][i] = (*sol->_Sol[solWIndex[K]]) (iWDof); // global to local solution
        SYSDOF[DIM * (nxDofs + nYDofs) + K * nWDofs + i] = pdeSys->GetSystemDof (solWIndex[K], solWPdeIndex[K], i, iel);  // global to global mapping between solution node and pdeSys dof
      }
    }

    if (volumeConstraint) {
      SYSDOF[sizeAll - 1u - areaConstraint ] = lambda1PdeDof;
    }
    if (areaConstraint) {
      SYSDOF[sizeAll - 1u ] = lambda2PdeDof;
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

   
    
    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solxType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phix_uv[dim]; // local test function first order partial derivatives

      const double *phiY;  // local test function

      const double *phiW;  // local test function
      const double *phiW_uv[dim]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // *** get gauss point weight, test function and test function partial derivatives ***
      phix = msh->_finiteElement[ielGeom][solxType]->GetPhi (ig);
      phix_uv[0] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDXi (ig); //derivative in u
      phix_uv[1] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDEta (ig); //derivative in v

      phiY = msh->_finiteElement[ielGeom][solYType]->GetPhi (ig);

      phiW = msh->_finiteElement[ielGeom][solWType]->GetPhi (ig);
      phiW_uv[0] = msh->_finiteElement[ielGeom][solWType]->GetDPhiDXi (ig);
      phiW_uv[1] = msh->_finiteElement[ielGeom][solWType]->GetDPhiDEta (ig);

      weight = msh->_finiteElement[ielGeom][solxType]->GetGaussWeight (ig);

      
      
      adept::adouble solYg[3] = {0., 0., 0.};
      adept::adouble solWg[3] = {0., 0., 0.};
      
      adept::adouble solxNewg[3] = {0., 0., 0.};
      double solxOldg[3] = {0., 0., 0.};
      
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


      adept::adouble g[dim][dim] = {{0., 0.}, {0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      adept::adouble detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];

      adept::adouble normal[DIM];
      normal[0] = normalSign * (solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) / sqrt (detg);
      normal[1] = normalSign * (solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);
      normal[2] = normalSign * (solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);
      
      adept::adouble YdotN = 0.;
      adept::adouble YdotY = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        YdotN += solYg[K] * normal[K];
        YdotY += solYg[K] * solYg[K];
      }
      double signYdotN = ( YdotN.value() >= 0.)? 1. : -1.;

      adept::adouble sumP1 = 0.;
      adept::adouble sumP2 = 0.;
      adept::adouble sumP3 = 0.;
      for (unsigned p = 0; p < 3; p++) {
        double signP = (P[p]%2u == 0) ? 1. : signYdotN;
        sumP1 += signP * ap[p] * P[p] * pow (YdotY, (P[p] - 2.)/2.);
        sumP2 += signP * ap[p] * (1. - P[p]) * pow (YdotY , P[p]/2.);
        //sumP2 += signP * (ap[p] - ap[p] * P[p]) * pow (YdotY , P[p]/2.);
        sumP3 += signP * ap[p] * pow (YdotY, P[p]/2.);
      }

      adept::adouble gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;

      adept::adouble Area = weight * sqrt (detg);

      adept::adouble Jir[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            Jir[i][J] += gi[i][k] * solx_uv[J][k];
          }
        }
      }

      adept::adouble solx_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solxNew_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solxOld_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solWNew_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solW_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solWOld_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};


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

      std::vector < adept::adouble > phiW_Xtan[DIM];
      std::vector < adept::adouble > phix_Xtan[DIM];

      for (unsigned J = 0; J < DIM; J++) {
        phix_Xtan[J].assign (nxDofs, 0.);
        for (unsigned inode  = 0; inode < nxDofs; inode++) {
          for (unsigned k = 0; k < dim; k++) {
            phix_Xtan[J][inode] += phix_uv[k][inode] * Jir[k][J];
          }
        }

        phiW_Xtan[J].assign (nWDofs, 0.);
        for (unsigned inode  = 0; inode < nWDofs; inode++) {
          for (unsigned k = 0; k < dim; k++) {
            phiW_Xtan[J][inode] += phiW_uv[k][inode] * Jir[k][J];
          }
        }
      }

      for (unsigned K = 0; K < DIM; K++) {
        // curvature equation
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;
          for (unsigned J = 0; J < DIM; J++) {
            term1 +=  solxNew_Xtan[K][J] * phix_Xtan[J][i];
          }
          aResx[K][i] += (solYg[K] * phix[i] + term1) * Area;
        }
        
        // W - sum_P aP P <Y,N>^{P-2} Y = 0
        for (unsigned i = 0; i < nWDofs; i++) {
          aResY[K][i] += (solWg[K] - sumP1 * solYg[K]) * phiY[i] * Area;
        }

        //P-Willmore equation
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
          aResW[K][i] += ( ( ( solLambda1 /*- YdotN * solLambda2*/ ) * normal[K] + (solxNewg[K] - solxOldg[K])  / dt) * phiW[i]
                           + solLambda2 * term1
                           - term0
                           + sumP2 * term1
                           - term2 * phiW_Xtan[K][i]
                           + term3
                         ) * Area;
        }
        // Lagrange Multiplier Euquation Dx . N =0
        if (volumeConstraint) {
          aResLambda1 += ( (solxNewg[K] - solxOldg[K]) * normal[K]) * Area;
        }
        if(areaConstraint){
          //aResLambda2 += ( -YdotN * (solxNewg[K] - solxOldg[K]) * normal[K]) * Area;
          
          adept::adouble term1t = 0.;
          for (unsigned J = 0; J < DIM; J++) {
            term1t +=  solx_Xtan[K][J] * (solxNew_Xtan[K][J] - solxOld_Xtan[K][J]) ;
          }
          aResLambda2 += term1t * Area;
          
        }
      }

      for (unsigned K = 0; K < DIM; K++) {
        surface += Area.value();
      }
      
      
      for (unsigned K = 0; K < DIM; K++) {
        volume += normalSign * (solxNewg[K].value()  * normal[K].value()) * Area.value();
      }
      energy += sumP3.value() * Area.value();

    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
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

    Jac.resize (sizeAll * sizeAll);

    // define the dependent variables
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
    if(areaConstraint){
      s.dependent (&aResLambda2, 1);
    }
    // define the dependent variables
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
    if(areaConstraint){
      s.independent (&solLambda2, 1);
    }
    
    // get the jacobian matrix (ordered by row)
    s.jacobian (&Jac[0], true);

    KK->add_matrix_blocked (Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();


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


  // ***************** END ASSEMBLY *******************
}



void AssembleInit (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem> ("Init");   // pointer to the linear implicit system named "Poisson"

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

  unsigned solxType;
  solxType = mlSol->GetSolutionType (solDxIndex[0]);  // get the finite element type for "U"

  std::vector < double > solx[DIM];  // surface coordinates

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  unsigned solYIndex[DIM];
  solYIndex[0] = mlSol->GetIndex ("Y1");   // get the position of "Y1" in the ml_sol object
  solYIndex[1] = mlSol->GetIndex ("Y2");   // get the position of "Y2" in the ml_sol object
  solYIndex[2] = mlSol->GetIndex ("Y3");   // get the position of "Y3" in the ml_sol object

  unsigned solYType;
  solYType = mlSol->GetSolutionType (solYIndex[0]);  // get the finite element type for "Y"

  unsigned solYPdeIndex[DIM];
  solYPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("Y1");   // get the position of "Y1" in the pdeSys object
  solYPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("Y2");   // get the position of "Y2" in the pdeSys object
  solYPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("Y3");   // get the position of "Y3" in the pdeSys object

  unsigned solWIndex[DIM];
  solWIndex[0] = mlSol->GetIndex ("W1");   // get the position of "W1" in the ml_sol object
  solWIndex[1] = mlSol->GetIndex ("W2");   // get the position of "W2" in the ml_sol object
  solWIndex[2] = mlSol->GetIndex ("W3");   // get the position of "W3" in the ml_sol object

  unsigned solWType;
  solWType = mlSol->GetSolutionType (solWIndex[0]);  // get the finite element type for "W"

  unsigned solWPdeIndex[DIM];
  solWPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("W1");   // get the position of "W1" in the pdeSys object
  solWPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("W2");   // get the position of "W2" in the pdeSys object
  solWPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("W3");   // get the position of "W3" in the pdeSys object

  std::vector < adept::adouble > solY[DIM]; // local Y solution
  std::vector < adept::adouble > solW[DIM]; // local W solution

  std::vector< int > SYSDOF; // local to global pdeSys dofs

  vector< double > Res; // local redidual vector
  std::vector< adept::adouble > aResY[3]; // local redidual vector
  std::vector< adept::adouble > aResW[3]; // local redidual vector

  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)

  KK->zero();  // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solxType);   // number of solution element dofs
    unsigned nYDofs  = msh->GetElementDofNumber (iel, solYType);   // number of solution element dofs
    unsigned nWDofs  = msh->GetElementDofNumber (iel, solWType);   // number of solution element dofs

    for (unsigned K = 0; K < DIM; K++) {
      solx[K].resize (nxDofs);
      solY[K].resize (nYDofs);
      solW[K].resize (nWDofs);
    }

    // resize local arrays
    SYSDOF.resize (DIM * (nYDofs + nWDofs));

    Res.resize (DIM * (nYDofs + nWDofs));       //resize

    for (unsigned K = 0; K < DIM; K++) {
      aResY[K].assign (nYDofs, 0.);  //resize and zet to zero
      aResW[K].assign (nWDofs, 0.);  //resize and zet to zero
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nxDofs; i++) {
      unsigned iDDof = msh->GetSolutionDof (i, iel, solxType); // global to local mapping between solution node and solution dof
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);
      for (unsigned K = 0; K < DIM; K++) {
        solx[K][i] = (*msh->_topology->_Sol[K]) (iXDof) + (*sol->_Sol[solDxIndex[K]]) (iDDof);
      }
    }


    // local storage of global mapping and solution
    for (unsigned i = 0; i < nYDofs; i++) {
      unsigned iYDof = msh->GetSolutionDof (i, iel, solYType); // global to local mapping between solution node and solution dof
      for (unsigned K = 0; K < DIM; K++) {
        solY[K][i] = (*sol->_Sol[solYIndex[K]]) (iYDof); // global to local solution
        SYSDOF[ K * nYDofs + i] = pdeSys->GetSystemDof (solYIndex[K], solYPdeIndex[K], i, iel); // global to global mapping between solution node and pdeSys dof
      }
    }


    // local storage of global mapping and solution
    for (unsigned i = 0; i < nWDofs; i++) {
      unsigned iWDof = msh->GetSolutionDof (i, iel, solWType); // global to local mapping between solution node and solution dof
      for (unsigned K = 0; K < DIM; K++) {
        solW[K][i] = (*sol->_Sol[solWIndex[K]]) (iWDof); // global to local solution
        SYSDOF[ DIM * nYDofs + K * nWDofs + i] = pdeSys->GetSystemDof (solWIndex[K], solWPdeIndex[K], i, iel); // global to global mapping between solution node and pdeSys dof
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solxType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phix_uv[dim]; // local test function first order partial derivatives

      const double *phiY;  // local test function
      const double *phiY_uv[dim]; // local test function first order partial derivatives

      const double *phiW;  // local test function

      double weight; // gauss point weight

      // *** get gauss point weight, test function and test function partial derivatives ***
      phix = msh->_finiteElement[ielGeom][solxType]->GetPhi (ig);
      phix_uv[0] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDXi (ig); //derivative in u
      phix_uv[1] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDEta (ig); //derivative in v

      phiY = msh->_finiteElement[ielGeom][solYType]->GetPhi (ig);
      phiY_uv[0] = msh->_finiteElement[ielGeom][solYType]->GetDPhiDXi (ig);
      phiY_uv[1] = msh->_finiteElement[ielGeom][solYType]->GetDPhiDEta (ig);

      phiW = msh->_finiteElement[ielGeom][solWType]->GetPhi (ig);

      weight = msh->_finiteElement[ielGeom][solxType]->GetGaussWeight (ig);

      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solY_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      double solxg[3] = {0., 0., 0.};
      adept::adouble solYg[3] = {0., 0., 0.};
      adept::adouble solWg[3] = {0., 0., 0.};


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
            solx_uv[K][j]    += phix_uv[j][i] * solx[K][i];
          }
        }
        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nYDofs; i++) {
            solY_uv[K][j] += phiY_uv[j][i] * solY[K][i];
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

      adept::adouble normal[DIM];
      normal[0] = normalSign * (solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) / sqrt (detg);
      normal[1] = normalSign * (solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);
      normal[2] = normalSign * (solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);

      adept::adouble YdotN = 0.;
      adept::adouble YdotY = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        YdotN += solYg[K] * normal[K];
        YdotY += solYg[K] * solYg[K];
      }
      double signYdotN = ( YdotN.value() >= 0.)? 1. : -1.;
      adept::adouble sumP1 = 0.;
      for (unsigned p = 0; p < 3; p++) {
        double signP = (P[p]%2u == 0) ? 1. : signYdotN;
        sumP1 += signP * ap[p] * P[p] * pow (YdotY, (P[p] - 2.)/2.);
      }


      double gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;

      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          double id = 0.;
          for (unsigned j2 = 0; j2 < dim; j2++) {
            id +=  g[i][j2] * gi[j2][j];
          }
          if (i == j && fabs (id - 1.) > 1.0e-10) std::cout << id << " error ";
          else if (i != j && fabs (id) > 1.0e-10) std::cout << id << " error ";
        }
      }

      double Area = weight * sqrt (detg);

      double Jir[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            Jir[i][J] += gi[i][k] * solx_uv[J][k];
          }
        }
      }

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


      std::vector < double > phix_Xtan[DIM];
      std::vector < double > phiY_Xtan[DIM];

      for (unsigned J = 0; J < DIM; J++) {
        phix_Xtan[J].assign (nxDofs, 0.);
        for (unsigned inode  = 0; inode < nxDofs; inode++) {
          for (unsigned k = 0; k < dim; k++) {
            phix_Xtan[J][inode] += phix_uv[k][inode] * Jir[k][J];
          }
        }

        phiY_Xtan[J].assign (nYDofs, 0.);
        for (unsigned inode  = 0; inode < nYDofs; inode++) {
          for (unsigned k = 0; k < dim; k++) {
            phiY_Xtan[J][inode] += phiY_uv[k][inode] * Jir[k][J];
          }
        }
      }

      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nYDofs; i++) {
          adept::adouble term1 = 0.;
          adept::adouble term2 = 0.;
          for (unsigned J = 0; J < DIM; J++) {
            term1 +=  solx_Xtan[K][J] * phiY_Xtan[J][i];
            term2 +=  solY_Xtan[K][J] * phiY_Xtan[J][i];
          }

          double eps = 0.005;
          //           std::cout << A <<" ";
          aResY[K][i] += (solYg[K] * phiY[i] + eps * term2 + term1) * Area;
        }

        for (unsigned i = 0; i < nWDofs; i++) {
          aResW[K][i] += (solWg[K] - sumP1 * solYg[K]) * phiW[i] * weight;
        }
      }
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store


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

    Jac.resize (DIM * (nYDofs + nWDofs) * DIM * (nYDofs + nWDofs));

    // define the dependent variables

    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResY[K][0], nYDofs);
    }

    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResW[K][0], nWDofs);
    }

    // define the dependent variables

    for (int K = 0; K < DIM; K++) {
      s.independent (&solY[K][0], nYDofs);
    }

    for (int K = 0; K < DIM; K++) {
      s.independent (&solW[K][0], nWDofs);
    }

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

void CopyDisplacement(MultiLevelSolution &mlSol,  const bool &forward) {
  
  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  
  Solution* solution  = mlSol.GetSolutionLevel ( level );
  Mesh* msh = mlSol._mlMesh->GetLevel ( level );
    
  unsigned DIM = 3u;
  vector < unsigned > solDxIndex ( DIM );
  solDxIndex[0] = mlSol.GetIndex ( "Dx1" ); 
  solDxIndex[1] = mlSol.GetIndex ( "Dx2" ); 
  solDxIndex[2] = mlSol.GetIndex ( "Dx3" ); 
  
  vector < unsigned > solNDxIndex ( DIM );
  solNDxIndex[0] = mlSol.GetIndex ( "nDx1" ); 
  solNDxIndex[1] = mlSol.GetIndex ( "nDx2" ); 
  solNDxIndex[2] = mlSol.GetIndex ( "nDx3" ); 
  
  if(forward){  
    for(unsigned i = 0; i < DIM; i++) {
      *(solution->_Sol[solNDxIndex[i]]) = *(solution->_Sol[solDxIndex[i]]);
    }
  }
  else {  
    for(unsigned i = 0; i < DIM; i++) {
      *(solution->_Sol[solDxIndex[i]]) = *(solution->_Sol[solNDxIndex[i]]);
    }
  }
  
  
}



void AssembleNProj (MultiLevelProblem& ml_prob) {
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
    aResL.assign(nLDofs, 0.);
    
    
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
            term1 +=  ( (solNDx_Xtan[K][J] - 0*solDxOld_Xtan[K][J]) +
                        (solNDx_Xtan[J][K] - 0*solDxOld_Xtan[J][K]) )* phix_Xtan[J][i];
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
    
    
    
    Jac.resize (  (DIM * nxDofs + nLDofs) * (DIM * nxDofs + nLDofs));
    
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



