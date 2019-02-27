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

const double P = 2;
using namespace femus;

const bool volumeConstraint = true;


void AssemblePWillmore (MultiLevelProblem&);

void AssembleInit (MultiLevelProblem&);

double GetTimeStep (const double time) {
  return 0.0000000000000005;
}

bool SetBoundaryCondition (const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true;
  
  if (!strcmp ("Dx1", SolName) || !strcmp ("Dx2", SolName) || !strcmp ("Dx3", SolName)) {
    value = 0.;
  }
  else if (!strcmp ("W1", SolName)) {
    value = -2. * x[0];
  }
  else if (!strcmp ("W2", SolName)) {
    value = -2. * x[1];
  }
  else if (!strcmp ("W3", SolName)) {
    value = 0.;
  }
  return dirichlet;
}

double InitalValueW1 (const std::vector < double >& x) {
  return -2. * pow (2., P - 2) * x[0];
}

double InitalValueW2 (const std::vector < double >& x) {
  return -2. * pow (2., P - 2) * x[1];
}

double InitalValueW3 (const std::vector < double >& x) {
  return -2. * pow (2., P - 2) * x[2];
}


double InitalValueX1 (const std::vector < double >& x) {
  return x[0];
}

double InitalValueX2 (const std::vector < double >& x) {
  return x[1];
}

double InitalValueX3 (const std::vector < double >& x) {
  return x[2];
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
  mlMsh.ReadCoarseMesh ("./input/ellipsoidRef3.neu", "seventh", scalingFactor);
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
  //   mlSol.AddSolution ("Dx1", LAGRANGE, SECOND, 2);
  //   mlSol.AddSolution ("Dx2", LAGRANGE, SECOND, 2);
  //   mlSol.AddSolution ("Dx3", LAGRANGE, SECOND, 2);
  
  mlSol.AddSolution("Dx1", LAGRANGE, FIRST, 2);
  mlSol.AddSolution("Dx2", LAGRANGE, FIRST, 2);
  mlSol.AddSolution("Dx3", LAGRANGE, FIRST, 2);
  //
  //   mlSol.AddSolution ("W1", LAGRANGE, SECOND, 2);
  //   mlSol.AddSolution ("W2", LAGRANGE, SECOND, 2);
  //   mlSol.AddSolution ("W3", LAGRANGE, SECOND, 2);
  
  mlSol.AddSolution ("Lambda", DISCONTINOUS_POLYNOMIAL, ZERO, 0);
  
  mlSol.AddSolution("W1", LAGRANGE, FIRST, 2);
  mlSol.AddSolution("W2", LAGRANGE, FIRST, 2);
  mlSol.AddSolution("W3", LAGRANGE, FIRST, 2);
  
  
  mlSol.Initialize ("All");
  mlSol.Initialize ("W1", InitalValueW1);
  mlSol.Initialize ("W2", InitalValueW2);
  mlSol.Initialize ("W3", InitalValueW3);
  
  //   mlSol.Initialize("Dx1", InitalValueX1);
  //   mlSol.Initialize("Dx2", InitalValueX2);
  //   mlSol.Initialize("Dx3", InitalValueX3);
  
  
  mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
  mlSol.GenerateBdc ("All");
  
  //   mlSol.FixSolutionAtOnePoint("Dx1");
  //   mlSol.FixSolutionAtOnePoint("Dx2");
  //   mlSol.FixSolutionAtOnePoint("Dx3");
  //   mlSol.FixSolutionAtOnePoint("W1");
  //   mlSol.FixSolutionAtOnePoint("W2");
  //   mlSol.FixSolutionAtOnePoint("W3");
  
  
  MultiLevelProblem mlProb (&mlSol);
  
  
  //BEGIN Get initial curbvature data
  NonLinearImplicitSystem& system0 = mlProb.add_system < NonLinearImplicitSystem > ("Init");
  
  // add solution "X", "Y", "Z" and "H" to the system
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
  //END Get initial curbvature data
  
  
  
  // add system Wilmore in mlProb as a Linear Implicit System
  TransientNonlinearImplicitSystem& system = mlProb.add_system < TransientNonlinearImplicitSystem > ("PWillmore");
  
  // add solution "X", "Y", "Z" and "H" to the system
  system.AddSolutionToSystemPDE ("Lambda");
  system.AddSolutionToSystemPDE ("Dx1");
  system.AddSolutionToSystemPDE ("Dx2");
  system.AddSolutionToSystemPDE ("Dx3");
  system.AddSolutionToSystemPDE ("W1");
  system.AddSolutionToSystemPDE ("W2");
  system.AddSolutionToSystemPDE ("W3");
  
  system.SetMaxNumberOfNonLinearIterations (20);
  system.SetNonLinearConvergenceTolerance (1.e-9);
  
  // attach the assembling function to system
  system.SetAssembleFunction (AssemblePWillmore);
  
  system.AttachGetTimeIntervalFunction (GetTimeStep);
  
  // initilaize and solve the system
  system.init();
  
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
  
  unsigned numberOfTimeSteps = 100000;
  for (unsigned time_step = 0; time_step < numberOfTimeSteps; time_step++) {
    
    system.CopySolutionToOldSolution();
    system.MGsolve();
    if ( (time_step + 1) % 1 == 0)
      mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "linear", variablesToBePrinted, (time_step + 1) / 1);
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
  
  double dt = GetTimeStep (0);
  
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
  
  unsigned solLambdaIndex;
  solLambdaIndex = mlSol->GetIndex ("Lambda");
  unsigned solLambdaType = mlSol->GetSolutionType (solLambdaIndex);
  
  unsigned solLambaPdeIndex;
  solLambaPdeIndex = mlPdeSys->GetSolPdeIndex ("Lambda");
  
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
  
  
  std::vector < adept::adouble > solLambda; // local lambda solution
  std::vector < adept::adouble > solW[DIM]; // local W solution
  std::vector < double > solWOld[DIM];  // surface coordinates
  std::vector< int > SYSDOF; // local to global pdeSys dofs
  
  vector< double > Res; // local redidual vector
  std::vector< adept::adouble > aResLambda;
  std::vector< adept::adouble > aResx[3]; // local redidual vector
  std::vector< adept::adouble > aResW[3]; // local redidual vector
  adept::adouble aResLambda0;
  
  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)
  
  //MatSetOption ( ( static_cast<PetscMatrix*> ( KK ) )->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE );
  
  KK->zero();  // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual
  
  unsigned lambda0SolDof;
  double solLambda0;
  unsigned lambda0PdeDof;
  if (iproc == 0) {
    lambda0SolDof = msh->GetSolutionDof (0, 0, solLambdaType); // global to local mapping between solution node and solution dof
    solLambda0 = (*sol->_Sol[solLambdaIndex]) (lambda0SolDof); // global to local solution
    lambda0PdeDof = pdeSys->GetSystemDof (solLambdaIndex, solLambaPdeIndex, 0, 0);
  }
  MPI_Bcast (&lambda0SolDof, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast (&solLambda0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&lambda0PdeDof, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  
  double volume = 0.;
  
  double energy = 0.;
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    
    bool iel0 = (iel == 0) ? true : false;
    
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solxType);   // number of solution element dofs
    unsigned nWDofs  = msh->GetElementDofNumber (iel, solWType);   // number of solution element dofs
    unsigned nLambdaDofs  = msh->GetElementDofNumber (iel, solLambdaType);   // number of solution element dofs
    
    for (unsigned K = 0; K < DIM; K++) {
      solx[K].resize (nxDofs);
      solxOld[K].resize (nxDofs);
      solW[K].resize (nWDofs);
      solWOld[K].resize (nWDofs);
    }
    solLambda.resize (nLambdaDofs);
    
    // resize local arrays
    SYSDOF.resize (DIM * (nxDofs + nWDofs) + nLambdaDofs + !iel0);
    
    Res.resize (DIM * (nxDofs + nWDofs) + nLambdaDofs + !iel0);       //resize
    
    for (unsigned K = 0; K < DIM; K++) {
      aResx[K].assign (nxDofs, 0.);  //resize and set to zero
      aResW[K].assign (nWDofs, 0.);  //resize and zet to zero
    }
    aResLambda.assign (nLambdaDofs, 0.);
    aResLambda0 = 0.;
    
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
    for (unsigned i = 0; i < nWDofs; i++) {
      unsigned iWDof = msh->GetSolutionDof (i, iel, solWType); // global to local mapping between solution node and solution dof
      for (unsigned K = 0; K < DIM; K++) {
        solWOld[K][i] = (*sol->_SolOld[solWIndex[K]]) (iWDof); // global to local solution
        solW[K][i] = (*sol->_Sol[solWIndex[K]]) (iWDof); // global to local solution
        SYSDOF[DIM * nxDofs + K * nWDofs + i] = pdeSys->GetSystemDof (solWIndex[K], solWPdeIndex[K], i, iel); // global to global mapping between solution node and pdeSys dof
      }
    }
    
    
    for (unsigned i = 0; i < nLambdaDofs; i++) {
      unsigned iLambdaDof = msh->GetSolutionDof (i, iel, solLambdaType); // global to local mapping between solution node and solution dof
      solLambda[i] = (*sol->_Sol[solLambdaIndex]) (iLambdaDof); // global to local solution
      SYSDOF[DIM * nxDofs + DIM * nWDofs + i] = pdeSys->GetSystemDof (solLambdaIndex, solLambaPdeIndex, i, iel); // global to global mapping between solution node and pdeSys dof
    }
    
    adept::adouble lambda0 = solLambda0;
    if (!iel0) SYSDOF[DIM * (nxDofs + nWDofs) + nLambdaDofs ] = lambda0PdeDof;
    
    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();
    if (iel != 0 || !volumeConstraint) aResLambda[0] = solLambda[0];
    
    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solxType]->GetGaussPointNumber(); ig++) {
      
      const double *phix;  // local test function
      const double *phix_uv[dim]; // local test function first order partial derivatives
      
      const double *phiW;  // local test function
      const double *phiW_uv[dim]; // local test function first order partial derivatives
      
      double weight; // gauss point weight
      
      // *** get gauss point weight, test function and test function partial derivatives ***
      phix = msh->_finiteElement[ielGeom][solxType]->GetPhi (ig);
      phix_uv[0] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDXi (ig); //derivative in u
      phix_uv[1] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDEta (ig); //derivative in v
      
      phiW = msh->_finiteElement[ielGeom][solWType]->GetPhi (ig);
      phiW_uv[0] = msh->_finiteElement[ielGeom][solWType]->GetDPhiDXi (ig);
      phiW_uv[1] = msh->_finiteElement[ielGeom][solWType]->GetDPhiDEta (ig);
      
      weight = msh->_finiteElement[ielGeom][solxType]->GetGaussWeight (ig);
      
      double sumPhi = 0.;
      double sumPhi_u = 0.;
      double sumPhi_v = 0.;
      for (unsigned i = 0; i < nxDofs; i++) {
        //std::cout << phix[i] << " " << phix_uv[0][i] << " " << phix_uv[1][i] <<std::endl;
        sumPhi += phix[i];
        sumPhi_u += phix_uv[0][i];
        sumPhi_v += phix_uv[1][i];
      }
      if (fabs (sumPhi - 1.) > 1.0e-12 || fabs (sumPhi_u) > 1.0e-12 || fabs (sumPhi_v) > 1.0e-12) {
        std::cout << sumPhi << " " << sumPhi_u << " " << sumPhi_v << std::endl;
      }
      
      adept::adouble solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solW_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solWg[3] = {0., 0., 0.};
      adept::adouble solxg[3] = {0., 0., 0.};
      
      double solxOld_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solWOld_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solWOldg[3] = {0., 0., 0.};
      double solxOldg[3] = {0., 0., 0.};
      
      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          solxg[K] += phix[i] * solx[K][i];
          solxOldg[K] += phix[i] * solxOld[K][i];
        }
        for (unsigned i = 0; i < nWDofs; i++) {
          solWg[K] += phiW[i] * solW[K][i];
          solWOldg[K] += phiW[i] * solWOld[K][i];
        }
        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * solx[K][i];
            solxOld_uv[K][j] += phix_uv[j][i] * solxOld[K][i];
          }
        }
        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nWDofs; i++) {
            solW_uv[K][j] += phiW_uv[j][i] * solW[K][i];
            solWOld_uv[K][j] += phiW_uv[j][i] * solWOld[K][i];
          }
        }
      }
      
      //       double solXOldnormPm2 = 0.;
      //       for (unsigned K = 0; K < DIM; K++) {
      //         solXOldnormPm2 += solxOldg[K] * solxOldg[K];
      //       }
      //       solXOldnormPm2 =  pow (solXOldnormPm2 , (P - 2.) / 2.);
      //       solXOldnormPm2 = 1.;
      
      adept::adouble normY = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        normY += solWg[K] * solWg[K];
      }
      normY = pow (normY , 1. / ( 2. * ( P - 1.) ) );
      //normY =  pow (normY , 1. / (2.* P));
      //std::cout << normY << " ";
      
      double g[dim][dim] = {{0., 0.}, {0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned K = 0; K < DIM; K++) {
            g[i][j] += solxOld_uv[K][i] * solxOld_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      
      double normal[DIM];
      normal[0] = (solxOld_uv[1][0] * solxOld_uv[2][1] - solxOld_uv[2][0] * solxOld_uv[1][1]) / sqrt (detg);
      normal[1] = (solxOld_uv[2][0] * solxOld_uv[0][1] - solxOld_uv[0][0] * solxOld_uv[2][1]) / sqrt (detg);;
      normal[2] = (solxOld_uv[0][0] * solxOld_uv[1][1] - solxOld_uv[1][0] * solxOld_uv[0][1]) / sqrt (detg);;
      
      std::cout << normal[0] <<" "<<normal[1]<<" " <<normal[2]<<std::endl;
      
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
          //           if(i == j && fabs(id-1.) > 1.0e-10) std::cout<<id << " error0 ";
          //           else if (i != j && fabs(id) > 1.0e-10) std::cout<< id <<" error0 ";
        }
      }
      
      
      double Area = weight * sqrt (detg);
      
      double Jir[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            Jir[i][J] += gi[i][k] * solxOld_uv[J][k];
          }
        }
      }
      
      adept::adouble solx_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      double solxOld_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solW_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      double solWOld_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      
      
      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            solx_Xtan[I][J] += solx_uv[I][k] * Jir[k][J];
            solxOld_Xtan[I][J] += solxOld_uv[I][k] * Jir[k][J];
            solW_Xtan[I][J] += solW_uv[I][k] * Jir[k][J];
            solWOld_Xtan[I][J] += solWOld_uv[I][k] * Jir[k][J];
          }
        }
      }

      std::vector < double > phiW_Xtan[DIM];
      std::vector < double > phix_Xtan[DIM];
      
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
      
      double D[3][3];
      double Pr[3][3];
      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          D[I][J] = phiW_Xtan[I][0] + phiW_Xtan[J][0];
          double deltaij = (I == J) ? 1. : 0.;
          Pr[I][J] = deltaij - normal[I] * normal[J];
        }
      }
      double D1[3][3];
      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          D1[I][J] = 0.;
          for (unsigned K = 0; K < DIM; K++) {
            D1[I][J] += solxOld_Xtan[I][K] * D[K][J]; 
          }
        }
      }
      double D2[3][3];
      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          D2[I][J] = 0.;
          for (unsigned K = 0; K < DIM; K++) {
            D2[I][J] += D1[I][K] * solxOld_Xtan[K][J]; 
          }
          //std::cout << D[I][J]<< " ";
        }
        //std::cout << std::endl;
      }
      //std::cout << std::endl;
      
      
      double DX[3][3];
      double DX2[3][3];
      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          DX[I][J] = 0.;  
          DX2[I][J] = 0.; 
          for (unsigned K = 0; K < DIM; K++) {
            DX[I][J] += D[I][K] * solWOld_Xtan[K][J];
            DX2[I][J] += D2[I][K] * solWOld_Xtan[K][J];
          }
        }
      }
      double DXW = 0;
      double DXW2 = 0;
      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          DXW += DX[I][J] * solxOld_Xtan[I][J];
          DXW2 += DX2[I][J] * solxOld_Xtan[I][J];
          //std::cout << solxOld_Xtan[I][J] -Pr[I][J] << " ";
        }
        //std::cout << std::endl;
      }
      //std::cout << std::endl;
      
      //std::cout << DXW << " " << DXW2 << "\t";
      
      double Pr2[3][3];
      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          Pr2[I][J] = 0.;  
          for (unsigned K = 0; K < DIM; K++) {
            Pr2[I][J] += Pr[I][K] * solWOld_Xtan[K][J];
          }
          std::cout << solWOld_Xtan[I][J] - Pr2[I][J] <<" ";
        }
        std::cout<<std::endl;
      }
      std::cout<<std::endl;
      
      exit(0);
      
      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;
          for (unsigned J = 0; J < DIM; J++) {
            term1 +=  solx_Xtan[K][J] * phix_Xtan[J][i];
          }
          aResW[K][i] += (solWg[K] * pow (normY, 2. - P)  * phix[i] + term1) * Area;
        }
        for (unsigned i = 0; i < nWDofs; i++) {
          adept::adouble term0 = 0.;
          adept::adouble term1 = 0.;
          adept::adouble term2 = 0.;
          adept::adouble term3 = 0.;
          
          adept::adouble term5 = 0.;
          
          for (unsigned J = 0; J < DIM; J++) {
            
            
            
            term0 +=  solW_Xtan[K][J] * phiW_Xtan[J][i];
            
            term1 +=  solx_Xtan[K][J] * phiW_Xtan[J][i];
            
            term2 +=  solW_Xtan[J][J];
            
            adept::adouble term4 = 0.;
            
            term5 += phiW_Xtan[J][i] * solW_Xtan[J][K];
            
            for (unsigned L = 0; L < DIM; L++) {
              term5 -= normal[K] * normal[J] * phiW_Xtan[L][i] * solW_Xtan[J][L];
              term4 += solx_Xtan[J][L] * solW_Xtan[K][L] + solx_Xtan[K][L] * solW_Xtan[J][L];
            }
            
            term3 += phiW_Xtan[J][i] * term4;
            
          }
          aResx[K][i] += (P * (volumeConstraint * lambda0 * normal[K] +
          (solxg[K] - solxOldg[K])  / dt) * phiW[i]
          //pow (normY, P - 2.) * (solxg[K] - solxOldg[K])  / dt) * phiW[i]
          - P * term0
          + (1. - P) * pow (normY , P) * term1
          - P * term2 * phiW_Xtan[K][i]
          //+ 0.5 * P * term3 + 0. * P * term3.value()
          //+ P * term5
          ) * Area;
          
          //std::cout << term5 << " " << -term0 + term3 << "\t";              
        }
        
        
      }
      for (unsigned K = 0; K < DIM; K++) {
        if (volumeConstraint) {
          if (iel != 0)aResLambda0 += ( (solxg[K] - solxOldg[K]) * normal[K]) * Area;
          else aResLambda[0] += ( (solxg[K] - solxOldg[K]) * normal[K]) * Area;
        }
        
        volume += (solxg[K].value()  * normal[K]) * Area;
      }
      
      energy += pow (normY.value(), P) * Area;
      
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
      for (int i = 0; i < nWDofs; i++) {
        Res[DIM * nxDofs + K * nWDofs + i] = -aResW[K][i].value();
      }
    }
    for (int i = 0; i < nLambdaDofs; i++) {
      Res[DIM * nxDofs + DIM * nWDofs + i] = -aResLambda[i].value();
    }
    if (!iel0) Res[DIM * (nxDofs + nWDofs) +  nLambdaDofs ] = - aResLambda0.value();
    
    
    RES->add_vector_blocked (Res, SYSDOF);
    
    Jac.resize ( (DIM * (nxDofs + nWDofs) + nLambdaDofs + !iel0) * (DIM * (nxDofs + nWDofs) + nLambdaDofs + !iel0));
    
    // define the dependent variables
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResx[K][0], nxDofs);
    }
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResW[K][0], nWDofs);
    }
    s.dependent (&aResLambda[0], nLambdaDofs);
    
    if (!iel0) s.dependent (&aResLambda0, 1);
    
    // define the dependent variables
    for (int K = 0; K < DIM; K++) {
      s.independent (&solx[K][0], nxDofs);
    }
    for (int K = 0; K < DIM; K++) {
      s.independent (&solW[K][0], nWDofs);
    }
    s.independent (&solLambda[0], nLambdaDofs);
    if (!iel0) s.independent (&lambda0, 1);
    
    // get the jacobian matrix (ordered by row)
    s.jacobian (&Jac[0], true);
    
    KK->add_matrix_blocked (Jac, SYSDOF, SYSDOF);
    
    s.clear_independents();
    s.clear_dependents();
    
  } //end element loop for each process
  
  RES->close();
  KK->close();
  
  
  double volumeAll;
  MPI_Reduce (&volume, &volumeAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  std::cout << "VOLUME = " << volumeAll << std::endl;
  
  
  double energyAll;
  MPI_Reduce (&energy, &energyAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  std::cout << "ENERGY = " << energyAll << std::endl;
  
  //VecView ( (static_cast<PetscVector*> (RES))->vec(),  PETSC_VIEWER_STDOUT_SELF);
  //MatView ( (static_cast<PetscMatrix*> (KK))->mat(), PETSC_VIEWER_STDOUT_SELF);
  
  //   PetscViewer    viewer;
  //   PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
  //   PetscObjectSetName ( (PetscObject) viewer, "PWilmore matrix");
  //   PetscViewerPushFormat (viewer, PETSC_VIEWER_DRAW_LG);
  //   MatView ( (static_cast<PetscMatrix*> (KK))->mat(), viewer);
  //   double a;
  //   std::cin >> a;
  
  
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
  
  std::vector< int > SYSDOF; // local to global pdeSys dofs
  
  vector< double > Res; // local redidual vector
  std::vector< adept::adouble > aResW[3]; // local redidual vector
  
  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)
  
  KK->zero();  // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual
  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solxType);   // number of solution element dofs
    unsigned nWDofs  = msh->GetElementDofNumber (iel, solWType);   // number of solution element dofs
    
    for (unsigned K = 0; K < DIM; K++) {
      solx[K].resize (nxDofs);
      solW[K].resize (nWDofs);
    }
    
    // resize local arrays
    SYSDOF.resize (DIM * nWDofs);
    
    Res.resize (DIM * nWDofs);     //resize
    
    for (unsigned K = 0; K < DIM; K++) {
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
    for (unsigned i = 0; i < nWDofs; i++) {
      unsigned iWDof = msh->GetSolutionDof (i, iel, solWType); // global to local mapping between solution node and solution dof
      for (unsigned K = 0; K < DIM; K++) {
        solW[K][i] = (*sol->_Sol[solWIndex[K]]) (iWDof); // global to local solution
        SYSDOF[ K * nWDofs + i] = pdeSys->GetSystemDof (solWIndex[K], solWPdeIndex[K], i, iel); // global to global mapping between solution node and pdeSys dof
      }
    }
    
    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();
    
    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solxType]->GetGaussPointNumber(); ig++) {
      
      const double *phix;  // local test function
      const double *phix_uv[dim]; // local test function first order partial derivatives
      
      const double *phiW;  // local test function
      const double *phiW_uv[dim]; // local test function first order partial derivatives
      
      double weight; // gauss point weight
      
      // *** get gauss point weight, test function and test function partial derivatives ***
      phix = msh->_finiteElement[ielGeom][solxType]->GetPhi (ig);
      phix_uv[0] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDXi (ig); //derivative in u
      phix_uv[1] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDEta (ig); //derivative in v
      
      phiW = msh->_finiteElement[ielGeom][solWType]->GetPhi (ig);
      phiW_uv[0] = msh->_finiteElement[ielGeom][solWType]->GetDPhiDXi (ig);
      phiW_uv[1] = msh->_finiteElement[ielGeom][solWType]->GetDPhiDEta (ig);
      
      weight = msh->_finiteElement[ielGeom][solxType]->GetGaussWeight (ig);
      
      double sumPhi = 0.;
      double sumPhi_u = 0.;
      double sumPhi_v = 0.;
      for (unsigned i = 0; i < nxDofs; i++) {
        //std::cout << phix[i] << " " << phix_uv[0][i] << " " << phix_uv[1][i] <<std::endl;
        sumPhi += phix[i];
        sumPhi_u += phix_uv[0][i];
        sumPhi_v += phix_uv[1][i];
      }
      if (fabs (sumPhi - 1.) > 1.0e-12 || fabs (sumPhi_u) > 1.0e-12 || fabs (sumPhi_v) > 1.0e-12) {
        std::cout << sumPhi << " " << sumPhi_u << " " << sumPhi_v << std::endl;
      }
      
      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solW_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solxg[3] = {0., 0., 0.};
      adept::adouble solWg[3] = {0., 0., 0.};
      
      
      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          solxg[K] += phix[i] * solx[K][i];
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
          for (unsigned i = 0; i < nWDofs; i++) {
            solW_uv[K][j] += phiW_uv[j][i] * solW[K][i];
          }
        }
      }
      
      adept::adouble solWnorm2 = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        solWnorm2 += solWg[K] * solWg[K];
      }
      
      adept::adouble normY = pow (solWnorm2 , 1. / ( (2.* P - 2.)));
      //adept::adouble normY = pow (solWnorm2 , 1 / (2.* P));
      
      double g[dim][dim] = {{0., 0.}, {0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      
      //       double normal[DIM];
      //       normal[0] = - (solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) / sqrt (detg);
      //       normal[1] = - (solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);;
      //       normal[2] = - (solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);;
      
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
      adept::adouble solW_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      
      
      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            solx_Xtan[I][J] += solx_uv[I][k] * Jir[k][J];
            solW_Xtan[I][J] += solW_uv[I][k] * Jir[k][J];
          }
        }
      }
      
      
      std::vector < double > phiW_Xtan[DIM];
      std::vector < double > phix_Xtan[DIM];
      
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
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;
          for (unsigned J = 0; J < DIM; J++) {
            term1 +=  solx_Xtan[K][J] * phiW_Xtan[J][i];
          }
          if (fabs (term1 - phix_Xtan[K][i]) > 1.0e-10) {
            std::cout << " error " << term1 << " " << phiW_Xtan[K][i];
            //abort();
          }
          //           std::cout << A <<" ";
          aResW[K][i] += (solWg[K] * pow (normY.value(), 2. - P) * phiW[i] + term1) * Area;
        }
      }
    } // end gauss point loop
    
    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    
    //copy the value of the adept::adoube aRes in double Res and store
    
    
    for (int K = 0; K < DIM; K++) {
      for (int i = 0; i < nWDofs; i++) {
        Res[ K * nWDofs + i] = -aResW[K][i].value();
      }
    }
    
    RES->add_vector_blocked (Res, SYSDOF);
    
    Jac.resize (DIM * nWDofs * DIM * nWDofs);
    
    // define the dependent variables
    
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResW[K][0], nWDofs);
    }
    
    // define the dependent variables
    
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


