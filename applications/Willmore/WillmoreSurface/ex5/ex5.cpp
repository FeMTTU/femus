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

const unsigned P = 3;
using namespace femus;

const bool volumeConstraint = true;


void AssemblePWillmore (MultiLevelProblem&);

void AssembleInit (MultiLevelProblem&);

double GetTimeStep (const double time) {
  return 0.0002;
}

bool SetBoundaryCondition (const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true;

  if (!strcmp ("Dx1", SolName) || !strcmp ("Dx2", SolName) || !strcmp ("Dx3", SolName)) {
    value = 0.;
  }
  else if (!strcmp ("Y1", SolName)) {
    value = -2. * x[0];
  }
  else if (!strcmp ("Y2", SolName)) {
    value = -2. * x[1];
  }
  else if (!strcmp ("Y3", SolName)) {
    value = 0.;
  }
  return dirichlet;
}

double InitalValueY1 (const std::vector < double >& x) {
  return -2. * pow (2., P - 2) * x[0];
}

double InitalValueY2 (const std::vector < double >& x) {
  return -2. * pow (2., P - 2) * x[1];
}

double InitalValueY3 (const std::vector < double >& x) {
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
  mlMsh.ReadCoarseMesh ("./input/ellipsoidRef2.neu", "seventh", scalingFactor);
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
  mlSol.AddSolution ("Dx1", LAGRANGE, SECOND, 2);
  mlSol.AddSolution ("Dx2", LAGRANGE, SECOND, 2);
  mlSol.AddSolution ("Dx3", LAGRANGE, SECOND, 2);

//   mlSol.AddSolution("Dx1", LAGRANGE, FIRST, 2);
//   mlSol.AddSolution("Dx2", LAGRANGE, FIRST, 2);
//   mlSol.AddSolution("Dx3", LAGRANGE, FIRST, 2);
//
  mlSol.AddSolution ("Y1", LAGRANGE, SECOND, 2);
  mlSol.AddSolution ("Y2", LAGRANGE, SECOND, 2);
  mlSol.AddSolution ("Y3", LAGRANGE, SECOND, 2);

  mlSol.AddSolution ("Lambda", DISCONTINOUS_POLYNOMIAL, ZERO, 0);

//   mlSol.AddSolution("Y1", LAGRANGE, FIRST, 2);
//   mlSol.AddSolution("Y2", LAGRANGE, FIRST, 2);
//   mlSol.AddSolution("Y3", LAGRANGE, FIRST, 2);


  mlSol.Initialize ("All");
  mlSol.Initialize ("Y1", InitalValueY1);
  mlSol.Initialize ("Y2", InitalValueY2);
  mlSol.Initialize ("Y3", InitalValueY3);

//   mlSol.Initialize("Dx1", InitalValueX1);
//   mlSol.Initialize("Dx2", InitalValueX2);
//   mlSol.Initialize("Dx3", InitalValueX3);


  mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
  mlSol.GenerateBdc ("All");

//   mlSol.FixSolutionAtOnePoint("Dx1");
//   mlSol.FixSolutionAtOnePoint("Dx2");
//   mlSol.FixSolutionAtOnePoint("Dx3");
//   mlSol.FixSolutionAtOnePoint("Y1");
//   mlSol.FixSolutionAtOnePoint("Y2");
//   mlSol.FixSolutionAtOnePoint("Y3");


  MultiLevelProblem mlProb (&mlSol);


  //BEGIN Get initial curbvature data
  NonLinearImplicitSystem& system0 = mlProb.add_system < NonLinearImplicitSystem > ("Init");

  // add solution "X", "Y", "Z" and "H" to the system
  system0.AddSolutionToSystemPDE ("Y1");
  system0.AddSolutionToSystemPDE ("Y2");
  system0.AddSolutionToSystemPDE ("Y3");

  system0.SetMaxNumberOfNonLinearIterations (10);
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
  system.AddSolutionToSystemPDE ("Dx1");
  system.AddSolutionToSystemPDE ("Dx2");
  system.AddSolutionToSystemPDE ("Dx3");
  system.AddSolutionToSystemPDE ("Y1");
  system.AddSolutionToSystemPDE ("Y2");
  system.AddSolutionToSystemPDE ("Y3");
  system.AddSolutionToSystemPDE ("Lambda");

  system.SetMaxNumberOfNonLinearIterations (1);
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
  mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

  unsigned numberOfTimeSteps = 100000;
  for (unsigned time_step = 0; time_step < numberOfTimeSteps; time_step++) {

    system.CopySolutionToOldSolution();
    system.MGsolve();
    if ( (time_step + 1) % 1 == 0)
      mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, (time_step + 1) / 1);
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


  std::vector < adept::adouble > solLambda; // local Y solution
  std::vector < adept::adouble > solY[DIM]; // local Y solution
  std::vector < double > solYOld[DIM];  // surface coordinates
  std::vector< int > SYSDOF; // local to global pdeSys dofs

  vector< double > Res; // local redidual vector
  std::vector< adept::adouble > aResLambda;
  std::vector< adept::adouble > aResx[3]; // local redidual vector
  std::vector< adept::adouble > aResY[3]; // local redidual vector
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
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    bool iel0 = (iel == 0) ? true : false;

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solxType);   // number of solution element dofs
    unsigned nYDofs  = msh->GetElementDofNumber (iel, solYType);   // number of solution element dofs
    unsigned nLambdaDofs  = msh->GetElementDofNumber (iel, solLambdaType);   // number of solution element dofs

    for (unsigned K = 0; K < DIM; K++) {
      solx[K].resize (nxDofs);
      solxOld[K].resize (nxDofs);
      solY[K].resize (nYDofs);
      solYOld[K].resize (nYDofs);
    }
    solLambda.resize (nLambdaDofs);

    // resize local arrays
    SYSDOF.resize (DIM * (nxDofs + nYDofs) + nLambdaDofs + !iel0);

    Res.resize (DIM * (nxDofs + nYDofs) + nLambdaDofs + !iel0);       //resize

    for (unsigned K = 0; K < DIM; K++) {
      aResx[K].assign (nxDofs, 0.);  //resize and set to zero
      aResY[K].assign (nYDofs, 0.);  //resize and zet to zero
    }
    aResLambda.assign (nLambdaDofs, 0.);
    aResLambda0 = 0;

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


    for (unsigned i = 0; i < nLambdaDofs; i++) {
      unsigned iLambdaDof = msh->GetSolutionDof (i, iel, solLambdaType); // global to local mapping between solution node and solution dof
      solLambda[i] = (*sol->_Sol[solLambdaIndex]) (iLambdaDof); // global to local solution
      SYSDOF[DIM * nxDofs + DIM * nYDofs + i] = pdeSys->GetSystemDof (solLambdaIndex, solLambaPdeIndex, i, iel); // global to global mapping between solution node and pdeSys dof
    }

    adept::adouble lambda0 = solLambda0;
    if (!iel0) SYSDOF[DIM * (nxDofs + nYDofs) + nLambdaDofs ] = lambda0PdeDof;

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();
    if (iel != 0 || !volumeConstraint) aResLambda[0] = solLambda[0];

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solxType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phix_uv[dim]; // local test function first order partial derivatives

      const double *phiY;  // local test function
      const double *phiY_uv[dim]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // *** get gauss point weight, test function and test function partial derivatives ***
      phix = msh->_finiteElement[ielGeom][solxType]->GetPhi (ig);
      phix_uv[0] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDXi (ig); //derivative in u
      phix_uv[1] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDEta (ig); //derivative in v

      phiY = msh->_finiteElement[ielGeom][solYType]->GetPhi (ig);
      phiY_uv[0] = msh->_finiteElement[ielGeom][solYType]->GetDPhiDXi (ig);
      phiY_uv[1] = msh->_finiteElement[ielGeom][solYType]->GetDPhiDEta (ig);

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
      adept::adouble solY_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solYg[3] = {0., 0., 0.};
      adept::adouble solxg[3] = {0., 0., 0.};

      double solxOld_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solYOldg[3] = {0., 0., 0.};
      double solxOldg[3] = {0., 0., 0.};

      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          solxg[K] += phix[i] * solx[K][i];
          solxOldg[K] += phix[i] * solxOld[K][i];
        }
        for (unsigned i = 0; i < nYDofs; i++) {
          solYg[K] += phiY[i] * solY[K][i];
          solYOldg[K] += phiY[i] * solYOld[K][i];
        }
        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * solx[K][i];
            solxOld_uv[K][j] += phix_uv[j][i] * solxOld[K][i];
          }
        }
        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nYDofs; i++) {
            solY_uv[K][j] += phiY_uv[j][i] * solY[K][i];
          }
        }
      }

      double solYOldnorm2 = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        solYOldnorm2 += solYOldg[K] * solYOldg[K];
      }

      double A = pow (solYOldnorm2 , 1. / ( (2.* P - 2)));
//       std::cout << solYnorm2 << " ";

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
      adept::adouble solY_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};


      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            solx_Xtan[I][J] += solx_uv[I][k] * Jir[k][J];
            solxOld_Xtan[I][J] += solxOld_uv[I][k] * Jir[k][J];
            solY_Xtan[I][J] += solY_uv[I][k] * Jir[k][J];
          }
        }
      }


      std::vector < double > phiY_Xtan[DIM];
      std::vector < double > phix_Xtan[DIM];

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
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;
          for (unsigned J = 0; J < DIM; J++) {
            term1 +=  solx_Xtan[K][J] * phix_Xtan[J][i];
          }
          aResY[K][i] += (solYg[K] * pow (A, 2. - P)  * phix[i] + term1) * Area;
        }
        for (unsigned i = 0; i < nYDofs; i++) {
          adept::adouble term0 = 0.;
          adept::adouble term1 = 0.;
          adept::adouble term2 = 0.;
          adept::adouble term3 = 0.;

          adept::adouble term5 = 0.;
          adept::adouble term6 = 0.;


          for (unsigned J = 0; J < DIM; J++) {

            term0 -=  P * solY_Xtan[K][J] * phiY_Xtan[J][i];

            term1 -= (P - 1.) * solx_Xtan[K][J] * phiY_Xtan[J][i];

            term2 -=  P * solY_Xtan[J][J];

            adept::adouble term4 = 0.;

            for (unsigned L = 0; L < DIM; L++) {
              term4 += solxOld_Xtan[J][L] * solY_Xtan[K][L] + solxOld_Xtan[K][L] * solY_Xtan[J][L];
            }

            term3 += P * phiY_Xtan[J][i] * term4;

          }
          aResx[K][i] += (P * ( volumeConstraint * lambda0 * normal[K] + (solxg[K] - solxOldg[K])  / dt) * phiY[i]
                          + term0
                          + pow (A , P) * term1.value()
                          + term2.value() * phiY_Xtan[K][i]
                          + term3.value()) * Area;
        }
      }
      for (unsigned K = 0; K < DIM; K++) {
        if ( volumeConstraint )  aResLambda0 += ( (solxg[K] - solxOldg[K]) * normal[K]) * Area;

        volume += (solxg[K].value()  * normal[K]) * Area;
      }

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
    for (int i = 0; i < nLambdaDofs; i++) {
      Res[DIM * nxDofs + DIM * nYDofs + i] = -aResLambda[i].value();
    }
    if (!iel0) Res[DIM * (nxDofs + nYDofs) +  nLambdaDofs ] = - aResLambda0.value();


    RES->add_vector_blocked (Res, SYSDOF);

    Jac.resize ( (DIM * (nxDofs + nYDofs) + nLambdaDofs + !iel0) * (DIM * (nxDofs + nYDofs) + nLambdaDofs + !iel0));

    // define the dependent variables
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResx[K][0], nxDofs);
    }
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResY[K][0], nYDofs);
    }
    s.dependent (&aResLambda[0], nLambdaDofs);

    if (!iel0) s.dependent (&aResLambda0, 1);

    // define the dependent variables
    for (int K = 0; K < DIM; K++) {
      s.independent (&solx[K][0], nxDofs);
    }
    for (int K = 0; K < DIM; K++) {
      s.independent (&solY[K][0], nYDofs);
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

  //VecView((static_cast<PetscVector*>(RES))->vec(),  PETSC_VIEWER_STDOUT_SELF );

  //abort();
// MatView((static_cast<PetscMatrix*>(KK))->mat(), PETSC_VIEWER_STDOUT_SELF );

//     PetscViewer    viewer;
//     PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,900,900,&viewer);
//     PetscObjectSetName((PetscObject)viewer,"FSI matrix");
//     PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
//     MatView((static_cast<PetscMatrix*>(KK))->mat(),viewer);
//     double a;
//     std::cin>>a;


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

  std::vector < adept::adouble > solY[DIM]; // local Y solution

  std::vector< int > SYSDOF; // local to global pdeSys dofs

  vector< double > Res; // local redidual vector
  std::vector< adept::adouble > aResY[3]; // local redidual vector

  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)

  KK->zero();  // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solxType);   // number of solution element dofs
    unsigned nYDofs  = msh->GetElementDofNumber (iel, solYType);   // number of solution element dofs

    for (unsigned K = 0; K < DIM; K++) {
      solx[K].resize (nxDofs);
      solY[K].resize (nYDofs);
    }

    // resize local arrays
    SYSDOF.resize (DIM * nYDofs);

    Res.resize (DIM * nYDofs);     //resize

    for (unsigned K = 0; K < DIM; K++) {
      aResY[K].assign (nYDofs, 0.);  //resize and zet to zero
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

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solxType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phix_uv[dim]; // local test function first order partial derivatives

      const double *phiY;  // local test function
      const double *phiY_uv[dim]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // *** get gauss point weight, test function and test function partial derivatives ***
      phix = msh->_finiteElement[ielGeom][solxType]->GetPhi (ig);
      phix_uv[0] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDXi (ig); //derivative in u
      phix_uv[1] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDEta (ig); //derivative in v

      phiY = msh->_finiteElement[ielGeom][solYType]->GetPhi (ig);
      phiY_uv[0] = msh->_finiteElement[ielGeom][solYType]->GetDPhiDXi (ig);
      phiY_uv[1] = msh->_finiteElement[ielGeom][solYType]->GetDPhiDEta (ig);

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
      adept::adouble solY_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solxg[3] = {0., 0., 0.};
      adept::adouble solYg[3] = {0., 0., 0.};


      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          solxg[K] += phix[i] * solx[K][i];
        }
        for (unsigned i = 0; i < nYDofs; i++) {
          solYg[K] += phiY[i] * solY[K][i];
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

      adept::adouble solYnorm2 = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        solYnorm2 += solYg[K] * solYg[K];
      }

      adept::adouble A = pow (solYnorm2 , 1 / ( (2.* P - 2)));


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
      normal[1] = (solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);;
      normal[2] = (solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);;

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


      std::vector < double > phiY_Xtan[DIM];
      std::vector < double > phix_Xtan[DIM];

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
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;
          for (unsigned J = 0; J < DIM; J++) {
            term1 +=  solx_Xtan[K][J] * phiY_Xtan[J][i];
          }
          if (fabs (term1 - phix_Xtan[K][i]) > 1.0e-10) {
            std::cout << " error " << term1 << " " << phiY_Xtan[K][i];
            //abort();
          }
//           std::cout << A <<" ";
          aResY[K][i] += (solYg[K] * pow (A.value(), 2. - P)  * phiY[i] + term1) * Area;
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

    RES->add_vector_blocked (Res, SYSDOF);

    Jac.resize (DIM * nYDofs * DIM * nYDofs);

    // define the dependent variables

    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResY[K][0], nYDofs);
    }

    // define the dependent variables

    for (int K = 0; K < DIM; K++) {
      s.independent (&solY[K][0], nYDofs);
    }

    // get the jacobian matrix (ordered by row)
    s.jacobian (&Jac[0], true);

    KK->add_matrix_blocked (Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();

  //VecView((static_cast<PetscVector*>(RES))->vec(),  PETSC_VIEWER_STDOUT_SELF );

  //abort();
  // MatView((static_cast<PetscMatrix*>(KK))->mat(), PETSC_VIEWER_STDOUT_SELF );

  //     PetscViewer    viewer;
  //     PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,900,900,&viewer);
  //     PetscObjectSetName((PetscObject)viewer,"FSI matrix");
  //     PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
  //     MatView((static_cast<PetscMatrix*>(KK))->mat(),viewer);
  //     double a;
  //     std::cin>>a;


  // ***************** END ASSEMBLY *******************
}


