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

using namespace femus;

void AssembleConformalMinimization (MultiLevelProblem&);  //stable and not bad
void AssembleShearMinimization (MultiLevelProblem&);  //vastly inferior
void UpdateScale (MultiLevelProblem& ml_prob, const double &scalingFactor) ;

// IBVs.  No boundary, and IVs set to sphere (just need something).
bool SetBoundaryCondition (const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {

  bool dirichlet = true;
  value = 0.;

//   if (!strcmp (solName, "Dx1")) {
//     if (1 == faceName ) {
//       dirichlet = false;
//     }
//     if (4 == faceName || 3 == faceName ) {
//       value = (0.5 + 0.4 * cos ( (x[1] - 0.5) * acos (-1.))) * (0.5 - x[0]);
//     }
//   }
//   else if (!strcmp (solName, "Dx2")) {
//     if (2 == faceName) {
//       dirichlet = false;
//     }
//   }


  if (!strcmp (solName, "Dx1")) {
    if (1 == faceName || 3 == faceName) {
      dirichlet = false;
    }
    if (4 == faceName) {
      value = 0.5 * sin (x[1] / 0.5 * acos (-1.));
    }
  }
  else if (!strcmp (solName, "Dx2")) {
    if (2 == faceName) {
      dirichlet = false;
    }
  }

  return dirichlet;
}

double InitalValueScale (const std::vector < double >& x) {
  return 1;
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

  //mlMsh.ReadCoarseMesh ("../input/squareTri.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh ("../input/square.neu", "seventh", scalingFactor);

  // Set number of mesh levels.
  unsigned numberOfUniformLevels = 5;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh (numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // Erase all the coarse mesh levels.
  mlMsh.EraseCoarseLevels (numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  const unsigned  dim = mlMsh.GetDimension();

  // Define the multilevel solution and attach the mlMsh object to it.
  MultiLevelSolution mlSol (&mlMsh);

  // Add variables X,Y,W to mlSol.
//   mlSol.AddSolution ("Dx1", LAGRANGE, FIRST, 0);
//   mlSol.AddSolution ("Dx2", LAGRANGE, FIRST, 0);
//   mlSol.AddSolution ("Dx3", LAGRANGE, FIRST, 0);

  mlSol.AddSolution ("Dx1", LAGRANGE, SECOND, 0);
  mlSol.AddSolution ("Dx2", LAGRANGE, SECOND, 0);
  mlSol.AddSolution ("Dx3", LAGRANGE, SECOND, 0);

  mlSol.AddSolution ("eScale", DISCONTINUOUS_POLYNOMIAL, ZERO, 0);
  mlSol.AddSolution ("eCounter", DISCONTINUOUS_POLYNOMIAL, ZERO, 0);
//   mlSol.AddSolution ("vScale", LAGRANGE, SECOND, 0);
//   mlSol.AddSolution ("vCounter", LAGRANGE, SECOND, 0);

  // Initialize the variables and attach boundary conditions.
  mlSol.Initialize ("All");
  mlSol.Initialize ("eScale", InitalValueScale);

  mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
  mlSol.GenerateBdc ("All");

  MultiLevelProblem mlProb (&mlSol);

  // Add system Conformal or Shear Minimization in mlProb.
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("conformal"); //for conformal

  // Add solutions newDX, Lambda1 to system.
  system.AddSolutionToSystemPDE ("Dx1");
  system.AddSolutionToSystemPDE ("Dx2");
  if (dim == 3) system.AddSolutionToSystemPDE ("Dx3");

  // Parameters for convergence and # of iterations.
  system.SetMaxNumberOfNonLinearIterations (100);
  system.SetNonLinearConvergenceTolerance (1.e-10);

  // Attach the assembling function to system and initialize.
  //system.SetAssembleFunction (AssembleShearMinimization);
  system.SetAssembleFunction (AssembleConformalMinimization);
  system.init();

  mlSol.SetWriter (VTK);
  std::vector<std::string> mov_vars;
  mov_vars.push_back ("Dx1");
  mov_vars.push_back ("Dx2");
  mlSol.GetWriter()->SetMovingMesh (mov_vars);

  // and this?
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back ("All");
  mlSol.GetWriter()->SetDebugOutput (true);
  //mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "linear", variablesToBePrinted, 0);
  mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

  system.MGsolve();

  //mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "linear", variablesToBePrinted, 1);
  mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 1);

  return 0;
}



// Building the Conformal Minimization system.
void AssembleConformalMinimization (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  UpdateScale (ml_prob, 1.1);

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  Extract pointers to the several objects that we are going to use.
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem> ("conformal");   // pointer to the linear implicit system named "Poisson"

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
  const unsigned  dim = msh->GetDimension();
  const unsigned  DIM = 3;

  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // Setting the reference elements to be equilateral triangles.
  std::vector < std::vector < double > > xT (2);
  xT[0].resize (7);
  xT[0][0] = -0.5;
  xT[0][1] = 0.5;
  xT[0][2] = 0.;
  xT[0][3] = 0.;
  xT[0][4] = 0.25;
  xT[0][5] = -0.25;
  xT[0][6] = 0.;

  xT[1].resize (7);
  xT[1][0] = 0.;
  xT[1][1] = 0.;
  xT[1][2] = sqrt (3.) / 2.;
  xT[1][3] = 0.;
  xT[1][4] = sqrt (3.) / 4.;
  xT[1][5] = sqrt (3.) / 4.;
  xT[1][6] = sqrt (3.) / 6.;;

  std::vector<double> phi_uv0;
  std::vector<double> phi_uv1;

  std::vector< double > stdVectorPhi;
  std::vector< double > stdVectorPhi_uv;

  // Extract positions of Dx in ml_sol object.
  std::vector < unsigned >  solDxIndex (DIM);
  solDxIndex[0] = mlSol->GetIndex ("Dx1");
  solDxIndex[1] = mlSol->GetIndex ("Dx2");
  solDxIndex[2] = mlSol->GetIndex ("Dx3");

  // Extract finite element type for the solution.
  unsigned solType;
  solType = mlSol->GetSolutionType (solDxIndex[0]);

  unsigned solScaleIndex = mlSol->GetIndex ("eScale");

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get the positions of Y in the pdeSys object.
  std::vector < unsigned > solDxPdeIndex (dim);
  solDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("Dx1");
  solDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("Dx2");
  if (dim == 3) solDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("Dx3");

  // Local solution vectors for Nx and NDx.
  std::vector < std::vector < adept::adouble > > solDx (DIM);
  std::vector < std::vector < adept::adouble > > solx (DIM);
  std::vector < std::vector < double > > solxHat (DIM);

  // Local-to-global pdeSys dofs.
  std::vector < int > SYSDOF;

  // Local residual vectors.
  vector< double > Res;
  std::vector < std::vector< adept::adouble > > aResDx (dim);

  // Local Jacobian matrix (ordered by column).
  vector < double > Jac;

  KK->zero();  // Zero all the entries of the Global Matrix
  RES->zero(); // Zero all the entries of the Global Residual

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    double scaleValue = (*sol->_Sol[solScaleIndex]) (iel);

    // Numer of solution element dofs.
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solType);
    //unsigned nLDofs  = msh->GetElementDofNumber (iel, solLType);

    // Resize local arrays.
    for (unsigned K = 0; K < DIM; K++) {
      solDx[K].resize (nxDofs);
      solx[K].resize (nxDofs);
      solxHat[K].resize (nxDofs);
    }

    // Resize local arrays
    SYSDOF.resize (dim * nxDofs);
    Res.resize (dim * nxDofs);

    for (unsigned k = 0; k < dim; k++) {
      aResDx[k].assign (nxDofs, 0.);
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nxDofs; i++) {
      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof (i, iel, solType);
      for (unsigned K = 0; K < DIM; K++) {
        solDx[K][i] = (*sol->_Sol[solDxIndex[K]]) (iDDof);
        // Global-to-global mapping between NDx solution node and pdeSys dof.
        if (K < dim) {
          SYSDOF[ K * nxDofs + i] = pdeSys->GetSystemDof (solDxIndex[K], solDxPdeIndex[K], i, iel);
        }
      }
    }

    // start a new recording of all the operations involving adept variables.
    s.new_recording();
    for (unsigned i = 0; i < nxDofs; i++) {
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);
      for (unsigned K = 0; K < DIM; K++) {
        solxHat[K][i] = (*msh->_topology->_Sol[K]) (iXDof);
        solx[K][i] = (*msh->_topology->_Sol[K]) (iXDof) + solDx[K][i];
      }
    }

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      //const double *phiL;  // local test function
      const double *phix_uv[dim]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // Get Gauss point weight, test function, and first order derivatives.
      if (ielGeom == QUAD) {
        phix = msh->_finiteElement[ielGeom][solType]->GetPhi (ig);
        //phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi (ig);

        phix_uv[0] = msh->_finiteElement[ielGeom][solType]->GetDPhiDXi (ig);
        phix_uv[1] = msh->_finiteElement[ielGeom][solType]->GetDPhiDEta (ig);

        weight = msh->_finiteElement[ielGeom][solType]->GetGaussWeight (ig);
      }

      // Special adjustments for triangles.
      else {
        msh->_finiteElement[ielGeom][solType]->Jacobian (xT, ig, weight, stdVectorPhi, stdVectorPhi_uv);
        phix = &stdVectorPhi[0];
        phi_uv0.resize (nxDofs);
        phi_uv1.resize (nxDofs);
        for (unsigned i = 0; i < nxDofs; i++) {
          phi_uv0[i] = stdVectorPhi_uv[i * dim];
          phi_uv1[i] = stdVectorPhi_uv[i * dim + 1];
        }
        phix_uv[0] = &phi_uv0[0];
        phix_uv[1] = &phi_uv1[0];
      }

      // Initialize and compute values of x, Dx, NDx, x_uv at the Gauss points.
      double solxHatg[DIM] = {0., 0., 0.};
      adept::adouble solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          solxHatg[K] += phix[i] * solxHat[K][i];
        }
        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]   += phix_uv[j][i] * solx[K][i];
          }
        }
      }

      // Compute the metric, metric determinant, and area element.
      std::vector < std::vector < adept::adouble > > g (dim);
      for (unsigned i = 0; i < dim; i++) g[i].assign (dim, 0.);

      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }

      adept::adouble detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      adept::adouble Area = weight * sqrt (detg);
      adept::adouble Area2 = weight;// Trick to give equal weight to each element.

      // Compute components of the unit normal N.
      adept::adouble normal[DIM];
//       normal[0] = (solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) / sqrt (detg);
//       normal[1] = (solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);
//       normal[2] = (solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);

      normal[0] = 0.;
      normal[1] = 0.;
      normal[2] = 1.;

      // Discretize the equation \delta CD = 0 on the basis d/du, d/dv.
      adept::adouble V[DIM];
      V[0] = solx_uv[0][1] - normal[1] * solx_uv[2][0] + normal[2] * solx_uv[1][0];
      V[1] = solx_uv[1][1] - normal[2] * solx_uv[0][0] + normal[0] * solx_uv[2][0];
      V[2] = solx_uv[2][1] - normal[0] * solx_uv[1][0] + normal[1] * solx_uv[0][0];

      adept::adouble W[DIM];
      W[0] = (solx_uv[0][0] + normal[1] * solx_uv[2][1] - normal[2] * solx_uv[1][1]);
      W[1] = (solx_uv[1][0] + normal[2] * solx_uv[0][1] - normal[0] * solx_uv[2][1]);
      W[2] = (solx_uv[2][0] + normal[0] * solx_uv[1][1] - normal[1] * solx_uv[0][1]);

      adept::adouble M[DIM][dim];
      M[0][0] = W[0] - normal[2] * V[1] + normal[1] * V[2];
      M[1][0] = W[1] - normal[0] * V[2] + normal[2] * V[0];
      M[2][0] = W[2] - normal[1] * V[0] + normal[0] * V[1];

      M[0][1] = V[0] + normal[2] * W[1] - normal[1] * W[2];
      M[1][1] = V[1] + normal[0] * W[2] - normal[2] * W[0];
      M[2][1] = V[2] + normal[1] * W[0] - normal[0] * W[1];

      // Implement the Conformal Minimization equations.
      for (unsigned k = 0; k < dim; k++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;
          for (unsigned j = 0; j < dim; j++) {
            term1 +=  M[k][j] * phix_uv[j][i];
          }
          // Conformal energy equation (with trick).
          aResDx[k][i] += term1 * Area2 * scaleValue;
        }
      }
    } // end GAUSS POINT LOOP

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    //copy the value of the adept::adoube aRes in double Res and store

    for (int k = 0; k < dim; k++) {
      for (int i = 0; i < nxDofs; i++) {
        Res[ k * nxDofs + i] = -aResDx[k][i].value();
      }
    }

    RES->add_vector_blocked (Res, SYSDOF);

    // Resize Jacobian.
    Jac.resize ( (dim * nxDofs) * (dim * nxDofs));

    // Define the dependent variables.
    for (int k = 0; k < dim; k++) {
      s.dependent (&aResDx[k][0], nxDofs);
    }

    // Define the independent variables.
    for (int k = 0; k < dim; k++) {
      s.independent (&solDx[k][0], nxDofs);
    }

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

  UpdateScale (ml_prob, 1.3);


  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem> ("conformal");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh *msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem *el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution *mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  SparseMatrix *KK = pdeSys->_KK;  // pointer to the global stiffness matrix object in pdeSys (level)
  NumericVector *RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension();

  std::vector <double> phi;  // local test function for velocity
  std::vector <adept::adouble> phi_x; // local test function first order partial derivatives
  adept::adouble weight; // gauss point weight

  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  std::vector < unsigned > solDxIndex (dim);
  solDxIndex[0] = mlSol->GetIndex ("Dx1"); // get the position of "DX" in the ml_sol object
  solDxIndex[1] = mlSol->GetIndex ("Dx2"); // get the position of "DY" in the ml_sol object
  if (dim == 3) solDxIndex[2] = mlSol->GetIndex ("Dx3"); // get the position of "DY" in the ml_sol object

  unsigned solType;
  solType = mlSol->GetSolutionType (solDxIndex[0]);  // get the finite element type for "U"

  unsigned solScaleIndex = mlSol->GetIndex ("eScale");

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector < unsigned > solDxPdeIndex (dim);
  solDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("Dx1");   // get the position of "Dx1" in the pdeSys object
  solDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("Dx2");   // get the position of "Dx2" in the pdeSys object
  if (dim == 3) solDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("Dx3");  // get the position of "Dx3" in the pdeSys object

  std::vector < std::vector < adept::adouble > > solDx (dim); // local Y solution
  std::vector < std::vector < adept::adouble > > x (dim);

  std::vector< int > SYSDOF; // local to global pdeSys dofs

  vector< double > Res; // local redidual vector
  std::vector< adept::adouble > aResDx[dim]; // local redidual vector

  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)

  KK->zero();  // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    double scaleValue = (*sol->_Sol[solScaleIndex]) (iel);

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solType);   // number of solution element dofs

    for (unsigned k = 0; k < dim; k++) {
      solDx[k].resize (nxDofs);
      x[k].resize (nxDofs);
    }

    // resize local arrays
    SYSDOF.resize (dim * nxDofs);
    Res.resize (dim * nxDofs);       //resize

    for (unsigned k = 0; k < dim; k++) {
      aResDx[k].assign (nxDofs, 0.);  //resize and zet to zero
    }


    // local storage of global mapping and solution
    for (unsigned i = 0; i < nxDofs; i++) {
      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof (i, iel, solType);
      for (unsigned k = 0; k < dim; k++) {
        solDx[k][i] = (*sol->_Sol[solDxIndex[k]]) (iDDof);
        // Global-to-global mapping between NDx solution node and pdeSys dof.
        SYSDOF[ k * nxDofs + i] = pdeSys->GetSystemDof (solDxIndex[k], solDxPdeIndex[k], i, iel);
      }
    }

    // start a new recording of all the operations involving adept variables.
    s.new_recording();
    for (unsigned i = 0; i < nxDofs; i++) {
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);
      for (unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k]) (iXDof);// + solDx[k][i];
      }
    }

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielGeom][solType]->Jacobian (x, ig, weight, phi, phi_x);
      weight = msh->_finiteElement[ielGeom][solType]->GetGaussWeight (ig);

      std::vector < std::vector < adept::adouble > > gradSolDx (dim);

      for (unsigned  k = 0; k < dim; k++) {
        gradSolDx[k].assign (dim, 0.);
      }

      for (unsigned i = 0; i < nxDofs; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
            gradSolDx[k][j] += solDx[k][i] * phi_x[i * dim + j];
          }
        }
      }

      for (unsigned i = 0; i < nxDofs; i++) {
        for (unsigned  k = 0; k < dim; k++) {
          adept::adouble term = 0.;
          for (unsigned j = 0; j < dim; j++) {
            term  +=  phi_x[i * dim + j] * (gradSolDx[k][j] + gradSolDx[j][k]);
          }
          aResDx[k][i] += term * scaleValue * weight;
        }
      }
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store


    for (int k = 0; k < dim; k++) {
      for (int i = 0; i < nxDofs; i++) {
        Res[ k * nxDofs + i] = -aResDx[k][i].value();
      }
    }

    RES->add_vector_blocked (Res, SYSDOF);

    Jac.resize ( (dim * nxDofs) * (dim * nxDofs));

    // define the dependent variables

    for (int k = 0; k < dim; k++) {
      s.dependent (&aResDx[k][0], nxDofs);
    }

    // define the dependent variables

    for (int k = 0; k < dim; k++) {
      s.independent (&solDx[k][0], nxDofs);
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

void UpdateScale (MultiLevelProblem& ml_prob, const double &elScalingFactor) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem> ("conformal");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh *msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem *el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution *mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  const unsigned  dim = msh->GetDimension();

  std::vector <double> phi;  // local test function for velocity
  std::vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned elScaleIndex = mlSol->GetIndex ("eScale");
  //unsigned vtScaleIndex = mlSol->GetIndex ("vScale");
  //unsigned vtScaleType = mlSol->GetSolutionType (vtScaleIndex);
  unsigned counterIndex = mlSol->GetIndex ("eCounter");


  //solution variable
  std::vector < unsigned > solDxIndex (dim);
  solDxIndex[0] = mlSol->GetIndex ("Dx1"); // get the position of "DX" in the ml_sol object
  solDxIndex[1] = mlSol->GetIndex ("Dx2"); // get the position of "DY" in the ml_sol object
  if (dim == 3) solDxIndex[2] = mlSol->GetIndex ("Dx3"); // get the position of "DY" in the ml_sol object
  unsigned solDxType;
  solDxType = mlSol->GetSolutionType (solDxIndex[0]);  // get the finite element type for "U"

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector < std::vector < double > > x (dim);

  //sol->_Sol[vtScaleIndex]->zero();
  //sol->_Sol[counterIndex]->zero();


  // element loop: each process loops only on the elements that owns
  for (unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    double scaleValue = (*sol->_Sol[elScaleIndex]) (iel);
    double counter = (*sol->_Sol[counterIndex]) (iel);

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solDxType);   // number of solution element dofs
    //unsigned vtDofs  = msh->GetElementDofNumber (iel, vtScaleType);   // number of solution element dofs

    for (unsigned k = 0; k < dim; k++) {
      x[k].resize (nxDofs);
    }

    for (unsigned i = 0; i < nxDofs; i++) {
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);
      unsigned iDxDof  = msh->GetSolutionDof (i, iel, solDxType);
      for (unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k]) (iXDof) + (*sol->_Sol[solDxIndex[k]]) (iDxDof);
      }
    }

    //double vtScalingFactor = 1.;
    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solDxType]->GetGaussPointNumber(); ig++) {
      msh->_finiteElement[ielGeom][solDxType]->Jacobian (x, ig, weight, phi, phi_x);

      weight /= msh->_finiteElement[ielGeom][solDxType]->GetGaussWeight (ig);

      if (weight < 0.) {
        std::cout << "warning element inversion at iel = " << iel << " " << scaleValue << std::endl;
        sol->_Sol[elScaleIndex]->set (iel, (elScalingFactor + 0.05 * counter) * scaleValue);
        //vtScalingFactor = elScalingFactor;
        sol->_Sol[counterIndex]->add (iel, 1.);
        break;
      }
    } // end Gauss point loop
//     for (unsigned i = 0; i < vtDofs; i++) {
//       unsigned iDof  = msh->GetSolutionDof (i, iel, vtScaleType);
//       sol->_Sol[vtScaleIndex]->add (iDof, vtScalingFactor * scaleValue);
//       sol->_Sol[counterIndex]->add (iDof, 1.);
//     }
  }
  sol->_Sol[elScaleIndex]->close();
  //sol->_Sol[vtScaleIndex]->close();
  sol->_Sol[counterIndex]->close();

//   for (unsigned i = msh->_dofOffset[vtScaleType][iproc]; i < msh->_dofOffset[vtScaleType][iproc + 1]; i++) {
//     double value = (*sol->_Sol[vtScaleIndex]) (i);
//     sol->_Sol[vtScaleIndex]->set (i, value / (*sol->_Sol[counterIndex]) (i));
//   }
//   sol->_Sol[vtScaleIndex]->close();
//
//   for (unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
//      double value = .99 * (*sol->_Sol[elScaleIndex]) (iel);
//      unsigned vtDofs  = msh->GetElementDofNumber (iel, vtScaleType);
//      for (unsigned i = 0; i < vtDofs; i++) {
//       unsigned iDof  = msh->GetSolutionDof (i, iel, vtScaleType);
//       value += 0.01 * (*sol->_Sol[vtScaleIndex]) (iDof) / vtDofs;
//     }
//     sol->_Sol[elScaleIndex]->set (iel, value);
//   }
//   sol->_Sol[elScaleIndex]->close();
//

  // ***************** END ASSEMBLY *******************
}

