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

void AssembleShearMinimization (MultiLevelProblem& ml_prob);
void AssembleConformalMinimization (MultiLevelProblem&);  //stable and not bad
void UpdateMesh (MultiLevelSolution& mlSol);

// IBVs.  No boundary, and IVs set to sphere (just need something).
bool SetBoundaryCondition (const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {

  bool dirichlet = true;
  value = 0.;

//   if (!strcmp (solName, "Dx1")) {
//     if (1 == faceName ) {
//       dirichlet = false;
//     }
//     if (4 == faceName || 3 == faceName ) {
//       value = (0.5 + 0.499 * cos ( (x[1] - 0.5) * acos (-1.))) * (0.5 - x[0]);
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


  if (!strcmp (solName, "U1")) {
    if (1 == faceName || 3 == faceName) {
      dirichlet = false;
    }
    if (4 == faceName) {
      value = 0.5 * sin (x[1] / 0.5 * acos (-1.));
    }
  }
  else if (!strcmp (solName, "U2")) {
    if (2 == faceName) {
      dirichlet = false;
    }
  }

//     if (!strcmp (solName, "U1")) {
//     if (1 == faceName ) {
//       dirichlet = false;
//     }
//     if (4 == faceName || 3 == faceName ) {
//       value = (0.5 + 0.499 * cos ( (x[1] - 0.5) * acos (-1.))) * (0.5 - x[0]);
//     }
//   }
//   else if (!strcmp (solName, "U2")) {
//     if (2 == faceName) {
//       dirichlet = false;
//     }
//   }


//   if (!strcmp (solName, "Dx1")) {
//     if (1 == faceName || 3 == faceName) {
//       dirichlet = false;
//     }
//     if (4 == faceName) {
//       value = 0;
//     }
//   }
//   else if (!strcmp (solName, "Dx2")) {
//     if (2 == faceName) {
//       dirichlet = false;
//     }
//   }

  return dirichlet;
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
  //mlMsh.ReadCoarseMesh ("../input/squareReg.neu", "seventh", scalingFactor);
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
  mlSol.AddSolution ("U1", LAGRANGE, SECOND, 0);
  mlSol.AddSolution ("U2", LAGRANGE, SECOND, 0);
  mlSol.AddSolution ("U3", LAGRANGE, SECOND, 0);

  mlSol.AddSolution ("Dx1", LAGRANGE, SECOND, 0);
  mlSol.AddSolution ("Dx2", LAGRANGE, SECOND, 0);
  mlSol.AddSolution ("Dx3", LAGRANGE, SECOND, 0);

  // Initialize the variables and attach boundary conditions.
  mlSol.Initialize ("All");

  mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
  mlSol.GenerateBdc ("All");

  MultiLevelProblem mlProb (&mlSol);

  LinearImplicitSystem& system0 = mlProb.add_system < LinearImplicitSystem > ("shear"); //for
  system0.AddSolutionToSystemPDE ("U1");
  system0.AddSolutionToSystemPDE ("U2");
  if (dim == 3) system0.AddSolutionToSystemPDE ("U3");

  system0.SetAssembleFunction (AssembleShearMinimization);
  system0.init();


  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back ("All");

  mlSol.SetWriter (VTK);
  mlSol.GetWriter()->SetDebugOutput (true);
  mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

  system0.MGsolve();

  //UpdateMesh (mlSol);
  std::vector<std::string> mov_vars1;
  mov_vars1.push_back ("U1");
  mov_vars1.push_back ("U2");
  mlSol.GetWriter()->SetMovingMesh (mov_vars1);
  

  mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 1);

  // Add system Conformal or Shear Minimization in mlProb.
  NonLinearImplicitSystem& system1 = mlProb.add_system < NonLinearImplicitSystem > ("conformal"); //for conformal

  // Add solutions newDX, Lambda1 to system.
  system1.AddSolutionToSystemPDE ("Dx1");
  system1.AddSolutionToSystemPDE ("Dx2");
  if (dim == 3) system1.AddSolutionToSystemPDE ("Dx3");

  // Parameters for convergence and # of iterations.
  system1.SetMaxNumberOfNonLinearIterations (1);
  system1.SetNonLinearConvergenceTolerance (1.e-10);

  system1.SetAssembleFunction (AssembleConformalMinimization);
  system1.init();

  for(unsigned j=0;j<10;j++)  {
  
    system1.MGsolve();

    mlSol.SetWriter (VTK);
    std::vector<std::string> mov_vars2;
    mov_vars2.push_back ("Dx1");
    mov_vars2.push_back ("Dx2");
    mlSol.GetWriter()->SetMovingMesh (mov_vars2);

    mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 2+j);
  }

  return 0;
}

unsigned counter = 0;

// Building the Conformal Minimization system.
void AssembleConformalMinimization (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled


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

  std::vector <double> phiHat;  // local test function for velocity
  std::vector <double> phiHat_x; // local test function first order partial derivatives
  double weightHat; // gauss point weight

  std::vector <double> phi1;  // local test function for velocity
  std::vector <adept::adouble> phi1_x; // local test function first order partial derivatives
  adept::adouble weight1; // gauss point weight


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
  xT[1][6] = sqrt (3.) / 6.;

  std::vector<double> phi_uv0;
  std::vector<double> phi_uv1;

  std::vector< double > stdVectorPhi;
  std::vector< double > stdVectorPhi_uv;

  // Extract positions of Dx in ml_sol object.
  std::vector < unsigned >  solDxIndex (DIM);
  solDxIndex[0] = mlSol->GetIndex ("Dx1");
  solDxIndex[1] = mlSol->GetIndex ("Dx2");
  solDxIndex[2] = mlSol->GetIndex ("Dx3");

  std::vector < unsigned >  solUIndex (DIM);
  solUIndex[0] = mlSol->GetIndex ("U1");
  solUIndex[1] = mlSol->GetIndex ("U2");
  solUIndex[2] = mlSol->GetIndex ("U3");

  //unsigned solLIndex = mlSol->GetIndex ("Lambda");

  // Extract finite element type for the solution.
  unsigned solType;
  solType = mlSol->GetSolutionType (solDxIndex[0]);

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get the positions of Y in the pdeSys object.
  std::vector < unsigned > solDxPdeIndex (dim);
  solDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("Dx1");
  solDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("Dx2");
  if (dim == 3) solDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("Dx3");

  // Local solution vectors for Nx and NDx.
  std::vector < std::vector < adept::adouble > > solDx (DIM);
  std::vector < std::vector < double > > solU (DIM);

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


    // Numer of solution element dofs.
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solType);

    // Resize local arrays.
    for (unsigned K = 0; K < DIM; K++) {
      solU[K].resize (nxDofs);
      solDx[K].resize (nxDofs);
      solx[K].resize (nxDofs);
      solxHat[K].resize (nxDofs);
    }
    //solL.resize (nLDofs);

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
        solU[K][i] = (*sol->_Sol[solUIndex[K]]) (iDDof);
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
        
        solxHat[K][i] = (*msh->_topology->_Sol[K]) (iXDof) + (counter != 0) * solDx[K][i].value();
        solx[K][i] = (*msh->_topology->_Sol[K]) (iXDof) + solDx[K][i];
      }
    }

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielGeom][solType]->Jacobian (solxHat, ig, weightHat, phiHat, phiHat_x);
      msh->_finiteElement[ielGeom][solType]->Jacobian (solx, ig, weight1, phi1, phi1_x);

      double xg = 0;
      
      std::vector < adept::adouble > solDxg (dim, 0.);
      std::vector < double > solUg (dim, 0.);
      std::vector < std::vector < adept::adouble > > gradSolDx (dim);
      std::vector < std::vector < double > > gradSolU (dim);
      for (unsigned  k = 0; k < dim; k++) {
        gradSolDx[k].assign (dim, 0.);
        gradSolU[k].assign (dim, 0.);
      }

      for (unsigned i = 0; i < nxDofs; i++) {
        xg += solxHat[0][i] * phiHat[i];
        for (unsigned k = 0; k < dim; k++) {
          solUg[k] += solU[k][i] * phiHat[i];
          solDxg[k] += solDx[k][i] * phiHat[i];
          for (unsigned  j = 0; j < dim; j++) {
            //if(counter == 0){  
              gradSolDx[k][j] += solDx[k][i] * phiHat_x[i * dim + j];
            //}
            //else{
              //gradSolDx[k][j] += solDx[k][i] * phi1_x[i * dim + j].value();
            //}
            gradSolU[k][j] += solU[k][i] * phiHat_x[i * dim + j];
          }
        }
      }

      double energy1 = 0;
      double energy2 = 0;
      for (unsigned k = 0; k < dim; k++) {
        energy1 += solUg[k] * solUg[k];
        for (unsigned j = 0; j < dim; j++) {
          energy2 += (gradSolU[k][j] + gradSolU[j][k]) * (gradSolU[k][j] + gradSolU[j][k]);
        }
      }
      double penalty = 0.8 * pow (10 * energy1, 1) + 0.2 * pow (0.1 * energy2, 4);

      const double *phi;  // local test function
      const double *phi_uv[dim]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // Get Gauss point weight, test function, and first order derivatives.
      if (ielGeom == QUAD) {
        phi = msh->_finiteElement[ielGeom][solType]->GetPhi (ig);
        phi_uv[0] = msh->_finiteElement[ielGeom][solType]->GetDPhiDXi (ig);
        phi_uv[1] = msh->_finiteElement[ielGeom][solType]->GetDPhiDEta (ig);

        weight = msh->_finiteElement[ielGeom][solType]->GetGaussWeight (ig);
      }

      // Special adjustments for triangles.
      else {
        msh->_finiteElement[ielGeom][solType]->Jacobian (xT, ig, weight, stdVectorPhi, stdVectorPhi_uv);
        phi = &stdVectorPhi[0];
        phi_uv0.resize (nxDofs);
        phi_uv1.resize (nxDofs);
        for (unsigned i = 0; i < nxDofs; i++) {
          phi_uv0[i] = stdVectorPhi_uv[i * dim];
          phi_uv1[i] = stdVectorPhi_uv[i * dim + 1];
        }
        phi_uv[0] = &phi_uv0[0];
        phi_uv[1] = &phi_uv1[0];
      }

      // Initialize and compute values of x, Dx, NDx, x_uv at the Gauss points.

      adept::adouble solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      for (unsigned K = 0; K < DIM; K++) {
        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]   += phi_uv[j][i] * solx[K][i];
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

      double normal[DIM];
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

      double dir[2] = {1, 0.05};


      adept::adouble F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      adept::adouble B[3][3];
      adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
      adept::adouble Cauchy[3][3];

      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          F[i][j] += gradSolDx[i][j];
        }
      }

      adept::adouble J_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                              - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          B[i][j] = 0.;

          for (int k = 0; k < 3; k++) {
            //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
            B[i][j] += F[i][k] * F[j][k];
          }
        }
      }

      adept::adouble I1_B = B[0][0] + B[1][1] + B[2][2];

      double ni = 0.3;
      double E = 1.;
      double lambda = (ni * E) / ( (1. + ni) * (1. - 2. * ni));
      double mu = E / (2. * (1. + ni));


      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          //  Cauchy[i][j] = mu * (B[i][j] - I1_B * Id2th[i][j] / 3.) / pow(J_hat, 5. / 3.)
          //                 + K * (J_hat - 1.) * Id2th[i][j];  //Generalized Neo-Hookean solid, in Allan-Bower's book, for rubbers with very limited compressibility and K >> mu

          Cauchy[i][j] = lambda * log (J_hat) / J_hat * Id2th[i][j] + mu / J_hat * (B[i][j] - Id2th[i][j]); //alternative formulation

        }
      }




      // Implement the Conformal Minimization equations.
      for (unsigned i = 0; i < nxDofs; i++) {

        adept::adouble CauchyDIR[3] = {0., 0., 0.};

        for (int k = 0.; k < dim; k++) {
          for (int j = 0.; j < dim; j++) {
            CauchyDIR[k] += phi1_x[i * dim + j] * Cauchy[k][j];
          }
        }

        for (int k = 0; k < dim; k++) {
          //aResDx[k][i] +=  CauchyDIR[k] * weight1;
        }

        for (unsigned k = 0; k < dim; k++) {
          adept::adouble term1 = 0.;
          for (unsigned j = 0; j < dim; j++) {
            term1 +=  M[k][j] * phi_uv[j][i];
          }

          adept::adouble term2 = 0.;
          if(k==0){
            if ( counter == 0 ){
              term2 += /*45 * pow((0.5 - xg), 4)*/ 100000 * gradSolDx[k][k] * phiHat_x[k + i*dim] * pow(weightHat, 2./dim);
            }
            else{
              term2 += /*45 * pow((0.5 - xg), 4)*/ 100000 * gradSolDx[k][k] * phi1_x[k + i*dim] * pow(weight1, 2./dim);
            }
            //term2 += /*45 * pow((0.5 - xg), 4) */ 0 * gradSolDx[k][1] * phiHat_x[1 + i*dim] * pow(weightHat, 2./dim);
          }
          aResDx[k][i] += 1 * (term1 + 0. * dir[k] * penalty * solDxg[k] * phi[i]) * Area2 +  term2;//1./weightHat;

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
    Jac.resize ( (dim * nxDofs /*+ nLDofs*/) * (dim * nxDofs /*+ nLDofs*/));

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

  counter++;

} // end AssembleConformalMinimization.

void AssembleShearMinimization (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< LinearImplicitSystem> ("shear");   // pointer to the linear implicit system named "Poisson"

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
  solDxIndex[0] = mlSol->GetIndex ("U1"); // get the position of "DX" in the ml_sol object
  solDxIndex[1] = mlSol->GetIndex ("U2"); // get the position of "DY" in the ml_sol object
  if (dim == 3) solDxIndex[2] = mlSol->GetIndex ("U3"); // get the position of "DY" in the ml_sol object

  unsigned solType;
  solType = mlSol->GetSolutionType (solDxIndex[0]);  // get the finite element type for "U"

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector < unsigned > solDxPdeIndex (dim);
  solDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("U1");   // get the position of "Dx1" in the pdeSys object
  solDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("U2");   // get the position of "Dx2" in the pdeSys object
  if (dim == 3) solDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("U3");  // get the position of "Dx3" in the pdeSys object

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
        x[k][i] = (*msh->_topology->_Sol[k]) (iXDof);
      }
    }

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielGeom][solType]->Jacobian (x, ig, weight, phi, phi_x);

      std::vector < std::vector < adept::adouble > > gradSolDx (dim);

      for (unsigned  k = 0; k < dim; k++) {
        gradSolDx[k].assign (dim, 0.);
      }

      for (unsigned i = 0; i < nxDofs; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
            gradSolDx[k][j] += (x[k][i] + solDx[k][i]) * phi_x[i * dim + j];
          }
        }
      }

      for (unsigned i = 0; i < nxDofs; i++) {
        for (unsigned  k = 0; k < dim; k++) {
          adept::adouble term = 0.;
          term  +=  phi_x[i * dim + k] * (gradSolDx[k][k]);
          aResDx[k][i] += term * weight;
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

void UpdateMesh (MultiLevelSolution& mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1u;

  Solution* sol  = mlSol.GetSolutionLevel (level);
  Mesh* msh = mlSol._mlMesh->GetLevel (level);
  elem* myel =  msh->el;

  const unsigned dim = msh->GetDimension();

  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  std::vector < unsigned > solUIndex (dim);
  solUIndex[0] = mlSol.GetIndex ("U1"); // get the position of "DX" in the ml_sol object
  solUIndex[1] = mlSol.GetIndex ("U2"); // get the position of "DY" in the ml_sol object
  if (dim == 3) solUIndex[2] = mlSol.GetIndex ("U3"); // get the position of "DY" in the ml_sol

  for (unsigned k = 0; k < dim; k++) {
    (*msh->_topology->_Sol[k]).add (*sol->_Sol[solUIndex[k]]);
  }
}

