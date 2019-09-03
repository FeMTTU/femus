/* This example is for quasi-conformal minimization */
/* mu controls the beltrami coefficient */

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
// Comment back in for working code
//const double mu[2] = {0.8, 0.};


void AssembleConformalMinimization (MultiLevelProblem&);  //stable and not bad

// IBVs.  No boundary, and IVs set to sphere (just need something).
bool SetBoundaryCondition (const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {

  bool dirichlet = true;
  value = 0.;

  // if (!strcmp (solName, "Dx1")) {
  //   if (1 == faceName ) {
  //     dirichlet = false;
  //   }
  //   if (4 == faceName || 3 == faceName ) {
  //     value = (0.5 + 0.4 * cos ( (x[1] - 0.5) * acos (-1.))) * (0.5 - x[0]);
  //   }
  // }
  // else if (!strcmp (solName, "Dx2")) {
  //   if (2 == faceName) {
  //     dirichlet = false;
  //   }
  // }


  if (!strcmp (solName, "Dx1")) {
    if (1 == faceName || 3 == faceName) {
      dirichlet = false;
    }
    if (4 == faceName) {
       //value = 0.04 * sin (4*(x[1] / 0.5 * acos (-1.)));
      value = 0.5 * sin ((x[1] / 0.5 * acos (-1.)));
      //dirichlet = false;
    }
  }
  else if (!strcmp (solName, "Dx2")) {
    if (2 == faceName) {
      dirichlet = false;
    }
  }

  // if (!strcmp (solName, "Dx2")) {
  //    if (2 == faceName || 4 == faceName) {
  //      dirichlet = false;
  //    }
  //    if (1 == faceName) {
  //      value = 0.5 * sin ((x[0] / 0.5 * acos (-1.)));
  //    }
  //  }
  //  else if (!strcmp (solName, "Dx1")) {
  //    if (3 == faceName) {
  //      dirichlet = false;
  //    }
  //  }

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

  mlSol.AddSolution ("Dx1", LAGRANGE, SECOND, 0);
  mlSol.AddSolution ("Dx2", LAGRANGE, SECOND, 0);

  // Initialize the variables and attach boundary conditions.
  mlSol.Initialize ("All");

  mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
  mlSol.GenerateBdc ("All");

  MultiLevelProblem mlProb (&mlSol);

  // Add system Conformal or Shear Minimization in mlProb.
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("conformal"); //for conformal

  // Add solutions newDX, Lambda1 to system.
  system.AddSolutionToSystemPDE ("Dx1");
  system.AddSolutionToSystemPDE ("Dx2");

  // Parameters for convergence and # of iterations.
  system.SetMaxNumberOfNonLinearIterations (20);
  system.SetNonLinearConvergenceTolerance (1.e-10);

  // Attach the assembling function to system and initialize.
  //system.SetAssembleFunction (AssembleShearMinimization);
  //system.SetAssembleFunction (AssembleConformalMinimization);
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

//   system.SetAssembleFunction (AssembleConformalMinimization);
//   system.MGsolve();

// mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 2);

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
  const unsigned  DIM = 2;

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

  // Extract finite element type for the solution.
  unsigned solType;
  solType = mlSol->GetSolutionType (solDxIndex[0]);

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get the positions of Y in the pdeSys object.
  std::vector < unsigned > solDxPdeIndex (dim);
  solDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("Dx1");
  solDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("Dx2");

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
      double solxHatg[DIM] = {0., 0.};
      adept::adouble solx_uv[2][2] = {{0., 0.}, {0., 0.}};
      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          solxHatg[K] += phix[i] * solxHat[K][i];
        }
        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j] += phix_uv[j][i] * solx[K][i];
          }
        }
        // for (unsigned i = 0; i < nxDofs; i++) {
        //   solx_z[K] +=
        //   solx_zBar[K] +=
        // }
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

      adept::adouble norm2Xz = (1./4.) * ( pow( (solx_uv[0][0] + solx_uv[1][1]), 2) + pow( (solx_uv[1][0] - solx_uv[0][1]), 2) );

      // Discretize the equation \delta CD = 0 on the basis d/du, d/dv.
      adept::adouble XzBarXz_Bar[DIM];

      // Comment out for working code
      

      XzBarXz_Bar[0] = (1./4.) * ( pow(solx_uv[0][0], 2) + pow(solx_uv[1][0], 2) - pow(solx_uv[0][1], 2) - pow(solx_uv[1][1], 2) );
      XzBarXz_Bar[1] = (1./2.) * ( solx_uv[0][0] * solx_uv[0][1] + solx_uv[1][0] * solx_uv[1][1] );

      // Comment out for working code
      
      adept::adouble mu[2]={0.,0.};
      if( counter == 0) mu[0] = 0.8;
      
      for (unsigned K = 0; K < DIM; K++) {
        if(counter > 0 && norm2Xz.value() > 0.){
          mu[K] += (1. / norm2Xz.value()) * XzBarXz_Bar[K].value();
        }
        //if(counter % 2 == 0) mu[K]*=1.01;
        //else mu[K]/=1.01;
      }
      
      //std::cout << mu[0] <<" "<< mu[1]<<" ";
      
      adept::adouble V[DIM];
      V[0] = (1 - mu[0]) * solx_uv[0][0] - (1 + mu[0]) * solx_uv[1][1] + mu[1] * (solx_uv[1][0] - solx_uv[0][1]);
      V[1] = (1 - mu[0]) * solx_uv[1][0] + (1 + mu[0]) * solx_uv[0][1] - mu[1] * (solx_uv[0][0] + solx_uv[1][1]);


      adept::adouble M[DIM][dim];

      M[0][0] = (1 - mu[0]) * V[0] - mu[1] * V[1];
      M[1][0] = (1 - mu[0]) * V[1] + mu[1] * V[0];
      //M[0][0] = (1 - mu1) * V[0] - mu2 * V[1];
      //M[1][0] = (1 - mu1) * V[1] + mu2 * V[0];

      M[0][1] = (1 + mu[0]) * V[1] - mu[1] * V[0];
      M[1][1]= -(1 + mu[0]) * V[0] - mu[1] * V[1];
      //M[0][1] = (1 + mu1) * V[1] - mu2 * V[0];
      //M[1][1]= -(1 + mu1) * V[0] - mu2 * V[1];


      // Implement the Conformal Minimization equations.
      for (unsigned k = 0; k < dim; k++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;
          for (unsigned j = 0; j < dim; j++) {
            term1 += 2 * M[k][j] * phix_uv[j][i];
          }
          // Conformal energy equation (with trick).
          aResDx[k][i] += term1 * Area2;
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
  
  counter++;

} // end AssembleConformalMinimization.

