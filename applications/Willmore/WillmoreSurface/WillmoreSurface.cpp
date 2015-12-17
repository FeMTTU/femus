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


using namespace femus;

double a = sqrt(2);

// Torus

bool SetBoundaryConditionTorus(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet

  double u = x[0];
  double v = x[1];

  if (!strcmp("X", SolName)) {
    value = (a + cos(u)) * cos(v);
  }
  else if (!strcmp("Y", SolName)) {
    value = (a + cos(u)) * sin(v);
  }
  else if (!strcmp("Z", SolName)) {
    value = sin(u);
  }
  else if (!strcmp("H", SolName)) {
    value = 0.5 * (1. + cos(u) / (a + cos(u)));
  }

  return dirichlet;
}

double InitalValueXTorus(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return (a + cos(u)) * cos(v);

}

double InitalValueYTorus(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return (a + cos(u)) * sin(v);

}

double InitalValueZTorus(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return sin(u);

}

double InitalValueHTorus(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return 0.5 * (1. + cos(u) / (a + cos(u)));

}


bool SetBoundaryConditionSphere(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet

  double u = x[0];
  double v = x[1];

  if (!strcmp("X", SolName)) {
    value = a * sin(v) * cos(u);
  }
  else if (!strcmp("Y", SolName)) {
    value = a * sin(v) * sin(u);
  }
  else if (!strcmp("Z", SolName)) {
    value = a * cos(v);
  }
  else if (!strcmp("H", SolName)) {
    value = 1. / a;
  }

  return dirichlet;
}

double InitalValueXSphere(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return a * sin(v) * cos(u);

}

double InitalValueYSphere(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return  a * sin(v) * sin(u);

}

double InitalValueZSphere(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return a * cos(v);

}

double InitalValueHSphere(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return 1. / a;

}




double SetVariableTimeStep(const double time) {
  double dt = 1.;
  return dt;
}


void AssembleWillmoreFlow_AD(MultiLevelProblem& ml_prob);

std::pair < double, double > GetErrorNorm(MultiLevelSolution* mlSol);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);


  // define multilevel mesh


  unsigned maxNumberOfMeshes;
  maxNumberOfMeshes = 1;

  vector < vector < double > > l2Norm;
  l2Norm.resize(maxNumberOfMeshes);

  vector < vector < double > > semiNorm;
  semiNorm.resize(maxNumberOfMeshes);

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

    std::ostringstream filename;

    filename << "./input/sphere.neu";

    MultiLevelMesh mlMsh;
    // read coarse level mesh and generate finers level meshes
    double scalingFactor = 2.;
    //mlMsh.ReadCoarseMesh("./input/circle_quad.neu","seventh", scalingFactor);
    mlMsh.ReadCoarseMesh(filename.str().c_str(), "seventh", scalingFactor);
    /* "seventh" is the order of accuracy that is used in the gauss integration scheme
    probably in the furure it is not going to be an argument of this function   */
    unsigned dim = mlMsh.GetDimension();

    unsigned numberOfUniformLevels = 4;
    unsigned numberOfSelectiveLevels = 0;
    mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

    // erase all the coarse mesh levels
    mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

    // print mesh info
    mlMsh.PrintInfo();

    FEOrder feOrder = SECOND;
    l2Norm[i].resize(1);
    semiNorm[i].resize(1);

    // define the multilevel solution and attach the mlMsh object to it
    MultiLevelSolution mlSol(&mlMsh);

    // add variables to mlSol
    mlSol.AddSolution("X", LAGRANGE, feOrder, 2);
    mlSol.AddSolution("Y", LAGRANGE, feOrder, 2);
    mlSol.AddSolution("Z", LAGRANGE, feOrder, 2);
    mlSol.AddSolution("H", LAGRANGE, feOrder, 2);

    mlSol.Initialize("X", InitalValueXTorus);
    mlSol.Initialize("Y", InitalValueYTorus);
    mlSol.Initialize("Z", InitalValueZTorus);
    mlSol.Initialize("H", InitalValueHTorus);
    // attach the boundary condition function and generate boundary data
    mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionTorus);

//     mlSol.Initialize("X", InitalValueXSphere);
//     mlSol.Initialize("Y", InitalValueYSphere);
//     mlSol.Initialize("Z", InitalValueZSphere);
//     mlSol.Initialize("H", InitalValueHSphere);
//         // attach the boundary condition function and generate boundary data
//     mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionSphere);


    mlSol.GenerateBdc("X", "Steady");
    mlSol.GenerateBdc("Y", "Steady");
    mlSol.GenerateBdc("Z", "Steady");
    mlSol.GenerateBdc("H", "Steady");

    // define the multilevel problem attach the mlSol object to it
    MultiLevelProblem mlProb(&mlSol);

    // add system Wilmore in mlProb as a Linear Implicit System
    TransientNonlinearImplicitSystem& system = mlProb.add_system < TransientNonlinearImplicitSystem > ("Willmore");

    // add solution "X", "Y", "Z" and "H" to the system
    system.AddSolutionToSystemPDE("X");
    system.AddSolutionToSystemPDE("Y");
    system.AddSolutionToSystemPDE("Z");
    system.AddSolutionToSystemPDE("H");

    system.SetMaxNumberOfNonLinearIterations(6);

    // attach the assembling function to system
    system.SetAssembleFunction(AssembleWillmoreFlow_AD);

    // initilaize and solve the system
    system.init();

    std::vector < std::string > variablesToBePrinted;
    variablesToBePrinted.push_back("All");

    std::vector < std::string > surfaceVariables;
    surfaceVariables.push_back("X");
    surfaceVariables.push_back("Y");
    surfaceVariables.push_back("Z");

    VTKWriter vtkIO(&mlSol);
    vtkIO.SetSurfaceVariables(surfaceVariables);

    GMVWriter gmvIO(&mlSol);
    gmvIO.SetSurfaceVariables(surfaceVariables);
    gmvIO.SetDebugOutput(true);

    // time loop parameter
    system.AttachGetTimeIntervalFunction(SetVariableTimeStep);
    const unsigned int n_timesteps = 1;

    vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);
    gmvIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

    system.SetMgType(V_CYCLE);

    for (unsigned time_step = 0; time_step < n_timesteps; time_step++) {

      if (time_step > 0)
        system.SetMgType(V_CYCLE);

      system.MGsolve();

      system.UpdateSolution();


      vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, time_step + 1);
      gmvIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, time_step + 1);
    }




    /*  system.MGsolve();

    //       std::pair< double , double > norm = GetErrorNorm(&mlSol);
    //       l2Norm[i][j]  = norm.first;
    //       semiNorm[i][j] = norm.second;
    //       // print solutions
      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back("All");

      std::vector < std::string > surfaceVariables;
      surfaceVariables.push_back("X");
      surfaceVariables.push_back("Y");
      surfaceVariables.push_back("Z");

      VTKWriter vtkIO(&mlSol);
      vtkIO.SetSurfaceVariables(surfaceVariables);
      vtkIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, i);

      GMVWriter gmvIO(&mlSol);
      gmvIO.SetSurfaceVariables(surfaceVariables);
      gmvIO.SetDebugOutput(true);
      gmvIO.Pwrite(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, i);
    */

  }


//     std::vector < std::string > variablesToBePrinted;
//     variablesToBePrinted.push_back("All");
//
//     std::vector < std::string > surfaceVariables;
//     surfaceVariables.push_back("X");
//     surfaceVariables.push_back("Y");
//     surfaceVariables.push_back("Z");
//
//     VTKWriter vtkIO(&mlSol);
//     vtkIO.SetSurfaceVariables(surfaceVariables);
//
//     GMVWriter gmvIO(&mlSol);
//     gmvIO.SetSurfaceVariables(surfaceVariables);
//     gmvIO.SetDebugOutput(true);
//
//     // time loop parameter
//     //system.AttachGetTimeIntervalFunction(SetVariableTimeStep);
//     const unsigned int n_timesteps = 500;
//
//     vtkIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);
//     gmvIO.Pwrite(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);
//
//     //system.SetMgType(F_CYCLE);
//
//     for (unsigned time_step = 0; time_step < n_timesteps; time_step++) {
//
//       //if( time_step > 0 )
//         //system.SetMgType(V_CYCLE);
//
//       system.MGsolve();
//
//       //system.UpdateSolution();
//
//
//       vtkIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, time_step+1);
//       gmvIO.Pwrite(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, time_step+1);
//     }


//   // print the seminorm of the error and the order of convergence between different levels
//   std::cout << std::endl;
//   std::cout << std::endl;
//   std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
//   std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";
//
//   for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
//     std::cout << i + 1 << "\t";
//     std::cout.precision(14);
//
//     for (unsigned j = 0; j < 3; j++) {
//       std::cout << l2Norm[i][j] << "\t";
//     }
//
//     std::cout << std::endl;
//
//     if (i < maxNumberOfMeshes - 1) {
//       std::cout.precision(3);
//       std::cout << "\t\t";
//
//       for (unsigned j = 0; j < 3; j++) {
//         std::cout << log(l2Norm[i][j] / l2Norm[i + 1][j]) / log(2.) << "\t\t\t";
//       }
//
//       std::cout << std::endl;
//     }
//
//   }
//
//   std::cout << std::endl;
//   std::cout << std::endl;
//   std::cout << "SEMINORM ERROR and ORDER OF CONVERGENCE:\n\n";
//   std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";
//
//   for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
//     std::cout << i + 1 << "\t";
//     std::cout.precision(14);
//
//     for (unsigned j = 0; j < 3; j++) {
//       std::cout << semiNorm[i][j] << "\t";
//     }
//
//     std::cout << std::endl;
//
//     if (i < maxNumberOfMeshes - 1) {
//       std::cout.precision(3);
//       std::cout << "\t\t";
//
//       for (unsigned j = 0; j < 3; j++) {
//         std::cout << log(semiNorm[i][j] / semiNorm[i + 1][j]) / log(2.) << "\t\t\t";
//       }
//
//       std::cout << std::endl;
//     }
//
//   }



  return 0;
}




/**
 * Given the non linear problem
 *
 *      \Delta^2 u  = f(x),
 *      u(\Gamma) = 0
 *      \Delta u(\Gamma) = 0
 *
 * in the unit box \Omega centered in the origin with boundary \Gamma, where
 *
 *                      f(x) = \Delta^2 u_e ,
 *                    u_e = \cos ( \pi * x ) * \cos( \pi * y ),
 *
 * the following function assembles the system:
 *
 *      \Delta u = v
 *      \Delta v = f(x) = 4. \pi^4 u_e
 *      u(\Gamma) = 0
 *      v(\Gamma) = 0
 *
 * using automatic differentiation
 **/

void AssembleWillmoreFlow_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  TransientNonlinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientNonlinearImplicitSystem> ("Willmore");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh        = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned solRIndex[3];
  solRIndex[0] = mlSol->GetIndex("X");    // get the position of "X" in the ml_sol object
  solRIndex[1] = mlSol->GetIndex("Y");    // get the position of "Y" in the ml_sol object
  solRIndex[2] = mlSol->GetIndex("Z");    // get the position of "Z" in the ml_sol object
  unsigned solRType[3];
  solRType[0] = mlSol->GetSolutionType(solRIndex[0]);   // get the finite element type for "R"
  solRType[1] = mlSol->GetSolutionType(solRIndex[1]);   // get the finite element type for "R"
  solRType[2] = mlSol->GetSolutionType(solRIndex[2]);   // get the finite element type for "R"

  unsigned solRPdeIndex[3];
  solRPdeIndex[0] = mlPdeSys->GetSolPdeIndex("X");    // get the position of "X" in the pdeSys object
  solRPdeIndex[1] = mlPdeSys->GetSolPdeIndex("Y");    // get the position of "Y" in the pdeSys object
  solRPdeIndex[2] = mlPdeSys->GetSolPdeIndex("Z");    // get the position of "Z" in the pdeSys object

  vector < adept::adouble >  solR[3]; // local solution
  vector < double > solR_old[3];


  unsigned solHIndex;
  solHIndex = mlSol->GetIndex("H");    // get the position of "H" in the ml_sol object
  unsigned solHType = mlSol->GetSolutionType(solHIndex);    // get the finite element type for "H"

  unsigned solHPdeIndex;
  solHPdeIndex = mlPdeSys->GetSolPdeIndex("H");    // get the position of "H" in the pdeSys object

  vector < adept::adouble >  solH; // local solution

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector< int > sysDof; // local to global pdeSys dofs
  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  vector< double > Res; // local redidual vector
  vector< adept::adouble > aResR[3]; // local redidual vector
  vector< adept::adouble > aResH; // local redidual vector


  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  for (int i = 0; i < 3; i++) {
    solR[i].reserve(maxSize);
    solR_old[i].reserve(maxSize);
  }

  solH.reserve(maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);

  sysDof.reserve(4 * maxSize);
  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  phi_xx.reserve(maxSize * dim2);

  Res.reserve(4 * maxSize);

  for (int i = 0; i < 3; i++) {
    aResR[i].reserve(maxSize);
  }

  aResH.reserve(maxSize);

  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)
  Jac.reserve(4 * maxSize * 4 * maxSize);


  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs  = msh->GetElementDofNumber(iel, solHType);    // number of solution element dofs
    unsigned nDofs2 = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    
    // resize local arrays
    sysDof.resize(4 * nDofs);

    for (int i = 0; i < 3; i++) {
      solR[i].resize(nDofs);
      solR_old[i].resize(nDofs);
    }

    solH.resize(nDofs);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofs2);
    }

    Res.resize(4 * nDofs);    //resize

    for (int i = 0; i < 3; i++) {
      aResR[i].resize(nDofs);    //resize
    }

    aResH.resize(nDofs);    //resize

    for (int i = 0; i < 3; i++) {
      std::fill(aResR[i].begin(), aResR[i].end(), 0);    //set aRes to zero
    }

    std::fill(aResH.begin(), aResH.end(), 0);    //set aRes to zero

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, solHType);    // global to global mapping between solution node and solution dof

      for (int k = 0; k < 3; k++) {
        solR[k][i] = (*sol->_Sol[solRIndex[k]])(solDof);      // global extraction and local storage for the solution
        solR_old[k][i] = (*sol->_SolOld[solRIndex[k]])(solDof);      // global extraction and local storage for the solution

      }

      solH[i] = (*sol->_Sol[solHIndex])(solDof);      // global extraction and local storage for the solution

      for (int k = 0; k < 3; k++) {
        sysDof[k * nDofs + i] = pdeSys->GetSystemDof(solRIndex[k], solRPdeIndex[k], i, iel);  // global to global mapping between solution node and pdeSys dof
      }

      sysDof[3 * nDofs + i] = pdeSys->GetSystemDof(solHIndex, solHPdeIndex, i, iel);  // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs2; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned idim = 0; idim < dim; idim++) {
        x[idim][i] = (*msh->_topology->_Sol[idim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }


    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();


    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solHType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solHType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      adept::adouble solRGauss[3];
      double solRGaussOld[3];
      adept::adouble solRGauss_x[3][2];
      adept::adouble solRGauss_xx[3][2][2];

      adept::adouble sol_x[2];
      sol_x[0] = sol_x[1] = 0.;

      for (int k = 0; k < 3; k++) {
        solRGauss[k] = 0.;
        solRGaussOld[k] = 0.;

        for (int i = 0; i < dim; i++) {
          solRGauss_x[k][i] = 0.;

          for (int j = 0; j < dim; j++) {
            solRGauss_xx[k][i][j] = 0.;
          }
        }
      }

      adept::adouble solHGauss = 0;
      adept::adouble solHGauss_x[2] = {0., 0.};

      for (unsigned i = 0; i < nDofs; i++) {
        for (int k = 0; k < 3; k++) {
          solRGauss[k] += phi[i] * solR[k][i];
          solRGaussOld[k] += phi[i] * solR_old[k][i];
        }

        solHGauss += phi[i] * solH[i];

        for (unsigned u = 0; u < dim; u++) {
          sol_x[u] += phi[i] * x[u][i];
        }

        for (unsigned u = 0; u < dim; u++) { // gradient
          for (int k = 0; k < 3; k++) {
            solRGauss_x[k][u] += phi_x[i * dim + u] * solR[k][i];
          }

          solHGauss_x[u] += phi_x[i * dim + u] * solH[i];
        }

        for (unsigned u = 0; u < dim; u++) {  // hessian
          for (unsigned v = 0; v < dim; v++) {
            unsigned uvindex = 0; //_uu

            if (u != v) uvindex = 2;  //_uv or _vu
            else if (u == 1) uvindex = 1;  //_vv

            for (int k = 0; k < 3; k++) {
              solRGauss_xx[k][u][v] += phi_xx[i * dim2 + uvindex] * solR[k][i];
            }
          }
        }
      }

      adept::adouble g[2][2];


      g[0][0] = g[0][1] = g[1][0] = g[1][1] = 0.;

      for (int k = 0; k < 3; k++) {
        for (int u = 0; u < dim; u++) {
          for (int v = 0; v < dim; v++) {
            g[u][v] += solRGauss_x[k][u] * solRGauss_x[k][v];
          }
        }
      }

      adept::adouble detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];

      adept::adouble  A = sqrt(detg);

      adept::adouble gI[2][2];

      gI[0][0] =  g[1][1] / detg;
      gI[0][1] = -g[0][1] / detg;
      gI[1][0] = -g[1][0] / detg;
      gI[1][1] =  g[0][0] / detg;

      adept::adouble N[3];

      N[0] = (solRGauss_x[1][0] * solRGauss_x[2][1] - solRGauss_x[1][1] * solRGauss_x[2][0]) / A;
      N[1] = (solRGauss_x[2][0] * solRGauss_x[0][1] - solRGauss_x[2][1] * solRGauss_x[0][0]) / A;
      N[2] = (solRGauss_x[0][0] * solRGauss_x[1][1] - solRGauss_x[0][1] * solRGauss_x[1][0]) / A;

      adept::adouble h[2][2];

      h[0][0] = h[0][1] = h[1][0] = h[1][1] = 0.;

      for (int k = 0; k < 3; k++) {
        for (int u = 0; u < dim; u++) {
          for (int v = 0; v < dim; v++) {
            h[u][v] += solRGauss_xx[k][u][v] * N[k];

          }
        }
      }

      //adept::adouble K = cos(sol_x[0])/(a+cos(sol_x[0]));//(h[0][0]*h[1][1]-h[0][1]*h[1][0])/detg;

      adept::adouble K = (h[0][0] * h[1][1] - h[0][1] * h[1][0]) / detg;

      adept::adouble H_exact = 0.5 * (1. + cos(sol_x[0]) / (a + cos(sol_x[0])));

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {

        for (int k = 0; k < 3; k++) {
          for (int u = 0; u < dim; u++) {
            adept::adouble AgIgradRgradPhi = 0;

            for (int v = 0; v < dim; v++) {
              AgIgradRgradPhi += A * gI[u][v].value() * solRGauss_x[k][v];
            }

            aResR[k][i] += AgIgradRgradPhi * phi_x[i * dim + u] * weight;
          }

          aResR[k][i] += 2.* A * solHGauss.value() * N[k] * phi[i] * weight;
        }


        for (int u = 0; u < dim; u++) {
          adept::adouble AgIgradHgradPhi = 0;

          for (int v = 0; v < dim; v++) {
            AgIgradHgradPhi += A * gI[u][v].value() * solHGauss_x[v];
          }

          aResH[i] -= AgIgradHgradPhi * phi_x[i * dim + u] * weight;
        }

        aResH[i] += A * (-0 * (solRGauss[0] - solRGaussOld[0]) * N[0].value()
                         - 0 * (solRGauss[1] - solRGaussOld[1]) * N[1].value()
                         - 0 * (solRGauss[2] - solRGaussOld[2]) * N[2].value()
                         + 2. * solHGauss * (solHGauss * solHGauss  - K.value())) * phi[i] * weight;

      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    for (int i = 0; i < nDofs; i++) {
      for (int k = 0; k < 3; k++) {
        Res[ k * nDofs + i] = -aResR[k][i].value();
      }

      Res[ 3 * nDofs + i] = -aResH[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);


    Jac.resize((4 * nDofs) * (4 * nDofs));

    // define the dependent variables
    for (int k = 0; k < 3; k++) {
      s.dependent(&aResR[k][0], nDofs);
    }

    s.dependent(&aResH[0], nDofs);

    // define the independent variables
    for (int k = 0; k < 3; k++) {
      s.independent(&solR[k][0], nDofs);
    }

    s.independent(&solH[0], nDofs);
    // get the jacobian matrix (ordered by row)
    s.jacobian(&Jac[0], true);

    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}

// functions post processing

double GetExactSolutionValueSphere(const std::vector < double >& x) {
  return 0;
};

void GetExactSolutionGradientSphere(const std::vector < double >& x, vector < double >& solGrad) {

};


double GetExactSolutionValueTorus(const std::vector < double >& x) {

  return 0;

};

void GetExactSolutionGradientTorus(const std::vector < double >& x, vector < double >& solGrad) {

};




std::pair < double, double > GetErrorNorm(MultiLevelSolution* mlSol) {
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*          msh          = mlSol->_mlMesh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)
  Solution*    sol        = mlSol->GetSolutionLevel(level);    // pointer to the solution (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  vector < double >  solu; // local solution

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  solu.reserve(maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  phi_xx.reserve(maxSize * dim2);

  double seminorm = 0.;
  double l2norm = 0.;

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofs2 = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    solu.resize(nDofs);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofs2);
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs2; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned idim = 0; idim < dim; idim++) {
        x[idim][i] = (*msh->_topology->_Sol[idim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double soluGauss = 0;
      vector < double > soluGauss_x(dim, 0.);
      vector < double > xGauss(dim, 0.);

      for (unsigned i = 0; i < nDofs; i++) {
        soluGauss += phi[i] * solu[i];

        for (unsigned idim = 0; idim < dim; idim++) {
          soluGauss_x[idim] += phi_x[i * dim + idim] * solu[i];
          xGauss[idim] += x[idim][i] * phi[i];
        }
      }

      double exactSol;
      vector <double> solGrad(dim);

//       if (simulation == 1) {
//         exactSol = GetExactSolutionValueSphere(xGauss);
//         GetExactSolutionGradientSphere(xGauss, solGrad);
//       } else if (simulation == 2) {
//         exactSol = GetExactSolutionValueTorus(xGauss);
//         GetExactSolutionGradientTorus(xGauss, solGrad);
//       }

      l2norm += (exactSol - soluGauss) * (exactSol - soluGauss) * weight;

      for (unsigned j = 0; j < dim ; j++) {
        seminorm   += ((soluGauss_x[j] - solGrad[j]) * (soluGauss_x[j] - solGrad[j])) * weight;
      }


    } // end gauss point loop
  } //end element loop for each process

  // add the norms of all processes
  NumericVector* norm_vec;
  norm_vec = NumericVector::build().release();
  norm_vec->init(msh->n_processors(), 1 , false, AUTOMATIC);

  norm_vec->set(iproc, l2norm);
  norm_vec->close();
  l2norm = norm_vec->l1_norm();

  norm_vec->set(iproc, seminorm);
  norm_vec->close();
  seminorm = norm_vec->l1_norm();

  delete norm_vec;

  std::pair < double, double > norm;
  norm.first  = sqrt(l2norm);
  norm.second = sqrt(seminorm);

  return norm;

}
