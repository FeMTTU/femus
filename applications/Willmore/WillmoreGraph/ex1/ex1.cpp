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
#include "adept.h"
#include <cstdlib>


using namespace femus;

int simulation = 1; // =1 sphere (default) = 2 torus

//Sphere

double thetaSphere = acos(-1.) / 6;

bool SetBoundaryConditionSphere(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet

  if (!strcmp("u", SolName)) {
    value = tan(thetaSphere);
  } else if (!strcmp("W", SolName)) {
    value = -1. / tan(thetaSphere);
  }

  return dirichlet;
}

double InitalValueUSphere(const std::vector < double >& x) {
  return tan(thetaSphere);
}

double InitalValueWSphere(const std::vector < double >& x) {
  return -1. / tan(thetaSphere);
}

// Torus

bool SetBoundaryConditionTorus(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet

  double theta = acos(-1.) / 6;
  double z = sin(theta);

  if (!strcmp("u", SolName)) {
    value = z;
  } else if (!strcmp("W", SolName)) {
    double theta1 = theta;
    double theta2 = acos(-1.) - theta1;
    double A = 1. / z;
    double H1 =  0.5 * (1. + cos(theta1) / (sqrt(2.) + cos(theta1)));
    double H2 =  0.5 * (1. + cos(theta2) / (sqrt(2.) + cos(theta2)));

    if (facename == 1) {
      value = - A * H1;
    } else if (facename == 2) {
      value = - A * H2;
    }
  }

  return dirichlet;
}

double InitalValueUTorus(const std::vector < double >& x) {
  double r = sqrt(x[0] * x[0] + x[1] * x[1]);
  double cosu = r - sqrt(2) ;
  return sqrt(1 - cosu * cosu);
  //return 0.5;
}

double InitalValueWTorus(const std::vector < double >& x) {

  double r = sqrt(x[0] * x[0] + x[1] * x[1]);
  double cosu = r - sqrt(2) ;
  double sinu = sqrt(1 - cosu * cosu);

  return - 0.5 / sinu * (1 + cosu / (sqrt(2.) + cosu));
  //return 0.25;
}


void AssembleWillmoreProblem_AD(MultiLevelProblem& ml_prob);

std::pair < double, double > GetErrorNorm(MultiLevelSolution* mlSol);

int main(int argc, char** args) {


  if (argc >= 2) {
    if (!strcmp("sphere", args[1]) || !strcmp("Sphere", args[1]) || !strcmp("SPHERE", args[1])) {
      simulation = 1;

      if (argc >= 3) {
        std::string str;
        std::stringstream ss;
        ss << args[2];
        ss >> str;
        int angle = atoi(str.c_str());
        thetaSphere = acos(-1.) / 180 * angle;
      }
    } else if (!strcmp("torus", args[1]) || !strcmp("Torus", args[1]) || !strcmp("TORUS", args[1])) simulation = 2;
    else {
      std::cout << "Wrong input, using default argument: simulation = 1 (Sphere)" << std::endl;
    }
  } else {
    std::cout << "No input argument, using default argument: simulation = 1 (Sphere)" << std::endl;
  }


  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);


  // define multilevel mesh


  unsigned maxNumberOfMeshes;
  maxNumberOfMeshes = 5;

  vector < vector < double > > l2Norm;
  l2Norm.resize(maxNumberOfMeshes);

  vector < vector < double > > semiNorm;
  semiNorm.resize(maxNumberOfMeshes);

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

    std::ostringstream filename;

    if (simulation == 1) {
      filename << "./input/circle_quad" << i << ".neu";
    } else if (simulation == 2) {
      filename << "./input/torus30_" << i << ".neu";
    }

    MultiLevelMesh mlMsh;
    // read coarse level mesh and generate finers level meshes
    double scalingFactor = 1.;
    //mlMsh.ReadCoarseMesh("./input/circle_quad.neu","seventh", scalingFactor);
    mlMsh.ReadCoarseMesh(filename.str().c_str(), "seventh", scalingFactor);
    /* "seventh" is the order of accuracy that is used in the gauss integration scheme
    probably in the furure it is not going to be an argument of this function   */
    unsigned dim = mlMsh.GetDimension();

    unsigned numberOfUniformLevels = 1;
    unsigned numberOfSelectiveLevels = 0;
    mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

    // erase all the coarse mesh levels
    mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

    // print mesh info
    mlMsh.PrintInfo();

    FEOrder feOrder[3] = {FIRST, SERENDIPITY, SECOND};
    l2Norm[i].resize(3);
    semiNorm[i].resize(3);

    for (unsigned j = 0; j < 3; j++) {   // loop on the FE Order
      // define the multilevel solution and attach the mlMsh object to it
      MultiLevelSolution mlSol(&mlMsh);

      // add variables to mlSol
      mlSol.AddSolution("u", LAGRANGE, feOrder[j]);
      mlSol.AddSolution("W", LAGRANGE, feOrder[j]);

      if (simulation == 1) {
        mlSol.Initialize("u", InitalValueUSphere);
        mlSol.Initialize("W", InitalValueWSphere);
        // attach the boundary condition function and generate boundary data
        mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionSphere);
        mlSol.GenerateBdc("u");
        mlSol.GenerateBdc("W");
      } else if (simulation == 2) {
        mlSol.Initialize("u", InitalValueUTorus);
        mlSol.Initialize("W", InitalValueWTorus);
        // attach the boundary condition function and generate boundary data
        mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionTorus);
        mlSol.GenerateBdc("u");
        mlSol.GenerateBdc("W");
      }

      // define the multilevel problem attach the mlSol object to it
      MultiLevelProblem mlProb(&mlSol);

      // add system Poisson in mlProb as a Linear Implicit System
      NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("Willmore");

      // add solution "u" to system
      system.AddSolutionToSystemPDE("u");
      system.AddSolutionToSystemPDE("W");

      // attach the assembling function to system
      system.SetAssembleFunction(AssembleWillmoreProblem_AD);

      // initilaize and solve the system
      system.init();
      system.MGsolve();

      std::pair< double , double > norm = GetErrorNorm(&mlSol);
      l2Norm[i][j]  = norm.first;
      semiNorm[i][j] = norm.second;
      // print solutions
      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back("All");

      VTKWriter vtkIO(&mlSol);
      vtkIO.SetGraphVariable("u");
      vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, i);

    }
  }

  // print the seminorm of the error and the order of convergence between different levels
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < 3; j++) {
      std::cout << l2Norm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < 3; j++) {
        std::cout << log(l2Norm[i][j] / l2Norm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "SEMINORM ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < 3; j++) {
      std::cout << semiNorm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < 3; j++) {
        std::cout << log(semiNorm[i][j] / semiNorm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }



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

void AssembleWillmoreProblem_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("Willmore");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solu; // local solution

  unsigned solWIndex;
  solWIndex = mlSol->GetIndex("W");    // get the position of "v" in the ml_sol object
  unsigned solWType = mlSol->GetSolutionType(solWIndex);    // get the finite element type for "v"

  unsigned solWPdeIndex;
  solWPdeIndex = mlPdeSys->GetSolPdeIndex("W");    // get the position of "v" in the pdeSys object

  vector < adept::adouble >  solW; // local solution



  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector< int > sysDof; // local to global pdeSys dofs
  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  vector< double > Res; // local redidual vector
  vector< adept::adouble > aResu; // local redidual vector
  vector< adept::adouble > aResW; // local redidual vector


  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  solu.reserve(maxSize);
  solW.reserve(maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);

  sysDof.reserve(2 * maxSize);
  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  phi_xx.reserve(maxSize * dim2);

  Res.reserve(2 * maxSize);
  aResu.reserve(maxSize);
  aResW.reserve(maxSize);

  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)
  Jac.reserve(4 * maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // element geometry type
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofs2 = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    // resize local arrays
    sysDof.resize(2 * nDofs);
    solu.resize(nDofs);
    solW.resize(nDofs);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofs2);
    }

    Res.resize(2 * nDofs);  //resize
    aResu.resize(nDofs);    //resize
    aResW.resize(nDofs);    //resize

    std::fill(aResu.begin(), aResu.end(), 0);    //set aRes to zero
    std::fill(aResW.begin(), aResW.end(), 0);    //set aRes to zero

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      solW[i] = (*sol->_Sol[solWIndex])(solDof);      // global extraction and local storage for the solution
      sysDof[i]       = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
      sysDof[nDofs + i] = pdeSys->GetSystemDof(solWIndex, solWPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
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
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      adept::adouble soluGauss = 0;
      vector < adept::adouble > soluGauss_x(dim, 0.);

      adept::adouble solWGauss = 0;
      vector < adept::adouble > solWGauss_x(dim, 0.);

      vector < double > xGauss(dim, 0.);

      for (unsigned i = 0; i < nDofs; i++) {
        soluGauss += phi[i] * solu[i];
        solWGauss += phi[i] * solW[i];

        for (unsigned idim = 0; idim < dim; idim++) {
          soluGauss_x[idim] += phi_x[i * dim + idim] * solu[i];
          solWGauss_x[idim] += phi_x[i * dim + idim] * solW[i];
          xGauss[idim] += x[idim][i] * phi[i];
        }
      }

      double c = 0.;
      double Id[2][2] = {{1., 0.}, {0., 1.}};
      adept::adouble A2 = 1.;
      vector < vector < adept::adouble> > B(dim);

      for (unsigned idim = 0; idim < dim; idim++) {
        B[idim].resize(dim);
        A2 += soluGauss_x[idim] * soluGauss_x[idim];
      }

      adept::adouble A = sqrt(A2);

      for (unsigned idim = 0; idim < dim; idim++) {
        for (unsigned jdim = 0; jdim < dim; jdim++) {
          B[idim][jdim] = Id[idim][jdim] - (soluGauss_x[idim] * soluGauss_x[jdim]) / A2;
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {

        adept::adouble nonLinearLaplaceU = 0.;
        adept::adouble nonLinearLaplaceW = 0.;

        for (unsigned idim = 0; idim < dim; idim++) {

          nonLinearLaplaceU +=  - 1. / A  * soluGauss_x[idim] * phi_x[i * dim + idim];

          nonLinearLaplaceW +=   -1. / A * ((B[idim][0] * solWGauss_x[0] +
                                             B[idim][1] * solWGauss_x[1])
                                             - (solWGauss * solWGauss / A2 + c) *
                                                soluGauss_x[idim]) * phi_x[i * dim + idim];

        }

        aResu[i] += (2.*solWGauss / A * phi[i] - nonLinearLaplaceU) * weight;
        aResW[i] += nonLinearLaplaceW * weight;

      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    for (int i = 0; i < nDofs; i++) {
      Res[i]       = -aResu[i].value();
      Res[nDofs + i] = -aResW[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    Jac.resize((2. * nDofs) *(2. * nDofs));
    // define the dependent variables
    s.dependent(&aResu[0], nDofs);
    s.dependent(&aResW[0], nDofs);

    // define the independent variables
    s.independent(&solu[0], nDofs);
    s.independent(&solW[0], nDofs);
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
  return sqrt(1. / (cos(thetaSphere) * cos(thetaSphere)) - x[0] * x[0] - x[1] * x[1]);
};

void GetExactSolutionGradientSphere(const std::vector < double >& x, vector < double >& solGrad) {
  double pi = acos(-1.);
  solGrad[0]  = -x[0] / sqrt(1. / (cos(thetaSphere) * cos(thetaSphere)) - x[0] * x[0] - x[1] * x[1]);
  solGrad[1]  = -x[1] / sqrt(1. / (cos(thetaSphere) * cos(thetaSphere)) - x[0] * x[0] - x[1] * x[1]);
};


double GetExactSolutionValueTorus(const std::vector < double >& x) {

  double r = sqrt(x[0] * x[0] + x[1] * x[1]);
  double cosu =  r - sqrt(2) ;
  return sqrt(1 - cosu * cosu);

};

void GetExactSolutionGradientTorus(const std::vector < double >& x, vector < double >& solGrad) {

  double r = sqrt(x[0] * x[0] + x[1] * x[1]);
  double cosu =  r - sqrt(2) ;
  double z = sqrt(1 - cosu * cosu);

  solGrad[0] =  - cosu / (r * z) * x[0];
  solGrad[1] =  - cosu / (r * z) * x[1];

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

    // element geometry type
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

      if (simulation == 1) {
        exactSol = GetExactSolutionValueSphere(xGauss);
        GetExactSolutionGradientSphere(xGauss, solGrad);
      } else if (simulation == 2) {
        exactSol = GetExactSolutionValueTorus(xGauss);
        GetExactSolutionGradientTorus(xGauss, solGrad);
      }

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
