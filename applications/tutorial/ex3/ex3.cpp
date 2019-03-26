/** tutorial/Ex3
 * This example shows how to set and solve the weak form of the nonlinear problem
 *                     -\Delta u + < u,u,u > \cdot \nabla u = f(x) \text{ on }\Omega,
 *            u=0 \text{ on } \Gamma,
 * on a box domain $\Omega$ with boundary $\Gamma$;
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


using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0;

  if (facename == 2)
    dirichlet = false;

  return dirichlet;
}

void AssembleNonlinearProblem(MultiLevelProblem& ml_prob);

void AssembleNonlinearProblem_AD(MultiLevelProblem& ml_prob);

std::pair < double, double > GetErrorNorm(MultiLevelSolution* mlSol);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);


  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/square_quad.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/square_tri.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/square_mixed.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_wedge.neu","seventh",scalingFactor);
  mlMsh.ReadCoarseMesh("./input/cube_tet.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_mixed.neu","seventh",scalingFactor);

  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
    probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();
  unsigned maxNumberOfMeshes;

  if (dim == 2) {
    maxNumberOfMeshes = 7;
  } else {
    maxNumberOfMeshes = 6;
  }

  vector < vector < double > > l2Norm;
  l2Norm.resize(maxNumberOfMeshes);

  vector < vector < double > > semiNorm;
  semiNorm.resize(maxNumberOfMeshes);

  for (unsigned i = 1; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

    unsigned numberOfUniformLevels = i ;
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
      mlSol.Initialize("All");

      // attach the boundary condition function and generate boundary data
      mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
      mlSol.GenerateBdc("u");

      // define the multilevel problem attach the mlSol object to it
      MultiLevelProblem mlProb(&mlSol);

      // add system Poisson in mlProb as a Linear Implicit System
      NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("Poisson");

      // add solution "u" to system
      system.AddSolutionToSystemPDE("u");

      // attach the assembling function to system
      system.SetAssembleFunction(AssembleNonlinearProblem_AD);

      // initilaize and solve the system
      system.init();
      system.SetOuterSolver(PREONLY);
      system.MGsolve();

      std::pair< double , double > norm = GetErrorNorm(&mlSol);
      l2Norm[i][j]  = norm.first;
      semiNorm[i][j] = norm.second;
      // print solutions
      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back("All");

      VTKWriter vtkIO(&mlSol);
      vtkIO.SetDebugOutput(true);
      vtkIO.Write(DEFAULT_OUTPUTDIR, "linear", variablesToBePrinted, i+j*10);
      vtkIO.Write(DEFAULT_OUTPUTDIR, "quadratic", variablesToBePrinted, i+j*10);
      vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, i+j*10);


//       GMVWriter gmvIO(&mlSol);
//       gmvIO.SetDebugOutput(true);
//       gmvIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted,i);
    }
  }

  // print the seminorm of the error and the order of convergence between different levels
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 1; i < maxNumberOfMeshes; i++) {
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

  for (unsigned i = 1; i < maxNumberOfMeshes; i++) {
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





}


double GetExactSolutionValue(const std::vector < double >& x) {
  double pi = acos(-1.);
  return cos(pi * x[0]) * cos(pi * x[1]);
};


void GetExactSolutionGradient(const std::vector < double >& x, vector < double >& solGrad) {
  double pi = acos(-1.);
  solGrad[0]  = -pi * sin(pi * x[0]) * cos(pi * x[1]);
  solGrad[1] = -pi * cos(pi * x[0]) * sin(pi * x[1]);
};


double GetExactSolutionLaplace(const std::vector < double >& x) {
  double pi = acos(-1.);
  return -2.*pi * pi * cos(pi * x[0]) * cos(pi * x[1]);       // - pi*pi*cos(pi*x[0])*cos(pi*x[1]);
};

/**
 * Given the non linear problem
 *
 *      - \Delta u + < u, u, u > \cdot \nabla u = f(x)
 * in the unit box centered in the origin with
 *                      f(x) = - \Delta u_e + < u_e, u_e, u_e > \cdot \nabla u_e
 *                    u_e = \cos ( \pi * x ) * \cos( \pi * y ),
 * the following function assembles (explicitely) the Jacobian matrix J(u^i) and the residual vector Res(u^i)
 * for the Newton iteration, i.e.
 *
 *                  J(u^i) w = Res(u^i) = f(x) - ( - \Delta u^i + < u^i, u^i, u^i > \cdot \nabla u^i ),
 *        u^{i+1} = u^i + w
 *        where
 *        J(u^i) w = - \Delta w  + < w , w , w >  \cdot \nabla u^i + < u^i , u^i , u^i >  \cdot \nabla w
 *
 **/

void AssembleNonlinearProblem(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem > ("Poisson");   // pointer to the linear implicit system named "Poisson"

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

  vector < double >  solu; // local solution

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector< int > sysDof; // local to global pdeSys dofs
  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  vector< double > Res; // local redidual vector
  vector< double > J; // local Jacobian matrix

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  solu.reserve(maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);

  sysDof.reserve(maxSize);
  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  phi_xx.reserve(maxSize * dim2);
  Res.reserve(maxSize);
  J.reserve(maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {


    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofs2 = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    // resize local arrays
    sysDof.resize(nDofs);
    solu.resize(nDofs);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofs2);
    }

    Res.resize(nDofs);    //resize
    std::fill(Res.begin(), Res.end(), 0);    //set Res to zero

    J.resize(nDofs * nDofs);    //resize
    std::fill(J.begin(), J.end(), 0);    //set K to zero

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      sysDof[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs2; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
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

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          soluGauss_x[jdim] += phi_x[i * dim + jdim] * solu[i];
          xGauss[jdim] += x[jdim][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {

        double nonLinearTerm = 0.;
        double mLaplace = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          mLaplace   +=  phi_x[i * dim + jdim] * soluGauss_x[jdim];
          nonLinearTerm += soluGauss * soluGauss_x[jdim] * phi[i];
        }

        double exactSolValue = GetExactSolutionValue(xGauss);
        vector < double > exactSolGrad(dim);
        GetExactSolutionGradient(xGauss , exactSolGrad);
        double exactSolLaplace = GetExactSolutionLaplace(xGauss);


        double f = (- exactSolLaplace + exactSolValue * (exactSolGrad[0] + exactSolGrad[1])) * phi[i] ;
        Res[i] += (f - (mLaplace + nonLinearTerm)) * weight;

        // *** phi_j loop ***
        for (unsigned j = 0; j < nDofs; j++) {
          mLaplace = 0.;
          nonLinearTerm = 0.;

          for (unsigned kdim = 0; kdim < dim; kdim++) {
            mLaplace += (phi_x[i * dim + kdim] * phi_x[j * dim + kdim]);
            nonLinearTerm +=   phi[i] * (phi[j] * soluGauss_x[kdim] + soluGauss * phi_x[j * dim + kdim]);
          }

          J[i * nDofs + j] += (mLaplace + nonLinearTerm) * weight;
        } // end phi_j loop
      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    RES->add_vector_blocked(Res, sysDof);

    KK->add_matrix_blocked(J, sysDof, sysDof);
  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}


/**
 * Given the non linear problem
 *
 *      - \Delta u + < u, u, u > \cdot \nabla u = f(x),
 *
 * in the unit box centered in the origin with
 *
 *                      f(x) = - \Delta u_e + < u_e, u_e, u_e > \cdot \nabla u_e,
 *                    u_e = \cos ( \pi * x ) * \cos( \pi * y ),
 *
 *the following function assembles the residual vector Res(u^i) and using automatic differentiation gets
 *the exact Jacobian matrix J(u^i) for the Newton iteration, i.e.
 *
 *                   J(u^i) w = Res(u^i) = f(x) - ( - \Delta u^i + < u^i, u^i, u^i > \cdot \nabla u^i ),
 *         u^{i+1} = u^i + w,
 *        where
 *        J(u^i) w = - \Delta w  + < w , w , w >  \cdot \nabla u^i + < u^i , u^i , u^i >  \cdot \nabla w.
 *
 **/

void AssembleNonlinearProblem_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use

  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem > ("Poisson");   // pointer to the linear implicit system named "Poisson"
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

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector< int > sysDof; // local to global pdeSys dofs
  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  vector< double > Res; // local redidual vector
  vector< adept::adouble > aRes; // local redidual vector


  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  solu.reserve(maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);

  sysDof.reserve(maxSize);
  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  phi_xx.reserve(maxSize * dim2);
  Res.reserve(maxSize);
  aRes.reserve(maxSize);

  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)
  Jac.reserve(maxSize * maxSize);
  vector< double > Jact; // local Jacobian matrix (ordered by raw, PETSC)
  Jact.reserve(maxSize * maxSize);



  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {


    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofs2 = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    // resize local arrays
    sysDof.resize(nDofs);
    solu.resize(nDofs);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofs2);
    }

    Res.resize(nDofs);    //resize
    std::fill(Res.begin(), Res.end(), 0);    //set Res to zero

    aRes.resize(nDofs);    //resize
    std::fill(aRes.begin(), aRes.end(), 0);    //set aRes to zero


    Jact.resize(nDofs * nDofs);
    Jac.resize(nDofs * nDofs);


    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      sysDof[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs2; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
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
      vector < double > xGauss(dim, 0.);

      for (unsigned i = 0; i < nDofs; i++) {
        soluGauss += phi[i] * solu[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          soluGauss_x[jdim] += phi_x[i * dim + jdim] * solu[i];
          xGauss[jdim] += x[jdim][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {

        adept::adouble nonLinearTerm = 0.;
        adept::adouble mLaplace = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          mLaplace   +=  phi_x[i * dim + jdim] * soluGauss_x[jdim];
          //nonLinearTerm += soluGauss * soluGauss_x[jdim] * phi[i];
        }

        double exactSolValue = GetExactSolutionValue(xGauss);
        vector < double > exactSolGrad(dim);
        GetExactSolutionGradient(xGauss , exactSolGrad);
        double exactSolLaplace = GetExactSolutionLaplace(xGauss);


        double f = (- exactSolLaplace + 0*exactSolValue * (exactSolGrad[0] + exactSolGrad[1])) * phi[i] ;
        aRes[i] += (f - (mLaplace + nonLinearTerm)) * weight;

      } // end phi_i loop
    } // end gauss point loop


    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    for (int i = 0; i < nDofs; i++) {
      Res[i] = aRes[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);



    // define the dependent variables
    s.dependent(&aRes[0], nDofs);

    // define the independent variables
    s.independent(&solu[0], nDofs);

    // get the jacobian matrix (ordered by column)
    s.jacobian(&Jac[0]);

    // get the jacobian matrix (ordered by raw, i.e. Jact=Jac^t)
    for (int inode = 0; inode < nDofs; inode++) {
      for (int jnode = 0; jnode < nDofs; jnode++) {
        Jact[inode * nDofs + jnode] = -Jac[jnode * nDofs + inode];
      }
    }

    //store Jact in the global matrix KK
    KK->add_matrix_blocked(Jact, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}





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

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
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

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          soluGauss_x[jdim] += phi_x[i * dim + jdim] * solu[i];
          xGauss[jdim] += x[jdim][i] * phi[i];
        }
      }

      vector <double> solGrad(dim);
      GetExactSolutionGradient(xGauss, solGrad);

      for (unsigned j = 0; j < dim ; j++) {
        seminorm   += ((soluGauss_x[j] - solGrad[j]) * (soluGauss_x[j] - solGrad[j])) * weight;
      }

      double exactSol = GetExactSolutionValue(xGauss);
      l2norm += (exactSol - soluGauss) * (exactSol - soluGauss) * weight;
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
