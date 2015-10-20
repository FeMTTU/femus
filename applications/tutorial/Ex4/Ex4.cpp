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
#include "adept.h"


using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0;
}

void AssembleBilaplaceProblem_AD(MultiLevelProblem& ml_prob);

std::pair < double, double > GetErrorNorm(MultiLevelSolution* mlSol);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);


  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  mlMsh.ReadCoarseMesh("./input/square_quad.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
    probably in the furure it is not going to be an argument of this function   */

  unsigned maxNumberOfMeshes = 5;

  vector < vector < double > > l2Norm;
  l2Norm.resize(maxNumberOfMeshes);

  vector < vector < double > > semiNorm;
  semiNorm.resize(maxNumberOfMeshes);

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

    unsigned numberOfUniformLevels = i + 1;
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

      std::cout << "level = " << i << " FEM = " << j << std::endl;

      // define the multilevel solution and attach the mlMsh object to it
      MultiLevelSolution mlSol(&mlMsh);

      // add variables to mlSol
      mlSol.AddSolution("u", LAGRANGE, feOrder[j]);
      mlSol.AddSolution("v", LAGRANGE, feOrder[j]);
      mlSol.Initialize("All");

      // attach the boundary condition function and generate boundary data
      mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
      mlSol.GenerateBdc("u");
      mlSol.GenerateBdc("v");

      // define the multilevel problem attach the mlSol object to it
      MultiLevelProblem mlProb(&mlSol);

      // add system Poisson in mlProb as a Linear Implicit System
      NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("Poisson");

      // add solution "u" to system
      system.AddSolutionToSystemPDE("u");
      system.AddSolutionToSystemPDE("v");

      // attach the assembling function to system
      system.SetAssembleFunction(AssembleBilaplaceProblem_AD);

      // initilaize and solve the system
      system.init();
      system.solve();

      std::pair< double , double > norm = GetErrorNorm(&mlSol);
      l2Norm[i][j]  = norm.first;
      semiNorm[i][j] = norm.second;
      // print solutions
      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back("All");

      VTKWriter vtkIO(&mlSol);
      vtkIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, i);

//       GMVWriter gmvIO(&mlSol);
//       gmvIO.SetDebugOutput(true);
//       gmvIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);
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

void AssembleBilaplaceProblem_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use

  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("Poisson");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const unsigned levelMax = mlPdeSys->GetLevelMax();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

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

  unsigned solvIndex;
  solvIndex = mlSol->GetIndex("v");    // get the position of "v" in the ml_sol object
  unsigned solvType = mlSol->GetSolutionType(solvIndex);    // get the finite element type for "v"

  unsigned solvPdeIndex;
  solvPdeIndex = mlPdeSys->GetSolPdeIndex("v");    // get the position of "v" in the pdeSys object

  vector < adept::adouble >  solv; // local solution



  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector< int > KKDof; // local to global pdeSys dofs
  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  vector< double > Res; // local redidual vector
  vector< adept::adouble > aResu; // local redidual vector
  vector< adept::adouble > aResv; // local redidual vector


  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  solu.reserve(maxSize);
  solv.reserve(maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);

  KKDof.reserve(2 * maxSize);
  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  phi_xx.reserve(maxSize * dim2);

  Res.reserve(2 * maxSize);
  aResu.reserve(maxSize);
  aResv.reserve(maxSize);

  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)
  Jac.reserve(4 * maxSize * maxSize);
  vector< double > Jact; // local Jacobian matrix (ordered by raw, PETSC)
  Jact.reserve(4 * maxSize * maxSize);


  if (assembleMatrix)
    KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->IS_Mts2Gmt_elem_offset[iproc]; iel < msh->IS_Mts2Gmt_elem_offset[iproc + 1]; iel++) {

    unsigned kel = msh->IS_Mts2Gmt_elem[iel]; // mapping between paralell dof and mesh dof
    short unsigned kelGeom = el->GetElementType(kel);    // element geometry type
    unsigned nDofs  = el->GetElementDofNumber(kel, soluType);    // number of solution element dofs
    unsigned nDofs2 = el->GetElementDofNumber(kel, xType);    // number of coordinate element dofs

    // resize local arrays
    KKDof.resize(2 * nDofs);
    solu.resize(nDofs);
    solv.resize(nDofs);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofs2);
    }

    Res.resize(2 * nDofs);    //resize
    aResu.resize(nDofs);    //resize
    aResv.resize(nDofs);    //resize

    std::fill(aResu.begin(), aResu.end(), 0);    //set aRes to zero
    std::fill(aResv.begin(), aResv.end(), 0);    //set aRes to zero

    if (assembleMatrix) {   //resize
      Jact.resize(4 * nDofs * nDofs);
      Jac.resize(4 * nDofs * nDofs);
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned iNode = el->GetMeshDof(kel, i, soluType);    // local to global solution node
      unsigned solDof = msh->GetMetisDof(iNode, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      solv[i] = (*sol->_Sol[solvIndex])(solDof);      // global extraction and local storage for the solution
      KKDof[i]       = pdeSys->GetKKDof(soluIndex, soluPdeIndex, iNode);    // global to global mapping between solution node and pdeSys dof
      KKDof[nDofs + i] = pdeSys->GetKKDof(solvIndex, solvPdeIndex, iNode);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs2; i++) {
      unsigned iNode = el->GetMeshDof(kel, i, xType);    // local to global coordinates node
      unsigned xDof  = msh->GetMetisDof(iNode, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_coordinate->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    if (level == levelMax || !el->GetRefinedElementIndex(kel)) {      // do not care about this if now (it is used for the AMR)
      // start a new recording of all the operations involving adept::adouble variables
      s.new_recording();

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][soluType]->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[kelGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

        // evaluate the solution, the solution derivatives and the coordinates in the gauss point
        adept::adouble soluGauss = 0;
        vector < adept::adouble > soluGauss_x(dim, 0.);

        adept::adouble solvGauss = 0;
        vector < adept::adouble > solvGauss_x(dim, 0.);

        vector < double > xGauss(dim, 0.);

        for (unsigned i = 0; i < nDofs; i++) {
          soluGauss += phi[i] * solu[i];
          solvGauss += phi[i] * solv[i];

          for (unsigned jdim = 0; jdim < dim; jdim++) {
            soluGauss_x[jdim] += phi_x[i * dim + jdim] * solu[i];
            solvGauss_x[jdim] += phi_x[i * dim + jdim] * solv[i];
            xGauss[jdim] += x[jdim][i] * phi[i];
          }
        }

        // *** phi_i loop ***
        for (unsigned i = 0; i < nDofs; i++) {

          adept::adouble Laplace_u = 0.;
          adept::adouble Laplace_v = 0.;

          for (unsigned jdim = 0; jdim < dim; jdim++) {
            Laplace_u   +=  - phi_x[i * dim + jdim] * soluGauss_x[jdim];
            Laplace_v   +=  - phi_x[i * dim + jdim] * solvGauss_x[jdim];
          }

          double exactSolValue = GetExactSolutionValue(xGauss);

          double pi = acos(-1.);

          aResv[i] += (solvGauss * phi[i] -  Laplace_u) * weight;
          aResu[i] += (4.*pi * pi * pi * pi * exactSolValue * phi[i] -  Laplace_v) * weight;

        } // end phi_i loop
      } // end gauss point loop
    } // endif single element not refined or fine grid loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    for (int i = 0; i < nDofs; i++) {
      Res[i]       = aResu[i].value();
      Res[nDofs + i] = aResv[i].value();
    }

    RES->add_vector_blocked(Res, KKDof);

    if (assembleMatrix) {

      // define the dependent variables
      s.dependent(&aResu[0], nDofs);
      s.dependent(&aResv[0], nDofs);

      // define the independent variables
      s.independent(&solu[0], nDofs);
      s.independent(&solv[0], nDofs);
      // get the jacobian matrix (ordered by column)
      s.jacobian(&Jac[0]);

      // get the jacobian matrix (ordered by raw, i.e. Jact=Jac^t)
      for (int inode = 0; inode < 2 * nDofs; inode++) {
        for (int jnode = 0; jnode < 2 * nDofs; jnode++) {
          Jact[inode * 2 * nDofs + jnode] = -Jac[jnode * 2 * nDofs + inode];
        }
      }

      //store Jact in the global matrix KK
      KK->add_matrix_blocked(Jact, KKDof, KKDof);

      s.clear_independents();
      s.clear_dependents();
    }
  } //end element loop for each process

  RES->close();

  if (assembleMatrix) KK->close();

  // ***************** END ASSEMBLY *******************
}





std::pair < double, double > GetErrorNorm(MultiLevelSolution* mlSol) {
  unsigned level = mlSol->_ml_msh->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*          msh          = mlSol->_ml_msh->GetLevel(level);    // pointer to the mesh (level) object
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
  for (int iel = msh->IS_Mts2Gmt_elem_offset[iproc]; iel < msh->IS_Mts2Gmt_elem_offset[iproc + 1]; iel++) {

    unsigned kel = msh->IS_Mts2Gmt_elem[iel]; // mapping between paralell dof and mesh dof
    short unsigned kelGeom = el->GetElementType(kel);    // element geometry type
    unsigned nDofs  = el->GetElementDofNumber(kel, soluType);    // number of solution element dofs
    unsigned nDofs2 = el->GetElementDofNumber(kel, xType);    // number of coordinate element dofs

    // resize local arrays
    solu.resize(nDofs);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofs2);
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned iNode = el->GetMeshDof(kel, i, soluType);    // local to global solution node
      unsigned solDof = msh->GetMetisDof(iNode, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs2; i++) {
      unsigned iNode = el->GetMeshDof(kel, i, xType);    // local to global coordinates node
      unsigned xDof  = msh->GetMetisDof(iNode, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_coordinate->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[kelGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

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
