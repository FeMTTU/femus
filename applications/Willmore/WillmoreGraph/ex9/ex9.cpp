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


using namespace femus;

//int simulation = 1; // =1 sphere (default) = 2 torus

unsigned P = 2;

//Sphere

double thetaSphere = acos (-1.) / 4;

bool SetBoundaryConditionSphere (const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet

  if (!strcmp ("u", SolName)) {
      value = 0;

    // if(x[1] >= 0){
    //   value = 0.1 * sin( 4. * acos(x[0]) );
    // }
    // else{
    //   value = 0.1 * sin( 4. * (2. * M_PI - acos(x[0])) );
    // }
  }
  else if (!strcmp ("H", SolName)) {
    //value = -1. / tan(thetaSphere);
    value = 0.;
  }
  else if (!strcmp ("W", SolName)) {
    //value = -1. / tan(thetaSphere);
    value =  0.;
  }
  return dirichlet;
}


double InitalValueHSphere (const std::vector < double >& x) {
  double r = sqrt( x[0] * x[0] + x[1] * x[1]);
  return sin(0.5 * M_PI * (1.-r));
}

double GetTimeStep (const double time) {
  return 0.005;
}


void AssembleWillmoreProblem_AD (MultiLevelProblem& ml_prob);

void AssembleInitialPWillmore (MultiLevelProblem& ml_prob);

void AssembleInitialM1 (MultiLevelProblem& ml_prob);

int main (int argc, char** args) {

  FemusInit mpinit (argc, args, MPI_COMM_WORLD);

  std::ostringstream filename;
  filename << "./input/circle_quad4.neu";


  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/circle_quad.neu","seventh", scalingFactor);
  mlMsh.ReadCoarseMesh (filename.str().c_str(), "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
   *    probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh (numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // print mesh info
  mlMsh.PrintInfo();
  MultiLevelSolution mlSol (&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution ("u", LAGRANGE, SECOND, 2);
  mlSol.AddSolution ("H", LAGRANGE, SECOND);
  mlSol.AddSolution ("W", LAGRANGE, SECOND);

  mlSol.Initialize ("All");
  mlSol.Initialize ("H", InitalValueHSphere);

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction (SetBoundaryConditionSphere);

  mlSol.GenerateBdc ("u");
  mlSol.GenerateBdc ("H");
  mlSol.GenerateBdc ("W");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb (&mlSol);

  LinearImplicitSystem& systemM1 = mlProb.add_system < LinearImplicitSystem > ("InitialM1");

  // add solution "u" to system
  systemM1.AddSolutionToSystemPDE ("u");

  // attach the assembling function to system
  systemM1.SetAssembleFunction (AssembleInitialM1);


  // initilaize and solve the system
  systemM1.init();

  systemM1.MGsolve();


  NonLinearImplicitSystem& system0 = mlProb.add_system < NonLinearImplicitSystem > ("Initial");

  // add solution "u" to system
  system0.AddSolutionToSystemPDE ("u");
  system0.AddSolutionToSystemPDE ("W");

  // attach the assembling function to system
  system0.SetAssembleFunction (AssembleInitialPWillmore);


  // initilaize and solve the system
  system0.init();

  system0.MGsolve();

  mlSol.SetWriter (VTK);

  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back ("All");


  mlSol.GetWriter()->SetGraphVariable("u");
  mlSol.GetWriter()->SetDebugOutput (true);
  mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);

  // add system Wilmore in mlProb as a Linear Implicit System
  TransientNonlinearImplicitSystem& system = mlProb.add_system < TransientNonlinearImplicitSystem > ("Willmore");

  // add solution "u" to system
  system.AddSolutionToSystemPDE ("u");
  system.AddSolutionToSystemPDE ("H");
  system.AddSolutionToSystemPDE ("W");

  // attach the assembling function to system
  system.SetAssembleFunction (AssembleWillmoreProblem_AD);
  system.AttachGetTimeIntervalFunction (GetTimeStep);

  // initilaize and solve the system
  system.init();


  unsigned numberOfTimeSteps = 300;
  for (unsigned time_step = 0; time_step < numberOfTimeSteps; time_step++) {

    system.CopySolutionToOldSolution();
    system.MGsolve();
    if ( (time_step + 1) % 1 == 0)
      mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, (time_step + 1) );
  }

  return 0;
}


void AssembleWillmoreProblem_AD (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  TransientNonlinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientNonlinearImplicitSystem> ("Willmore");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh        = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex ("u");   // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType (soluIndex);   // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex ("u");   // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solu; // local solution
  vector < double >  soluOld; // local solution

  unsigned solHIndex;
  solHIndex = mlSol->GetIndex ("H");   // get the position of "v" in the ml_sol object
  unsigned solHType = mlSol->GetSolutionType (solHIndex);   // get the finite element type for "v"

  unsigned solHPdeIndex;
  solHPdeIndex = mlPdeSys->GetSolPdeIndex ("H");   // get the position of "v" in the pdeSys object


  vector < adept::adouble >  solH; // local solution


  unsigned solWIndex;
  solWIndex = mlSol->GetIndex ("W");   // get the position of "u" in the ml_sol object
  unsigned solWType = mlSol->GetSolutionType (solWIndex);   // get the finite element type for "u"

  unsigned solWPdeIndex;
  solWPdeIndex = mlPdeSys->GetSolPdeIndex ("W");   // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solW; // local solution


  vector < vector < double > > x (dim);   // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE BIQUADRATIC)

  vector< int > sysDof; // local to global pdeSys dofs
  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives

  double weight; // gauss point weight

  vector< double > Res; // local redidual vector
  vector< adept::adouble > aResu; // local redidual vector
  vector< adept::adouble > aResH; // local redidual vector
  vector< adept::adouble > aResW; // local redidual vector


  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27
  solu.reserve (maxSize);
  soluOld.reserve (maxSize);
  solH.reserve (maxSize);
  solW.reserve (maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve (maxSize);

  sysDof.reserve (3 * maxSize);
  phi.reserve (maxSize);
  phi_x.reserve (maxSize * dim);
  unsigned dim2 = (3 * (dim - 1) + ! (dim - 1));       // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)


  Res.reserve (3 * maxSize);
  aResu.reserve (maxSize);
  aResH.reserve (maxSize);
  aResW.reserve (maxSize);

  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)
  Jac.reserve (9 * maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  double dt = GetTimeStep (0);

  double HPintegralLocal = 0.;

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // element geometry type
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofs  = msh->GetElementDofNumber (iel, soluType);   // number of solution element dofs
    unsigned nDofs2 = msh->GetElementDofNumber (iel, xType);   // number of coordinate element dofs
    // resize local arrays
    sysDof.resize (3 * nDofs);
    solu.resize (nDofs);
    soluOld.resize(nDofs);
    solH.resize (nDofs);
    solW.resize (nDofs);

    for (int i = 0; i < dim; i++) {
      x[i].resize (nDofs2);
    }

    Res.resize (3 * nDofs); //resize
    aResu.assign (nDofs, 0.);  //resize
    aResH.assign (nDofs, 0.);   //resize
    aResW.assign (nDofs, 0.);  //resize

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, soluType);   // local to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex]) (solDof);     // global extraction and local storage for the solution
      soluOld[i] = (*sol->_SolOld[soluIndex]) (solDof);     // global extraction and local storage for the solution
      solH[i] = (*sol->_Sol[solHIndex]) (solDof);     // global extraction and local storage for the solution
      solW[i] = (*sol->_Sol[solWIndex]) (solDof);     // global extraction and local storage for the solution
      sysDof[i]       = pdeSys->GetSystemDof (soluIndex, soluPdeIndex, i, iel);   // local to global mapping between solution node and pdeSys dof
      sysDof[nDofs + i] = pdeSys->GetSystemDof (solHIndex, solHPdeIndex, i, iel);   // local to global mapping between solution node and pdeSys dof
      sysDof[2 * nDofs + i] = pdeSys->GetSystemDof (solWIndex, solWPdeIndex, i, iel);   // local to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs2; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // local to global mapping between coordinates node and coordinate dof

      for (unsigned idim = 0; idim < dim; idim++) {
        x[idim][i] = (*msh->_topology->_Sol[idim]) (xDof);     // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian (x, ig, weight, phi, phi_x);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      //adept::adouble soluGauss = 0;

      vector < adept::adouble > u_x (dim, 0.);

      adept::adouble W = 0;
      vector < adept::adouble > W_x (dim, 0.);

      adept::adouble H = 0;

      adept::adouble U = 0;
      double UOld = 0;

      for (unsigned i = 0; i < nDofs; i++) {
        U += phi[i] * solu[i];
        UOld += phi[i] * soluOld[i];

        H += phi[i] * solH[i];
        W += phi[i] * solW[i];
        for (unsigned idim = 0; idim < dim; idim++) {
          u_x[idim] += phi_x[i * dim + idim] * solu[i];
          W_x[idim] += phi_x[i * dim + idim] * solW[i];
        }
      }

      adept::adouble HPm1 = H;
      for (unsigned i = 0; i < P - 2; i++) {
        HPm1 *= H;
      }

      adept::adouble u_xNorm2 = 0;
      for (unsigned idim = 0; idim < dim; idim++) {
        u_xNorm2 += u_x[idim] * u_x[idim];
      }

      adept::adouble A = sqrt (1. + u_xNorm2);
      adept::adouble A2 = A * A;

      HPintegralLocal += pow( H.value(), P) * A.value() * weight;

      double Id[2][2] = {{1., 0.}, {0., 1.}};
      vector < vector < adept::adouble> > B (dim);

      for (unsigned idim = 0; idim < dim; idim++) {
        B[idim].resize (dim);
        for (unsigned jdim = 0; jdim < dim; jdim++) {
          B[idim][jdim] = Id[idim][jdim] - (u_x[idim] * u_x[jdim]) / A2;
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {

        adept::adouble nonLinearLaplaceU = 0.;
        adept::adouble nonLinearLaplaceW = 0.;


        for (unsigned idim = 0; idim < dim; idim++) {

          nonLinearLaplaceU +=  - 1. / A  * u_x[idim] * phi_x[i * dim + idim];

          nonLinearLaplaceW +=  - (P / (2.* A) * (B[idim][0] * W_x[0] + B[idim][1] * W_x[1])
                                    - (W * H) / A2 * u_x[idim]
                                   ) * phi_x[i * dim + idim];

        }

        aResu[i] += (2. * H * phi[i] - nonLinearLaplaceU) * weight;
        aResH[i] += (W - A * HPm1) * phi[i] * weight;
        aResW[i] += ( (U-UOld)/ (A * dt) * phi[i] + nonLinearLaplaceW ) * weight;

      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    for (int i = 0; i < nDofs; i++) {
      Res[i]       = -aResu[i].value();
      Res[nDofs + i] = -aResH[i].value();
      Res[2 * nDofs + i] = -aResW[i].value();
    }

    RES->add_vector_blocked (Res, sysDof);

    Jac.resize ( (3 * nDofs) * (3 * nDofs));
    // define the dependent variables
    s.dependent (&aResu[0], nDofs);
    s.dependent (&aResH[0], nDofs);
    s.dependent (&aResW[0], nDofs);

    // define the independent variables
    s.independent (&solu[0], nDofs);
    s.independent (&solH[0], nDofs);
    s.independent (&solW[0], nDofs);
    // get the jacobian matrix (ordered by row)
    s.jacobian (&Jac[0], true);

    KK->add_matrix_blocked (Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();
  } //end element loop for each process

  RES->close();
  KK->close();

  double HPIntegral = 0.;
  MPI_Allreduce ( &HPintegralLocal, &HPIntegral, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  std::cout.precision(14);
  std::cout << "int_S H^p dS = " << HPIntegral << std::endl;


  // ***************** END ASSEMBLY *******************
}


void AssembleInitialPWillmore (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("Initial");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh        = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex ("u");   // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType (soluIndex);   // get the finite element type for "u"
  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex ("u");   // get the position of "u" in the pdeSys object
  vector < adept::adouble >  solu; // local solution


  unsigned solHIndex;
  solHIndex = mlSol->GetIndex ("H");   // get the position of "v" in the ml_sol object
  unsigned solHType = mlSol->GetSolutionType (solHIndex);   // get the finite element type for "v"
  vector < double >  solH; // local solution


  unsigned solWIndex;
  solWIndex = mlSol->GetIndex ("W");   // get the position of "u" in the ml_sol object
  unsigned solWType = mlSol->GetSolutionType (solWIndex);   // get the finite element type for "u"
  unsigned solWPdeIndex;
  solWPdeIndex = mlPdeSys->GetSolPdeIndex ("W");   // get the position of "u" in the pdeSys object
  vector < adept::adouble >  solW; // local solution


  vector < vector < double > > x (dim);   // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE BIQUADRATIC)

  vector< int > sysDof; // local to global pdeSys dofs
  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives

  double weight; // gauss point weight

  vector< double > Res; // local redidual vector
  vector< adept::adouble > aResu; // local redidual vector
  vector< adept::adouble > aResW; // local redidual vector


  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27
  solu.reserve (maxSize);
  solH.reserve (maxSize);
  solW.reserve (maxSize);

  for (unsigned i = 0; i < dim; i++){
    x[i].reserve (maxSize);
  }

  sysDof.reserve (2 * maxSize);
  phi.reserve (maxSize);
  phi_x.reserve (maxSize * dim);

  Res.reserve (2 * maxSize);
  aResu.reserve (maxSize);
  aResW.reserve (maxSize);

  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)
  Jac.reserve (4 * maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // element geometry type
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofs  = msh->GetElementDofNumber (iel, soluType);   // number of solution element dofs
    unsigned nDofs2 = msh->GetElementDofNumber (iel, xType);   // number of coordinate element dofs
    // resize local arrays
    sysDof.resize (2 * nDofs);
    solu.resize (nDofs);
    solH.resize (nDofs);
    solW.resize (nDofs);

    for (int i = 0; i < dim; i++) {
      x[i].resize (nDofs2);
    }

    Res.resize (2 * nDofs); //resize
    aResu.assign (nDofs, 0.);  //resize
    aResW.assign (nDofs, 0.);  //resize

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, soluType);   // local to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex]) (solDof);     // global extraction and local storage for the solution
      solH[i] = (*sol->_Sol[solHIndex]) (solDof);     // global extraction and local storage for the solution
      solW[i] = (*sol->_Sol[solWIndex]) (solDof);     // global extraction and local storage for the solution
      sysDof[i]       = pdeSys->GetSystemDof (soluIndex, soluPdeIndex, i, iel);   // local to global mapping between solution node and pdeSys dof
      sysDof[ nDofs + i] = pdeSys->GetSystemDof (solWIndex, solWPdeIndex, i, iel);   // local to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs2; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // local to global mapping between coordinates node and coordinate dof

      for (unsigned idim = 0; idim < dim; idim++) {
        x[idim][i] = (*msh->_topology->_Sol[idim]) (xDof);     // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian (x, ig, weight, phi, phi_x);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      //adept::adouble soluGauss = 0;

      vector < adept::adouble > u_x (dim, 0.);

      adept::adouble W = 0;

      double H = 0;

      for (unsigned i = 0; i < nDofs; i++) {
        H += phi[i] * solH[i];
        W += phi[i] * solW[i];
        for (unsigned idim = 0; idim < dim; idim++) {
          u_x[idim] += phi_x[i * dim + idim] * solu[i];
        }
      }

      double HPm1 = H;
      for (unsigned i = 0; i < P - 2; i++) {
        HPm1 *= H;
      }

      adept::adouble u_xNorm2 = 0;
      for (unsigned idim = 0; idim < dim; idim++) {
        u_xNorm2 += u_x[idim] * u_x[idim];
      }

      adept::adouble A = sqrt (1. + u_xNorm2);
      adept::adouble A2 = A * A;

      double Id[2][2] = {{1., 0.}, {0., 1.}};
      vector < vector < adept::adouble> > B (dim);

      for (unsigned idim = 0; idim < dim; idim++) {
        B[idim].resize (dim);
        for (unsigned jdim = 0; jdim < dim; jdim++) {
          B[idim][jdim] = Id[idim][jdim] - (u_x[idim] * u_x[jdim]) / A2;
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {

        adept::adouble nonLinearLaplaceU = 0.;



        for (unsigned idim = 0; idim < dim; idim++) {

          nonLinearLaplaceU +=  - 1. / A  * u_x[idim] * phi_x[i * dim + idim];

        }

        aResu[i] += (2. * H * phi[i] - nonLinearLaplaceU) * weight;
        aResW[i] += (W - A * HPm1) * phi[i] * weight;


      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    for (int i = 0; i < nDofs; i++) {
      Res[i]       = -aResu[i].value();
      Res[nDofs + i] = -aResW[i].value();
    }

    RES->add_vector_blocked (Res, sysDof);

    Jac.resize ( (2 * nDofs) * (2 * nDofs));
    // define the dependent variables
    s.dependent (&aResu[0], nDofs);
    s.dependent (&aResW[0], nDofs);

    // define the independent variables
    s.independent (&solu[0], nDofs);
    s.independent (&solW[0], nDofs);
    // get the jacobian matrix (ordered by row)
    s.jacobian (&Jac[0], true);

    KK->add_matrix_blocked (Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();
  } //end element loop for each process

  RES->close();
  KK->close();


  // ***************** END ASSEMBLY *******************
}


void AssembleInitialM1 (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("InitialM1");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh        = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex ("u");   // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType (soluIndex);   // get the finite element type for "u"
  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex ("u");   // get the position of "u" in the pdeSys object
  vector < adept::adouble >  solu; // local solution

  vector < vector < double > > x (dim);   // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE BIQUADRATIC)

  vector< int > sysDof; // local to global pdeSys dofs
  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives

  double weight; // gauss point weight

  vector< double > Res; // local redidual vector
  vector< adept::adouble > aResu; // local redidual vector


  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27
  solu.reserve (maxSize);


  for (unsigned i = 0; i < dim; i++){
    x[i].reserve (maxSize);
  }

  sysDof.reserve (maxSize);
  phi.reserve (maxSize);
  phi_x.reserve (maxSize * dim);

  Res.reserve ( maxSize);
  aResu.reserve (maxSize);

  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)
  Jac.reserve (maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // element geometry type
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofs  = msh->GetElementDofNumber (iel, soluType);   // number of solution element dofs
    unsigned nDofs2 = msh->GetElementDofNumber (iel, xType);   // number of coordinate element dofs
    // resize local arrays
    sysDof.resize (nDofs);
    solu.resize (nDofs);

    for (int i = 0; i < dim; i++) {
      x[i].resize (nDofs2);
    }

    Res.resize (nDofs); //resize
    aResu.assign (nDofs, 0.);  //resize

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, soluType);   // local to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex]) (solDof);     // global extraction and local storage for the solution
      sysDof[i]       = pdeSys->GetSystemDof (soluIndex, soluPdeIndex, i, iel);   // local to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs2; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // local to global mapping between coordinates node and coordinate dof

      for (unsigned idim = 0; idim < dim; idim++) {
        x[idim][i] = (*msh->_topology->_Sol[idim]) (xDof);     // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian (x, ig, weight, phi, phi_x);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      //adept::adouble soluGauss = 0;

      vector < adept::adouble > u_x (dim, 0.);

      for (unsigned i = 0; i < nDofs; i++) {
        for (unsigned idim = 0; idim < dim; idim++) {
          u_x[idim] += phi_x[i * dim + idim] * solu[i];
        }
      }




      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {

        adept::adouble LinearLaplaceU = 0.;



        for (unsigned idim = 0; idim < dim; idim++) {

          LinearLaplaceU +=  - u_x[idim] * phi_x[i * dim + idim];

        }

        aResu[i] += (LinearLaplaceU) * weight;



      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    for (int i = 0; i < nDofs; i++) {
      Res[i]       = -aResu[i].value();
    }

    RES->add_vector_blocked (Res, sysDof);

    Jac.resize ( nDofs * nDofs);
    // define the dependent variables
    s.dependent (&aResu[0], nDofs);


    // define the independent variables
    s.independent (&solu[0], nDofs);

    // get the jacobian matrix (ordered by row)
    s.jacobian (&Jac[0], true);

    KK->add_matrix_blocked (Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();
  } //end element loop for each process

  RES->close();
  KK->close();


  // ***************** END ASSEMBLY *******************
}
