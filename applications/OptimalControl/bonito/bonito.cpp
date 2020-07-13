/** tutorial/Ex2
 * This example shows how to set and solve the weak form of the Poisson problem
 *                    $$ \Delta u = \Delta u_exact \text{ on }\Omega, $$
 *                    $$ u=0 \text{ on } \Gamma, $$
 * on a square domain $\Omega$ with boundary $\Gamma$;
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
#include "LinearImplicitSystem.hpp"
#include "adept.h"


using namespace femus;

#define N_UNIFORM_LEVELS  5
#define N_ERASED_LEVELS   4
#define S_FRAC            0.25

#define q_step            1.
// #define N              10

#define EX_1              -1.
#define EX_2               1.
#define EY_1              -1.
#define EY_2               1.


bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0;

  return dirichlet;
}

void BuildU(MultiLevelSolution& mlSol);

void AssemblePoissonProblem(MultiLevelProblem& ml_prob);

int n_sys;

int main(int argc, char** args) {
  
  const std::string fe_quad_rule_1 = "seventh";

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
  unsigned numberOfUniformLevels = N_UNIFORM_LEVELS;

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;
//   mlMsh.ReadCoarseMesh("./input/square_quad.neu", fe_quad_rule_1, scalingFactor);
//   //mlMsh.ReadCoarseMesh("./input/cube_tet.neu", "seventh", scalingFactor);
//   /* "seventh" is the order of accuracy that is used in the gauss integration scheme
//     probably in furure it is not going to be an argument of this function   */
  
//   mlMsh.GenerateCoarseBoxMesh(2, 0, 0, EX_1, EX_2, 0., 0., 0., 0., EDGE3, fe_quad_rule_1.c_str());
//   mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);
  
  mlMsh.GenerateCoarseBoxMesh(2, 2, 0, EX_1, EX_2, EY_1, EY_2, 0., 0., QUAD9, fe_quad_rule_1.c_str());
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);
  
  unsigned dim = mlMsh.GetDimension();

  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("u", LAGRANGE, SECOND);
  mlSol.Initialize("All");
  
  int N_plus = round( pow( M_PI, 2 ) / ( 4 * (1 - S_FRAC) * pow( q_step, 2 ) ) );
  int N_minus = round( pow( M_PI, 2 ) / ( 4  * S_FRAC * pow( q_step, 2 )) ) ;

  for (int i = - N_minus; i < N_plus + 1; i++) {
//   for (int i = - N; i < N + 1; i++) {
    char solName[10];
    sprintf (solName, "w%d", i);
    mlSol.AddSolution (solName, LAGRANGE, SECOND);
  }
 
  // attach the boundary condition function and generate boundary data
  
  mlSol.Initialize ("All");
  mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
  mlSol.GenerateBdc ("All");
  
  

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  
  for (int i = - N_minus; i < N_plus + 1; i++) {
//   for (int i = - N; i < N + 1; i++) {
    char sysName[20];
    sprintf (sysName, "Poisson%d", i);
    LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > (sysName);
    char solName[20];
    sprintf (solName, "w%d", i);
    system.AddSolutionToSystemPDE(solName);
    system.SetAssembleFunction(AssemblePoissonProblem);
     
    system.init();
    n_sys = i;
    system.SetOuterSolver(PREONLY);
    system.MGsolve();
  }

  BuildU(mlSol);
  
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
 
  //vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);
  return 0;
}


void AssemblePoissonProblem(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled
  // call the adept stack object

  adept::Stack& s = FemusInit::_adeptStack;
  //  extract pointers to the several objects that we are going to use

  char sysName[20];
  sprintf (sysName, "Poisson%d", n_sys);
  char solName[10];
  sprintf (solName, "w%d", n_sys);  
  
  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> (sysName);   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex(solName);    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex(solName);    // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solu; // local solution
  solu.reserve(maxSize);

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  phi_xx.reserve(maxSize * dim2);

  vector< adept::adouble > aRes; // local redidual vector
  aRes.reserve(maxSize);

  vector< int > l2GMap; // local to global mapping
  l2GMap.reserve(maxSize);
  vector< double > Res; // local redidual vector
  Res.reserve(maxSize);
  vector < double > Jac;
  Jac.reserve(maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
     
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    l2GMap.resize(nDofu);
    solu.resize(nDofu);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    aRes.resize(nDofu);    //resize
    std::fill(aRes.begin(), aRes.end(), 0);    //set aRes to zero

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      l2GMap[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
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
      adept::adouble solu_gss = 0;
      vector < adept::adouble > gradSolu_gss(dim, 0.);
      vector < double > x_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i]; // We dont use this one for this problem.

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
          x_gss[jdim] += x[jdim][i] * phi[i];
        }
      }
      
//       double q_step = 1 / ( sqrt( N ) );

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofu; i++) {

        adept::adouble laplace = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          laplace   +=  phi_x[i * dim + jdim] * gradSolu_gss[jdim];
          x_gss[jdim] += x[jdim][i] * phi[i];
        }

        aRes[i] += ( + exp( 2 * S_FRAC * q_step * n_sys ) /** sin(2 * acos(0.0) * x[0][i]) * sin(2 * acos(0.0) * x[1][i]) */* phi[i] - solu_gss * phi[i] - exp( 2 * q_step * n_sys ) * laplace) * weight ;

      } // end phi_i loop
      
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    Res.resize(nDofu);    //resize

    for (int i = 0; i < nDofu; i++) {
      Res[i] = - aRes[i].value();
    }

    RES->add_vector_blocked(Res, l2GMap);

    // define the dependent variables
    s.dependent(&aRes[0], nDofu);

    // define the independent variables
    s.independent(&solu[0], nDofu);

    // get the jacobian matrix (ordered by row major )
    Jac.resize(nDofu * nDofu);    //resize
    s.jacobian(&Jac[0], true);

    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac, l2GMap, l2GMap);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}

void BuildU(MultiLevelSolution& mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel (level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel (level);
  unsigned iproc  = msh->processor_id();
  
//   double q_step = 1. / sqrt( N );
  
  int N_plus = round( pow( M_PI, 2 ) / ( 4. * (1 - S_FRAC) * pow( q_step, 2 ) ) );
  int N_minus = round( pow( M_PI, 2 ) / ( 4. * S_FRAC * pow( q_step, 2 )) ) ;

  std::vector< unsigned > wIndex(N_minus + N_plus + 1);
//   std::vector< unsigned > wIndex(2 * N + 1);
  
  for (int i = - N_minus; i < N_plus + 1; i++) {
//   for (int i = - N; i < N + 1; i++) {
    char solName[10];
    sprintf (solName, "w%d", i);
    wIndex[i + N_minus] = mlSol.GetIndex (solName);
//     wIndex[i + N] = mlSol.GetIndex (solName);
  }
  unsigned uIndex = mlSol.GetIndex ("u");
  
  unsigned solType = mlSol.GetSolutionType (uIndex);

  sol->_Sol[uIndex]->zero();
  
  for(unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    double value = 0.;
    
    double Cs =  2 * sin( M_PI * S_FRAC) / M_PI ;
    
//     for (unsigned k = 0; k < N_minus + N_plus + 1; k++) {
// //     for (int k = 0; k < 2 * N + 1; k++) {
//       double weight = exp( 2 * S_FRAC * q_step * ( k - N_minus) );
// //       double weight = exp( /*2 **/ S_FRAC * q_step * ( k - N) );
//       
//       value += Cs * q_step * weight * (*sol->_Sol[wIndex[k]])(i);
//       std::cout.precision(14);
//       if(i == 32) std::cout<<k-N_minus<< " " << Cs * q_step * weight << "  " <<
//         (*sol->_Sol[wIndex[k]])(i) << "  " <<
//         Cs * q_step * weight * (*sol->_Sol[wIndex[k]])(i) << "  " << value <<"\n";
//       
// //       if(k == 0 || k == 1) std::cout<< k<< "  " << i<< "  " <<  (*sol->_Sol[wIndex[k]])(i) << "\n" ;
//       
//     }
        for (int j = -N_minus ; j < N_plus + 1; j++) {
//     for (int k = 0; k < 2 * N + 1; k++) {
      double weight = 1.;//exp( 2 * S_FRAC * q_step * j );
//       double weight = exp( /*2 **/ S_FRAC * q_step * ( k - N) );
      
      value += Cs * q_step * weight * (*sol->_Sol[wIndex[j+N_minus]])(i);
      std::cout.precision(14);
//       if(i == 32) std::cout<< j << " " << Cs * q_step * weight << "  " <<
//         (*sol->_Sol[wIndex[j+N_minus]])(i) << "  " <<
//         Cs * q_step * weight * (*sol->_Sol[wIndex[j+N_minus]])(i) << "  " << value <<"\n";
      
//       if(k == 0 || k == 1) std::cout<< k<< "  " << i<< "  " <<  (*sol->_Sol[wIndex[k]])(i) << "\n" ;
      
    }
    sol->_Sol[uIndex]->set(i,value);
  }
  
  sol->_Sol[uIndex]->close();

}







