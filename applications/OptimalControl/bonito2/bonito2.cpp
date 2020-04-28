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

#include "PetscMatrix.hpp"



using namespace femus;

#define N_UNIFORM_LEVELS  6
#define N_ERASED_LEVELS   5
#define S_FRAC            0.5

#define q_step            .5
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

  mlMsh.GenerateCoarseBoxMesh(2, 0, 0, EX_1, EX_2, 0., 0., 0., 0., EDGE3, fe_quad_rule_1.c_str());
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);

//   mlMsh.GenerateCoarseBoxMesh(2, 2, 0, EX_1, EX_2, EY_1, EY_2, 0., 0., QUAD9, fe_quad_rule_1.c_str());
//   mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);

  unsigned dim = mlMsh.GetDimension();

  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);



  int N_plus = round(pow(M_PI, 2) / (4 * (1 - S_FRAC) * pow(q_step, 2)));
  int N_minus = round(pow(M_PI, 2) / (4  * S_FRAC * pow(q_step, 2))) ;

  for(int i = - N_minus; i < N_plus + 1; i++) {
//   for (int i = - N; i < N + 1; i++) {
    char solName[10];
    sprintf(solName, "w%d", i);
    mlSol.AddSolution(solName, LAGRANGE, SECOND);
  }

  // add variables to mlSol
  mlSol.AddSolution("u", LAGRANGE, SECOND);
  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data

  mlSol.Initialize("All");
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");



  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System

  LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("Poisson");
  char solName[20];
  for(int i = - N_minus; i < N_plus + 1; i++) {
    sprintf(solName, "w%d", i);
    system.AddSolutionToSystemPDE(solName);
  }
  system.AddSolutionToSystemPDE("u");
  system.SetAssembleFunction(AssemblePoissonProblem);


  unsigned N = N_minus + N_plus + 1;
  FieldSplitTree **field1;
  field1 = new FieldSplitTree * [N + 1];

  std::vector < FieldSplitTree *> split1;
  split1.reserve(N + 1);

  for(int i = - N_minus; i < N_plus + 1; i++) {
    sprintf(solName, "w%d", i);

    std::vector < unsigned > fieldw(1);
    fieldw[0] = system.GetSolPdeIndex(solName);

    std::vector < unsigned > solutionTypew(1);
    solutionTypew[0] = mlSol.GetSolutionType(solName);

    field1[i + N_minus] = new FieldSplitTree(PREONLY, MLU_PRECOND, fieldw, solutionTypew, solName);

    split1.push_back(field1[i + N_minus]);
  }
  
  {
    std::vector < unsigned > fieldu(1);
    fieldu[0] = system.GetSolPdeIndex("u");

    std::vector < unsigned > solutionTypeu(1);
    solutionTypeu[0] = mlSol.GetSolutionType("u");

    field1[N] = new FieldSplitTree(PREONLY, MLU_PRECOND, fieldu, solutionTypeu, "u");

    split1.push_back(field1[N]);
  }
  
  
//   FieldSplitTree FieldSplit1(PREONLY, FIELDSPLIT_ADDITIVE_PRECOND, split1, "FS1");
// 
//   FieldSplit1.PrintFieldSplitTree();
// 
// 
//   FieldSplitTree **field2;
//   field2 = new FieldSplitTree * [1];
// 
//   std::vector < FieldSplitTree *> split2;
//   split2.reserve(1);
// 
//   {
//     std::vector < unsigned > fieldu(1);
//     fieldu[0] = system.GetSolPdeIndex("u");
// 
//     std::vector < unsigned > solutionTypeu(1);
//     solutionTypeu[0] = mlSol.GetSolutionType("u");
// 
//     field2[0] = new FieldSplitTree(PREONLY, MLU_PRECOND, fieldu, solutionTypeu, "u");
// 
//     split2.push_back(field2[0]);
//   }
//   FieldSplitTree FieldSplit2(PREONLY, FIELDSPLIT_ADDITIVE_PRECOND, split2, "FS2");
// 
//   FieldSplit2.PrintFieldSplitTree();
//   
//  
//   std::vector < FieldSplitTree *> FSAll(2);
//   FSAll.reserve(2);
//   FSAll.push_back(&FieldSplit1);
//   FSAll.push_back(&FieldSplit2);


  FieldSplitTree uSG(PREONLY, FIELDSPLIT_MULTIPLICATIVE_PRECOND, split1, "uSG");

  uSG.PrintFieldSplitTree();

  //return 1;

  system.SetOuterSolver(RICHARDSON);
  system.SetLinearEquationSolverType(FEMuS_FIELDSPLIT, INCLUDE_COARSE_LEVEL_TRUE);

  system.init();

  system.SetFieldSplitTree(&uSG);

  system.MGsolve();


  //BuildU(mlSol);

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


  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Poisson");   // pointer to the linear implicit system named "Poisson"
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

  int N_plus = round(pow(M_PI, 2) / (4 * (1 - S_FRAC) * pow(q_step, 2)));
  int N_minus = round(pow(M_PI, 2) / (4  * S_FRAC * pow(q_step, 2))) ;

  unsigned N = N_plus + N_minus + 1;

  std::vector<unsigned> solwIndex(N);
  std::vector<unsigned> solwPdeIndex(N);
  for(int i = - N_minus; i < N_plus + 1; i++) {
    char solName[10];
    sprintf(solName, "w%d", i);
    solwIndex[i + N_minus] = mlSol->GetIndex(solName);
    solwPdeIndex[i + N_minus] = mlPdeSys->GetSolPdeIndex(solName);
    std::cout << i + N_minus << " " << solName << " " << solwIndex[i + N_minus] << " " <<  solwPdeIndex[i + N_minus] << std::endl;
  }

  //solution variable
  unsigned soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluPdeIndex = mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object

  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  std::vector < std::vector < adept::adouble > >  solw(N); // local solution
  std::vector < adept::adouble >  solu; // local solution

  std::vector < std::vector < adept::adouble > >  aResw(N); // local solution
  std::vector< adept::adouble > aResu; // local redidual vector

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight


  vector< unsigned > l2GMap; // local to global mapping
  vector< double > Res; // local redidual vector
  vector< double > Jac;

  KK->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofAll  = (N + 1) * nDofu;

    // resize local arrays
    l2GMap.resize(nDofAll);

    for(unsigned i = 0; i < N; i++) {
      solw[i].resize(nDofu);
      aResw[i].assign(nDofu, 0.);
    }
    solu.resize(nDofu);
    aResu.assign(nDofu, 0.);

    for(int i = 0; i < dim; i++) {
      x[i].resize(nDofu);
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      for(unsigned j = 0; j < N; j++) {
        solw[j][i] = (*sol->_Sol[solwIndex[j]])(solDof);      // global extraction and local storage for the solution
        l2GMap[j * nDofu + i] = pdeSys->GetSystemDof(solwIndex[j], solwPdeIndex[j], i, iel);    // global to global mapping between solution node and pdeSys dof
      }

      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      l2GMap[N * nDofu + i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofu; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
      }
    }


    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      adept::adouble solu_gss = 0;
      std::vector < adept::adouble > solw_gss(N, 0);
      std::vector < std::vector < adept::adouble > > gradSolw_gss(N);

      for(unsigned j = 0; j < N; j++) gradSolw_gss[j].assign(dim, 0.);

      for(unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i]; // We dont use this one for this problem.
        for(unsigned j = 0; j < N; j++) {
          solw_gss[j] += phi[i] * solw[j][i];
          for(unsigned k = 0; k < dim; k++) {
            gradSolw_gss[j][k] += phi_x[i * dim + k] * solw[j][i];
          }
        }
      }
      double Cs =  2 * sin(M_PI * S_FRAC) / M_PI ;
      // *** phi_i loop ***
      for(unsigned i = 0; i < nDofu; i++) {
        aResu[i] -= solu[i] * phi[i] * weight ;
        for(int j = 0; j < N; j++) {
          adept::adouble laplace = 0.;
          for(unsigned k = 0; k < dim; k++) {
            laplace   +=  phi_x[i * dim + k] * gradSolw_gss[j][k];
          }
          aResw[j][i] += (+ exp(2 * S_FRAC * q_step * (j - N_minus)) * phi[i]
                          - solw_gss[j] * phi[i]
                          - exp(2 * q_step * (j - N_minus)) *   laplace
                         ) * weight ;
          aResu[i] += Cs * q_step * solw[j][i] * phi[i] * weight ;
        }
      } // end phi_i loop

    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    Res.resize(nDofAll);    //resize

    for(unsigned j = 0; j < N; j++) {
      for(unsigned i = 0; i < nDofu; i++) {
        Res[j * nDofu + i] = -aResw[j][i].value();
      }
    }
    for(int i = 0; i < nDofu; i++) {
      Res[N * nDofu + i] = -aResu[i].value();
    }
    RES->add_vector_blocked(Res, l2GMap);

    // define the dependent and independent variables
    for(unsigned j = 0; j < N; j++) {
      s.dependent(&aResw[j][0], nDofu);
    }
    s.dependent(&aResu[0], nDofu);

    for(unsigned j = 0; j < N; j++) {
      s.independent(&solw[j][0], nDofu);
    }

    s.independent(&solu[0], nDofu);

    // get the jacobian matrix (ordered by row major )
    Jac.resize(nDofAll * nDofAll);    //resize
    s.jacobian(&Jac[0], true);

    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac, l2GMap, l2GMap);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();

  KK->close();

//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//   PetscObjectSetName((PetscObject)viewer, "FSI matrix");
//   PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(), viewer);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(),  PETSC_VIEWER_STDOUT_WORLD);
//   double a;
//   std::cin >> a;

  // ***************** END ASSEMBLY *******************
}

void BuildU(MultiLevelSolution& mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel(level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel(level);
  unsigned iproc  = msh->processor_id();

//   double q_step = 1. / sqrt( N );

  int N_plus = round(pow(M_PI, 2) / (4. * (1 - S_FRAC) * pow(q_step, 2)));
  int N_minus = round(pow(M_PI, 2) / (4. * S_FRAC * pow(q_step, 2))) ;

  std::vector< unsigned > wIndex(N_minus + N_plus + 1);
//   std::vector< unsigned > wIndex(2 * N + 1);

  for(int i = - N_minus; i < N_plus + 1; i++) {
//   for (int i = - N; i < N + 1; i++) {
    char solName[10];
    sprintf(solName, "w%d", i);
    wIndex[i + N_minus] = mlSol.GetIndex(solName);
//     wIndex[i + N] = mlSol.GetIndex (solName);
  }
  unsigned uIndex = mlSol.GetIndex("u");

  unsigned solType = mlSol.GetSolutionType(uIndex);

  sol->_Sol[uIndex]->zero();

  for(unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    double value = 0.;

    double Cs =  2 * sin(M_PI * S_FRAC) / M_PI ;

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
    for(int j = -N_minus ; j < N_plus + 1; j++) {
//     for (int k = 0; k < 2 * N + 1; k++) {
      double weight = 1.;//exp( 2 * S_FRAC * q_step * j );
//       double weight = exp( /*2 **/ S_FRAC * q_step * ( k - N) );

      value += Cs * q_step * weight * (*sol->_Sol[wIndex[j + N_minus]])(i);
      std::cout.precision(14);
      if(i == 32) std::cout << j << " " << Cs * q_step * weight << "  " <<
                              (*sol->_Sol[wIndex[j + N_minus]])(i) << "  " <<
                              Cs * q_step * weight * (*sol->_Sol[wIndex[j + N_minus]])(i) << "  " << value << "\n";

//       if(k == 0 || k == 1) std::cout<< k<< "  " << i<< "  " <<  (*sol->_Sol[wIndex[k]])(i) << "\n" ;

    }
    sol->_Sol[uIndex]->set(i, value);
  }

  sol->_Sol[uIndex]->close();

}








