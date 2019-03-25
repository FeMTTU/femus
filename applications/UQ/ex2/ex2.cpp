
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"

#include "NumericVector.hpp"
#include "adept.h"

#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

#include "slepceps.h"

#include "../include/sfem_assembly.hpp"

using namespace femus;


bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;
  return dirichlet;
}

void GetEigenPair(MultiLevelProblem& ml_prob, const int& numberOfEigPairs, std::vector < std::pair<double, double> >& eigenvalues);

void GetQuantityOfInterest(MultiLevelProblem& ml_prob, std::vector < double >&  QoI, const unsigned& m, const double& domainMeasure);

void GetStochasticData(std::vector <double>& QoI);

void PlotStochasticData();

//BEGIN stochastic data
double L = 0.1 ; // correlation length of the covariance function
double domainMeasure = 1.; //measure of the domain
unsigned totMoments = 6;
std::vector <double> moments(totMoments, 0.); //initialization
std::vector <double> momentsStandardized(totMoments, 0.); //initialization
std::vector <double> cumulants(totMoments, 0.); //initialization
std::vector <double> cumulantsStandardized(totMoments, 0.); //initialization
double meanQoI = 0.; //initialization
double varianceQoI = 0.; //initialization
double stdDeviationQoI = 0.; //initialization
double startPoint = - 3.8;
double endPoint = 3.8;
unsigned M = 10000; //number of samples for the Monte Carlo
double deltat;
int pdfHistogramSize;
//END

unsigned numberOfUniformLevels = 4;

int main(int argc, char** argv) {


  //BEGIN eigenvalue problem instances
  PetscErrorCode ierr;
  ierr = SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

//   numberOfEigPairs = 5; //number of eigenpairs desired

  eigenvalues.resize(numberOfEigPairs); //this is where we store the eigenvalues

  //END


  //BEGIN deterministic FEM instances

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.ReadCoarseMesh("../input/square.neu", "fifth", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);

  unsigned dim = mlMsh.GetDimension();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("u", LAGRANGE, SECOND, 2);

  for(unsigned i = 0; i < numberOfEigPairs; i++) {
    char name[10];
    sprintf(name, "egnf%d", i);
    mlSol.AddSolution(name, LAGRANGE, SECOND, 0, false);
  }

  mlSol.Initialize("All");

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc("All");

  MultiLevelProblem ml_prob(&mlSol);

  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("UQ");
  system.AddSolutionToSystemPDE("u");

  // ******* System FEM Assembly *******
  system.SetAssembleFunction(AssembleUQSys);
  system.SetMaxNumberOfLinearIterations(1);
  //system.SetAssembleFunction(AssembleFEM);
  // ******* set MG-Solver *******
  system.SetMgType(V_CYCLE);

  system.SetAbsoluteLinearConvergenceTolerance(1.e-50);
  //   system.SetNonLinearConvergenceTolerance(1.e-9);
//   system.SetMaxNumberOfNonLinearIterations(20);

  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);

  // ******* Set Preconditioner *******
  system.SetLinearEquationSolverType(FEMuS_DEFAULT);

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids(GMRES);

  system.SetPreconditionerFineGrids(ILU_PRECOND);

  system.SetTolerances(1.e-20, 1.e-20, 1.e+50, 100);
  //END

  //BEGIN stochastic index set test
//   std::vector < std::vector <unsigned> > Jp;
//   ComputeIndexSetJp(Jp, 3);
//   for(unsigned i=0; i<Jp.size();i++){
//     std::cout << "Jp[" << i << "]=( ";
//     for(unsigned j=0; j<numberOfEigPairs; j++){
//      std::cout << Jp[i][j] << " ";
//    }
//    std::cout << ")" << std::endl;
//   }
//   std::vector < std::vector < std::vector < double > > >  G;
//   EvaluateStochasticMassMatrices(5, 3, G);

  //END stochastic index set test


  GetEigenPair(ml_prob, numberOfEigPairs, eigenvalues); //solve the generalized eigenvalue problem and compute the eigenpairs

  for(int i = 0; i < numberOfEigPairs; i++) {
    std::cout << eigenvalues[i].first << " " << eigenvalues[i].second << std::endl;
  }



  std::vector <double> QoI(M, 0.);


  for(unsigned m = 0; m < M; m++) {

    std::cout << " --------------------------------------------------- m = " << m << " ---------------------------------------------------  " << std::endl;

    system.MGsolve();

    GetQuantityOfInterest(ml_prob, QoI, m, domainMeasure);

  }

//   for(unsigned m = 0; m < M; m++) {
//     std::cout << "QoI[" << m << "] = " << QoI[m] << std::endl;
//   }

  GetStochasticData(QoI);

  PlotStochasticData();

  // ******* Print solution *******
  mlSol.SetWriter(VTK);
  std::vector<std::string> print_vars;
  print_vars.push_back("All");
  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);

  //ierr = SlepcFinalize();
  //CHKERRQ(ierr);

  return 0;

} //end main

void GetEigenPair(MultiLevelProblem& ml_prob, const int& numberOfEigPairs, std::vector < std::pair<double, double> >& eigenvalues) {
//void GetEigenPair(MultiLevelProblem & ml_prob, Mat &CCSLEPc, Mat &MMSLEPc) {

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("UQ");   // pointer to the linear implicit system named "Poisson"

  unsigned level = numberOfUniformLevels - 1;

  double varianceInput = stdDeviationInput * stdDeviationInput;

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             MM = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned solType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector < vector < double > > x1(dim);    // local coordinates
  vector < vector < double > > x2(dim);    // local coordinates
  for(unsigned k = 0; k < dim; k++) {
    x1[k].reserve(maxSize);
    x2[k].reserve(maxSize);
  }

  vector <double> phi_x; // local test function first order partial derivatives

  phi_x.reserve(maxSize * dim);

  vector< int > l2GMap1; // local to global mapping
  vector< int > l2GMap2; // local to global mapping
  l2GMap1.reserve(maxSize);
  l2GMap2.reserve(maxSize);

  vector < double > MMlocal;
  MMlocal.reserve(maxSize * maxSize);

  vector < double > CClocal;
  CClocal.reserve(maxSize * maxSize);

  MM->zero(); // Set to zero all the entries of the Global Matrix

  int MM_size = msh->_dofOffset[solType][nprocs];
  int MM_local_size = msh->_dofOffset[solType][iproc + 1] - msh->_dofOffset[solType][iproc];

  SparseMatrix* CC;
  CC = SparseMatrix::build().release();
  CC->init(MM_size, MM_size, MM_local_size, MM_local_size, MM_local_size, MM_size - MM_local_size);
  CC->zero();

  for(int kproc = 0; kproc < nprocs; kproc++) {
    for(int jel = msh->_elementOffset[kproc]; jel < msh->_elementOffset[kproc + 1]; jel++) {

      short unsigned ielGeom2;
      unsigned nDof2;
      unsigned nDofx2;

      if(iproc == kproc) {
        ielGeom2 = msh->GetElementType(jel);
        nDof2  = msh->GetElementDofNumber(jel, solType);    // number of solution element dofs
        nDofx2 = msh->GetElementDofNumber(jel, xType);    // number of coordinate element dofs
      }

      MPI_Bcast(&ielGeom2, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD);
      MPI_Bcast(&nDof2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      MPI_Bcast(&nDofx2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);

      // resize local arrays
      l2GMap2.resize(nDof2);

      for(int k = 0; k < dim; k++) {
        x2[k].resize(nDofx2);
      }

      // local storage of global mapping and solution
      if(iproc == kproc) {
        for(unsigned j = 0; j < nDof2; j++) {
          l2GMap2[j] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, j, jel);  // global to global mapping between solution node and pdeSys dof
        }
      }
      MPI_Bcast(&l2GMap2[0], nDof2, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);

      // local storage of coordinates
      if(iproc == kproc) {
        for(unsigned j = 0; j < nDofx2; j++) {
          unsigned xDof  = msh->GetSolutionDof(j, jel, xType);  // global to global mapping between coordinates node and coordinate dof
          for(unsigned k = 0; k < dim; k++) {
            x2[k][j] = (*msh->_topology->_Sol[k])(xDof);  // global extraction and local storage for the element coordinates
          }
        }
      }
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& x2[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }

      unsigned jgNumber = msh->_finiteElement[ielGeom2][solType]->GetGaussPointNumber();
      vector < vector < double > > xg2(jgNumber);
      vector <double> weight2(jgNumber);
      vector < vector <double> > phi2(jgNumber);  // local test function

      for(unsigned jg = 0; jg < jgNumber; jg++) {
        msh->_finiteElement[ielGeom2][solType]->Jacobian(x2, jg, weight2[jg], phi2[jg], phi_x);

        xg2[jg].assign(dim, 0.);

        for(unsigned j = 0; j < nDof2; j++) {
          for(unsigned k = 0; k < dim; k++) {
            xg2[jg][k] += x2[k][j] * phi2[jg][j];
          }
        }
      }

      // element loop: each process loops only on the elements that owns
      for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

        short unsigned ielGeom1 = msh->GetElementType(iel);
        unsigned nDof1  = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
        unsigned nDofx1 = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

        // resize local arrays
        l2GMap1.resize(nDof1);
        //std::vector<bool>bdcDirichlet(nDof1);

        for(int k = 0; k < dim; k++) {
          x1[k].resize(nDofx1);
        }

        // local storage of global mapping and solution
        for(unsigned i = 0; i < nDof1; i++) {
          l2GMap1[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
          //unsigned solDof = msh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
          //bdcDirichlet[i] = ( (*sol->_Bdc[soluIndex])(solDof) < 1.5)? false:false;
        }

        // local storage of coordinates
        for(unsigned i = 0; i < nDofx1; i++) {
          unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
          for(unsigned k = 0; k < dim; k++) {
            x1[k][i] = (*msh->_topology->_Sol[k])(xDof);  // global extraction and local storage for the element coordinates
          }
        }

        if(iel == jel) MMlocal.assign(nDof1 * nDof1, 0.);  //resize
        CClocal.assign(nDof1 * nDof2, 0.);   //resize

        // *** Gauss point loop ***
        unsigned igNumber = msh->_finiteElement[ielGeom1][solType]->GetGaussPointNumber();
        double weight1;
        vector <double> phi1;  // local test function
        for(unsigned ig = 0; ig < igNumber; ig++) {

          msh->_finiteElement[ielGeom1][solType]->Jacobian(x1, ig, weight1, phi1, phi_x);

          // evaluate the solution, the solution derivatives and the coordinates in the gauss point
          vector < double > xg1(dim, 0.);
          for(unsigned i = 0; i < nDof1; i++) {
            for(unsigned k = 0; k < dim; k++) {
              xg1[k] += x1[k][i] * phi1[i];
            }
          }

          if(iel == jel) {
            for(unsigned i = 0; i < nDof1; i++) {
              for(unsigned i1 = 0; i1 < nDof1; i1++) {
                MMlocal[ i * nDof1 + i1 ] += phi1[i] * phi1[i1] * weight1;
              }
            }
          }

          for(unsigned jg = 0; jg < jgNumber; jg++) {
            double dist = 0.;
            for(unsigned k = 0; k < dim; k++) {
              dist += fabs(xg1[k] - xg2[jg][k]);
            }
            double C = varianceInput * exp(- dist / L);
            for(unsigned i = 0; i < nDof1; i++) {
              for(unsigned j = 0; j < nDof2; j++) {
                CClocal[i * nDof2 + j] += weight1 * phi1[i] * C * phi2[jg][j] * weight2[jg];
              }//endl j loop
            } //endl i loop
          } //endl jg loop
        } //endl ig loop
        if(iel == jel) MM->add_matrix_blocked(MMlocal, l2GMap1, l2GMap1);
        CC->add_matrix_blocked(CClocal, l2GMap1, l2GMap2);
      } // end iel loop
    } //end jel loop
  } //end kproc loop

  MM->close();
  CC->close();

  //BEGIN solve the eigenvalue problem

  int ierr;
  EPS eps;
  PetscInt convergedSolns, numberOfIterations;

  ierr = EPSCreate(PETSC_COMM_WORLD, &eps);
  CHKERRABORT(MPI_COMM_WORLD, ierr);
  ierr = EPSSetOperators(eps, (static_cast<PetscMatrix*>(CC))->mat(), (static_cast<PetscMatrix*>(MM))->mat());
  CHKERRABORT(MPI_COMM_WORLD, ierr);
  ierr = EPSSetFromOptions(eps);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  //ierr = EPSSetDimensions(eps, numberOfEigPairs, 8 * numberOfEigPairs, 600);
  ierr = EPSSetDimensions(eps, numberOfEigPairs, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRABORT(MPI_COMM_WORLD, ierr);
  ierr = EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  //ierr = EPSSetTolerances(eps,1.0e-10,1000);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  ierr = EPSSolve(eps);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  //ierr = EPSView(eps, PETSC_VIEWER_STDOUT_SELF);

  std::cout << " -----------------------------------------------------------------" << std::endl;

  ierr = EPSGetConverged(eps, &convergedSolns);
  CHKERRABORT(MPI_COMM_WORLD, ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of converged eigenpairs: %D\n\n", convergedSolns);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  if(convergedSolns > 0) {

    for(unsigned i = 0; i < numberOfEigPairs; i++) {

      char name[10];
      sprintf(name, "egnf%d", i);
      soluIndex = mlSol->GetIndex(name);    // get the position of "u" in the ml_sol object

      // Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and ki (imaginary part)

      ierr = EPSGetEigenpair(eps, i, &eigenvalues[i].first, &eigenvalues[i].second, (static_cast<PetscVector*>(sol->_Sol[soluIndex]))->vec(), NULL);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

    }
  }

  ierr = EPSDestroy(&eps);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  delete CC;

  //BEGIN OLD
//    std::vector <unsigned> eigfIndex(numberOfEigPairs);
//   char name[10];
//   for(unsigned i = 0; i < numberOfEigPairs; i++) {
//     sprintf(name, "egnf%d", i);
//     eigfIndex[i] = mlSol->GetIndex(name);    // get the position of "u" in the ml_sol object
//   }
//
//   std::vector < double > local_integral(numberOfEigPairs, 0.);
//   std::vector < double > local_norm2(numberOfEigPairs, 0.);
//
//   vector < vector < double > > eigenFunction(numberOfEigPairs); // local solution
//
//   for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
//
//     short unsigned ielGeom = msh->GetElementType(iel);
//     unsigned nDofu  = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
//     unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
//
//     // resize local arrays
//     for(unsigned i = 0; i < numberOfEigPairs; i++) {
//       eigenFunction[i].resize(nDofu);
//     }
//
//     for(int i = 0; i < dim; i++) {
//       x1[i].resize(nDofx);
//     }
//
//     // local storage of global mapping and solution
//     for(unsigned i = 0; i < nDofu; i++) {
//       unsigned solDof = msh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
//       for(unsigned j = 0; j < numberOfEigPairs; j++) {
//         eigenFunction[j][i] = (*sol->_Sol[eigfIndex[j]])(solDof);
//       }
//     }
//
//     // local storage of coordinates
//     for(unsigned i = 0; i < nDofx; i++) {
//       unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
//       for(unsigned jdim = 0; jdim < dim; jdim++) {
//         x1[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
//       }
//     }
//     double weight;
//     vector <double> phi;  // local test function
//     // *** Gauss point loop ***
//     for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
//       // *** get gauss point weight, test function and test function partial derivatives ***
//       msh->_finiteElement[ielGeom][solType]->Jacobian(x1, ig, weight, phi, phi_x, *nullDoublePointer);
//       for(unsigned j = 0; j < numberOfEigPairs; j++) {
//         double eigenFunction_gss = 0.;
//         for(unsigned i = 0; i < nDofu; i++) {
//           eigenFunction_gss += phi[i] * eigenFunction[j][i];
//         }
//         local_integral[j] += eigenFunction_gss * weight;
//         local_norm2[j] += eigenFunction_gss * eigenFunction_gss * weight;
//       }
//     }
//   }
//   for(unsigned j = 0; j < numberOfEigPairs; j++) {
//     double integral = 0.;
//     MPI_Allreduce(&local_integral[j], &integral, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//     double sign = (integral >= 0) ? 1 : -1;
//     double norm2 = 0.;
//     MPI_Allreduce(&local_norm2[j], &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//     double inorm = /*0.01 * sign */ 1. / sqrt(norm2);
//     std::cout << "BBBBBBBBBBBBBBBBBB  " << inorm << std::endl;
//     sol->_Sol[eigfIndex[j]]->scale(inorm);
//   }
//
  //END OLD



  //BEGIN GRAM SCHMIDT ORTHONORMALIZATION

  std::vector <unsigned> eigfIndex(numberOfEigPairs);
  char name[10];
  for(unsigned i = 0; i < numberOfEigPairs; i++) {
    sprintf(name, "egnf%d", i);
    eigfIndex[i] = mlSol->GetIndex(name);    // get the position of "u" in the ml_sol object
  }

  vector < double >  eigenFunction(numberOfEigPairs); // local solution
  vector < double >  eigenFunctionOld(numberOfEigPairs); // local solution

  std::vector < std::vector < double > > coeffsGS_local(numberOfEigPairs);
  std::vector < std::vector < double > > coeffsGS_global(numberOfEigPairs);
  for(unsigned i = 0; i < numberOfEigPairs; i++) {
    coeffsGS_local[i].assign(numberOfEigPairs, 0.);
    coeffsGS_global[i].assign(numberOfEigPairs, 0.);
  }

  for(unsigned iGS = 0; iGS < numberOfEigPairs; iGS++) {

    if(iGS > 0) {

      for(unsigned jGS = 0; jGS < iGS; jGS++) {

        //BEGIN COMPUTE coeffsGS LOCAL

        for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

          short unsigned ielGeom = msh->GetElementType(iel);
          unsigned nDofu  = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
          unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

          eigenFunction.resize(nDofu);
          eigenFunctionOld.resize(nDofu);

          for(int i = 0; i < dim; i++) {
            x1[i].resize(nDofx);
          }

          // local storage of global mapping and solution
          for(unsigned i = 0; i < nDofu; i++) {
            unsigned solDof = msh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
            eigenFunction[i] = (*sol->_Sol[eigfIndex[iGS]])(solDof);
            eigenFunctionOld[i] = (*sol->_Sol[eigfIndex[jGS]])(solDof);
          }

          // local storage of coordinates
          for(unsigned i = 0; i < nDofx; i++) {
            unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
            for(unsigned jdim = 0; jdim < dim; jdim++) {
              x1[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
            }
          }
          double weight;
          vector <double> phi;  // local test function
          // *** Gauss point loop ***
          for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
            // *** get gauss point weight, test function and test function partial derivatives ***
            msh->_finiteElement[ielGeom][solType]->Jacobian(x1, ig, weight, phi, phi_x);
            double eigenFunction_gss = 0.;
            double eigenFunction_gss_old = 0.;
            for(unsigned i = 0; i < nDofu; i++) {
              eigenFunction_gss += phi[i] * eigenFunction[i];
              eigenFunction_gss_old += phi[i] * eigenFunctionOld[i];
            }
            coeffsGS_local[iGS][jGS] += eigenFunction_gss * eigenFunction_gss_old * weight;
          }
        }

        //END COMPUTE coeffsGS LOCAL

        MPI_Allreduce(&coeffsGS_local[iGS][jGS], &coeffsGS_global[iGS][jGS], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      }

      for(unsigned idof = msh->_dofOffset[solType][iproc]; idof < msh->_dofOffset[solType][iproc + 1]; idof++) {
        double sum = 0.;
        for(unsigned jGS = 0; jGS < iGS; jGS++) {
          sum += coeffsGS_global[iGS][jGS] * (*sol->_Sol[eigfIndex[jGS]])(idof);
        }
        double valueToSet = (*sol->_Sol[eigfIndex[iGS]])(idof) - sum;
        sol->_Sol[eigfIndex[iGS]]->set(idof, valueToSet);
      }

    }

    sol->_Sol[eigfIndex[iGS]]->close();
    
    double local_norm2 = 0.;
    for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofu  = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
      unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

      eigenFunction.resize(nDofu);

      for(int i = 0; i < dim; i++) {
        x1[i].resize(nDofx);
      }

      // local storage of global mapping and solution
      for(unsigned i = 0; i < nDofu; i++) {
        unsigned solDof = msh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
        eigenFunction[i] = (*sol->_Sol[eigfIndex[iGS]])(solDof);
      }

      // local storage of coordinates
      for(unsigned i = 0; i < nDofx; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
        for(unsigned jdim = 0; jdim < dim; jdim++) {
          x1[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
        }
      }
      double weight;
      vector <double> phi;  // local test function
      // *** Gauss point loop ***
      for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[ielGeom][solType]->Jacobian(x1, ig, weight, phi, phi_x);
        double eigenFunction_gss = 0.;
        for(unsigned i = 0; i < nDofu; i++) {
          eigenFunction_gss += phi[i] * eigenFunction[i];
        }
        local_norm2 += eigenFunction_gss * eigenFunction_gss * weight;
      }
    }

    double norm2 = 0.;
    MPI_Allreduce(&local_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double norm = sqrt(norm2);
    std::cout << "norm = " << norm << std::endl;
    sol->_Sol[eigfIndex[iGS]]->scale(1. / norm);

    sol->_Sol[eigfIndex[iGS]]->close();

  }

  //END GRAM SCHMIDT ORTHONORMALIZATION

  //BEGIN GRAM SCHMIDT CHECK
  vector < double >  eigenFunctionCheck(numberOfEigPairs); // local solution
  vector < double >  eigenFunctionOldCheck(numberOfEigPairs); // local solution


  for(unsigned i1 = 0; i1 < numberOfEigPairs; i1++) {
    for(unsigned j1 = 0; j1 < numberOfEigPairs; j1++) {

      double integral = 0.;
      for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

        short unsigned ielGeom = msh->GetElementType(iel);
        unsigned nDofu  = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
        unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

        eigenFunctionCheck.resize(nDofu);
        eigenFunctionOldCheck.resize(nDofu);

        for(int i = 0; i < dim; i++) {
          x1[i].resize(nDofx);
        }

        // local storage of global mapping and solution
        for(unsigned i = 0; i < nDofu; i++) {
          unsigned solDof = msh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
          eigenFunctionCheck[i] = (*sol->_Sol[eigfIndex[i1]])(solDof);
          eigenFunctionOldCheck[i] = (*sol->_Sol[eigfIndex[j1]])(solDof);
        }

        // local storage of coordinates
        for(unsigned i = 0; i < nDofx; i++) {
          unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
          for(unsigned jdim = 0; jdim < dim; jdim++) {
            x1[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
          }
        }
        double weight;
        vector <double> phi;  // local test function
        // *** Gauss point loop ***
        for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
          // *** get gauss point weight, test function and test function partial derivatives ***
          msh->_finiteElement[ielGeom][solType]->Jacobian(x1, ig, weight, phi, phi_x);
          double eigenFunction_gss = 0.;
          double eigenFunction_gss_old = 0.;
          for(unsigned i = 0; i < nDofu; i++) {
            eigenFunction_gss += phi[i] * eigenFunctionCheck[i];
            eigenFunction_gss_old += phi[i] * eigenFunctionOldCheck[i];
          }
          integral += eigenFunction_gss * eigenFunction_gss_old * weight;
        }
      }

      double globalIntegral = 0.;
      MPI_Allreduce(&integral, &globalIntegral, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      std::cout << "i = " << i1 << " , " << "j = " << j1 << " , " << "integral = " << globalIntegral << std::endl;
    }
  }

  //END GRAM SCHMIDT CHECK

  // ***************** END ASSEMBLY *******************
}


void GetQuantityOfInterest(MultiLevelProblem& ml_prob, std::vector < double >&  QoI, const unsigned& m, const double& domainMeasure) {

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("UQ");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

//   unsigned soluPdeIndex;
//   soluPdeIndex = mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object

  vector < double >  solu; // local solution
  solu.reserve(maxSize);

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for(unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);

  double quantityOfInterest = 0.;

  // element loop: each process loops only on the elements that owns

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    solu.resize(nDofu);

    for(int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for(unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }


    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x);

      double solu_gss = 0.;
      for(unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i];
      }
//       quantityOfInterest += solu_gss * solu_gss *  weight ; // this is the integral of the square.
      quantityOfInterest +=  solu_gss *  weight / domainMeasure; // this is the spatial average over the domain.

    } // end gauss point loop

  } //end element loop for each process

  QoI[m] = 0.;
  MPI_Allreduce(&quantityOfInterest, &QoI[m], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //QoI[m] = quantityOfInterest;

//   // add the quantityOfInterest of all processes
//   NumericVector* norm_vec;
//   norm_vec = NumericVector::build().release();
//   norm_vec->init (msh->n_processors(), 1 , false, AUTOMATIC);
//
//   norm_vec->set (iproc, quantityOfInterest);
//   norm_vec->close();
//
//   quantityOfInterest = norm_vec->l1_norm();
//
//   delete norm_vec;
//
//   QoI[m] = quantityOfInterest;

}

void GetStochasticData(std::vector <double>& QoI) {

  //let's standardize the quantity of interest after finding moments and standard deviation

  if(totMoments <= 0) {

    std::cout << "ERROR: total number of moments has to be a positive integer" << std::endl;

  }

  else {

    //BEGIN computation of the mean of QoI
    meanQoI = 0.;
    unsigned meanCounter = 0;
    for(unsigned m = 0; m < M; m++) {
      if(QoI[m] < endPoint && QoI[m] > startPoint) {
        meanQoI += QoI[m];
        meanCounter++;
      }
    }
    meanQoI /= meanCounter;
    //END


    //BEGIN computation of the variance and standard deviation of QoI
    varianceQoI = 0;
    unsigned varianceCounter = 0;
    for(unsigned m = 0; m < M; m++) {
      if(QoI[m] < endPoint && QoI[m] > startPoint) {
        varianceQoI += (QoI[m] - meanQoI) * (QoI[m] - meanQoI);
        varianceCounter++;
      }
    }
    varianceQoI /= varianceCounter;

    stdDeviationQoI = sqrt(varianceQoI);
    //END

    pdfHistogramSize = static_cast <int>(floor(1. + 3.3 * log(M)));
    std::vector <double> pdfHistogram(pdfHistogramSize, 0.);
    double lengthOfTheInterval = fabs(endPoint - startPoint);
    deltat = lengthOfTheInterval / pdfHistogramSize;

    std::vector < double > QoIStandardized(M, 0.);
    //BEGIN standardization of QoI before computing the moments
    for(unsigned m = 0; m < M; m++) {
      QoIStandardized[m] = (QoI[m] - meanQoI) / stdDeviationQoI ;
//      std::cout<< "standardized QoI " << QoIStandardized[m] << std::endl;

      //BEGIN estimation of the PDF
      bool sampleCaptured = false;
      for(unsigned i = 0; i < pdfHistogramSize; i++) {
        double leftBound = startPoint + i*deltat;
        double rightBound = startPoint + (i+1)*deltat;
        if(leftBound <=  QoIStandardized[m] && QoIStandardized[m] < rightBound) {
          pdfHistogram[i]++;
//           std::cout << "leftBound = " << leftBound << " " << "rightBound = " << rightBound << " " << " standardized QoI = " << QoIStandardized[m] << std::endl;
          sampleCaptured = true;
          break;
        }
      }
      if(sampleCaptured == false) {
        std::cout << "WARNING: sample " << QoIStandardized[m] << "is not in any interval" << std::endl;
      }
      //END
    }
    //END

    double pdfIntegral = 0;
    for(unsigned i = 0; i < pdfHistogramSize; i++) {
      pdfIntegral += pdfHistogram[i] * deltat;
    }


    //BEGIN histogram check
    double checkHistogram = 0;
    for(unsigned i = 0; i < pdfHistogramSize; i++) {
      double point = (startPoint + i * deltat + startPoint + (i+1) * deltat) * 0.5;
      double pdfCheck = pdfHistogram[i] / M;
      pdfHistogram[i] /= pdfIntegral;
      std::cout << point << "  " << pdfHistogram[i]  << std::endl;
      //       std::cout << "{" << point << "," << pdfHistogram[i]  << "}," << std::endl;
      checkHistogram += pdfCheck;
    }
    std::cout << "checkHistogram = " << checkHistogram << std::endl;
    //END



    //BEGIN computation of the raw moments
    for(unsigned p = 0; p < totMoments; p++) {
      moments[p] = 0.;
      momentsStandardized[p] = 0.;
      unsigned momentsCounter = 0;
      for(unsigned m = 0; m < M; m++) {
        if(QoI[m] < endPoint && QoI[m] > startPoint) {
          moments[p] += pow(QoI[m], p + 1);
          momentsStandardized[p] += pow(QoIStandardized[m], p + 1);
          momentsCounter++;
        }
      }
      moments[p] /= momentsCounter;
      momentsStandardized[p] /= momentsCounter;
    }
    //END

    cumulants[0] = moments[0];
    cumulantsStandardized[0] = momentsStandardized[0];

    if(totMoments > 1) {
      cumulants[1] = moments[1] - moments[0] * moments[0];

      cumulantsStandardized[1] = momentsStandardized[1] - momentsStandardized[0] * momentsStandardized[0];
//       std::cout.precision(14);
//       std::cout << "AAAAAAAAAAAAAAA" << cumulants[1] << std::endl;
      if(totMoments > 2) {
        cumulants[2] = moments[2] - 3. * moments[1] * moments[0] + 2. * pow(moments[0], 3);

        cumulantsStandardized[2] = momentsStandardized[2] - 3. * momentsStandardized[1] * momentsStandardized[0] + 2. * pow(momentsStandardized[0], 3);
        if(totMoments > 3) {
          cumulants[3] = moments[3] - 4. * moments[2] * moments[0] - 3. * moments[1] * moments[1] + 12. * moments[1] * moments[0] * moments[0] - 6. * pow(moments[0], 4);

          cumulantsStandardized[3] = momentsStandardized[3] - 4. * momentsStandardized[2] * momentsStandardized[0] - 3. * momentsStandardized[1] * momentsStandardized[1]
                                     + 12. * momentsStandardized[1] * momentsStandardized[0] * momentsStandardized[0] - 6. * pow(momentsStandardized[0], 4);
          if(totMoments > 4) {
            cumulants[4] = moments[4] - 5. * moments[3] * moments[0] - 10. * moments[2] * moments[1] + 20. * moments[2] * moments[0] * moments[0]
                           + 30. * moments[1] * moments[1] * moments[0] - 60. * moments[1] * pow(moments[0], 3) + 24. * pow(moments[0], 5);

            cumulantsStandardized[4] = momentsStandardized[4] - 5. * momentsStandardized[3] * momentsStandardized[0] - 10. * momentsStandardized[2] * momentsStandardized[1]
                                       + 20. * momentsStandardized[2] * momentsStandardized[0] * momentsStandardized[0]
                                       + 30. * momentsStandardized[1] * momentsStandardized[1] * momentsStandardized[0]
                                       - 60. * momentsStandardized[1] * pow(momentsStandardized[0], 3) + 24. * pow(momentsStandardized[0], 5);
            if(totMoments > 5) {
              cumulants[5] = moments[5] - 6. * moments[4] * moments[0] - 15. * moments[3] * moments[1] + 30. * moments[3] * moments[0] * moments[0]
                             - 10. * moments[2] * moments[2] + 120. * moments[2] * moments[1] * moments[0] - 120. * moments[2] * pow(moments[0], 3)
                             + 30. * pow(moments[1], 3) - 270. * pow(moments[1], 2) * pow(moments[0], 2) + 360. * moments[1] * pow(moments[0], 4) - 120. * pow(moments[0], 6);

              cumulantsStandardized[5] = momentsStandardized[5] - 6. * momentsStandardized[4] * momentsStandardized[0] - 15. * momentsStandardized[3] * momentsStandardized[1]
                                         + 30. * momentsStandardized[3] * momentsStandardized[0] * momentsStandardized[0]
                                         - 10. * momentsStandardized[2] * momentsStandardized[2] + 120. * momentsStandardized[2] * momentsStandardized[1] * momentsStandardized[0]
                                         - 120. * momentsStandardized[2] * pow(momentsStandardized[0], 3) + 30. * pow(momentsStandardized[1], 3)
                                         - 270. * pow(momentsStandardized[1], 2) * pow(momentsStandardized[0], 2) + 360. * momentsStandardized[1] * pow(momentsStandardized[0], 4)
                                         - 120. * pow(momentsStandardized[0], 6);
            }
          }
        }
      }
    }
  }
}


void PlotStochasticData() {

  std::cout.precision(14);
  std::cout << " the number of MC samples is " << M << std::endl;
  std::cout << " the mean is " << meanQoI << std::endl;
  std::cout << " the standard deviation is " << sqrt(varianceQoI) << std::endl;
  std::cout << " the variance is " << varianceQoI << std::endl;

  std::cout << "Standardized Moments" << std::endl;
  for(unsigned p = 0; p < totMoments; p++) {
    std::cout << " & " << momentsStandardized[p] << "  ";
  }
  std::cout << std::endl;
  std::cout << "Standardized Cumulants" << std::endl;
  for(unsigned p = 0; p < totMoments; p++) {
    std::cout << " & " << cumulantsStandardized[p] << "  ";
  }
  std::cout << std::endl;
  std::cout << " Moments " << std::endl;
  for(unsigned p = 0; p < totMoments; p++) {
    std::cout << " & " << moments[p] << "  ";
  }
  std::cout << std::endl;
  std::cout << " Cumulants " << std::endl;
  for(unsigned p = 0; p < totMoments; p++) {
    std::cout << " & " << cumulants[p] << "  ";
  }

  std::cout << std::endl;
  std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;

  double edgeworth1Term = 0.;
  double edgeworth2Terms = 0.;
  double edgeworth3Terms = 0.;
  double edgeworth4Terms = 0.;
  double edgeworth5Terms = 0.;
  double edgeworth6Terms = 0.;

  double generalizedGC1Term = 0.;
  double generalizedGC2Terms = 0.;
  double generalizedGC3Terms = 0.;
  double generalizedGC4Terms = 0.;
  double generalizedGC5Terms = 0.;
  double generalizedGC6Terms = 0.;

  double lambda3 = 0.;
  double lambda4 = 0.;
  double lambda5 = 0.;
  double lambda6 = 0.;

  double d1gaussian;
  double d2gaussian;
  double d3gaussian;
  double d4gaussian;
  double d5gaussian;
  double d6gaussian;
  double d7gaussian;
  double d8gaussian;
  double d9gaussian;
  double d10gaussian;
  double d12gaussian;

//   cumulants[0] = 0; //decomment for nonStdGaussian

  //BEGIN GRAM CHARLIER PRINT
  std::cout << " ------------------------- GRAM CHARLIER ------------------------- " << std::endl;
  for(unsigned i = 0; i < pdfHistogramSize; i++) {
    double t = startPoint + i * deltat;
    std::cout << t << " ";
//     double t = x - meanQoI; //decomment for nonStdGaussian
    double gaussian = 1. / (sqrt(2 * acos(- 1))) * exp(- 0.5 * (t * t)) ;
    std::cout << gaussian << " ";

    d1gaussian = (- 1.) * gaussian * t ;

    generalizedGC1Term = gaussian - cumulantsStandardized[0] * d1gaussian;

    std::cout << generalizedGC1Term << " ";

    if(totMoments > 1) {

      d2gaussian = (1.) * gaussian * (t * t - 1.) ;
      d3gaussian = (- 1.) * gaussian * (t * t * t - 3. * t) ;

      generalizedGC2Terms = generalizedGC1Term + 0.5 * ((cumulantsStandardized[1] - 1.) + pow(cumulantsStandardized[0], 2)) * d2gaussian ;

      std::cout << generalizedGC2Terms << " ";

      if(totMoments > 2) {

        d4gaussian = (1.) * gaussian * (t * t * t * t - 6. * t * t + 3.) ;
        d6gaussian = (1.) * gaussian * (pow(t, 6) - 15 * pow(t, 4) + 45 * t * t - 15);

        generalizedGC3Terms = generalizedGC2Terms - 1. / 6 * (cumulantsStandardized[2] + 3 * (cumulantsStandardized[1] - 1.) * cumulantsStandardized[0]
                              + pow(cumulantsStandardized[0], 3)) * d3gaussian;

        std::cout << generalizedGC3Terms << " ";

        if(totMoments > 3) {

          d5gaussian = (- 1.) * gaussian * (pow(t, 5) - 10. * t * t * t + 15. * t);
          d7gaussian = (- 1.) * gaussian * (pow(t, 7) - 21. * pow(t, 5) + 105. * t * t * t -  105. * t) ;
          d9gaussian = (- 1.) * gaussian * (pow(t, 9) - 36. * pow(t, 7) + 378. * pow(t, 5) - 1260. * t * t * t + 945. * t) ;

          generalizedGC4Terms = generalizedGC3Terms + 1. / 24 * (cumulantsStandardized[3] + 4. * cumulantsStandardized[2] * cumulantsStandardized[0]
                                + 3. * pow((cumulantsStandardized[1] - 1.), 2) + 6. * (cumulantsStandardized[1] - 1.) * pow(cumulantsStandardized[0], 2)
                                + pow(cumulantsStandardized[0], 4)) * d4gaussian;

          std::cout << generalizedGC4Terms << " ";

          if(totMoments > 4) {

            generalizedGC5Terms = generalizedGC4Terms - 1. / 120 * (cumulantsStandardized[4] + 5. * cumulantsStandardized[3] * cumulantsStandardized[0]
                                  + 10. * cumulantsStandardized[2] * (cumulantsStandardized[1] - 1.) + 10. * cumulantsStandardized[2] * pow(cumulantsStandardized[0], 2)
                                  + 15. * pow((cumulantsStandardized[1] - 1.), 2) * cumulantsStandardized[0] + 10. * (cumulantsStandardized[1] - 1.) * pow(cumulantsStandardized[0], 3)
                                  + pow(cumulantsStandardized[0], 5)) * d5gaussian;

            std::cout << generalizedGC5Terms << " ";

            if(totMoments > 5) {

              generalizedGC6Terms = generalizedGC5Terms + 1. / 720 * (cumulantsStandardized[5] + 6. * cumulantsStandardized[4] * cumulantsStandardized[0]
                                    + 15. * cumulantsStandardized[3] * (cumulantsStandardized[1] - 1.) + 15. * cumulantsStandardized[3] * pow(cumulantsStandardized[0], 2)
                                    + 10. * pow(cumulantsStandardized[2], 2) + 60. * cumulantsStandardized[2] * (cumulantsStandardized[1] - 1.) * cumulantsStandardized[0]
                                    + 20. * cumulantsStandardized[2] * pow(cumulantsStandardized[0], 3) + 15. * pow((cumulantsStandardized[1] - 1.), 3)
                                    + 45. * pow((cumulantsStandardized[1] - 1.), 2) * pow(cumulantsStandardized[0], 2) + 15. * (cumulantsStandardized[1] - 1.) * pow(cumulantsStandardized[0], 4)
                                    +  pow(cumulantsStandardized[0], 6)) * d6gaussian;

              std::cout << generalizedGC6Terms << " \n ";

            }

          }
        }
      }
    }
  }

  //BEGIN EDGEWORTH PRINT
  std::cout << " ------------------------- EDGEWORTH ------------------------- " << std::endl;
  for(unsigned i = 0; i < pdfHistogramSize; i++) {
    double t = startPoint + i * deltat;
    std::cout << t << " ";
//     double t = x - meanQoI; //decomment for nonStdGaussian
    double gaussian = 1. / (sqrt(2 * acos(- 1))) * exp(- 0.5 * (t * t)) ;
    std::cout << gaussian << " ";

    d3gaussian = (- 1.) * gaussian * (t * t * t - 3. * t) ;

    lambda3 = cumulants[2] / pow(stdDeviationQoI, 3);

    edgeworth1Term = gaussian - lambda3 / 6. * d3gaussian;

    std::cout << edgeworth1Term << " ";

    if(totMoments > 1) {

      d4gaussian = (1.) * gaussian * (t * t * t * t - 6. * t * t + 3.) ;
      d6gaussian = (1.) * gaussian * (pow(t, 6) - 15 * pow(t, 4) + 45 * t * t - 15);

      lambda4 = cumulants[3] / pow(stdDeviationQoI, 4);

      edgeworth2Terms = edgeworth1Term + lambda4 / 24. * d4gaussian + lambda3 * lambda3 / 72. * d6gaussian;

      std::cout << edgeworth2Terms << " ";

      if(totMoments > 2) {

        d5gaussian = (- 1.) * gaussian * (pow(t, 5) - 10. * t * t * t + 15. * t);
        d7gaussian = (- 1.) * gaussian * (pow(t, 7) - 21. * pow(t, 5) + 105. * t * t * t -  105. * t) ;
        d9gaussian = (- 1.) * gaussian * (pow(t, 9) - 36. * pow(t, 7) + 378. * pow(t, 5) - 1260. * t * t * t + 945. * t) ;

        lambda5 = cumulants[4] / pow(stdDeviationQoI, 5);

        edgeworth3Terms = edgeworth2Terms - lambda5 / 120. * d5gaussian + lambda3 * lambda4 / 144. * d7gaussian + pow(lambda3, 3) / 1296. * d9gaussian;

        std::cout << edgeworth3Terms << " ";

      }
      if(totMoments > 3) {

        d8gaussian = (1.) * gaussian * (1. / 16. * (1680. - 6720. * t * t + 3360. * pow(t, 4) - 448. * pow(t, 6) + 16. * pow(t, 8)));
        d10gaussian = (1.) * gaussian * (1. / 32. * (- 30240. + 151200. *  t * t - 100800. * pow(t, 4) + 20160. * pow(t, 6) - 1440. * pow(t, 8)
                                         + 32. * pow(t, 10)));
        d12gaussian = (1.) * gaussian * (1. / 64. * (665280. - 3991680. * t * t + 3326400. * pow(t, 4) - 887040. * pow(t, 6) + 95040. * pow(t, 8)
                                         - 4224. * pow(t, 10) + 64. * pow(t, 12)));

        lambda6 = cumulants[5] / pow(stdDeviationQoI, 6);

        edgeworth4Terms = edgeworth3Terms + 1. / 720. * lambda6 * d6gaussian + (1. / 1152. * lambda4 * lambda4 + 1. / 720. * lambda3 * lambda5) * d8gaussian
                          + 1. / 1728. * lambda3 * lambda3 * lambda4 * d10gaussian + 1. / 31104. * pow(lambda3, 4) * d12gaussian;

        std::cout << edgeworth4Terms << " \n ";

      }
    }
  }
  //END

}





