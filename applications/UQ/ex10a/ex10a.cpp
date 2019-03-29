
#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "LinearImplicitSystem.hpp"
#include "Marker.hpp"
#include "NumericVector.hpp"
#include "adept.h"

#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

#include "slepceps.h"

#include "../include/sgfem_assembly_uq.hpp"

using namespace femus;

bool SetBoundaryCondition (const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;
  return dirichlet;
}

void GetEigenPair (MultiLevelProblem& ml_prob, const int& numberOfEigPairs, std::vector < std::pair<double, double> >& eigenvalues);
//
void GetCoefficientsForQuantityOfInterest (MultiLevelProblem& ml_prob, std::vector <double > &  alphas, const double& domainMeasure);
//
void GetMomentsAndCumulants (std::vector <double>& alphas);
//
void GetQoIStandardizedSamples (std::vector< double >& alphas, std::vector< std::vector <double > > & sgmQoIStandardized, std::vector< std::vector <double > > &                               sgmQoIStandardizedFinest, const unsigned & dimCoarseBox);
//
void GetHistogramAndKDE (std::vector< std::vector < double > > & sgmQoIStandardized, std::vector< std::vector <double > > & sgmQoIStandardizedFinest, MultiLevelProblem& ml_prob, MultiLevelProblem& ml_probFinest);
//
void PlotGCandEDExpansion();

void GetKDEIntegral (MultiLevelProblem& ml_prob);

void GetAverageL2Error (std::vector< std::vector <double > > & sgmQoIStandardized, MultiLevelProblem& ml_prob, MultiLevelProblem& ml_probFinest);

//BEGIN stochastic data for the PDE solution
double domainMeasure = 1.; //measure of the domain
unsigned totMoments = 6;
std::vector <double> moments (totMoments, 0.);   //initialization
std::vector <double> momentsStandardized (totMoments, 0.);   //initialization
std::vector <double> cumulants (totMoments, 0.);   //initialization
std::vector <double> cumulantsStandardized (totMoments, 0.);   //initialization
double meanQoI = 0.; //initialization
double varianceQoI = 0.; //initialization
double stdDeviationQoI = 0.; //initialization
double L = 0.1 ; // correlation length of the covariance function
unsigned kOrder = 7; //for order tests
unsigned numberOfSamples = 1000; //for MC sampling of the QoI
unsigned nxCoarseBox;
double xMinCoarseBox = - 5.; //-5.5 for Gaussian, -5. for SGM (average) with Gaussian KL, -3. for SGM (integral),  -1.5 for uniform, not KL: -0.6 for avg and int
double xMaxCoarseBox = 3.;  //5.5 for Gaussian, 3. for SGM (average) with Gaussian KL,  5.5 for SGM (integral), 1.5 for uniform, not KL: 0.6 for avg and 0.8 for int
unsigned nyCoarseBox;
double yMinCoarseBox = - 5.;
double yMaxCoarseBox = 3.;
unsigned nzCoarseBox;
double zMinCoarseBox = - 5.;
double zMaxCoarseBox = 3.;

unsigned numberOfSamplesFinest = 100000000; //10^6 for spatial average, 10^7 for "integral" of the square, 10^7 for SGM with random variable (not KL)
unsigned kOrderFinest = 8;
// unsigned nxCoarseBoxFinest = static_cast<unsigned> ( floor ( 1. + 3.3 * log ( numberOfSamplesFinest ) ) ); //for spatial average
// unsigned nxCoarseBoxFinest = static_cast<unsigned> ( floor ( 1. + 2. * log2 ( numberOfSamplesFinest ) ) ); //for integral of the square
unsigned nxCoarseBoxFinest = static_cast<unsigned> ( pow(2,kOrderFinest) );
unsigned nyCoarseBoxFinest = nxCoarseBoxFinest;
unsigned nzCoarseBoxFinest = nxCoarseBoxFinest;

bool histoFinest = true; //for SGM must be true
bool histoErr = false; //true only if the histogram error is to be calculated, for analytic sampling
double bLaplace = 1.5;
double muLaplace = 0.;
//END

unsigned numberOfUniformLevels = 5; //refinement for the PDE mesh

int main (int argc, char** argv) {

  PetscErrorCode ierr;
  ierr = SlepcInitialize (&argc, &argv, PETSC_NULL, PETSC_NULL);

  // init Petsc-MPI communicator
  FemusInit mpinit (argc, argv, MPI_COMM_WORLD);

  myuq.SetOutput (false);

  //BEGIN Add solutions to sol vector (not to PDE)
  eigenvalues.resize (numberOfEigPairs);   //this is where we store the eigenvalues

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.ReadCoarseMesh ("../input/square.neu", "fifth", scalingFactor);
  mlMsh.RefineMesh (numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);

//   unsigned dim = mlMsh.GetDimension();

  MultiLevelSolution mlSol (&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution ("u", LAGRANGE, SECOND, 2);

  for (unsigned i = 0; i < numberOfEigPairs; i++) {
    char name[10];
    sprintf (name, "egnf%d", i);
    mlSol.AddSolution (name, LAGRANGE, SECOND, 0, false);
  }

  const std::vector < std::vector <unsigned> > &Jp = myuq.GetIndexSet (pIndex, numberOfEigPairs);

  for (unsigned i = 0; i < Jp.size(); i++) {
    char name[10];
    sprintf (name, "uSG%d", i);
    mlSol.AddSolution (name, LAGRANGE, SECOND, 2);
  }

  mlSol.Initialize ("All");

  mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc ("All");
  //END


  //BEGIN define FEM system for eig problem
  MultiLevelProblem ml_prob (&mlSol);
  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("UQ");
  system.AddSolutionToSystemPDE ("u");

  // ******* System FEM Assembly *******
//   system.SetAssembleFunction(AssembleUQSys);
  system.SetMaxNumberOfLinearIterations (1);
  //system.SetAssembleFunction(AssembleFEM);
  // ******* set MG-Solver *******
  system.SetMgType (V_CYCLE);

  system.SetAbsoluteLinearConvergenceTolerance (1.e-50);
  system.SetNumberPreSmoothingStep (1);
  system.SetNumberPostSmoothingStep (1);

  // ******* Set Preconditioner *******
  system.SetLinearEquationSolverType (FEMuS_DEFAULT);

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids (GMRES);

  system.SetPreconditionerFineGrids (ILU_PRECOND);

  system.SetTolerances (1.e-20, 1.e-20, 1.e+50, 100);
  //END

  //BEGIN define SGM system for stochastic PDE
  MultiLevelProblem ml_probSG (&mlSol);

  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& systemSG = ml_probSG.add_system < LinearImplicitSystem > ("SG");


  for (unsigned i = 0; i < Jp.size(); i++) {
    char name[10];
    sprintf (name, "uSG%d", i);
    systemSG.AddSolutionToSystemPDE (name);
  }

  FieldSplitTree **FielduSGi;

  FielduSGi = new FieldSplitTree * [Jp.size()];

  std::vector < FieldSplitTree *> FSAll;
  FSAll.reserve (Jp.size());


  //BEGIN buid fieldSplitTree (only for FieldSplitPreconditioner)
  for (unsigned i = 0; i < Jp.size(); i++) {
    char name[10];
    sprintf (name, "uSG%d", i);
    std::vector < unsigned > fielduSGi (1);
    fielduSGi[0] = systemSG.GetSolPdeIndex (name);

    std::vector < unsigned > solutionTypeuSGi (1);
    solutionTypeuSGi[0] = mlSol.GetSolutionType (name);

    FielduSGi[i] = new FieldSplitTree (PREONLY, ILU_PRECOND, fielduSGi, solutionTypeuSGi, name);

    FSAll.push_back (FielduSGi[i]);
  }

  FieldSplitTree uSG (PREONLY, FIELDSPLIT_PRECOND, FSAll, "uSG");
  //uSG.SetRichardsonScaleFactor(1.);
  //END buid fieldSplitTree
   systemSG.SetLinearEquationSolverType (FEMuS_FIELDSPLIT);
  // ******* System FEM Assembly *******
  systemSG.SetAssembleFunction (AssembleSysSG);
  systemSG.SetMaxNumberOfLinearIterations (1);
  // ******* set MG-Solver *******
  systemSG.SetMgType (V_CYCLE);

  systemSG.SetAbsoluteLinearConvergenceTolerance (1.e-50);

  systemSG.SetNumberPreSmoothingStep (1);
  systemSG.SetNumberPostSmoothingStep (1);

  // ******* Set Preconditioner *******
  // systemSG.SetLinearEquationSolverType (FEMuS_DEFAULT);

  systemSG.init();

  // ******* Set Smoother *******
  //systemSG.SetSolverFineGrids (GMRES);
  systemSG.SetSolverFineGrids (RICHARDSON);
  
  systemSG.SetFieldSplitTree (&uSG);

  systemSG.SetPreconditionerFineGrids (ILU_PRECOND);

  systemSG.SetTolerances (1.e-20, 1.e-20, 1.e+50, 100);
  //END



//BEGIN solve eigenproblem to get the KL functions
  GetEigenPair (ml_prob, numberOfEigPairs, eigenvalues);   //solve the generalized eigenvalue problem and compute the eigenpairs

  for (int i = 0; i < numberOfEigPairs; i++) {
    std::cout << eigenvalues[i].first << " " << eigenvalues[i].second << std::endl;
  }

//END


//BEGIN solve SGM system
  systemSG.MGsolve();
//END


  //BEGIN post processing
  std::vector <double> alphas;
  GetCoefficientsForQuantityOfInterest (ml_probSG, alphas, domainMeasure);   //gets alpha for the QoI

  GetMomentsAndCumulants (alphas);   //computes moments and cumulants

  //PlotGCandEDExpansion();

  // ******* Print solution *******
  mlSol.SetWriter (VTK);
  std::vector<std::string> print_vars;
  print_vars.push_back ("All");
  //mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);
  //END

/*
  //BEGIN Define the instances of the problem for HISTOGRAM and KDE
  MultiLevelMesh mlMshHisto;
  MultiLevelMesh mlMshHistoFinest;

//     nxCoarseBox = static_cast<unsigned> ( floor ( 1. + 3.3 * log ( numberOfSamples ) ) );
//     nxCoarseBox = static_cast<unsigned> ( floor ( 1. + 2. * log2 ( numberOfSamples ) ) );
    nxCoarseBox = static_cast<unsigned> ( pow(2,kOrder) );
    nyCoarseBox = nxCoarseBox;
    nzCoarseBox = nxCoarseBox;

    mlMshHisto.GenerateCoarseBoxMesh ( nxCoarseBox, 0, 0, xMinCoarseBox, xMaxCoarseBox, 0., 0., 0., 0., EDGE3, "seventh" ); //for 1D
//     mlMshHisto.GenerateCoarseBoxMesh ( nxCoarseBox, nyCoarseBox, 0, xMinCoarseBox, xMaxCoarseBox, yMinCoarseBox, yMaxCoarseBox, 0., 0., QUAD9, "seventh" ); //for 2D
//     mlMshHisto.GenerateCoarseBoxMesh ( nxCoarseBox, nyCoarseBox, nzCoarseBox, xMinCoarseBox, xMaxCoarseBox, yMinCoarseBox, yMaxCoarseBox, zMinCoarseBox, zMaxCoarseBox, HEX27, "seventh" ); //for 3D

    mlMshHistoFinest.GenerateCoarseBoxMesh ( nxCoarseBoxFinest, 0, 0, xMinCoarseBox, xMaxCoarseBox, 0., 0., 0., 0., EDGE3, "seventh" ); //for 1D
//     mlMshHistoFinest.GenerateCoarseBoxMesh ( nxCoarseBoxFinest, nyCoarseBoxFinest, 0, xMinCoarseBox, xMaxCoarseBox, yMinCoarseBox, yMaxCoarseBox, 0., 0., QUAD9, "seventh" ); //for 2D
//     mlMshHistoFinest.GenerateCoarseBoxMesh ( nxCoarseBoxFinest, nyCoarseBoxFinest, nzCoarseBoxFinest, xMinCoarseBox, xMaxCoarseBox, yMinCoarseBox, yMaxCoarseBox, zMinCoarseBox, zMaxCoarseBox, HEX27, "seventh" ); //for 3D

  mlMshHisto.PrintInfo();

  unsigned dimCoarseBox = mlMshHisto.GetDimension();

  std::vector< std::vector <double > > sgmQoIStandardized;
  std::vector< std::vector <double > > sgmQoIStandardizedFinest;
  GetQoIStandardizedSamples (alphas, sgmQoIStandardized, sgmQoIStandardizedFinest, dimCoarseBox);


  MultiLevelSolution mlSolHisto (&mlMshHisto);
  MultiLevelSolution mlSolHistoFinest (&mlMshHistoFinest);

  mlSolHisto.AddSolution ("HISTO", DISCONTINUOUS_POLYNOMIAL, ZERO);
  mlSolHisto.AddSolution ("PROPOSED", LAGRANGE, FIRST);

  mlSolHistoFinest.AddSolution ("HISTO_F", DISCONTINUOUS_POLYNOMIAL, ZERO);

  mlSolHisto.Initialize ("All");

  mlSolHistoFinest.Initialize ("All");

  MultiLevelProblem ml_probHisto (&mlSolHisto);

  MultiLevelProblem ml_probHistoFinest (&mlSolHistoFinest);

  clock_t start_time = clock();

  GetHistogramAndKDE (sgmQoIStandardized, sgmQoIStandardizedFinest, ml_probHisto, ml_probHistoFinest);

  std::cout << std::endl << " RANNA in: " << std::setw (11) << std::setprecision (6) << std::fixed
            << static_cast<double> ( (clock() - start_time)) / CLOCKS_PER_SEC << " s" << std::endl;

//   GetKDEIntegral(ml_probHisto);

  GetAverageL2Error (sgmQoIStandardized, ml_probHisto, ml_probHistoFinest);

  mlSolHisto.SetWriter (VTK);
  std::vector<std::string> print_vars_2;
  print_vars_2.push_back ("All");
  //mlSolHisto.GetWriter()->SetDebugOutput(true);
  mlSolHisto.GetWriter()->Write (DEFAULT_OUTPUTDIR, "histo_and_proposed", print_vars_2, 0);


  mlSolHistoFinest.SetWriter (VTK);
  std::vector<std::string> print_vars_3;
  print_vars_3.push_back ("All");
  //mlSolHisto.GetWriter()->SetDebugOutput(true);
  mlSolHistoFinest.GetWriter()->Write (DEFAULT_OUTPUTDIR, "histo_finer", print_vars_3, 0);
  //END
  */
  
  for (unsigned i = 0; i < Jp.size(); i++) {
    delete FielduSGi[i];
  }
  delete [] FielduSGi;

  return 0;

} //end main

void GetEigenPair (MultiLevelProblem& ml_prob, const int& numberOfEigPairs, std::vector < std::pair<double, double> >& eigenvalues) {
//void GetEigenPair(MultiLevelProblem & ml_prob, Mat &CCSLEPc, Mat &MMSLEPc) {

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("UQ");   // pointer to the linear implicit system named "Poisson"

  unsigned level = numberOfUniformLevels - 1;

  double varianceInput = stdDeviationInput * stdDeviationInput;

  Mesh*                    msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             MM = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex ("u");   // get the position of "u" in the ml_sol object
  unsigned solType = mlSol->GetSolutionType (soluIndex);   // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex ("u");   // get the position of "u" in the pdeSys object

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector < vector < double > > x1 (dim);   // local coordinates
  vector < vector < double > > x2 (dim);   // local coordinates

  for (unsigned k = 0; k < dim; k++) {
    x1[k].reserve (maxSize);
    x2[k].reserve (maxSize);
  }

  vector <double> phi_x; // local test function first order partial derivatives

  phi_x.reserve (maxSize * dim);

  vector< int > l2GMap1; // local to global mapping
  vector< int > l2GMap2; // local to global mapping
  l2GMap1.reserve (maxSize);
  l2GMap2.reserve (maxSize);

  vector < double > MMlocal;
  MMlocal.reserve (maxSize * maxSize);

  vector < double > CClocal;
  CClocal.reserve (maxSize * maxSize);

  MM->zero(); // Set to zero all the entries of the Global Matrix

  int MM_size = msh->_dofOffset[solType][nprocs];
  int MM_local_size = msh->_dofOffset[solType][iproc + 1] - msh->_dofOffset[solType][iproc];

  SparseMatrix* CC;
  CC = SparseMatrix::build().release();
  CC->init (MM_size, MM_size, MM_local_size, MM_local_size, MM_local_size, MM_size - MM_local_size);
  CC->zero();

  for (int kproc = 0; kproc < nprocs; kproc++) {
    for (int jel = msh->_elementOffset[kproc]; jel < msh->_elementOffset[kproc + 1]; jel++) {

      short unsigned ielGeom2;
      unsigned nDof2;
      unsigned nDofx2;

      if (iproc == kproc) {
        ielGeom2 = msh->GetElementType (jel);
        nDof2  = msh->GetElementDofNumber (jel, solType);   // number of solution element dofs
        nDofx2 = msh->GetElementDofNumber (jel, xType);   // number of coordinate element dofs
      }

      MPI_Bcast (&ielGeom2, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD);
      MPI_Bcast (&nDof2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      MPI_Bcast (&nDofx2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);

      // resize local arrays
      l2GMap2.resize (nDof2);

      for (int k = 0; k < dim; k++) {
        x2[k].resize (nDofx2);
      }

      // local storage of global mapping and solution
      if (iproc == kproc) {
        for (unsigned j = 0; j < nDof2; j++) {
          l2GMap2[j] = pdeSys->GetSystemDof (soluIndex, soluPdeIndex, j, jel);   // global to global mapping between solution node and pdeSys dof
        }
      }

      MPI_Bcast (&l2GMap2[0], nDof2, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);

      // local storage of coordinates
      if (iproc == kproc) {
        for (unsigned j = 0; j < nDofx2; j++) {
          unsigned xDof  = msh->GetSolutionDof (j, jel, xType);   // global to global mapping between coordinates node and coordinate dof

          for (unsigned k = 0; k < dim; k++) {
            x2[k][j] = (*msh->_topology->_Sol[k]) (xDof);     // global extraction and local storage for the element coordinates
          }
        }
      }

      for (unsigned k = 0; k < dim; k++) {
        MPI_Bcast (& x2[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }

      unsigned jgNumber = msh->_finiteElement[ielGeom2][solType]->GetGaussPointNumber();
      vector < vector < double > > xg2 (jgNumber);
      vector <double> weight2 (jgNumber);
      vector < vector <double> > phi2 (jgNumber);   // local test function

      for (unsigned jg = 0; jg < jgNumber; jg++) {
        msh->_finiteElement[ielGeom2][solType]->Jacobian (x2, jg, weight2[jg], phi2[jg], phi_x);

        xg2[jg].assign (dim, 0.);

        for (unsigned j = 0; j < nDof2; j++) {
          for (unsigned k = 0; k < dim; k++) {
            xg2[jg][k] += x2[k][j] * phi2[jg][j];
          }
        }
      }

      // element loop: each process loops only on the elements that owns
      for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

        short unsigned ielGeom1 = msh->GetElementType (iel);
        unsigned nDof1  = msh->GetElementDofNumber (iel, solType);   // number of solution element dofs
        unsigned nDofx1 = msh->GetElementDofNumber (iel, xType);   // number of coordinate element dofs

        // resize local arrays
        l2GMap1.resize (nDof1);
        //std::vector<bool>bdcDirichlet(nDof1);

        for (int k = 0; k < dim; k++) {
          x1[k].resize (nDofx1);
        }

        // local storage of global mapping and solution
        for (unsigned i = 0; i < nDof1; i++) {
          l2GMap1[i] = pdeSys->GetSystemDof (soluIndex, soluPdeIndex, i, iel);   // global to global mapping between solution node and pdeSys dof
          //unsigned solDof = msh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
          //bdcDirichlet[i] = ( (*sol->_Bdc[soluIndex])(solDof) < 1.5)? false:false;
        }

        // local storage of coordinates
        for (unsigned i = 0; i < nDofx1; i++) {
          unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // global to global mapping between coordinates node and coordinate dof

          for (unsigned k = 0; k < dim; k++) {
            x1[k][i] = (*msh->_topology->_Sol[k]) (xDof);     // global extraction and local storage for the element coordinates
          }
        }

        if (iel == jel) {
          MMlocal.assign (nDof1 * nDof1, 0.);      //resize
        }

        CClocal.assign (nDof1 * nDof2, 0.);   //resize

        // *** Gauss point loop ***
        unsigned igNumber = msh->_finiteElement[ielGeom1][solType]->GetGaussPointNumber();
        double weight1;
        vector <double> phi1;  // local test function

        for (unsigned ig = 0; ig < igNumber; ig++) {

          msh->_finiteElement[ielGeom1][solType]->Jacobian (x1, ig, weight1, phi1, phi_x);

          // evaluate the solution, the solution derivatives and the coordinates in the gauss point
          vector < double > xg1 (dim, 0.);

          for (unsigned i = 0; i < nDof1; i++) {
            for (unsigned k = 0; k < dim; k++) {
              xg1[k] += x1[k][i] * phi1[i];
            }
          }

          if (iel == jel) {
            for (unsigned i = 0; i < nDof1; i++) {
              for (unsigned i1 = 0; i1 < nDof1; i1++) {
                MMlocal[ i * nDof1 + i1 ] += phi1[i] * phi1[i1] * weight1;
              }
            }
          }

          for (unsigned jg = 0; jg < jgNumber; jg++) {
            double dist = 0.;

            for (unsigned k = 0; k < dim; k++) {
              dist += fabs (xg1[k] - xg2[jg][k]);
            }

            double C = varianceInput * exp (- dist / L);

            for (unsigned i = 0; i < nDof1; i++) {
              for (unsigned j = 0; j < nDof2; j++) {
                CClocal[i * nDof2 + j] += weight1 * phi1[i] * C * phi2[jg][j] * weight2[jg];
              }//endl j loop
            } //endl i loop
          } //endl jg loop
        } //endl ig loop

        if (iel == jel) {
          MM->add_matrix_blocked (MMlocal, l2GMap1, l2GMap1);
        }

        CC->add_matrix_blocked (CClocal, l2GMap1, l2GMap2);
      } // end iel loop
    } //end jel loop
  } //end kproc loop

  MM->close();
  CC->close();

  //BEGIN solve the eigenvalue problem

  int ierr;
  EPS eps;
  PetscInt convergedSolns, numberOfIterations;

  ierr = EPSCreate (PETSC_COMM_WORLD, &eps);
  CHKERRABORT (MPI_COMM_WORLD, ierr);
  ierr = EPSSetOperators (eps, (static_cast<PetscMatrix*> (CC))->mat(), (static_cast<PetscMatrix*> (MM))->mat());
  CHKERRABORT (MPI_COMM_WORLD, ierr);
  ierr = EPSSetFromOptions (eps);
  CHKERRABORT (MPI_COMM_WORLD, ierr);

  //ierr = EPSSetDimensions(eps, numberOfEigPairs, 8 * numberOfEigPairs, 600);
  ierr = EPSSetDimensions (eps, numberOfEigPairs, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRABORT (MPI_COMM_WORLD, ierr);
  ierr = EPSSetWhichEigenpairs (eps, EPS_LARGEST_MAGNITUDE);
  CHKERRABORT (MPI_COMM_WORLD, ierr);

  //ierr = EPSSetTolerances(eps,1.0e-10,1000);
  CHKERRABORT (MPI_COMM_WORLD, ierr);

  ierr = EPSSolve (eps);
  CHKERRABORT (MPI_COMM_WORLD, ierr);

  //ierr = EPSView(eps, PETSC_VIEWER_STDOUT_SELF);

  std::cout << " -----------------------------------------------------------------" << std::endl;

  ierr = EPSGetConverged (eps, &convergedSolns);
  CHKERRABORT (MPI_COMM_WORLD, ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, " Number of converged eigenpairs: %D\n\n", convergedSolns);
  CHKERRABORT (MPI_COMM_WORLD, ierr);

  if (convergedSolns > 0) {

    for (unsigned i = 0; i < numberOfEigPairs; i++) {

      char name[10];
      sprintf (name, "egnf%d", i);
      soluIndex = mlSol->GetIndex (name);   // get the position of "u" in the ml_sol object

      // Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and ki (imaginary part)

      ierr = EPSGetEigenpair (eps, i, &eigenvalues[i].first, &eigenvalues[i].second, (static_cast<PetscVector*> (sol->_Sol[soluIndex]))->vec(), NULL);
      CHKERRABORT (MPI_COMM_WORLD, ierr);

    }
  }

  ierr = EPSDestroy (&eps);
  CHKERRABORT (MPI_COMM_WORLD, ierr);

  delete CC;

  //BEGIN GRAM SCHMIDT ORTHONORMALIZATION

  std::vector <unsigned> eigfIndex (numberOfEigPairs);
  char name[10];

  for (unsigned i = 0; i < numberOfEigPairs; i++) {
    sprintf (name, "egnf%d", i);
    eigfIndex[i] = mlSol->GetIndex (name);   // get the position of "u" in the ml_sol object
  }

  vector < double >  eigenFunction (numberOfEigPairs);   // local solution
  vector < double >  eigenFunctionOld (numberOfEigPairs);   // local solution

  std::vector < std::vector < double > > coeffsGS_local (numberOfEigPairs);
  std::vector < std::vector < double > > coeffsGS_global (numberOfEigPairs);

  for (unsigned i = 0; i < numberOfEigPairs; i++) {
    coeffsGS_local[i].assign (numberOfEigPairs, 0.);
    coeffsGS_global[i].assign (numberOfEigPairs, 0.);
  }

  for (unsigned iGS = 0; iGS < numberOfEigPairs; iGS++) {

    if (iGS > 0) {

      for (unsigned jGS = 0; jGS < iGS; jGS++) {

        //BEGIN COMPUTE coeffsGS LOCAL

        for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

          short unsigned ielGeom = msh->GetElementType (iel);
          unsigned nDofu  = msh->GetElementDofNumber (iel, solType);   // number of solution element dofs
          unsigned nDofx = msh->GetElementDofNumber (iel, xType);   // number of coordinate element dofs

          eigenFunction.resize (nDofu);
          eigenFunctionOld.resize (nDofu);

          for (int i = 0; i < dim; i++) {
            x1[i].resize (nDofx);
          }

          // local storage of global mapping and solution
          for (unsigned i = 0; i < nDofu; i++) {
            unsigned solDof = msh->GetSolutionDof (i, iel, solType);   // global to global mapping between solution node and solution dof
            eigenFunction[i] = (*sol->_Sol[eigfIndex[iGS]]) (solDof);
            eigenFunctionOld[i] = (*sol->_Sol[eigfIndex[jGS]]) (solDof);
          }

          // local storage of coordinates
          for (unsigned i = 0; i < nDofx; i++) {
            unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // global to global mapping between coordinates node and coordinate dof

            for (unsigned jdim = 0; jdim < dim; jdim++) {
              x1[jdim][i] = (*msh->_topology->_Sol[jdim]) (xDof);     // global extraction and local storage for the element coordinates
            }
          }

          double weight;
          vector <double> phi;  // local test function

          // *** Gauss point loop ***
          for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
            // *** get gauss point weight, test function and test function partial derivatives ***
            msh->_finiteElement[ielGeom][solType]->Jacobian (x1, ig, weight, phi, phi_x);
            double eigenFunction_gss = 0.;
            double eigenFunction_gss_old = 0.;

            for (unsigned i = 0; i < nDofu; i++) {
              eigenFunction_gss += phi[i] * eigenFunction[i];
              eigenFunction_gss_old += phi[i] * eigenFunctionOld[i];
            }

            coeffsGS_local[iGS][jGS] += eigenFunction_gss * eigenFunction_gss_old * weight;
          }
        }

        //END COMPUTE coeffsGS LOCAL

        MPI_Allreduce (&coeffsGS_local[iGS][jGS], &coeffsGS_global[iGS][jGS], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      }

      for (unsigned idof = msh->_dofOffset[solType][iproc]; idof < msh->_dofOffset[solType][iproc + 1]; idof++) {
        double sum = 0.;

        for (unsigned jGS = 0; jGS < iGS; jGS++) {
          sum += coeffsGS_global[iGS][jGS] * (*sol->_Sol[eigfIndex[jGS]]) (idof);
        }

        double valueToSet = (*sol->_Sol[eigfIndex[iGS]]) (idof) - sum;
        sol->_Sol[eigfIndex[iGS]]->set (idof, valueToSet);
      }

    }

    sol->_Sol[eigfIndex[iGS]]->close();

    double local_norm2 = 0.;

    for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

      short unsigned ielGeom = msh->GetElementType (iel);
      unsigned nDofu  = msh->GetElementDofNumber (iel, solType);   // number of solution element dofs
      unsigned nDofx = msh->GetElementDofNumber (iel, xType);   // number of coordinate element dofs

      eigenFunction.resize (nDofu);

      for (int i = 0; i < dim; i++) {
        x1[i].resize (nDofx);
      }

      // local storage of global mapping and solution
      for (unsigned i = 0; i < nDofu; i++) {
        unsigned solDof = msh->GetSolutionDof (i, iel, solType);   // global to global mapping between solution node and solution dof
        eigenFunction[i] = (*sol->_Sol[eigfIndex[iGS]]) (solDof);
      }

      // local storage of coordinates
      for (unsigned i = 0; i < nDofx; i++) {
        unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // global to global mapping between coordinates node and coordinate dof

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          x1[jdim][i] = (*msh->_topology->_Sol[jdim]) (xDof);     // global extraction and local storage for the element coordinates
        }
      }

      double weight;
      vector <double> phi;  // local test function

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[ielGeom][solType]->Jacobian (x1, ig, weight, phi, phi_x);
        double eigenFunction_gss = 0.;

        for (unsigned i = 0; i < nDofu; i++) {
          eigenFunction_gss += phi[i] * eigenFunction[i];
        }

        local_norm2 += eigenFunction_gss * eigenFunction_gss * weight;
      }
    }

    double norm2 = 0.;
    MPI_Allreduce (&local_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double norm = sqrt (norm2);
    std::cout << "norm = " << norm << std::endl;
    sol->_Sol[eigfIndex[iGS]]->scale (1. / norm);

    sol->_Sol[eigfIndex[iGS]]->close();

  }

  //END GRAM SCHMIDT ORTHONORMALIZATION

  //BEGIN GRAM SCHMIDT CHECK
  vector < double >  eigenFunctionCheck (numberOfEigPairs);   // local solution
  vector < double >  eigenFunctionOldCheck (numberOfEigPairs);   // local solution


  for (unsigned i1 = 0; i1 < numberOfEigPairs; i1++) {
    for (unsigned j1 = 0; j1 < numberOfEigPairs; j1++) {

      double integral = 0.;

      for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

        short unsigned ielGeom = msh->GetElementType (iel);
        unsigned nDofu  = msh->GetElementDofNumber (iel, solType);   // number of solution element dofs
        unsigned nDofx = msh->GetElementDofNumber (iel, xType);   // number of coordinate element dofs

        eigenFunctionCheck.resize (nDofu);
        eigenFunctionOldCheck.resize (nDofu);

        for (int i = 0; i < dim; i++) {
          x1[i].resize (nDofx);
        }

        // local storage of global mapping and solution
        for (unsigned i = 0; i < nDofu; i++) {
          unsigned solDof = msh->GetSolutionDof (i, iel, solType);   // global to global mapping between solution node and solution dof
          eigenFunctionCheck[i] = (*sol->_Sol[eigfIndex[i1]]) (solDof);
          eigenFunctionOldCheck[i] = (*sol->_Sol[eigfIndex[j1]]) (solDof);
        }

        // local storage of coordinates
        for (unsigned i = 0; i < nDofx; i++) {
          unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // global to global mapping between coordinates node and coordinate dof

          for (unsigned jdim = 0; jdim < dim; jdim++) {
            x1[jdim][i] = (*msh->_topology->_Sol[jdim]) (xDof);     // global extraction and local storage for the element coordinates
          }
        }

        double weight;
        vector <double> phi;  // local test function

        //  *** Gauss point loop ***
        for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
          //    *** get gauss point weight, test function and test function partial derivatives ***
          msh->_finiteElement[ielGeom][solType]->Jacobian (x1, ig, weight, phi, phi_x);
          double eigenFunction_gss = 0.;
          double eigenFunction_gss_old = 0.;

          for (unsigned i = 0; i < nDofu; i++) {
            eigenFunction_gss += phi[i] * eigenFunctionCheck[i];
            eigenFunction_gss_old += phi[i] * eigenFunctionOldCheck[i];
          }

          integral += eigenFunction_gss * eigenFunction_gss_old * weight;
        }
      }

      double globalIntegral = 0.;
      MPI_Allreduce (&integral, &globalIntegral, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      std::cout << "i = " << i1 << " , " << "j = " << j1 << " , " << "integral = " << globalIntegral << std::endl;
    }
  }

  //END GRAM SCHMIDT CHECK

}
//
//
void GetCoefficientsForQuantityOfInterest (MultiLevelProblem& ml_prob, std::vector <double > &  alphas, const double& domainMeasure) {

  //  extract pointers to the several objects that we are going to use

  //uq &myuq = FemusInit::_uqHermite;

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("SG");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh* msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  const std::vector < std::vector <unsigned> > &Jp = myuq.GetIndexSet (pIndex, numberOfEigPairs);

  std::vector <double > alphasTemp (Jp.size(), 0.);
  alphas.resize (Jp.size());

  //solution Index
  std::vector <unsigned> soluIndex (Jp.size());

  for (unsigned i = 0; i < Jp.size(); i++) {
    char name[10];
    sprintf (name, "uSG%d", i);
    soluIndex[i] = mlSol->GetIndex (name);   // get the position of "u" in the ml_sol object
  }

  unsigned soluType = mlSol->GetSolutionType (soluIndex[0]);

  vector < vector < double > >  solu (Jp.size());   // local solution

  vector < vector < double > > x (dim);   // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  phi.reserve (maxSize);
  phi_x.reserve (maxSize * dim);

  // element loop: each process loops only on the elements that owns

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofu  = msh->GetElementDofNumber (iel, soluType);   // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber (iel, xType);   // number of coordinate element dofs

    for (unsigned j = 0; j < Jp.size(); j++) {
      solu[j].resize (nDofu);
    }

    for (int i = 0; i < dim; i++) {
      x[i].resize (nDofx);
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, soluType);   // global to global mapping between solution node and solution dof

      for (unsigned j = 0; j < Jp.size(); j++) {
        solu[j][i] = (*sol->_Sol[soluIndex[j]]) (solDof);     // global extraction and local storage for the solution
      }
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim]) (xDof);     // global extraction and local storage for the element coordinates
      }
    }


    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian (x, ig, weight, phi, phi_x);

      for (unsigned j = 0; j < Jp.size(); j++) {
        double solu_gss = 0.;

        for (unsigned i = 0; i < nDofu; i++) {
          solu_gss += phi[i] * solu[j][i];
        }

//          alphasTemp[j] += solu_gss * solu_gss * weight ; // this is similar to the integral of the square.
        alphasTemp[j] +=  solu_gss *  weight / domainMeasure; // this is the spatial average over the domain.
      }
    } // end gauss point loop

  } //end element loop for each process


  for (unsigned j = 0; j < Jp.size(); j++) {
    alphas[j] = 0.;
    MPI_Allreduce (&alphasTemp[j], &alphas[j], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

//     for ( unsigned j = 0; j < Jp.size(); j++ ) {
//         std::cout << " alpha = " << std::setprecision ( 14 ) << alphas[j] << std::endl;
//     }

}
//
void GetMomentsAndCumulants (std::vector <double>& alphas) {

  //let's standardize the quantity of interest after finding moments and standard deviation

  //uq &myuq = FemusInit::_uqHermite;

  if (totMoments <= 0) {

    std::cout << "ERROR: total number of moments has to be a positive integer" << std::endl;

  }

  else {

    unsigned desiredQuadraturePoints = static_cast<double> (ceil ( (totMoments * pIndex + 1) * 0.5));

    unsigned numberOfQuadraturePoints = (desiredQuadraturePoints <= 16) ? desiredQuadraturePoints : 16;

    if (desiredQuadraturePoints > 16) {
      std::cout <<
                "------------------------------- WARNING: less quadrature points than needed were employed in function GetMomentsAndCumulants -------------------------------"
                << std::endl;
      std::cout << " Needed : " << desiredQuadraturePoints << " , " << " Used : " << 16 << std::endl;
    }

    const std::vector < std::vector <unsigned> > &Tp = myuq.GetTensorProductSet (numberOfQuadraturePoints, numberOfEigPairs);

    const std::vector < std::vector <unsigned> > &Jp = myuq.GetIndexSet (pIndex, numberOfEigPairs);

    const std::vector < std::vector < double > >  & multivariatePoly = myuq.GetMultivariatePolynomial (numberOfQuadraturePoints, pIndex, numberOfEigPairs);
    const std::vector < double > & multivariateQuadratureWeights = myuq.GetMultivariateWeights (numberOfQuadraturePoints, pIndex, numberOfEigPairs);


    //BEGIN computation of the raw moments
    for (unsigned p = 0; p < totMoments; p++) {
      moments[p] = 0.;

      for (unsigned j = 0; j < Tp.size(); j++) {
        double integrandFunction = 0.;

        for (unsigned i = 0; i < Jp.size(); i++) {
          integrandFunction += multivariatePoly[i][j] * alphas[i];
        }

        integrandFunction = pow (integrandFunction, p + 1);
        moments[p] += multivariateQuadratureWeights[j] * integrandFunction;
      }
    }

    //END

    //BEGIN computation of the mean of QoI (AKA first moment)
    meanQoI = moments[0];
    //END


    //BEGIN computation of the variance and standard deviation of QoI
//     varianceQoI = moments[1] - meanQoI * meanQoI;

    varianceQoI = 0;

    for (unsigned j = 0; j < Tp.size(); j++) {
      double integrandFunctionVariance = 0.;

      for (unsigned i = 0; i < Jp.size(); i++) {
        integrandFunctionVariance += multivariatePoly[i][j] * alphas[i];
      }

      integrandFunctionVariance = pow (integrandFunctionVariance - meanQoI, 2);
      varianceQoI += multivariateQuadratureWeights[j] * integrandFunctionVariance;
    }

    stdDeviationQoI = sqrt (varianceQoI);
    //END


    //BEGIN computation of the raw moments of the standardized variable
    for (unsigned p = 0; p < totMoments; p++) {
      momentsStandardized[p] = 0.;

      for (unsigned j = 0; j < Tp.size(); j++) {
        double integrandFunction = 0.;

        for (unsigned i = 0; i < Jp.size(); i++) {
          integrandFunction += multivariatePoly[i][j] * alphas[i];
        }

        integrandFunction = (integrandFunction - meanQoI) / stdDeviationQoI;   //standardization of the QoI
        integrandFunction = pow (integrandFunction, p + 1);
        momentsStandardized[p] += multivariateQuadratureWeights[j] * integrandFunction;
      }
    }

    //END


    //BEGIN computation of the CUMULANTS
    cumulants[0] = moments[0];
    cumulantsStandardized[0] = momentsStandardized[0];

    if (totMoments > 1) {
      cumulants[1] = moments[1] - moments[0] * moments[0];

      cumulantsStandardized[1] = momentsStandardized[1] - momentsStandardized[0] * momentsStandardized[0];

//       std::cout.precision(14);
//       std::cout << "AAAAAAAAAAAAAAA" << cumulants[1] << std::endl;
      if (totMoments > 2) {
        cumulants[2] = moments[2] - 3. * moments[1] * moments[0] + 2. * pow (moments[0], 3);

        cumulantsStandardized[2] = momentsStandardized[2] - 3. * momentsStandardized[1] * momentsStandardized[0] + 2. * pow (momentsStandardized[0], 3);

        if (totMoments > 3) {
          cumulants[3] = moments[3] - 4. * moments[2] * moments[0] - 3. * moments[1] * moments[1] + 12. * moments[1] * moments[0] * moments[0] - 6. * pow (moments[0], 4);

          cumulantsStandardized[3] = momentsStandardized[3] - 4. * momentsStandardized[2] * momentsStandardized[0] - 3. * momentsStandardized[1] * momentsStandardized[1]
                                     + 12. * momentsStandardized[1] * momentsStandardized[0] * momentsStandardized[0] - 6. * pow (momentsStandardized[0], 4);

          if (totMoments > 4) {
            cumulants[4] = moments[4] - 5. * moments[3] * moments[0] - 10. * moments[2] * moments[1] + 20. * moments[2] * moments[0] * moments[0]
                           + 30. * moments[1] * moments[1] * moments[0] - 60. * moments[1] * pow (moments[0], 3) + 24. * pow (moments[0], 5);

            cumulantsStandardized[4] = momentsStandardized[4] - 5. * momentsStandardized[3] * momentsStandardized[0] - 10. * momentsStandardized[2] * momentsStandardized[1]
                                       + 20. * momentsStandardized[2] * momentsStandardized[0] * momentsStandardized[0]
                                       + 30. * momentsStandardized[1] * momentsStandardized[1] * momentsStandardized[0]
                                       - 60. * momentsStandardized[1] * pow (momentsStandardized[0], 3) + 24. * pow (momentsStandardized[0], 5);

            if (totMoments > 5) {
              cumulants[5] = moments[5] - 6. * moments[4] * moments[0] - 15. * moments[3] * moments[1] + 30. * moments[3] * moments[0] * moments[0]
                             - 10. * moments[2] * moments[2] + 120. * moments[2] * moments[1] * moments[0] - 120. * moments[2] * pow (moments[0], 3)
                             + 30. * pow (moments[1], 3) - 270. * pow (moments[1], 2) * pow (moments[0], 2) + 360. * moments[1] * pow (moments[0], 4) - 120. * pow (moments[0], 6);

              cumulantsStandardized[5] = momentsStandardized[5] - 6. * momentsStandardized[4] * momentsStandardized[0] - 15. * momentsStandardized[3] * momentsStandardized[1]
                                         + 30. * momentsStandardized[3] * momentsStandardized[0] * momentsStandardized[0]
                                         - 10. * momentsStandardized[2] * momentsStandardized[2] + 120. * momentsStandardized[2] * momentsStandardized[1] * momentsStandardized[0]
                                         - 120. * momentsStandardized[2] * pow (momentsStandardized[0], 3) + 30. * pow (momentsStandardized[1], 3)
                                         - 270. * pow (momentsStandardized[1], 2) * pow (momentsStandardized[0], 2) + 360. * momentsStandardized[1] * pow (momentsStandardized[0], 4)
                                         - 120. * pow (momentsStandardized[0], 6);
            }
          }
        }
      }
    }

    //END


    //BEGIN Plot moments and cumulants

    std::cout.precision (14);
    std::cout << " the mean is " << meanQoI << std::endl;
    std::cout << " the standard deviation is " << stdDeviationQoI << std::endl;
    std::cout << " the variance is " << varianceQoI << std::endl;

    std::cout << "Standardized Moments" << std::endl;

    for (unsigned p = 0; p < totMoments; p++) {
      std::cout << " & " << momentsStandardized[p] << "  ";
    }

    std::cout << std::endl;
    std::cout << "Standardized Cumulants" << std::endl;

    for (unsigned p = 0; p < totMoments; p++) {
      std::cout << " & " << cumulantsStandardized[p] << "  ";
    }

    std::cout << std::endl;
    std::cout << " Moments " << std::endl;

    for (unsigned p = 0; p < totMoments; p++) {
      std::cout << " & " << moments[p] << "  ";
    }

    std::cout << std::endl;
    std::cout << " Cumulants " << std::endl;

    for (unsigned p = 0; p < totMoments; p++) {
      std::cout << " & " << cumulants[p] << "  ";
    }

    std::cout << std::endl;
    std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;

    //END

  }
}
//
//
void PlotGCandEDExpansion() {

  std::cout.precision (14);
  std::cout << " the mean is " << meanQoI << std::endl;
  std::cout << " the standard deviation is " << stdDeviationQoI << std::endl;
  std::cout << " the variance is " << varianceQoI << std::endl;

  std::cout << "Standardized Moments" << std::endl;

  for (unsigned p = 0; p < totMoments; p++) {
    std::cout << " & " << momentsStandardized[p] << "  ";
  }

  std::cout << std::endl;
  std::cout << "Standardized Cumulants" << std::endl;

  for (unsigned p = 0; p < totMoments; p++) {
    std::cout << " & " << cumulantsStandardized[p] << "  ";
  }

  std::cout << std::endl;
  std::cout << " Moments " << std::endl;

  for (unsigned p = 0; p < totMoments; p++) {
    std::cout << " & " << moments[p] << "  ";
  }

  std::cout << std::endl;
  std::cout << " Cumulants " << std::endl;

  for (unsigned p = 0; p < totMoments; p++) {
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

  double t = - 5.5;
  double dt = (11.) / 300.;

//   cumulants[0] = 0; //decomment for nonStdGaussian

  //BEGIN GRAM CHARLIER PRINT
  std::cout << " ------------------------- GRAM CHARLIER ------------------------- " << std::endl;

  for (unsigned i = 0; i <= 300; i++) {
    std::cout << t << " ";
//     double t = x - meanQoI; //decomment for nonStdGaussian
    double gaussian = 1. / (sqrt (2 * acos (- 1))) * exp (- 0.5 * (t * t)) ;
    std::cout << gaussian << " ";

    d1gaussian = (- 1.) * gaussian * t ;

    generalizedGC1Term = gaussian - cumulantsStandardized[0] * d1gaussian;

    std::cout << generalizedGC1Term << " ";

    if (totMoments > 1) {

      d2gaussian = (1.) * gaussian * (t * t - 1.) ;
      d3gaussian = (- 1.) * gaussian * (t * t * t - 3. * t) ;

      generalizedGC2Terms = generalizedGC1Term + 0.5 * ( (cumulantsStandardized[1] - 1.) + pow (cumulantsStandardized[0], 2)) * d2gaussian ;

      std::cout << generalizedGC2Terms << " ";

      if (totMoments > 2) {

        d4gaussian = (1.) * gaussian * (t * t * t * t - 6. * t * t + 3.) ;
        d6gaussian = (1.) * gaussian * (pow (t, 6) - 15 * pow (t, 4) + 45 * t * t - 15);

        generalizedGC3Terms = generalizedGC2Terms - 1. / 6 * (cumulantsStandardized[2] + 3 * (cumulantsStandardized[1] - 1.) * cumulantsStandardized[0]
                                                              + pow (cumulantsStandardized[0], 3)) * d3gaussian;

        std::cout << generalizedGC3Terms << " ";

        if (totMoments > 3) {

          d5gaussian = (- 1.) * gaussian * (pow (t, 5) - 10. * t * t * t + 15. * t);
          d7gaussian = (- 1.) * gaussian * (pow (t, 7) - 21. * pow (t, 5) + 105. * t * t * t -  105. * t) ;
          d9gaussian = (- 1.) * gaussian * (pow (t, 9) - 36. * pow (t, 7) + 378. * pow (t, 5) - 1260. * t * t * t + 945. * t) ;

          generalizedGC4Terms = generalizedGC3Terms + 1. / 24 * (cumulantsStandardized[3] + 4. * cumulantsStandardized[2] * cumulantsStandardized[0]
                                                                 + 3. * pow ( (cumulantsStandardized[1] - 1.), 2) + 6. * (cumulantsStandardized[1] - 1.) * pow (cumulantsStandardized[0], 2)
                                                                 + pow (cumulantsStandardized[0], 4)) * d4gaussian;

          std::cout << generalizedGC4Terms << " ";

          if (totMoments > 4) {

            generalizedGC5Terms = generalizedGC4Terms - 1. / 120 * (cumulantsStandardized[4] + 5. * cumulantsStandardized[3] * cumulantsStandardized[0]
                                                                    + 10. * cumulantsStandardized[2] * (cumulantsStandardized[1] - 1.) + 10. * cumulantsStandardized[2] * pow (cumulantsStandardized[0], 2)
                                                                    + 15. * pow ( (cumulantsStandardized[1] - 1.), 2) * cumulantsStandardized[0] + 10. * (cumulantsStandardized[1] - 1.) * pow (cumulantsStandardized[0], 3)
                                                                    + pow (cumulantsStandardized[0], 5)) * d5gaussian;

            std::cout << generalizedGC5Terms << " ";

            if (totMoments > 5) {

              generalizedGC6Terms = generalizedGC5Terms + 1. / 720 * (cumulantsStandardized[5] + 6. * cumulantsStandardized[4] * cumulantsStandardized[0]
                                                                      + 15. * cumulantsStandardized[3] * (cumulantsStandardized[1] - 1.) + 15. * cumulantsStandardized[3] * pow (cumulantsStandardized[0], 2)
                                                                      + 10. * pow (cumulantsStandardized[2], 2) + 60. * cumulantsStandardized[2] * (cumulantsStandardized[1] - 1.) * cumulantsStandardized[0]
                                                                      + 20. * cumulantsStandardized[2] * pow (cumulantsStandardized[0], 3) + 15. * pow ( (cumulantsStandardized[1] - 1.), 3)
                                                                      + 45. * pow ( (cumulantsStandardized[1] - 1.), 2) * pow (cumulantsStandardized[0], 2) + 15. * (cumulantsStandardized[1] - 1.) * pow (cumulantsStandardized[0], 4)
                                                                      +  pow (cumulantsStandardized[0], 6)) * d6gaussian;

              std::cout << generalizedGC6Terms << " \n ";

            }

          }
        }
      }
    }

    t += dt;
  }

  t = -  5.5;
  dt = (11.) / 300.;

  //BEGIN EDGEWORTH PRINT
  std::cout << " ------------------------- EDGEWORTH ------------------------- " << std::endl;

  for (unsigned i = 0; i <= 300; i++) {
    std::cout << t << " ";
//     double t = x - meanQoI; //decomment for nonStdGaussian
    double gaussian = 1. / (sqrt (2 * acos (- 1))) * exp (- 0.5 * (t * t)) ;
    std::cout << gaussian << " ";

    d3gaussian = (- 1.) * gaussian * (t * t * t - 3. * t) ;

    lambda3 = cumulants[2] / pow (stdDeviationQoI, 3);

    edgeworth1Term = gaussian - lambda3 / 6. * d3gaussian;

    std::cout << edgeworth1Term << " ";

    if (totMoments > 1) {

      d4gaussian = (1.) * gaussian * (t * t * t * t - 6. * t * t + 3.) ;
      d6gaussian = (1.) * gaussian * (pow (t, 6) - 15 * pow (t, 4) + 45 * t * t - 15);

      lambda4 = cumulants[3] / pow (stdDeviationQoI, 4);

      edgeworth2Terms = edgeworth1Term + lambda4 / 24. * d4gaussian + lambda3 * lambda3 / 72. * d6gaussian;

      std::cout << edgeworth2Terms << " ";

      if (totMoments > 2) {

        d5gaussian = (- 1.) * gaussian * (pow (t, 5) - 10. * t * t * t + 15. * t);
        d7gaussian = (- 1.) * gaussian * (pow (t, 7) - 21. * pow (t, 5) + 105. * t * t * t -  105. * t) ;
        d9gaussian = (- 1.) * gaussian * (pow (t, 9) - 36. * pow (t, 7) + 378. * pow (t, 5) - 1260. * t * t * t + 945. * t) ;

        lambda5 = cumulants[4] / pow (stdDeviationQoI, 5);

        edgeworth3Terms = edgeworth2Terms - lambda5 / 120. * d5gaussian + lambda3 * lambda4 / 144. * d7gaussian + pow (lambda3, 3) / 1296. * d9gaussian;

        std::cout << edgeworth3Terms << " ";

      }

      if (totMoments > 3) {

        d8gaussian = (1.) * gaussian * (1. / 16. * (1680. - 6720. * t * t + 3360. * pow (t, 4) - 448. * pow (t, 6) + 16. * pow (t, 8)));
        d10gaussian = (1.) * gaussian * (1. / 32. * (- 30240. + 151200. *  t * t - 100800. * pow (t, 4) + 20160. * pow (t, 6) - 1440. * pow (t, 8)
                                                     + 32. * pow (t, 10)));
        d12gaussian = (1.) * gaussian * (1. / 64. * (665280. - 3991680. * t * t + 3326400. * pow (t, 4) - 887040. * pow (t, 6) + 95040. * pow (t, 8)
                                                     - 4224. * pow (t, 10) + 64. * pow (t, 12)));

        lambda6 = cumulants[5] / pow (stdDeviationQoI, 6);

        edgeworth4Terms = edgeworth3Terms + 1. / 720. * lambda6 * d6gaussian + (1. / 1152. * lambda4 * lambda4 + 1. / 720. * lambda3 * lambda5) * d8gaussian
                          + 1. / 1728. * lambda3 * lambda3 * lambda4 * d10gaussian + 1. / 31104. * pow (lambda3, 4) * d12gaussian;

        std::cout << edgeworth4Terms << " \n ";

      }
    }

    t += dt;
  }

  //END


}


void GetQoIStandardizedSamples (std::vector< double >& alphas, std::vector< std::vector <double > > & sgmQoIStandardized, std::vector< std::vector <double > > &                               sgmQoIStandardizedFinest, const unsigned & dimCoarseBox) {

  //uq &myuq = FemusInit::_uqHermite;
  const std::vector < std::vector <unsigned> > &Jp = myuq.GetIndexSet (pIndex, numberOfEigPairs);

  //FOR STD GAUSSIAN SAMPLING
  boost::mt19937 rng;
  boost::normal_distribution<> nd (0., 1.);
  boost::variate_generator < boost::mt19937&,
        boost::normal_distribution<> > var_nor (rng, nd);

  //FOR UNIFORM DISTRIBUTION
  boost::mt19937 rng1; // I don't seed it on purpouse (it's not relevant)
  boost::random::uniform_real_distribution<> un (- 1., 1.);
  boost::variate_generator < boost::mt19937&, boost::random::uniform_real_distribution<> > var_unif (rng1, un);

  //FOR LAPLACE DISTRIBUTION
  boost::mt19937 rng2; // I don't seed it on purpouse (it's not relevant)
  boost::random::uniform_real_distribution<> un1 (- 0.5, 0.49999999999);
  boost::variate_generator < boost::mt19937&,
        boost::random::uniform_real_distribution<> > var_unif1 (rng2, un1);


  sgmQoIStandardized.resize (numberOfSamples);

  if (histoFinest) {
    sgmQoIStandardizedFinest.resize (numberOfSamplesFinest);
  }

  for (unsigned m = 0; m < sgmQoIStandardized.size(); m++) {

    sgmQoIStandardized[m].resize (dimCoarseBox);

    for (unsigned idim = 0; idim < dimCoarseBox; idim++) {

      double sgmQoI = 0;

      std::vector<double> samplePoints (numberOfEigPairs, 0.);

      for (unsigned k = 0; k < numberOfEigPairs; k++) {

        if (myuq.GetQuadratureType() == UQ_HERMITE) {
          samplePoints[k] = var_nor();
        }

        else if (myuq.GetQuadratureType() == UQ_LEGENDRE) {
          samplePoints[k] = var_unif();
        }

      }

      const std::vector < std::vector <double> > &polyHistogram =
        myuq.GetPolynomialHistogram (pIndex, samplePoints, numberOfEigPairs);

      std::vector<double> MultivariatePolyHistogram (Jp.size(), 1.);

      for (unsigned i = 0; i < Jp.size(); i++) {
        for (unsigned k = 0; k < numberOfEigPairs; k++) {
          MultivariatePolyHistogram[i] *= polyHistogram[Jp[i][k]][k];
        }

        sgmQoI += alphas[i] * MultivariatePolyHistogram[i]; //TODO with QoIs that are different from each other, alphas[i] will be alphas[idim][i]
      }

            sgmQoIStandardized[m][idim] = ( sgmQoI - meanQoI ) / stdDeviationQoI; //TODO with QoIs that are different from each other, meanQoI and stdDeviationQoI will depend on idim

//             double normalSample = var_nor();
//             sgmQoIStandardized[m][idim] = normalSample;

//             double uniformSample = var_unif();
//             sgmQoIStandardized[m][idim] = uniformSample;

      //mixed input
//             if ( idim == 0 ) {
//                 double U = var_unif1();
//                 double signU = 0.;
//
//                 if ( U < 0 ) {
//                     signU = - 1.;
//                 }
//
//                 else if ( U > 0 ) {
//                     signU = 1.;
//                 }
//
//                 sgmQoIStandardized[m][idim] = muLaplace - bLaplace * signU * log ( 1. - 2. * fabs ( U ) ) ;
//             }
//
//             else if ( idim == 1 ) {
//                 double normalSample = var_nor();
//                 sgmQoIStandardized[m][idim] = normalSample;
//             }


    }
  }

  //END

  myuq.ClearPolynomialHistogram ();

  if (histoFinest) {

    for (unsigned m = 0; m < sgmQoIStandardizedFinest.size(); m++) {

      sgmQoIStandardizedFinest[m].resize (dimCoarseBox);

      for (unsigned idim = 0; idim < dimCoarseBox; idim++) {

        double sgmQoI = 0;

        std::vector<double> samplePoints (numberOfEigPairs, 0.);

        for (unsigned k = 0; k < numberOfEigPairs; k++) {
          if (myuq.GetQuadratureType() == UQ_HERMITE) {
            samplePoints[k] = var_nor();
          }

          else if (myuq.GetQuadratureType() == UQ_LEGENDRE) {
            samplePoints[k] = var_unif();
          }
        }

        const std::vector < std::vector <double> > &polyHistogram =
          myuq.GetPolynomialHistogram (pIndex, samplePoints, numberOfEigPairs);

        std::vector<double> MultivariatePolyHistogram (Jp.size(), 1.);

        for (unsigned i = 0; i < Jp.size(); i++) {
          for (unsigned k = 0; k < numberOfEigPairs; k++) {
            MultivariatePolyHistogram[i] *= polyHistogram[Jp[i][k]][k];
          }

          sgmQoI += alphas[i] * MultivariatePolyHistogram[i]; //TODO with QoIs that are different from each other, alphas[i] will be alphas[idim][i]
        }

        sgmQoIStandardizedFinest[m][idim] = (sgmQoI - meanQoI) / stdDeviationQoI;   //TODO with QoIs that are different from each other, meanQoI and stdDeviationQoI will depend on idim

      }
    }

    //END

    myuq.ClearPolynomialHistogram ();


  }


}

void GetHistogramAndKDE (std::vector< std::vector <double > > & sgmQoIStandardized, std::vector< std::vector <double > > & sgmQoIStandardizedFinest, MultiLevelProblem& ml_prob, MultiLevelProblem& ml_probFinest) {

  unsigned level = 0.;

  Mesh*                    msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  char name[10];
  sprintf (name, "HISTO");
  double solIndexHISTO = mlSol->GetIndex (name);   // get the position of "Ti" in the sol object
  sprintf (name, "PROPOSED");
  double solIndexKDE = mlSol->GetIndex (name);   // get the position of "Ti" in the sol object

  unsigned solTypeHISTO = mlSol->GetSolutionType (solIndexHISTO);
  unsigned solTypeKDE = mlSol->GetSolutionType (solIndexKDE);

  vector < double > phi;
  vector < double> gradphi;
  double weight;

  double dx = (xMaxCoarseBox - xMinCoarseBox) / nxCoarseBox;   //mesh size assuming a coarse box is used
  double dy = (dim > 1) ? (yMaxCoarseBox - yMinCoarseBox) / nyCoarseBox : 1. ;
  double dz = (dim > 2) ? (zMaxCoarseBox - zMinCoarseBox) / nzCoarseBox : 1. ;

  double measure = numberOfSamples * dx * dy * dz;

  for (unsigned m = 0; m < numberOfSamples; m++) {

    if (dim == 1) {

      for (int iel = sol->GetMesh()->_elementOffset[iproc]; iel < sol->GetMesh()->_elementOffset[iproc + 1]; iel ++) {

        unsigned  xLeftDof = sol->GetMesh()->GetSolutionDof (0, iel, 2);
        unsigned  xRightDof = sol->GetMesh()->GetSolutionDof (1, iel, 2);

        double xLeft = (*sol->GetMesh()->_topology->_Sol[0]) (xLeftDof);
        double xRight = (*sol->GetMesh()->_topology->_Sol[0]) (xRightDof);

        if (sgmQoIStandardized[m][0] > xLeft && sgmQoIStandardized[m][0] <= xRight) {

          //BEGIN write HISTO solution
          double histoValue = 1. / measure;
          sol->_Sol[solIndexHISTO]->add (iel, histoValue);
          //END

          //BEGIN write KDE solution
          short unsigned ielType = msh->GetElementType (iel);
          unsigned nDofsKDE = msh->GetElementDofNumber (iel, solTypeKDE);

          std::vector < std::vector < double> > vx (1);
          vx[0].resize (nDofsKDE);
          vx[0][0] = xLeft;
          vx[0][1] = xRight;

          std::vector < double> sampleLocal (1, 0.);
          sampleLocal[0] = - 1. + 2. * (sgmQoIStandardized[m][0] - xLeft) / (xRight - xLeft);

          msh->_finiteElement[ielType][solTypeKDE]->Jacobian (vx, sampleLocal, weight, phi, gradphi);

          for (unsigned inode = 0; inode < nDofsKDE; inode++) {
            unsigned globalDof = msh->GetSolutionDof (inode, iel, solTypeKDE);
            double KDEvalue = phi[inode] / (numberOfSamples * dx);   //note: numberOfSamples * h is an approximation of the area under the histogram, aka the integral
            sol->_Sol[solIndexKDE]->add (globalDof, KDEvalue);
          }

          //END

          break;
        }
      }

    }

    else {

      Marker marker (sgmQoIStandardized[m], 0., VOLUME, mlSol->GetLevel (level), 2, true);
      unsigned iel = marker.GetMarkerElement();
      std::vector<double> sampleLocal;
      marker.GetMarkerLocalCoordinates (sampleLocal);

      if (iel >= sol->GetMesh()->_elementOffset[iproc]  &&  iel < sol->GetMesh()->_elementOffset[iproc + 1]) {

        //BEGIN write HISTO solution
        double histoValue = 1. / measure;
        sol->_Sol[solIndexHISTO]->add (iel, histoValue);
        //END

        //BEGIN write KDE solution
        short unsigned ielType = msh->GetElementType (iel);
        unsigned nDofsKDE = msh->GetElementDofNumber (iel, solTypeKDE);

        std::vector < std::vector < double> > vx (dim);

        for (int idim = 0; idim < dim; idim++) {
          vx[idim].resize (nDofsKDE);
        }

        for (unsigned inode = 0; inode < nDofsKDE; inode++) {
          unsigned idofVx = msh->GetSolutionDof (inode, iel, 2);

          for (int jdim = 0; jdim < dim; jdim++) {
            vx[jdim][inode] = (*msh->_topology->_Sol[jdim]) (idofVx);
          }
        }

        msh->_finiteElement[ielType][solTypeKDE]->Jacobian (vx, sampleLocal, weight, phi, gradphi);

        for (unsigned inode = 0; inode < nDofsKDE; inode++) {
          unsigned globalDof = msh->GetSolutionDof (inode, iel, solTypeKDE);
          double KDEvalue = phi[inode] / (numberOfSamples * dx * dy * dz);
          sol->_Sol[solIndexKDE]->add (globalDof, KDEvalue);
        }

        //END

      }

    }

    sol->_Sol[solIndexHISTO]->close();
    sol->_Sol[solIndexKDE]->close();
  }



//     double integralLocal = 0;
//     double integral;
//     for(unsigned i =  msh->_dofOffset[solTypeHISTO][iproc]; i <  msh->_dofOffset[solTypeHISTO][iproc + 1]; i++) {
//         integralLocal += (*sol->_Sol[solIndexHISTO])(i) * dx * dy * dz; // this is assuming the mesh is a coarse box
//     }
//
//     MPI_Allreduce(&integralLocal, &integral, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//
//     std::cout << "integral 1 = " << integral << " integral 2 = " << numberOfSamples * dx *dy *dz << std::endl;
//
//     for(unsigned i =  msh->_dofOffset[solTypeHISTO][iproc]; i <  msh->_dofOffset[solTypeHISTO][iproc + 1]; i++) {
//         double valueHISTO = (*sol->_Sol[solIndexHISTO])(i);
//         sol->_Sol[solIndexHISTO]->set(i, valueHISTO / integral);
//     }
//     sol->_Sol[solIndexHISTO]->close();


  if (histoFinest) {

//BEGIN computation of histo finest

    Mesh*                    mshFinest = ml_probFinest._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
    MultiLevelSolution*    mlSolFinest = ml_probFinest._ml_sol;  // pointer to the multilevel solution object
    Solution*                solFinest = ml_probFinest._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object
    unsigned    iprocFinest = mshFinest->processor_id(); // get the process_id (for parallel computation)

    sprintf (name, "HISTO_F");
    double solIndexHISTOF = mlSolFinest->GetIndex (name);   // get the position of "Ti" in the sol object


    double dxF = (xMaxCoarseBox - xMinCoarseBox) / nxCoarseBoxFinest;   //mesh size assuming a coarse box is used
    double dyF = (dim > 1) ? (yMaxCoarseBox - yMinCoarseBox) / nyCoarseBoxFinest : 1. ;
    double dzF = (dim > 2) ? (zMaxCoarseBox - zMinCoarseBox) / nzCoarseBoxFinest : 1. ;

    double measureFinest = numberOfSamplesFinest * dxF * dyF * dzF;


    for (unsigned m = 0; m < numberOfSamplesFinest; m++) {


      if (dim == 1) {

        //BEGIN write finest histogram solution
        for (int iel = solFinest->GetMesh()->_elementOffset[iprocFinest]; iel < solFinest->GetMesh()->_elementOffset[iprocFinest + 1]; iel ++) {

          unsigned  xLeftDof = solFinest->GetMesh()->GetSolutionDof (0, iel, 2);
          unsigned  xRightDof = solFinest->GetMesh()->GetSolutionDof (1, iel, 2);

          double xLeft = (*solFinest->GetMesh()->_topology->_Sol[0]) (xLeftDof);
          double xRight = (*solFinest->GetMesh()->_topology->_Sol[0]) (xRightDof);

          if (sgmQoIStandardizedFinest[m][0] > xLeft && sgmQoIStandardizedFinest[m][0] <= xRight) {
            double histoValue = 1. / measureFinest;
            solFinest->_Sol[solIndexHISTOF]->add (iel, histoValue);

            break;
          }
        }

        //END
      }

      else {

        //BEGIN write finest histogram solution
        Marker marker2 (sgmQoIStandardizedFinest[m], 0., VOLUME, mlSolFinest->GetLevel (level), 2, true);
        unsigned iel2 = marker2.GetMarkerElement();

        if (iel2 >= solFinest->GetMesh()->_elementOffset[iprocFinest]  &&  iel2 < solFinest->GetMesh()->_elementOffset[iprocFinest + 1]) {
          double histoValue = 1. / measureFinest;
          solFinest->_Sol[solIndexHISTOF]->add (iel2, histoValue);

        }

        //END

      }

      solFinest->_Sol[solIndexHISTOF]->close();
    }

//END
  }

}


void GetKDEIntegral (MultiLevelProblem& ml_prob) {

  unsigned level = 0.;

  Mesh*                    msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  char name[10];
  sprintf (name, "PROPOSED");
  double solIndexKDE = mlSol->GetIndex (name);   // get the position of "Ti" in the sol object

  unsigned solTypeKDE = mlSol->GetSolutionType (solIndexKDE);

  std::vector <double> solKDELocal;
  vector < vector < double > > xLocal (dim);   // local coordinates
  vector < double >  xGauss (dim);   //

  vector < double > phi;
  vector < double> gradphi;
  double weight;

  double local_integral = 0.;

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofu  = msh->GetElementDofNumber (iel, solTypeKDE);   // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber (iel, 2);   // number of coordinate element dofs

    solKDELocal.resize (nDofu);

    for (int i = 0; i < dim; i++) {
      xLocal[i].resize (nDofx);
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, solTypeKDE);   // global to global mapping between solution node and solution dof
      solKDELocal[i] = (*sol->_Sol[solIndexKDE]) (solDof);
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, 2);   // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        xLocal[jdim][i] = (*msh->_topology->_Sol[jdim]) (xDof);     // global extraction and local storage for the element coordinates
      }
    }

    double weight;
    vector <double> phi;  // local test function

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solTypeKDE]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solTypeKDE]->Jacobian (xLocal, ig, weight, phi, gradphi);
      double solKDEGauss = 0.;

      for (unsigned i = 0; i < nDofu; i++) {
        solKDEGauss += phi[i] * solKDELocal[i];

        for (unsigned j = 0; j < dim; j++) {
          xGauss[j] += xLocal[j][i] * phi[i];
        }
      }

      local_integral += solKDEGauss * weight;
    }
  }

  double integral = 0.;
  MPI_Allreduce (&local_integral, &integral, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  std::cout << "Proposed method integral = " << integral << std::endl;


}



void GetAverageL2Error (std::vector< std::vector <double > > & sgmQoIStandardized, MultiLevelProblem& ml_prob, MultiLevelProblem& ml_probFinest) {

  unsigned level = 0.;

  Mesh*                    msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  Mesh*                    mshFinest = ml_probFinest._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  MultiLevelSolution*    mlSolFinest = ml_probFinest._ml_sol;  // pointer to the multilevel solution object
  Solution*                solFinest = ml_probFinest._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object
  unsigned    iprocFinest = mshFinest->processor_id(); // get the process_id (for parallel computation)

  char name[10];
  sprintf (name, "PROPOSED");
  double solIndexKDE = mlSol->GetIndex (name);
  unsigned solTypeKDE = mlSol->GetSolutionType (solIndexKDE);

  sprintf (name, "HISTO");
  double solIndexHISTOError = mlSol->GetIndex (name);

  sprintf (name, "HISTO_F");
  double solIndexHISTO = mlSolFinest->GetIndex (name);

  std::vector <double> solKdeLocal;

  vector < double > phi;
  vector < double> gradphi;
  double weight;
  double PI = acos (-1.);

  double aL2ELocal = 0;
  double aL2ELocalHisto = 0;
  double aL2ELocalFinest = 0;

  for (unsigned m = 0; m < numberOfSamples; m++) {

    double solKDESample = 0.;
    double solHISTOSampleError = 0.;
    double solHISTOFSample = 0.;

    if (dim == 1) {
      for (int iel = sol->GetMesh()->_elementOffset[iproc]; iel < sol->GetMesh()->_elementOffset[iproc + 1]; iel ++) {

        unsigned  xLeftDof = sol->GetMesh()->GetSolutionDof (0, iel, 2);
        unsigned  xRightDof = sol->GetMesh()->GetSolutionDof (1, iel, 2);

        double xLeft = (*sol->GetMesh()->_topology->_Sol[0]) (xLeftDof);
        double xRight = (*sol->GetMesh()->_topology->_Sol[0]) (xRightDof);

        if (sgmQoIStandardized[m][0] > xLeft && sgmQoIStandardized[m][0] <= xRight) {

          if (histoErr) solHISTOSampleError = (*sol->_Sol[solIndexHISTOError]) (iel);

          //BEGIN evaluate the KDE at the sample
          short unsigned ielType = msh->GetElementType (iel);
          unsigned nDofsKDE = msh->GetElementDofNumber (iel, solTypeKDE);
          solKdeLocal.resize (nDofsKDE);

          std::vector < std::vector < double> > vx (1);
          vx[0].resize (nDofsKDE);
          vx[0][0] = xLeft;
          vx[0][1] = xRight;

          std::vector < double> sampleLocal (1, 0.);
          sampleLocal[0] = - 1. + 2. * (sgmQoIStandardized[m][0] - xLeft) / (xRight - xLeft);

          msh->_finiteElement[ielType][solTypeKDE]->Jacobian (vx, sampleLocal, weight, phi, gradphi);

          for (unsigned inode = 0; inode < nDofsKDE; inode++) {
            unsigned globalDof = msh->GetSolutionDof (inode, iel, solTypeKDE);
            solKdeLocal[inode] = (*sol->_Sol[solIndexKDE]) (globalDof);
          }

          for (unsigned inode = 0; inode < nDofsKDE; inode++) {
            solKDESample += solKdeLocal[inode] * phi[inode];
          }

          //END

                    double diffPHI = 0.5 * ( 1. + erf ( 5.5 / sqrt ( 2 ) ) ) - 0.5 * ( 1. + erf ( -5.5 / sqrt ( 2 ) ) );
                    double stdGaussian = ( exp ( - sgmQoIStandardized[m][0] * sgmQoIStandardized[m][0] * 0.5 ) / sqrt ( 2 * PI ) ) / diffPHI;

                    aL2ELocal += ( solKDESample - stdGaussian ) * ( solKDESample - stdGaussian );


//                     double uniform = ( fabs ( sgmQoIStandardized[m][0] ) <= 1. ) ? 0.5 : 0. ;
// 
//                     aL2ELocal += ( solKDESample - uniform ) * ( solKDESample - uniform );

//                     if ( histoErr ) aL2ELocalHisto += ( solHISTOSampleError - uniform ) * ( solHISTOSampleError - uniform );

          break;
        }
      }

      if (histoFinest) {

        //BEGIN evaluate the finest histogram at the sample
        for (int iel = solFinest->GetMesh()->_elementOffset[iprocFinest]; iel < solFinest->GetMesh()->_elementOffset[iprocFinest + 1]; iel ++) {
          unsigned  xLeftDof = solFinest->GetMesh()->GetSolutionDof (0, iel, 2);
          unsigned  xRightDof = solFinest->GetMesh()->GetSolutionDof (1, iel, 2);

          double xLeft = (*solFinest->GetMesh()->_topology->_Sol[0]) (xLeftDof);
          double xRight = (*solFinest->GetMesh()->_topology->_Sol[0]) (xRightDof);

          if (sgmQoIStandardized[m][0] > xLeft && sgmQoIStandardized[m][0] <= xRight) {
            solHISTOFSample = (*solFinest->_Sol[solIndexHISTO]) (iel);

            break;
          }
        }

        //END

        aL2ELocalFinest += (solKDESample - solHISTOFSample) * (solKDESample - solHISTOFSample);
      }

    }

    else {

      Marker marker (sgmQoIStandardized[m], 0., VOLUME, mlSol->GetLevel (level), 2, true);
      unsigned iel = marker.GetMarkerElement();
      std::vector<double> sampleLocal;
      marker.GetMarkerLocalCoordinates (sampleLocal);

      if (iel >= sol->GetMesh()->_elementOffset[iproc]  &&  iel < sol->GetMesh()->_elementOffset[iproc + 1]) {

        if (histoErr) solHISTOSampleError = (*sol->_Sol[solIndexHISTOError]) (iel);

        //BEGIN evaluate the KDE at the sample
        short unsigned ielType = msh->GetElementType (iel);
        unsigned nDofsKDE = msh->GetElementDofNumber (iel, solTypeKDE);
        solKdeLocal.resize (nDofsKDE);

        std::vector < std::vector < double> > vx (dim);

        for (int idim = 0; idim < dim; idim++) {
          vx[idim].resize (nDofsKDE);
        }

        for (unsigned inode = 0; inode < nDofsKDE; inode++) {
          unsigned idofVx = msh->GetSolutionDof (inode, iel, 2);

          for (int jdim = 0; jdim < dim; jdim++) {
            vx[jdim][inode] = (*msh->_topology->_Sol[jdim]) (idofVx);
          }
        }

        msh->_finiteElement[ielType][solTypeKDE]->Jacobian (vx, sampleLocal, weight, phi, gradphi);

        for (unsigned inode = 0; inode < nDofsKDE; inode++) {
          unsigned globalDof = msh->GetSolutionDof (inode, iel, solTypeKDE);
          solKdeLocal[inode] = (*sol->_Sol[solIndexKDE]) (globalDof);
        }

        for (unsigned inode = 0; inode < nDofsKDE; inode++) {
          solKDESample += solKdeLocal[inode] * phi[inode];
        }

                double dotProduct = 0.;

                for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                    dotProduct += sgmQoIStandardized[m][jdim] * sgmQoIStandardized[m][jdim];
                }


                double diffPHI = 0.5 * ( 1. + erf ( 5.5 / sqrt ( 2 ) ) ) - 0.5 * ( 1. + erf ( -5.5 / sqrt ( 2 ) ) );
                double stdGaussian = ( exp ( - dotProduct * 0.5 ) / ( 2 * PI ) ) / ( diffPHI * diffPHI );


                if ( dim == 3 ) stdGaussian /= ( sqrt ( 2 * PI ) * diffPHI );



                aL2ELocal += ( solKDESample - stdGaussian ) * ( solKDESample - stdGaussian );


//                     double uniform = (dim == 2) ? 0.25 : 0.125;
//
//                     for(unsigned kdim = 0; kdim < dim; kdim++) {
//                         if(fabs(sgmQoIStandardized[m][kdim]) > 1.) uniform = 0.;
//                     }
// // //
//                     aL2ELocal += (solKDESample - uniform) * (solKDESample - uniform);

//                 double diffPHI = 0.5 * ( 1. + erf ( 5.5 / sqrt ( 2 ) ) ) - 0.5 * ( 1. + erf ( -5.5 / sqrt ( 2 ) ) );
//                 double laplaceDist = ( 1. / ( 2. * bLaplace ) ) * exp ( - fabs ( sgmQoIStandardized[m][0] - muLaplace ) / bLaplace ) / ( 0.974438 );
//                 double stdGaussian = ( exp ( - sgmQoIStandardized[m][1] * sgmQoIStandardized[m][1] * 0.5 ) / sqrt ( 2 * PI ) ) / diffPHI;
// 
//                 double jointPDF = laplaceDist * stdGaussian;
// 
//                 aL2ELocal += ( solKDESample - jointPDF ) * ( solKDESample - jointPDF );

//                 if ( histoErr ) aL2ELocalHisto += ( solHISTOSampleError - jointPDF ) * ( solHISTOSampleError - jointPDF );

        //END

      }

      if (histoFinest) {

        //BEGIN evaluate the finest histogram at the sample
        Marker marker2 (sgmQoIStandardized[m], 0., VOLUME, mlSolFinest->GetLevel (level), 2, true);
        unsigned iel2 = marker2.GetMarkerElement();

        if (iel2 >= solFinest->GetMesh()->_elementOffset[iprocFinest]  &&  iel2 < solFinest->GetMesh()->_elementOffset[iprocFinest + 1]) {
          solHISTOFSample = (*solFinest->_Sol[solIndexHISTO]) (iel2);

        }

        //END

        aL2ELocalFinest += (solKDESample - solHISTOFSample) * (solKDESample - solHISTOFSample);
      }

    }
  }

  double aL2E = 0.;
  double aL2EHisto = 0.;
  double aL2EFinest = 0.;


  if (histoFinest == true) MPI_Allreduce (&aL2ELocalFinest, &aL2EFinest, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  else {
    MPI_Allreduce (&aL2ELocal, &aL2E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (&aL2ELocalHisto, &aL2EHisto, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

  aL2E = sqrt (aL2E) / numberOfSamples;

  aL2EHisto = sqrt (aL2EHisto) / numberOfSamples;

  aL2EFinest = sqrt (aL2EFinest) / numberOfSamples;

  if (histoFinest == true) std::cout << "Average L2 Error Finest = " << std::setprecision (11) << aL2EFinest << std::endl;

  else {
    std::cout << "Average L2 Error = " << std::setprecision (11) << aL2E << std::endl;
    std::cout << "Average L2 Error Histogram = " << std::setprecision (11) << aL2EHisto << std::endl;
  }

}






