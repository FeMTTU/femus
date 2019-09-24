/** \file Ex6.cpp
 *  \brief This example shows how to set and solve the weak form
 *   of the Navier-Stokes Equation
 *
 *  \f{eqnarray*}
 *  && \mathbf{V} \cdot \nabla \mathbf{V} - \nabla \cdot \nu (\nabla \mathbf{V} +(\nabla \mathbf{V})^T)
 *  +\nabla P = 0 \\
 *  && \nabla \cdot \mathbf{V} = 0
 *  \f}
 *  in a unit box domain (in 2D and 3D) with given vertical velocity 1 on
 *  the left boundary and walls elsewhere.
 *  \author Eugenio Aulisa
 */

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "ImplicitRungeKuttaSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"


using namespace femus;

double GetSolutionFluxes (MultiLevelSolution& mlSol);

double GetTimeStep (const double time) {
  double dt = .02;
  return dt;
}


bool SetBoundaryCondition (const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;

  double L = 5;
  double xl = -2.5;

  if (!strcmp (SolName, "U")) { // strcmp compares two string in lexiographic sense. why !strcmp??? not false ==true?
    if (facename == 2 || facename == 1) {
      if (x[0] <= -2.5 + 1.0e-10) {
        value = sin (M_PI * time);
      }
    }
    else if (facename == 3) {
      dirichlet = false;
    }
  }
  else if (!strcmp (SolName, "V")) {
    if (facename == 3) {
      dirichlet = false;
    }
  }
  else if (!strcmp (SolName, "W")) {
    if (facename == 3) {
      dirichlet = false;
    }
  }
  else if (!strcmp (SolName, "P")) {
    dirichlet = false;
  }
  else if (!strcmp (SolName, "DX")) { // strcmp compares two string in lexiographic sense. why !strcmp??? not false ==true?
    value = (1. - cos (M_PI * time)) / M_PI * (L - (x[0] - xl)) / L;
  }
  
  return dirichlet;
}


void AssembleBoussinesqAppoximation_AD (MultiLevelProblem& ml_prob);   //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );


int main (int argc, char** args) {



  // init Petsc-MPI communicator
  FemusInit mpinit (argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh ("./input/rectangle2.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 2;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh (numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels (numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol (&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution ("U", LAGRANGE, SECOND, 2);
  mlSol.AddSolution ("V", LAGRANGE, SECOND,  2);
  if (dim == 3) mlSol.AddSolution ("W", LAGRANGE, SECOND, 2);

  mlSol.AddSolution ("DX", LAGRANGE, SECOND, 2);
  mlSol.AddSolution ("DY", LAGRANGE, SECOND,  2);
  if (dim == 3) mlSol.AddSolution ("DZ", LAGRANGE, SECOND, 2);

  mlSol.AddSolution ("P",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  //  Taylor-hood
  //  mlSol.AddSolution("U", LAGRANGE, SERENDIPITY);
  //  mlSol.AddSolution("V", LAGRANGE, SERENDIPITY);
  //  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SERENDIPITY);
  //  mlSol.AddSolution("P", LAGRANGE, FIRST);


  mlSol.Initialize ("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
  //mlSol.FixSolutionAtOnePoint("P");


  mlSol.GenerateBdc ("U", "Time_dependent");
  mlSol.GenerateBdc ("V", "Steady");
  if (dim == 3) mlSol.GenerateBdc ("W", "Steady");
  
  mlSol.GenerateBdc ("DX", "Time_dependent");
  mlSol.GenerateBdc ("DY", "Steady");
  if (dim == 3) mlSol.GenerateBdc ("DZ", "Steady");
  
  mlSol.GenerateBdc ("P", "Steady");


  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb (&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  ImplicitRungeKuttaNonlinearImplicitSystem & system = mlProb.add_system < ImplicitRungeKuttaNonlinearImplicitSystem > ("NS");

  //system.SetImplicitRungeKuttaScheme (CROUZEIX2);
  system.SetImplicitRungeKuttaScheme (DIRK3);

  std::vector < std::string > solVName (3), solDName (3);
  solVName[0] = "U";
  solVName[1] = "V";
  solVName[2] = "W";
  solDName[0] = "DX";
  solDName[1] = "DY";
  solDName[2] = "DZ";
  
  for(unsigned k=0; k < dim; k++){
    system.AddSolutionToSystemPDE ( solVName[k].c_str());
  }
  for(unsigned k=0; k < dim; k++){
    system.AddSolutionToSystemPDE ( solDName[k].c_str());
  }
  system.AddSolutionToSystemPDE ("P");
  system.SetRKVariableType ("P", false);
    
  
  std::vector < std::vector < std::string > > solVkName(dim), solDkName(dim);
  for(unsigned k = 0; k<dim; k++){ 
    solVkName[k] = system.GetSolkiNames( solVName[k].c_str() );
    solDkName[k] = system.GetSolkiNames( solDName[k].c_str() );
  }
  const std::vector < std::string > solPkName = system.GetSolkiNames ("P");
  
  unsigned RK = system.GetRungeKuttaStages();
    
  FieldSplitTree **VDP;
  VDP = new FieldSplitTree * [RK];
  std::vector < FieldSplitTree *> VDPAll;
  VDPAll.reserve(RK);
  
  std::vector < std::vector < unsigned > > fieldVDP(RK);
  std::vector < std::vector < unsigned > > solutionTypeVDP(RK);
 
  for(unsigned i = 0; i < RK; i++){
    fieldVDP[i].resize(2 * dim + 1);
    solutionTypeVDP[i].resize(2 * dim + 1);
    for(unsigned k = 0; k < dim; k++){
      fieldVDP[i][k] = system.GetSolPdeIndex(solVkName[k][i].c_str());
      solutionTypeVDP[i][k] = mlSol.GetSolutionType(solVkName[k][i].c_str());
    }
    for(unsigned k = 0; k < dim; k++){
      fieldVDP[i][dim + k] = system.GetSolPdeIndex(solDkName[k][i].c_str());
      solutionTypeVDP[i][dim + k] = mlSol.GetSolutionType(solDkName[k][i].c_str());
    }
    fieldVDP[i][2 * dim] = system.GetSolPdeIndex(solPkName[i].c_str());
    solutionTypeVDP[i][2 * dim ] = mlSol.GetSolutionType(solPkName[i].c_str());

    char name[10];
    sprintf (name, "VDPi%d", i);
    VDP[i]= new FieldSplitTree (PREONLY, MLU_PRECOND, fieldVDP[i], solutionTypeVDP[i], name);  
    VDPAll.push_back(VDP[i]);  
  }
  
  FieldSplitTree FS(PREONLY, FIELDSPLIT_MULTIPLICATIVE_PRECOND, VDPAll, "RK");
  FS.SetRichardsonScaleFactor(1.);

  system.SetLinearEquationSolverType(FEMuS_FIELDSPLIT, INCLUDE_COARSE_LEVEL_TRUE);   // Field-Split preconditioned
  
  // attach the assembling function to system
  system.SetAssembleFunction (AssembleBoussinesqAppoximation_AD);
 
  // initilaize and solve the system
  system.init();
 
  system.SetSolverCoarseGrid(RICHARDSON);
  system.SetRichardsonScaleFactor(1.);

  //system.SetSolverCoarseGrid(PREONLY);
  
  system.SetOuterSolver(PREONLY);
  
  system.SetTolerances(1.e-10, 1.e-10, 1.e+50, 10, 10); //GMRES tolerances
  
  system.AttachGetTimeIntervalFunction (GetTimeStep);
  const unsigned int n_timesteps = 2;

  system.SetFieldSplitTree(&FS);
  
  // ******* Print solution *******
  mlSol.SetWriter (VTK);
  mlSol.GetWriter()->SetDebugOutput (false);

  std::vector<std::string> print_vars;
  print_vars.push_back ("All");
  //mlSol.GetWriter()->SetDebugOutput (true);

  std::vector<std::string> mov_vars;
  mov_vars.push_back ("DX");
  mov_vars.push_back ("DY");
  if ( dim==3 ) mov_vars.push_back ("DY");
  mlSol.GetWriter()->SetMovingMesh (mov_vars);

  mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);

  for (unsigned time_step = 0; time_step < n_timesteps; time_step++) {

    system.CopySolutionToOldSolution();

    system.MGsolve();

    std::cout << "Fluxes = " << GetSolutionFluxes (mlSol) << "    " << sin (M_PI * system.GetTime()) << std::endl;

    mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", print_vars, time_step + 1);
  }

  for (unsigned i = 0; i < RK; i++) {
    delete VDP[i];
  }
  delete [] VDP;
  
  return 0;
}


void AssembleBoussinesqAppoximation_AD (MultiLevelProblem& ml_prob) {

  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  ImplicitRungeKuttaNonlinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<ImplicitRungeKuttaNonlinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27

  double dt = GetTimeStep (0.);

  unsigned RK = mlPdeSys->GetRungeKuttaStages();

  std::vector < std::string > solVName (3), solDName (3);
  solVName[0] = "U";
  solVName[1] = "V";
  solVName[2] = "W";
  solDName[0] = "DX";
  solDName[1] = "DY";
  solDName[2] = "DZ";
  std::string  solPName = "P";

  //solution variable
  vector < unsigned > solVIndex (dim),solDIndex (dim);
  for(unsigned k = 0;k<dim; k++){
    solVIndex[k] = mlSol->GetIndex (solVName[k].c_str());
    solDIndex[k] = mlSol->GetIndex (solDName[k].c_str());
  }
  unsigned solPIndex = mlSol->GetIndex ("P");   // get the position of "P" in the ml_sol object

  unsigned solType = mlSol->GetSolutionType (solVIndex[0]);   // get the finite element type for "u"
  unsigned solPType = mlSol->GetSolutionType (solPIndex);   // get the finite element type for "u"

  std::vector < std::vector < std::string > > solVkName(dim),solDkName(dim);
  for(unsigned k = 0;k<dim; k++){ 
    solVkName[k] = mlPdeSys->GetSolkiNames( solVName[k].c_str() );
    solDkName[k] = mlPdeSys->GetSolkiNames( solDName[k].c_str() );
  }
  const std::vector < std::string > solPkName = mlPdeSys->GetSolkiNames ("P");

  std::vector < std::vector < unsigned > > solVkIndex (dim), solDkIndex (dim); 
  std::vector < std::vector < unsigned > > solVkPdeIndex (dim), solDkPdeIndex (dim);
  std::vector < unsigned > solPkIndex (RK), solPkPdeIndex (RK);
  for (unsigned k = 0; k < dim; k++) {
    solVkIndex[k].resize (RK);
    solDkIndex[k].resize (RK);
    solVkPdeIndex[k].resize (RK);
    solDkPdeIndex[k].resize (RK);
  }

  for (unsigned jj = 0; jj < RK; jj++) {
    for(unsigned k = 0; k < dim ;k++){
      solVkIndex[k][jj] = mlSol->GetIndex (solVkName[k][jj].c_str());
      solDkIndex[k][jj] = mlSol->GetIndex (solDkName[k][jj].c_str());
      solVkPdeIndex[k][jj] = mlPdeSys->GetSolPdeIndex (solVkName[k][jj].c_str());
      solDkPdeIndex[k][jj] = mlPdeSys->GetSolPdeIndex (solDkName[k][jj].c_str());
    }
    solPkIndex[jj] = mlSol->GetIndex (solPkName[jj].c_str());
    solPkPdeIndex[jj] = mlPdeSys->GetSolPdeIndex (solPkName[jj].c_str());
  }
  //solution variable

  std::vector < std::vector < std::vector < adept::adouble > > >  solV (dim), solD (dim);   // local solution
  std::vector < std::vector < std::vector < adept::adouble > > >  solVk (dim), solDk (dim);   // local solution
  std::vector < std::vector < double > >   solVOld (dim), solDOld (dim);   // local solution
  std::vector < std::vector< std::vector < adept::adouble > > > aResV (dim), aResD (dim);   // local redidual vector
  std::vector < std::vector < adept::adouble > > solP (RK), aResP (RK); // local solution
  for (unsigned k = 0; k < dim; k++) {
    solV[k].resize (RK);
    solVk[k].resize (RK);
    solD[k].resize (RK);
    solDk[k].resize (RK);
    aResV[k].resize (RK);
    aResD[k].resize (RK);
  }

  std::vector < std::vector < adept::adouble > > x (dim);   // local coordinates
  std::vector < std::vector < double > > xHat (dim);   // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phi;  // local test function for velocity
  std::vector <adept::adouble> phi_x; // local test function first order partial derivatives

  double* phiP; // local test function for the pressure
  adept::adouble weight; // gauss point weight

  std::vector< int > sysDof; // local to global pdeSys dofs
  std::vector< double > Res; // local redidual vector
  std::vector < double > Jac;
  
  std::vector < adept::adouble > solV_gss (dim), solVk_gss (dim), solDk_gss (dim), NSV (dim), DISP (dim);
  std::vector < std::vector < adept::adouble > > gradSolV_gss (dim), /*gradSolVk_gss (dim),*/ gradSolDk_gss (dim);
  
  RES->zero(); // Set to zero all the entries of the Global Residual vector
  KK->zero(); // Set to zero all the entries of the Global Matrix
  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);

    unsigned nDofs = msh->GetElementDofNumber (iel, solType);   // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber (iel, solPType);   // number of solution element dofs

    unsigned DStart = dim * nDofs * RK;
    unsigned PStart = 2 * DStart;
    
    unsigned nDofsVDP = PStart + nDofsP * RK;
    // resize local arrays
    sysDof.resize (nDofsVDP);

    for (unsigned jj = 0; jj < RK; jj++) {
      for (unsigned  k = 0; k < dim; k++) {
        solV[k][jj].resize (nDofs);
        solVk[k][jj].resize (nDofs);
        solD[k][jj].resize (nDofs);
        solDk[k][jj].resize (nDofs);
        aResV[k][jj].assign (nDofs, 0.);
        aResD[k][jj].assign (nDofs, 0.);
      }
      solP[jj].resize (nDofsP);
      aResP[jj].assign (nDofsP, 0.);
    }

    for (unsigned  k = 0; k < dim; k++) {
      solVOld[k].resize (nDofs);
      solDOld[k].resize (nDofs);
      xHat[k].resize (nDofs);
      x[k].resize (nDofs);
    }
    
    // local storage of global mapping and solution
   
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned idof = msh->GetSolutionDof (i, iel, solType);   // local to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solVOld[k][i] = (*sol->_SolOld[solVIndex[k]]) (idof);     // global extraction and local storage for the solution
        solDOld[k][i] = (*sol->_SolOld[solDIndex[k]]) (idof);     // global extraction and local storage for the solution
        for (unsigned jj = 0; jj < RK; jj++) {
          solVk[k][jj][i] = (*sol->_Sol[solVkIndex[k][jj]]) (idof);     // global extraction and local storage for the solution
          sysDof[k * RK * nDofs + jj * nDofs + i] = pdeSys->GetSystemDof (solVkIndex[k][jj], solVkPdeIndex[k][jj], i, iel);
          // stacking of the Uk1, Uk2, ..., UkRK,Vk1,Vk2,...,VkRK
          solDk[k][jj][i] = (*sol->_Sol[solDkIndex[k][jj]]) (idof);
          sysDof[ DStart + k * RK * nDofs + jj * nDofs + i] = pdeSys->GetSystemDof (solDkIndex[k][jj], solDkPdeIndex[k][jj], i, iel);
          // stacking of the DXk1, DXk2, ..., DXkRK,DYk1,DYk2,...,DYkRK after all the velocities (starting from dStart)
        }
      }
    }

    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof (i, iel, solPType);   // local to global mapping between solution node and solution dof
      //solPOld[i] = (*sol->_SolOld[solPIndex]) (solPDof);

      for (unsigned jj = 0; jj < RK; jj++) {
        solP[jj][i] = (*sol->_Sol[solPkIndex[jj]]) (solPDof);     // global extraction and local storage for the solution
        sysDof[ PStart + jj * nDofsP + i ] = pdeSys->GetSystemDof (solPkIndex[jj], solPkPdeIndex[jj], i, iel);
        // global to global mapping betweensolution node and pdeSys dof
      }
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned coordXDof  = msh->GetSolutionDof (i, iel, coordXType);   // local to global mapping between coordinates node and coordinate dof
      for (unsigned k = 0; k < dim; k++) {
        xHat[k][i] = (*msh->_topology->_Sol[k]) (coordXDof);     // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    for (unsigned  k = 0; k < dim; k++) {
      mlPdeSys->GetIntermediateSolutions (solVOld[k], solVk[k], solV[k]);
    }

    for (unsigned jj = 0; jj < RK; jj++) {

      // local storage of coordinates
      for (unsigned i = 0; i < nDofs; i++) {
        unsigned coordXDof  = msh->GetSolutionDof (i, iel, coordXType);   // local to global mapping between coordinates node and coordinate dof
        for (unsigned k = 0; k < dim; k++) {
          x[k][i] = xHat[k][i] + solD[k][jj][i];      // global extraction and local storage for the element coordinates
        }
      }

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[ielGeom][solType]->Jacobian (x, ig, weight, phi, phi_x);
        
        phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi (ig);

        solV_gss.assign(dim, 0);
        solVk_gss.assign(dim, 0);
        solDk_gss.assign(dim, 0);
       
        for (unsigned  k = 0; k < dim; k++) {
          gradSolV_gss[k].assign (dim, 0.);
          //gradSolVk_gss[k].assign (dim, 0.);
          gradSolDk_gss[k].assign (dim, 0.);
        }

        for (unsigned i = 0; i < nDofs; i++) {
          for (unsigned  k = 0; k < dim; k++) {
            solV_gss[k] += solV[k][jj][i] * phi[i];
            solVk_gss[k] += solVk[k][jj][i] * phi[i];
            solDk_gss[k] += solDk[k][jj][i] * phi[i];
          }
          for (unsigned j = 0; j < dim; j++) {
            for (unsigned  k = 0; k < dim; k++) {
              gradSolV_gss[k][j] += solV[k][jj][i] * phi_x[i * dim + j];
              //gradSolVk_gss[k][j] += solVk[k][jj][i] * phi_x[i * dim + j];
              gradSolDk_gss[k][j] += solDk[k][jj][i] * phi_x[i * dim + j];
            }
          }
        }

        adept::adouble solP_gss = 0;
        for (unsigned i = 0; i < nDofsP; i++) {
          solP_gss += phiP[i] * solP[jj][i];
        }

        double nu = 1.;

        // *** phiV_i loop ***
        for (unsigned i = 0; i < nDofs; i++) {
          
          
          NSV.assign (dim, 0.);
          DISP.assign (dim, 0.);

          for (unsigned  k = 0; k < dim; k++) { //momentum equation in k
            for (unsigned j = 0; j < dim; j++) { // second index j in each equation
              NSV[k]   +=  nu * phi_x[i * dim + j] * (gradSolV_gss[k][j] + gradSolV_gss[j][k]);   // weak laplace
              NSV[k]   +=  phi[i] * (solV_gss[j] - solDk_gss[j]) * gradSolV_gss[k][j];      // non-linear term
              DISP[k]  +=  phi_x[i * dim + j] * gradSolDk_gss[k][j]; // laplace
            }
            NSV[k] += -solP_gss * phi_x[i * dim + k]; // pressure gradient
          }
          for (unsigned  k = 0; k < dim; k++) {
            aResV[k][jj][i] += (- phi[i] * solVk_gss[k] - NSV[k]) * weight;
            aResD[k][jj][i] += (- DISP[k]) * weight;
          }
        } // end phiV_i loop

        // *** phiP_i loop ***
        for (unsigned i = 0; i < nDofsP; i++) {
          for (int k = 0; k < dim; k++) {
            aResP[jj][i] += - (gradSolV_gss[k][k]) * phiP[i]  * weight;
          }
        } // end phiP_i loop
      } // end gauss point loop
    }
    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize (nDofsVDP);   //resize


    for (unsigned  k = 0; k < dim; k++) {
      for (unsigned  jj = 0; jj < RK; jj++) {
        for (unsigned i = 0; i < nDofs; i++) {
          Res[ k * RK * nDofs + jj * nDofs + i ] = -aResV[k][jj][i].value();
          Res[ DStart + k * RK * nDofs + jj * nDofs + i ] = -aResD[k][jj][i].value();
        }
      }
    }

    for (unsigned  jj = 0; jj < RK; jj++) {
      for (int i = 0; i < nDofsP; i++) {
        Res[ PStart + jj * nDofsP + i] = -aResP[jj][i].value();
      }
    }

    RES->add_vector_blocked (Res, sysDof);

    //Extract and store the Jacobian

    Jac.resize (nDofsVDP * nDofsVDP);
    // define the dependent variables

    for (unsigned  k = 0; k < dim; k++) {
      for (unsigned  jj = 0; jj < RK; jj++) {
        s.dependent (&aResV[k][jj][0], nDofs);
      }
    }
    for (unsigned  k = 0; k < dim; k++) {
      for (unsigned  jj = 0; jj < RK; jj++) {
        s.dependent (&aResD[k][jj][0], nDofs);
      }
    }
    for (unsigned  jj = 0; jj < RK; jj++) {
      s.dependent (&aResP[jj][0], nDofsP);
    }

    // define the independent variables
    for (unsigned  k = 0; k < dim; k++) {
      for (unsigned  jj = 0; jj < RK; jj++) {
        s.independent (&solVk[k][jj][0], nDofs);
      }
    }
    for (unsigned  k = 0; k < dim; k++) {
      for (unsigned  jj = 0; jj < RK; jj++) {
        s.independent (&solDk[k][jj][0], nDofs);
      }
    }
    for (unsigned  jj = 0; jj < RK; jj++) {
      s.independent (&solP[jj][0], nDofsP);
    }

    // get the and store jacobian matrix (row-major)
    s.jacobian (&Jac[0], true); // This is rowwise order.
    KK->add_matrix_blocked (Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();
// ***************** END ASSEMBLY *******************
}





double GetSolutionFluxes (MultiLevelSolution& mlSol) {

  int  iproc, nprocs;
  MPI_Comm_rank (MPI_COMM_WORLD, &iproc);
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs);

  MyVector<double> pFlux (1, 0);
  pFlux.stack();

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* solution  = mlSol.GetSolutionLevel (level);
  Mesh* msh = mlSol._mlMesh->GetLevel (level);
  elem* myel =  msh->el;

  const unsigned dim = msh->GetDimension();
  const unsigned max_size = static_cast< unsigned > (ceil (pow (3, dim)));

  vector< vector < double> >  sol (dim);
  vector< vector < double> > x (dim);

  const char varname[6][3] = {"U", "V", "W", "DX", "DY", "DZ"};
  vector <unsigned> indVar (2 * dim);
  unsigned solType;

  for (unsigned ivar = 0; ivar < dim; ivar++) {
    for (unsigned k = 0; k < 2; k++) {
      indVar[ivar + k * dim] = mlSol.GetIndex (&varname[ivar + k * 3][0]);
    }
  }
  solType = mlSol.GetSolutionType (&varname[0][0]);


  std::vector < double > phi;
  std::vector < double > gradphi;
  //std::vector< double > xx(dim, 0.);
  double weight;

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    // loop on faces
    for (unsigned jface = 0; jface < msh->GetElementFaceNumber (iel); jface++) {
      int faceNumber = myel->GetBoundaryIndex (iel, jface);
      // look for boundary faces
      if (faceNumber == 3) {

        vector < double> normal (dim, 0);

        unsigned nve = msh->GetElementFaceDofNumber (iel, jface, solType);
        const unsigned felt = msh->GetElementFaceType (iel, jface);

        for (unsigned d = 0; d < dim; d++) {
          x[d].resize (nve);
          sol[d].resize (nve);
        }

        for (unsigned i = 0; i < nve; i++) {
          unsigned int ilocal = msh->GetLocalFaceVertexIndex (iel, jface, i);
          unsigned idof = msh->GetSolutionDof (ilocal, iel, 2);
          for (unsigned d = 0; d < dim; d++) {
            x[d][i] = (*msh->_topology->_Sol[d]) (idof) + (*solution->_Sol[indVar[d + dim]]) (idof);;
            sol[d][i] = (*solution->_Sol[indVar[d]]) (idof);;
          }
        }

        for (unsigned igs = 0; igs < msh->_finiteElement[felt][solType]->GetGaussPointNumber(); igs++) {
          msh->_finiteElement[felt][solType]->JacobianSur (x, igs, weight, phi, gradphi, normal);
          double value;
          for (unsigned i = 0; i < nve; i++) {
            value = 0.;
            for (unsigned d = 0; d < dim; d++) {
              value += normal[d] * sol[d][i];
            }
            value *= phi[i];
            pFlux[iproc] += value * weight;
          }
        }
      }
    }
  }

  double flux = 0.;

  for (int j = 0; j < nprocs; j++) {
    pFlux.broadcast (j);
    flux += pFlux[j];
    pFlux.clearBroadcast();
  }

  return flux;
}



