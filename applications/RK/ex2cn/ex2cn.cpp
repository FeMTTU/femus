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
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"


using namespace femus;

double GetSolutionFluxes(MultiLevelSolution& mlSol);

double GetTimeStep(const double time) {
  double dt = .02;
  return dt;
}


bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;
  
  double L = 5;
  double xl = -2.5;
  
  if (!strcmp(SolName, "U")) { // strcmp compares two string in lexiographic sense. why !strcmp??? not false ==true?
    if (facename == 2 || facename == 1 ) {
      if(x[0] <= -2.5 + 1.0e-10){  
        value = sin(M_PI * time);
      }
    }
    else if (facename == 3) {
      dirichlet = false;  
    }
  } 
  else if (!strcmp(SolName, "V")) {
    if (facename == 3) {
      dirichlet = false;  
    }
  }
  else if (!strcmp(SolName, "W")) {
    if (facename == 3) {  
      dirichlet = false;  
    }
  } 
  else if (!strcmp(SolName, "P")) {
    dirichlet = false;
  }
  else if (!strcmp(SolName, "DX")) { // strcmp compares two string in lexiographic sense. why !strcmp??? not false ==true?
    value = ( 1. - cos(M_PI * time) ) / M_PI * (L - ( x[0] - xl ) ) / L;   
  }
  return dirichlet;
}


void AssembleBoussinesqAppoximation_AD(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );


int main(int argc, char** args) {



  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh ( "./input/rectangle2.neu", "seventh", scalingFactor );
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

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND , 2);
  mlSol.AddSolution("V", LAGRANGE, SECOND,  2);
  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND, 2);
  
  mlSol.AddSolution("DX", LAGRANGE, SECOND , 2);
  mlSol.AddSolution("DY", LAGRANGE, SECOND,  2);
  if (dim == 3) mlSol.AddSolution("DZ", LAGRANGE, SECOND, 2);
    
  mlSol.AddSolution("P",  DISCONTINUOUS_POLYNOMIAL, FIRST);
  
  //  Taylor-hood
  //  mlSol.AddSolution("U", LAGRANGE, SERENDIPITY);
  //  mlSol.AddSolution("V", LAGRANGE, SERENDIPITY);
  //  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SERENDIPITY);
  //  mlSol.AddSolution("P", LAGRANGE, FIRST);
 
  
  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  //mlSol.FixSolutionAtOnePoint("P");
  
  
  mlSol.GenerateBdc("U", "Time_dependent");
  mlSol.GenerateBdc("V", "Steady");
  
  mlSol.GenerateBdc("DX", "Time_dependent");
  mlSol.GenerateBdc("DY", "Steady");
  
  mlSol.GenerateBdc("P", "Steady");
  
  
  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  TransientNonlinearImplicitSystem& system = mlProb.add_system < TransientNonlinearImplicitSystem > ("NS");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  
  system.AddSolutionToSystemPDE("DX");
  system.AddSolutionToSystemPDE("DY");

  if (dim == 3) system.AddSolutionToSystemPDE("W");

  system.AddSolutionToSystemPDE("P");

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleBoussinesqAppoximation_AD);

  // initilaize and solve the system
  system.init();
  
  system.AttachGetTimeIntervalFunction(GetTimeStep);
  const unsigned int n_timesteps = 100;
  
  // ******* Print solution *******
  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->SetDebugOutput(false);

  std::vector<std::string> print_vars;
  print_vars.push_back("All");
  mlSol.GetWriter()->SetDebugOutput(true);
  
  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  mlSol.GetWriter()->SetMovingMesh(mov_vars);
  
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR,"biquadratic",print_vars, 0);
 
  for (unsigned time_step = 0; time_step < n_timesteps; time_step++) {

    system.CopySolutionToOldSolution();

    system.MGsolve();
    
    std::cout<<"Fluxes = " << GetSolutionFluxes(mlSol) << "    " << sin( M_PI * system.GetTime() ) << std::endl;

    mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR,"biquadratic",print_vars, time_step+1);
  }
  
  return 0;
}


void AssembleBoussinesqAppoximation_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  TransientNonlinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientNonlinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
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

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  double dt = GetTimeStep(0.);
  
  //solution variable
  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object
  if (dim == 3) solVIndex[2] = mlSol->GetIndex("W");      // get the position of "V" in the ml_sol object
  
  //solution variable
  vector < unsigned > solDIndex(dim);
  solDIndex[0] = mlSol->GetIndex("DX");    // get the position of "U" in the ml_sol object
  solDIndex[1] = mlSol->GetIndex("DY");    // get the position of "V" in the ml_sol object
  if (dim == 3) solDIndex[2] = mlSol->GetIndex("DZ");      // get the position of "V" in the ml_sol object
 
  unsigned solType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"
  
  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object
  if (dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  vector < unsigned > solDPdeIndex(dim);
  solDPdeIndex[0] = mlPdeSys->GetSolPdeIndex("DX");    // get the position of "U" in the pdeSys object
  solDPdeIndex[1] = mlPdeSys->GetSolPdeIndex("DY");    // get the position of "V" in the pdeSys object
  if (dim == 3) solDPdeIndex[2] = mlPdeSys->GetSolPdeIndex("DZ");
  
  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  vector < vector < adept::adouble > >  solV(dim);    // local solution
  vector < vector < double > >  solVOld(dim);    // local solution
  
  vector < vector < adept::adouble > >  solD(dim);    // local solution
  vector < vector < double > >  solDOld(dim);    // local solution
  
  vector < adept::adouble >  solP; // local solution

  vector< vector < adept::adouble > > aResV(dim);    // local redidual vector
  vector< vector < adept::adouble > > aResD(dim);    // local redidual vector
  vector< adept::adouble > aResP; // local redidual vector

  vector < vector < adept::adouble > > x(dim);    // local coordinates
  vector < vector < double > > xHat(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
    solVOld[k].reserve(maxSize);
    solD[k].reserve(maxSize);
    solDOld[k].reserve(maxSize);
    aResV[k].reserve(maxSize);
    aResD[k].reserve(maxSize);
    xHat[k].reserve(maxSize);
    x[k].reserve(maxSize);
  }

  solP.reserve(maxSize);
  aResP.reserve(maxSize);

  
  vector <double> phiHat;  // local test function for velocity
  vector <double> phiHat_x; // local test function first order partial derivatives

  vector <double> phi;  // local test function for velocity
  vector <adept::adouble> phi_x; // local test function first order partial derivatives

  phiHat.reserve(maxSize);
  phiHat_x.reserve(maxSize * dim);
  
  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  
  double* phiP; // local test function for the pressure
  adept::adouble weight; // gauss point weight
  double weightHat; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve((2 * dim + 1) * maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve((2 * dim + 1) * maxSize);

  vector < double > Jac;
  Jac.reserve((2 * dim + 1) * maxSize * (2 * dim + 1) * maxSize);

  RES->zero(); // Set to zero all the entries of the Global Residual vector
  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
        
    unsigned nDofsVDP = 2 * dim * nDofs + nDofsP;
    // resize local arrays
    sysDof.resize(nDofsVDP);

    for (unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofs);
      solVOld[k].resize(nDofs);
      solD[k].resize(nDofs);
      solDOld[k].resize(nDofs);
      xHat[k].resize(nDofs);
      x[k].resize(nDofs);
    }
    solP.resize(nDofsP);

    for (unsigned  k = 0; k < dim; k++) {
      aResV[k].assign(nDofs,0.);
      aResD[k].assign(nDofs,0.);
    }
    aResP.assign(nDofsP,0.);

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solType);    // local to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(idof);      // global extraction and local storage for the solution
        solVOld[k][i] = (*sol->_SolOld[solVIndex[k]])(idof);      // global extraction and local storage for the solution
        
        solD[k][i] = (*sol->_Sol[solDIndex[k]])(idof);      // global extraction and local storage for the solution
        solDOld[k][i] = (*sol->_SolOld[solDIndex[k]])(idof);      // global extraction and local storage for the solution
        
        sysDof[k * nDofs + i] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
        sysDof[(dim + k) * nDofs + i] = pdeSys->GetSystemDof(solDIndex[k], solDPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // local to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      sysDof[ 2 * dim * nDofs + i ] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for (unsigned k = 0; k < dim; k++) {
        xHat[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();
    
    
     // local storage of coordinates
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for (unsigned k = 0; k < dim; k++) {
        x[k][i] = xHat[k][i] + 0.5 * (solDOld[k][i] + solD[k][i]) ;      // global extraction and local storage for the element coordinates
      }
    }
    
    
//     // *** Face Gauss point loop (boundary Integral) ***
//     for ( unsigned jface = 0; jface < msh->GetElementFaceNumber ( iel ); jface++ ) {
//       int faceIndex = el->GetBoundaryIndex(iel, jface);
//       // look for boundary faces
//       if ( faceIndex == 3 ) {  
//         const unsigned faceGeom = msh->GetElementFaceType ( iel, jface );
//         unsigned faceDofs = msh->GetElementFaceDofNumber (iel, jface, solVType);
//                     
//         vector  < vector  <  double> > faceCoordinates ( dim ); // A matrix holding the face coordinates rowwise.
//         for ( int k = 0; k < dim; k++ ) {
//           faceCoordinates[k].resize (faceDofs);
//         }
//         for ( unsigned i = 0; i < faceDofs; i++ ) {
//           unsigned inode = msh->GetLocalFaceVertexIndex ( iel, jface, i ); // face-to-element local node mapping.
//           for ( unsigned k = 0; k < dim; k++ ) {
//             faceCoordinates[k][i] =  coordX[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
//           }
//         }
//         for ( unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solVType]->GetGaussPointNumber(); ig++ ) { 
//             // We call the method GetGaussPointNumber from the object finiteElement in the mesh object msh. 
//           vector < double> normal;
//           msh->_finiteElement[faceGeom][solVType]->JacobianSur ( faceCoordinates, ig, weight, phiV, phiV_x, normal );
//             
//           vector< double > xg(dim,0.);
//           for ( unsigned i = 0; i < faceDofs; i++ ) {
//             for( unsigned k=0; k<dim; k++){
//               xg[k] += phiV[i] * faceCoordinates[k][i]; // xg(ig)= \sum_{i=0}^faceDofs phi[i](xig) facecoordinates[i]     
//             }
//           }
//           double tau; // a(u)*grad_u\cdot normal
//           SetBoundaryCondition( xg, "P", tau, faceIndex, 0. ); // return tau
//           // *** phi_i loop ***
//           for ( unsigned i = 0; i < faceDofs; i++ ) {
//             unsigned inode = msh->GetLocalFaceVertexIndex ( iel, jface, i );
//             for(unsigned k=0;k<dim;k++){
//               aResV[k][inode] +=  -phiV[i] * tau * normal[k] * weight;
//             }
//           }        
//         }
//       }
//     }   

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solType]->Jacobian(x, ig, weight, phi, phi_x);
      msh->_finiteElement[ielGeom][solType]->Jacobian(xHat, ig, weightHat, phiHat, phiHat_x);
      
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

      vector < adept::adouble > solV_gss(dim, 0);
      vector < double > solVOld_gss(dim, 0);
      vector < vector < adept::adouble > > gradSolV_gss(dim);
      vector < vector < adept::adouble > > gradSolVOld_gss(dim);
      
      vector < adept::adouble > solD_gss(dim, 0);
      vector < double > solDOld_gss(dim, 0);
      vector < vector < adept::adouble > > gradSolDHat_gss(dim);

      

      for (unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].assign(dim,0.);
        gradSolVOld_gss[k].assign(dim,0.);
        gradSolDHat_gss[k].assign(dim,0.);
      }

      for (unsigned i = 0; i < nDofs; i++) {
        for (unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += solV[k][i] * phi[i];
          solVOld_gss[k] += solVOld[k][i] * phi[i];
          solD_gss[k] += solD[k][i] * phi[i];
          solDOld_gss[k] += solDOld[k][i] * phi[i];
        }
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += solV[k][i] * phi_x[i * dim + j];
            gradSolVOld_gss[k][j] += solVOld[k][i] * phi_x[i * dim + j];
            gradSolDHat_gss[k][j] += solD[k][i] * phiHat_x[i * dim + j];
          }
        }
      }

      adept::adouble solP_gss = 0;
      for (unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }

      double nu = 1.;

      // *** phiV_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {
        vector < adept::adouble > NSV(dim, 0.);
        vector < adept::adouble > DISP(dim, 0.);
        
        for (unsigned  k = 0; k < dim; k++) { //momentum equation in k 
          for (unsigned j = 0; j < dim; j++) { // second index j in each equation
            NSV[k]   +=  nu * phi_x[i * dim + j] * 0.5 * ( (gradSolV_gss[k][j] + gradSolV_gss[j][k]) + (gradSolVOld_gss[k][j] + gradSolVOld_gss[j][k]) ); // laplace
            NSV[k]   +=  phi[i] * ( 0.5 * (solV_gss[j] + solVOld_gss[j]) - (solD_gss[j] - solDOld_gss[j]) / dt ) * 0.5 * ( gradSolV_gss[k][j] + gradSolVOld_gss[k][j]); // non-linear term
            DISP[k]  +=  phiHat_x[i * dim + j] * (gradSolDHat_gss[k][j] + gradSolDHat_gss[j][k]); // laplace
          }
          NSV[k] += -solP_gss * phi_x[i * dim + k]; // pressure gradient
        }
        for (unsigned  k = 0; k < dim; k++) {
          aResV[k][i] += ( - phi[i] * (solV_gss[k] + solVOld_gss[k]) / dt - NSV[k] ) * weight;
          aResD[k][i] += ( - DISP[k] ) * weightHat;
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for (unsigned i = 0; i < nDofsP; i++) {
        for (int k = 0; k < dim; k++) {
          aResP[i] += - (gradSolV_gss[k][k]) * phiP[i]  * weight;
        }
      } // end phiP_i loop

    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize(nDofsVDP);    //resize

    for (int i = 0; i < nDofs; i++) {
      for (unsigned  k = 0; k < dim; k++) {
        Res[ k * nDofs + i ] = -aResV[k][i].value();
        Res[ (dim + k) * nDofs + i ] = -aResD[k][i].value();
      }
    }

    for (int i = 0; i < nDofsP; i++) {
      Res[2 * dim * nDofs + i ] = -aResP[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extract and store the Jacobian

    Jac.resize(nDofsVDP * nDofsVDP);
    // define the dependent variables

    for (unsigned  k = 0; k < dim; k++) {
      s.dependent(&aResV[k][0], nDofs);
    }
    for (unsigned  k = 0; k < dim; k++) {
      s.dependent(&aResD[k][0], nDofs);
    }
    
    s.dependent(&aResP[0], nDofsP);

    // define the independent variables
    for (unsigned  k = 0; k < dim; k++) {
      s.independent(&solV[k][0], nDofs);
    }
    for (unsigned  k = 0; k < dim; k++) {
      s.independent(&solD[k][0], nDofs);
    }
    s.independent(&solP[0], nDofsP);

    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0] , true); // This is rowwise order.
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();
  // ***************** END ASSEMBLY *******************
}





double GetSolutionFluxes(MultiLevelSolution& mlSol){

  int  iproc, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  MyVector<double> pFlux(1,0);
  pFlux.stack();
  
  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* solution  = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);
  elem* myel =  msh->el;
  
  const unsigned dim = msh->GetDimension();
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

  vector< vector < double> >  sol(dim);
  vector< vector < double> > x(dim);
 
  const char varname[6][3] = {"U", "V", "W","DX", "DY", "DZ"};
  vector <unsigned> indVar(2 * dim);
  unsigned solType;

  for (unsigned ivar = 0; ivar < dim; ivar++) {
    for (unsigned k = 0; k < 2; k++) {
      indVar[ivar + k * dim] = mlSol.GetIndex(&varname[ivar + k * 3][0]);
    }
  }
  solType = mlSol.GetSolutionType(&varname[0][0]);
    
  
   std::vector < double > phi;
   std::vector < double > gradphi;
   //std::vector< double > xx(dim, 0.);
   double weight;
  
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    // loop on faces
    for (unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
      int faceNumber = myel->GetBoundaryIndex(iel, jface);
      // look for boundary faces
      if ( faceNumber == 3 ) {
          
        vector < double> normal(dim, 0);  
       
        unsigned nve = msh->GetElementFaceDofNumber(iel, jface, solType);
        const unsigned felt = msh->GetElementFaceType(iel, jface);
	
        for (unsigned d = 0; d < dim; d++) {
          x[d].resize(nve);
	      sol[d].resize(nve);
	    }
	
        for (unsigned i = 0; i < nve; i++) {
          unsigned int ilocal = msh->GetLocalFaceVertexIndex(iel, jface, i);
          unsigned idof = msh->GetSolutionDof(ilocal, iel, 2);
          for (unsigned d = 0; d < dim; d++) {
            x[d][i] = (*msh->_topology->_Sol[d])(idof) + (*solution->_Sol[indVar[d+dim]])(idof);;
	        sol[d][i] = (*solution->_Sol[indVar[d]])(idof);;
          }
        }

        for (unsigned igs = 0; igs < msh->_finiteElement[felt][solType]->GetGaussPointNumber(); igs++) {
          msh->_finiteElement[felt][solType]->JacobianSur(x, igs, weight, phi, gradphi, normal);
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
  
  for(int j = 0; j < nprocs; j++) {
    pFlux.broadcast(j);
    flux += pFlux[j]; 
    pFlux.clearBroadcast();
  } 
  
  return flux;
}



