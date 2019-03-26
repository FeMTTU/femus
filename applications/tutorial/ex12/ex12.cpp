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
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"


using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;
  
  if (!strcmp(SolName, "U")) { // strcmp compares two string in lexiographic sense. why !strcmp??? not false ==true?
    if (facename == 2) {
      value = -(x[1]-0.5) * (x[1]+0.5);
    }
    else if (facename == 3) {
      dirichlet = false;  
    }
  } 
  else if (!strcmp(SolName, "V")) {
     if (facename == 2) {
      value = 0.;
    }
    else if (facename == 3) {
      dirichlet = false;  
    }
  }
  else if (!strcmp(SolName, "W")) {
    value = 0.;
  } 
  else if (!strcmp(SolName, "P")) {
    dirichlet = false;
    value = 0.;
    if(facename == 3){
      value = 10.;
    }
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
  mlMsh.ReadCoarseMesh ( "./input/rectangle.neu", "seventh", scalingFactor );
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 6;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);
  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND);
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
  
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("NS");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");

  if (dim == 3) system.AddSolutionToSystemPDE("W");

  system.AddSolutionToSystemPDE("P");

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleBoussinesqAppoximation_AD);

  // initilaize and solve the system
  system.init();
  system.SetOuterSolver(PREONLY);
  system.MGsolve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

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
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
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

  //solution variable
  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object

  if (dim == 3) solVIndex[2] = mlSol->GetIndex("W");      // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"
 
  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object

  if (dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  vector < vector < adept::adouble > >  solV(dim);    // local solution
  vector < adept::adouble >  solP; // local solution

  vector< vector < adept::adouble > > aResV(dim);    // local redidual vector
  vector< adept::adouble > aResP; // local redidual vector

  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
    aResV[k].reserve(maxSize);
    coordX[k].reserve(maxSize);
  }

  solP.reserve(maxSize);
  aResP.reserve(maxSize);


  vector <double> phiV;  // local test function for velocity
  vector <double> phiV_x; // local test function first order partial derivatives

  phiV.reserve(maxSize);
  phiV_x.reserve(maxSize * dim);
  
  double* phiP; // local test function for the pressure
  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve((dim + 1) * maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve((dim + 1) * maxSize);

  vector < double > Jac;
  Jac.reserve((dim + 1) * maxSize * (dim + 1) * maxSize);

  RES->zero(); // Set to zero all the entries of the Global Residual vector
  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
        
    unsigned nDofsVP = dim * nDofsV + nDofsP;
    // resize local arrays
    sysDof.resize(nDofsVP);

    for (unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsV);
    }
    solP.resize(nDofsP);

    for (unsigned  k = 0; k < dim; k++) {
      aResV[k].assign(nDofsV,0.);
    }
    aResP.assign(nDofsP,0.);

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // local to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        sysDof[k * nDofsV + i] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // local to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      sysDof[dim * nDofsV + i ] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // local to global mapping between coordinates node and coordinate dof
      for (unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();
    
    // *** Face Gauss point loop (boundary Integral) ***
    for ( unsigned jface = 0; jface < msh->GetElementFaceNumber ( iel ); jface++ ) {
      int faceIndex = el->GetBoundaryIndex(iel, jface);
      // look for boundary faces
      if ( faceIndex == 3 ) {  
        const unsigned faceGeom = msh->GetElementFaceType ( iel, jface );
        unsigned faceDofs = msh->GetElementFaceDofNumber (iel, jface, solVType);
                    
        vector  < vector  <  double> > faceCoordinates ( dim ); // A matrix holding the face coordinates rowwise.
        for ( int k = 0; k < dim; k++ ) {
          faceCoordinates[k].resize (faceDofs);
        }
        for ( unsigned i = 0; i < faceDofs; i++ ) {
          unsigned inode = msh->GetLocalFaceVertexIndex ( iel, jface, i ); // face-to-element local node mapping.
          for ( unsigned k = 0; k < dim; k++ ) {
            faceCoordinates[k][i] =  coordX[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
          }
        }
        for ( unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solVType]->GetGaussPointNumber(); ig++ ) { 
            // We call the method GetGaussPointNumber from the object finiteElement in the mesh object msh. 
          vector < double> normal;
          msh->_finiteElement[faceGeom][solVType]->JacobianSur ( faceCoordinates, ig, weight, phiV, phiV_x, normal );
            
          vector< double > xg(dim,0.);
          for ( unsigned i = 0; i < faceDofs; i++ ) {
            for( unsigned k=0; k<dim; k++){
              xg[k] += phiV[i] * faceCoordinates[k][i]; // xg(ig)= \sum_{i=0}^faceDofs phi[i](xig) facecoordinates[i]     
            }
          }
          double tau; // a(u)*grad_u\cdot normal
          SetBoundaryCondition( xg, "P", tau, faceIndex, 0. ); // return tau
          // *** phi_i loop ***
          for ( unsigned i = 0; i < faceDofs; i++ ) {
            unsigned inode = msh->GetLocalFaceVertexIndex ( iel, jface, i );
            for(unsigned k=0;k<dim;k++){
              aResV[k][inode] +=  -phiV[i] * tau * normal[k] * weight;
            }
          }        
        }
      }
    }   

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x);
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

      vector < adept::adouble > solV_gss(dim, 0);
      vector < vector < adept::adouble > > gradSolV_gss(dim);

      for (unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].assign(dim,0.);
      }

      for (unsigned i = 0; i < nDofsV; i++) {
        for (unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += solV[k][i] * phiV[i];
        }
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += solV[k][i] * phiV_x[i * dim + j];
          }
        }
      }

      adept::adouble solP_gss = 0;
      for (unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }

      double nu = 1.;

      // *** phiV_i loop ***
      for (unsigned i = 0; i < nDofsV; i++) {
        vector < adept::adouble > NSV(dim, 0.);
        
        for (unsigned  k = 0; k < dim; k++) { //momentum equation in k 
          for (unsigned j = 0; j < dim; j++) { // second index j in each equation
            NSV[k]   +=  nu * phiV_x[i * dim + j] * (gradSolV_gss[k][j] + gradSolV_gss[j][k]); // laplace
            NSV[k]   +=  phiV[i] * (solV_gss[j] * gradSolV_gss[k][j]); // non-linear term
          }
          NSV[k] += -solP_gss * phiV_x[i * dim + k]; // pressure gradient
        }
        for (unsigned  k = 0; k < dim; k++) {
          aResV[k][i] += - NSV[k] * weight;
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
    Res.resize(nDofsVP);    //resize

    for (int i = 0; i < nDofsV; i++) {
      for (unsigned  k = 0; k < dim; k++) {
        Res[ k * nDofsV + i ] = -aResV[k][i].value();
      }
    }

    for (int i = 0; i < nDofsP; i++) {
      Res[ dim * nDofsV + i ] = -aResP[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extract and store the Jacobian

    Jac.resize(nDofsVP * nDofsVP);
    // define the dependent variables

    for (unsigned  k = 0; k < dim; k++) {
      s.dependent(&aResV[k][0], nDofsV);
    }
    s.dependent(&aResP[0], nDofsP);

    // define the independent variables
    for (unsigned  k = 0; k < dim; k++) {
      s.independent(&solV[k][0], nDofsV);
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




