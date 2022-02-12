#ifndef __femus_include_CFDSteadyStateAssembly_hpp__
#define __femus_include_CFDSteadyStateAssembly_hpp__
#endif

//#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "MultiLevelSolution.hpp"
#include "adept.h"


namespace femus {
    
    // ***************** BEGIN ASSEMBLY *******************

  void CFDSteadyStateAssembly(MultiLevelProblem& mlProb)
{
  //  mlProb is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem



  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &mlProb.get_system<NonLinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*           msh         = mlProb._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*           el          = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*   mlSol   = mlProb._ml_sol;  // pointer to the multilevel solution object
  Solution*   sol         = mlProb._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level];  // pointer to the equation (level) object


  bool assembleMatrix = mlPdeSys->GetAssembleMatrix();
  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;
  if (assembleMatrix) s.continue_recording();
  else s.pause_recording();


  SparseMatrix*   KK          = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*  RES         = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

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


  vector <double> phiV;  // local test function
  vector <double> gradPhiV; // local test function first order partial derivatives
  vector <double> nablaPhiV; // local test function second order partial derivatives

  phiV.reserve(maxSize);
  gradPhiV.reserve(maxSize * dim);
  nablaPhiV.reserve(maxSize * dim2);

  double* phiP;
  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve((dim + 1) *maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve((dim + 1) *maxSize);

  vector < double > Jac;
  Jac.reserve((dim + 1) *maxSize * (dim + 1) *maxSize);

  if (assembleMatrix) KK->zero();

  //BEGIN element loop for each process
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType(iel);

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs

    unsigned nDofsTVP = dim * nDofsV + nDofsP;
    // resize local arrays
    sysDof.resize(nDofsTVP);

    for (unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
    }
    solP.resize(nDofsP);

    for (unsigned  k = 0; k < dim; k++) {
      aResV[k].assign(nDofsV, 0.);   //resize
    }
    aResP.assign(nDofsP, 0.);   //resize

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof
      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        sysDof[i + k * nDofsV] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // global to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      sysDof[i + dim * nDofsV] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof
      for (unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    //BEGIN a new recording of all the operations involving adept::adouble variables
    if (assembleMatrix) s.new_recording();

    //BEGIN Gauss point loop
    for (unsigned ig = 0; ig < msh->_finiteElement[ielType][solVType]->GetGaussPointNumber(); ig++) {
      msh->_finiteElement[ielType][solVType]->Jacobian(coordX, ig, weight, phiV, gradPhiV, nablaPhiV);
      phiP = msh->_finiteElement[ielType][solPType]->GetPhi(ig);

      // evaluate the solution, the solution derivatives and the coordinates in the Gauss point

      vector < adept::adouble > solVig(dim, 0);
      vector < vector < adept::adouble > > gradSolVig(dim);

      for (unsigned  k = 0; k < dim; k++) {
        gradSolVig[k].assign(dim, 0);
      }

      for (unsigned i = 0; i < nDofsV; i++) {
        for (unsigned  k = 0; k < dim; k++) {
          solVig[k] += phiV[i] * solV[k][i];
        }

        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
            gradSolVig[k][j] += gradPhiV[i * dim + j] * solV[k][i];
          }
        }
      }

      adept::adouble solPig = 0;

      for (unsigned i = 0; i < nDofsP; i++) {
        solPig += phiP[i] * solP[i];
      }

      //BEGIN phiV loop (momentum)
      //double nu = 3.5 * 1.0e-6;
      double nu = 0.1;
      for (unsigned i = 0; i < nDofsV; i++) {
        vector < adept::adouble > NSV(dim, 0.);

        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
            NSV[k]   +=  nu * gradPhiV[i * dim + j] * (gradSolVig[k][j] + gradSolVig[j][k]);
            NSV[k]   +=  phiV[i] * (solVig[j] * gradSolVig[k][j]);
          }
        }

        for (unsigned  k = 0; k < dim; k++) {
          NSV[k] += -solPig * gradPhiV[i * dim + k];
        }

        for (unsigned  k = 0; k < dim; k++) {
          aResV[k][i] += - NSV[k] * weight;
        }
      }
      //END phiV loop

      //BEGIN phiP loop (continuity)
      for (unsigned i = 0; i < nDofsP; i++) {
        for (int k = 0; k < dim; k++) {
          aResP[i] += - (gradSolVig[k][k]) * phiP[i]  * weight;
        }
      }
      //END phiP loop

    }
    //END Gauss point loop

    //BEGIN Extract and store the residual
    Res.resize(nDofsTVP);    //resize
    for (int i = 0; i < nDofsV; i++) {
      for (unsigned  k = 0; k < dim; k++) {
        Res[ i + k * nDofsV ] = -aResV[k][i].value();
      }
    }
    for (int i = 0; i < nDofsP; i++) {
      Res[ i + dim * nDofsV ] = -aResP[i].value();
    }
    RES->add_vector_blocked(Res, sysDof);
    //END Extract and store the residual

    //BEGIN Extract and store the jacobian
    if (assembleMatrix) {
      Jac.resize(nDofsTVP * nDofsTVP);
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
      s.jacobian(&Jac[0] , true);
      KK->add_matrix_blocked(Jac, sysDof, sysDof);

      s.clear_independents();
      s.clear_dependents();
    }
    //END extract and store the jacobian
  }
  //END element loop for each process

  RES->close();
  if (assembleMatrix) KK->close();

  // ***************** END ASSEMBLY *******************
}



}
