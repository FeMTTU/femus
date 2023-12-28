#ifndef __femus_equations_navier_stokes_hpp__
#define __femus_equations_navier_stokes_hpp__





#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "LinearEquationSolver.hpp"
#include "SparseMatrix.hpp"
#include "NumericVector.hpp"

#include "adept.h"


namespace femus {
    
    

void AssembleNavierStokes_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->GetMeshElements();  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
 
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27


  bool assembleMatrix = mlPdeSys->GetAssembleMatrix();
  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;
  if( assembleMatrix ) s.continue_recording();
  else s.pause_recording();



  //solution variable
  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object

  if (dim == 3) solVIndex[2] = mlSol->GetIndex("W");      // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"
 
  const unsigned solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  const unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object

  if (dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  const unsigned solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  std::vector < std::vector < adept::adouble > >  solV(dim);    // local solution
  std::vector < adept::adouble >  solP; // local solution

  std::vector < std::vector < adept::adouble > > aResV(dim);    // local redidual vector
  std::vector < adept::adouble > aResP; // local redidual vector

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
    aResV[k].reserve(maxSize);
    coordX[k].reserve(maxSize);
  }

  solP.reserve(maxSize);
  aResP.reserve(maxSize);


  std::vector <double> phiV;  // local test function for velocity
  std::vector <double> phiV_x; // local test function first order partial derivatives

  phiV.reserve(maxSize);
  phiV_x.reserve(maxSize * dim);
  
  double* phiP; // local test function for the pressure
  double weight; // gauss point weight

  std::vector < int > sysDof; // local to global pdeSys dofs
  sysDof.reserve((dim + 1) * maxSize);

  std::vector < double > Res; // local redidual vector
  Res.reserve((dim + 1) * maxSize);

  std::vector < double > Jac;
  Jac.reserve((dim + 1) * maxSize * (dim + 1) * maxSize);

  
  // // // //pressure stabilization ================
  // // // std::vector < double > Mass_p;
  // // // std::vector < adept::adouble > aResMassP;
  // // // aResMassP.reserve(maxSize);
  // // // Mass_p.reserve(maxSize*maxSize);
  
  
  
  
  RES->zero(); // Set to zero all the entries of the Global Residual vector
  if(assembleMatrix) { KK->zero(); } // Set to zero all the entries of the Global Matrix

  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs
        
    unsigned nDofsVP = dim * nDofsV + nDofsP;
    // resize local arrays
    sysDof.resize(nDofsVP);

    for (unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
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
    for (unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof

            for(unsigned k=0;k<dim;k++){
        coordX[k][i] = (*msh->GetTopology()->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
        }
      }
      



    // // // //pressure stabilization ================
    // // // aResMassP.resize(nDofsP);    //resize
    // // // std::fill(aResMassP.begin(), aResMassP.end(), 0);    //set aRes to zero
      
      
      
      // start a new recording of all the operations involving adept::adouble variables
     if(assembleMatrix) {  s.new_recording();  }



    /// *** Boundary Integral for the "P" part   @todo add the "normal gradient" part  - BEGIN ************
    for ( unsigned jface = 0; jface < msh->GetElementFaceNumber ( iel ); jface++ ) {

      
       const int boundary_index = msh->GetMeshElements()->GetFaceElementIndex(iel, jface);
       
       if ( boundary_index < 0) { //I am on the boundary

      
      const int faceIndex = el->GetBoundaryIndex(iel, jface);


          
      //Here the integral is \int_{\partial \Omega}  p n \cdot \vec{\phi},  where \vec{\phi} is the velocity test function

     
         
      // compute the face coordinates - BEGIN
        const unsigned faceGeom = msh->GetElementFaceType ( iel, jface );
        unsigned faceDofs = msh->GetElementFaceDofNumber (iel, jface, solVType);
                    
        std::vector < std::vector <  double> > faceCoordinates ( dim ); // A matrix holding the face coordinates rowwise.
        for ( int k = 0; k < dim; k++ ) {
          faceCoordinates[k].resize (faceDofs);
        }
        for ( unsigned i = 0; i < faceDofs; i++ ) {
          unsigned inode = msh->GetLocalFaceVertexIndex ( iel, jface, i ); // face-to-element local node mapping.
          for ( unsigned k = 0; k < dim; k++ ) {
            faceCoordinates[k][i] =  coordX[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
          }
        }
      // compute the face coordinates - END

      // compute the face coordinates, center - BEGIN
       std::vector < double > faceCoordinates_center(dim, 0.);
          for ( unsigned d = 0; d < faceCoordinates_center.size(); d++ ) {
            
                    for ( unsigned i = 0; i < faceDofs; i++ ) {
                        faceCoordinates_center[d] +=  faceCoordinates[d][i];
                    }
                    
            faceCoordinates_center[d] /= faceDofs;
                    
          }
      // compute the face coordinates, center - END
      
      // retrieve the dirichlet flag for each component - BEGIN
       double tau;
       
       std::vector <std::string > velocity_vec(dim);
      
      velocity_vec[0] = "U";
      velocity_vec[1] = "V";
      if (dim == 3) { velocity_vec[2] = "W"; }

      std::vector< bool > velocity_vec_is_dirichlet(dim);
     for (unsigned d = 0; d < velocity_vec_is_dirichlet.size(); d++) {
     
          velocity_vec_is_dirichlet[d] =  mlSol->GetBdcFunctionMLProb()(& ml_prob, faceCoordinates_center, velocity_vec[d].c_str(), tau, faceIndex, 0.);                     
     }
      // retrieve the dirichlet flag for each component - END

      
      /// compute the "element normal" by using the first quadrature point as auxiliary (@todo we assume non-curved boundary elements) - BEGIN
       const unsigned quadrature_point_placeholder_for_computing_element_normal = 0;
       std::vector < double > normal;
       msh->_finiteElement[faceGeom][solVType]->JacobianSur ( faceCoordinates, quadrature_point_placeholder_for_computing_element_normal, weight, phiV, phiV_x, normal );
     
      /// compute the "element normal" by using the first quadrature point as auxiliary (@todo we assume non-curved boundary elements) - END

      /// using the normal, pick the normal component of the velocity (@todo only works for faces parallel to the main axes) - BEGIN
      unsigned u_dot_n_component = 0;
     for (unsigned d = 0; d < dim; d++) {
       if ( fabs(normal[d]) >= 1.e-4) u_dot_n_component = d;
     }
       
      /// using the normal, pick the normal component of the velocity (@todo only works for faces parallel to the main axes) - END
       
        
      // look for boundary faces
      if ( velocity_vec_is_dirichlet[u_dot_n_component] == false  ) {
        
        
        for ( unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solVType]->GetGaussPointNumber(); ig++ ) { 
            // We call the method GetGaussPointNumber from the object finiteElement in the mesh object msh. 
          
          std::fill(normal.begin(), normal.end(), 0.);
          msh->_finiteElement[faceGeom][solVType]->JacobianSur ( faceCoordinates, ig, weight, phiV, phiV_x, normal );
            
          std::vector < double > xg(dim,0.);
          for ( unsigned i = 0; i < faceDofs; i++ ) {
            for( unsigned k=0; k<dim; k++){
              xg[k] += phiV[i] * faceCoordinates[k][i]; // xg(ig)= \sum_{i=0}^faceDofs phi[i](xig) facecoordinates[i]     
            }
          }

          mlSol->GetBdcFunctionMLProb()(& ml_prob, xg, "P", tau, faceIndex, 0. ); // return tau
          // *** phi_i loop ***
          for ( unsigned i = 0; i < faceDofs; i++ ) {
            unsigned inode = msh->GetLocalFaceVertexIndex ( iel, jface, i );
            for(unsigned k=0;k<dim;k++){
              aResV[k][inode] +=  phiV[i] * tau * normal[k] * weight;
            }
          }        
        }
      }
      
     } //end "is on boundary"
     
   }   
    /// *** Boundary Integral for the "P" part   @todo add the "normal gradient" part- END ************
    

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x);
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

        // evaluate the solution, the solution derivatives and the coordinates in the gauss point

      std::vector < adept::adouble > solV_gss(dim, 0);
      std::vector < std::vector < adept::adouble > > gradSolV_gss(dim);

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
        std::vector < adept::adouble > NSV(dim, 0.);
        
        for (unsigned  k = 0; k < dim; k++) { //momentum equation in k 
          for (unsigned j = 0; j < dim; j++) { // second index j in each equation
            NSV[k]   +=  nu * phiV_x[i * dim + j] * (gradSolV_gss[k][j] /* + gradSolV_gss[j][k] */ ); ///@todo see difference btw "Grad U + (Grad U)^T" and "Lapl"
            NSV[k]   +=  phiV[i] * (solV_gss[j] * gradSolV_gss[k][j]); // non-linear term
        }
          NSV[k] += -solP_gss * phiV_x[i * dim + k]; // pressure gradient
          }

          for (unsigned  k = 0; k < dim; k++) {
          aResV[k][i] += NSV[k] * weight;
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
        for (unsigned i = 0; i < nDofsP; i++) {
        for (int k = 0; k < dim; k++) {
          aResP[i] += - (gradSolV_gss[k][k]) * phiP[i]  * weight;
        }
      } // end phiP_i loop
      

        // // // //pressure stabilization ================
        // // // for (unsigned i = 0; i < nDofsP; i++) {
        // // //     aResMassP[i] += phiP[i] * solP_gss * weight;
        // // // }
        
      
      
      
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

    if(assembleMatrix){
      //Extarct and store the Jacobian

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
    s.jacobian(&Jac[0] , true);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

    
      // // // //pressure stabilization ================
      // // // Mass_p.resize(nDofsP * nDofsP);
      // // // s.dependent(&aResMassP[0], nDofsP);
      // // // s.independent(&solP[0], nDofsP);
      // // // s.jacobian(&Mass_p[0] , true);
      // // // KK->add_matrix_blocked(Mass_p, sysDof_P, sysDof_P);
      // // // s.clear_independents();
      // // // s.clear_dependents();
    
    }
    
    
  } //end element loop for each process

  RES->close();

  if(assembleMatrix)
  KK->close();

  // ***************** END ASSEMBLY *******************
}




}



#endif
