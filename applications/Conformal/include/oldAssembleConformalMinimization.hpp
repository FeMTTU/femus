

// Building the Conformal Minimization system.
void AssembleConformalMinimizationOld(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  Extract pointers to the several objects that we are going to use.
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem> ("conformal");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Pointers to the mesh (level) object and elem object in mesh (level).
  Mesh *msh = ml_prob._ml_msh->GetLevel(level);
  elem *el = msh->GetMeshElements();

  // Pointers to the multilevel solution, solution (level) and equation (level).
  MultiLevelSolution *mlSol = ml_prob._ml_sol;

  if(counter > 0 && !stopIterate) {
    UpdateMu(*mlSol);
  }

  Solution *sol = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level];

  // Pointers to global stiffness matrix and residual vector in pdeSys (level).
  SparseMatrix *KK = pdeSys->_KK;
  NumericVector *RES = pdeSys->_RES;

  // Convenience variables to keep track of the dimension.
  const unsigned  dim = 2;
  const unsigned  DIM = 3;
  //const unsigned  DIM = msh->GetDimension();;

  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // std::vector < double > phi;  // local test function for velocity
  // std::vector < double > phi_x; // local test function first order partial derivatives
  // double weight; // gauss point weight

  //Setting the reference elements to be equilateral triangles.
  std::vector < std::vector < double > > xT(2);
  xT[0].resize(7);
  xT[0][0] = -0.5;
  xT[0][1] = 0.5;
  xT[0][2] = 0.;
  xT[0][3] = 0.;
  xT[0][4] = 0.25;
  xT[0][5] = -0.25;
  xT[0][6] = 0.;

  xT[1].resize(7);
  xT[1][0] = 0.;
  xT[1][1] = 0.;
  xT[1][2] = sqrt(3.) / 2.;
  xT[1][3] = 0.;
  xT[1][4] = sqrt(3.) / 4.;
  xT[1][5] = sqrt(3.) / 4.;
  xT[1][6] = sqrt(3.) / 6.;
//
  std::vector< double > phi_uv0;
  std::vector< double > phi_uv1;

  std::vector< double > stdVectorPhi;
  std::vector< double > stdVectorPhi_uv;

  // Extract positions of Dx in ml_sol object.
  std::vector < unsigned >  solDxIndex(DIM);
  solDxIndex[0] = mlSol->GetIndex("Dx1");
  solDxIndex[1] = mlSol->GetIndex("Dx2");
  solDxIndex[2] = mlSol->GetIndex("Dx3");

  std::vector < unsigned > solMuIndex(dim);
  solMuIndex[0] = mlSol->GetIndex("mu1");
  solMuIndex[1] = mlSol->GetIndex("mu2");
  unsigned solType1 = mlSol->GetSolutionType(solMuIndex[0]);

  // Extract finite element type for the solution.
  unsigned solType;
  solType = mlSol->GetSolutionType(solDxIndex[0]);

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get the positions of Y in the pdeSys object.
  std::vector < unsigned > solDxPdeIndex(DIM);
  solDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex("Dx1");
  solDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex("Dx2");
  solDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex("Dx3");

  // Get the position of "Lambda1" in the ml_sol object.
  unsigned solLIndex;
  solLIndex = mlSol->GetIndex("Lambda1");

  // Get the finite element type for "Lambda1".
  unsigned solLType;
  solLType = mlSol->GetSolutionType(solLIndex);

  // Get the position of "Lambda1" in the pdeSys object.
  unsigned solLPdeIndex;
  solLPdeIndex = mlPdeSys->GetSolPdeIndex("Lambda1");

  // Local Lambda1 solution.
  std::vector < adept::adouble > solL;

  // Local solution vectors for Nx and NDx.
  std::vector < std::vector < adept::adouble > > solDx(DIM);
  std::vector < std::vector < adept::adouble > > solNx(DIM);
  std::vector < std::vector < double > > xHat(DIM);

  std::vector < std::vector < double > > solMu(dim);

  // Local-to-global pdeSys dofs.
  std::vector < int > SYSDOF;

  // Local residual vectors.
  vector< double > Res;
  std::vector< adept::adouble > aResL;
  std::vector < std::vector< adept::adouble > > aResDx(DIM);

  // Local Jacobian matrix (ordered by column).
  vector < double > Jac;

  KK->zero();  // Zero all the entries of the Global Matrix
  RES->zero(); // Zero all the entries of the Global Residual

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for(int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    // Numer of solution element dofs.
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nxDofs  = msh->GetElementDofNumber(iel, solType);
    unsigned nDofs1  = msh->GetElementDofNumber(iel, solType1);
    unsigned nLDofs  = msh->GetElementDofNumber(iel, solLType);

    // Resize local arrays.
    for(unsigned K = 0; K < DIM; K++) {
      solDx[K].resize(nxDofs);
      solNx[K].resize(nxDofs);
      xHat[K].resize(nxDofs);
      // solMu[K].resize(nDofs1);
    }
    for(unsigned k = 0; k < dim; k++) {
      solMu[k].resize(nDofs1);
    }
    solL.resize(nLDofs);

    // Resize local arrays
    unsigned sizeAll = DIM * nxDofs + nLDofs;

    SYSDOF.resize(sizeAll);
    Res.assign(sizeAll, 0.);
    Jac.assign(sizeAll * sizeAll, 0.);

    // // Resize local arrays
    //   SYSDOF.resize(DIM * nxDofs);
    //   Res.resize(DIM * nxDofs);

    for(unsigned K = 0; K < DIM; K++) {
      aResDx[K].assign(nxDofs, 0.);
    }
    aResL.resize(nLDofs);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nxDofs; i++) {
      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof(i, iel, solType);
      for(unsigned K = 0; K < DIM; K++) {
        solDx[K][i] = (*sol->_Sol[solDxIndex[K]])(iDDof);
        // Global-to-global mapping between NDx solution node and pdeSys dof.
        SYSDOF[ K * nxDofs + i] = pdeSys->GetSystemDof(solDxIndex[K], solDxPdeIndex[K], i, iel);
      }
    }

    // Local storage of global mapping and solution.
    for(unsigned i = 0; i < nLDofs; i++) {

      // Global-to-local mapping between Lambda solution node and solution dof.
      unsigned iLDof = msh->GetSolutionDof(i, iel, solLType);
      solL[i] = (*sol->_Sol[solLIndex])(iLDof);

      // Global-to-global mapping between Lambda solution node and pdeSys dof.
      SYSDOF[DIM * nxDofs + i] =
        pdeSys->GetSystemDof(solLIndex, solLPdeIndex, i, iel);
    }

    for(unsigned i = 0; i < nDofs1; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solType1);
      for(unsigned k = 0; k < dim; k++) {
        solMu[k][i] = (*sol->_Sol[solMuIndex[k]])(iDof);
      }
    }


    // start a new recording of all the operations involving adept variables.
    s.new_recording();

    for(unsigned i = 0; i < nxDofs; i++) {
      unsigned iXDof  = msh->GetSolutionDof(i, iel, xType);
      for(unsigned K = 0; K < DIM; K++) {
        xHat[K][i] = (*msh->GetTopology()->_Sol[K])(iXDof);
        solNx[K][i] = xHat[K][i] + solDx[K][i];
      }
    }

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phi1;  // local test function
      const double *phiL;  // local test function
      const double *phix_uv[dim]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // Get Gauss point weight, test function, and first order derivatives.
      if(ielGeom == QUAD) {
        phix = msh->_finiteElement[ielGeom][solType]->GetPhi(ig);
        phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi(ig);

        phix_uv[0] = msh->_finiteElement[ielGeom][solType]->GetDPhiDXi(ig);
        phix_uv[1] = msh->_finiteElement[ielGeom][solType]->GetDPhiDEta(ig);

        weight = msh->_finiteElement[ielGeom][solType]->GetGaussWeight(ig);
      }

      // Special adjustments for triangles.
      else {
        msh->_finiteElement[ielGeom][solType]->Jacobian(xT, ig, weight, stdVectorPhi, stdVectorPhi_uv);

        phix = &stdVectorPhi[0];
        phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi(ig);

        phi_uv0.resize(nxDofs);
        phi_uv1.resize(nxDofs);


        for(unsigned i = 0; i < nxDofs; i++) {
          phi_uv0[i] = stdVectorPhi_uv[i * dim];
          phi_uv1[i] = stdVectorPhi_uv[i * dim + 1];
        }

        phix_uv[0] = &phi_uv0[0];
        phix_uv[1] = &phi_uv1[0];

      }

      //msh->_finiteElement[ielGeom][solType]->Jacobian(xHat, ig, weight, phi, phi_x);

      // std::vector < std::vector < adept::adouble > > gradSolx(dim);
      //
      // for(unsigned  k = 0; k < dim; k++) {
      //   gradSolx[k].assign(dim, 0.);
      // }
      //
      // for(unsigned i = 0; i < nxDofs; i++) {
      //   for(unsigned j = 0; j < dim; j++) {
      //     for(unsigned  k = 0; k < dim; k++) {
      //       gradSolx[k][j] += solx[k][i] * phi_x[i * dim + j];
      //     }
      //   }
      // }

      phi1 = msh->_finiteElement[ielGeom][solType1]->GetPhi(ig);

      // Initialize and compute values of x, Dx, NDx, x_uv at the Gauss points.
      double xHatg[3] = {0., 0., 0. };
      adept::adouble solDxg[3] = {0., 0., 0.};
      adept::adouble solNxg[3] = {0. , 0. , 0. };
      double xHat_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solNx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for(unsigned K = 0; K < DIM; K++) {
        for(unsigned i = 0; i < nxDofs; i++) {
          xHatg[K] += phix[i] * xHat[K][i];
          solNxg[K] += phix[i] * solNx[K][i];
          solDxg[K] += phix[i] * solDx[K][i];
        }
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            xHat_uv[K][j] += phix_uv[j][i] * xHat[K][i];
            solNx_uv[K][j] += phix_uv[j][i] * solNx[K][i];
          }
        }

      }

      adept::adouble solLg = 0.;
      for(unsigned i = 0; i < nLDofs; i++) {
        solLg += phiL[i] * solL[i];
      }

      // Compute the metric, metric determinant, and area element.
      std::vector < std::vector < double > > g(dim);
      for(unsigned i = 0; i < dim; i++) g[i].assign(dim, 0.);

      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned K = 0; K < DIM; K++) {
            g[i][j] += xHat_uv[K][i] * xHat_uv[K][j];
          }
        }
      }
      adept::adouble detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      adept::adouble Area = weight * sqrt(detg);
      adept::adouble Area2 = weight;// Trick to give equal weight to each element.

      // Computing the unit normal vector N.
      adept::adouble normal[DIM];
      normal[0] = (xHat_uv[1][0] * xHat_uv[2][1]
                   - xHat_uv[2][0] * xHat_uv[1][1]) / sqrt(detg);
      normal[1] = (xHat_uv[2][0] * xHat_uv[0][1]
                   - xHat_uv[0][0] * xHat_uv[2][1]) / sqrt(detg);
      normal[2] = (xHat_uv[0][0] * xHat_uv[1][1]
                   - xHat_uv[1][0] * xHat_uv[0][1]) / sqrt(detg);

      adept::adouble normalMSqrtDetg[DIM];
      normalMSqrtDetg[0] = (xHat_uv[1][0] * xHat_uv[2][1] - xHat_uv[2][0] * xHat_uv[1][1]);
      normalMSqrtDetg[1] = (xHat_uv[2][0] * xHat_uv[0][1] - xHat_uv[0][0] * xHat_uv[2][1]);
      normalMSqrtDetg[2] = (xHat_uv[0][0] * xHat_uv[1][1] - xHat_uv[1][0] * xHat_uv[0][1]);

//       // Comment out for working code

      double mu[2] = {0., 0.};

      for(unsigned i = 0; i < nDofs1; i++) {
        for(unsigned k = 0; k < 2; k++) {
          mu[k] += phi1[i] * solMu[k][i];
        }
      }
      double muPlus = pow(1 + mu[0] * mu[0], 2) + mu[1] * mu[1];
      double muMinus = pow(1 - mu[0] * mu[0], 2) + mu[1] * mu[1];
      double muSqr = 2 * (mu[0] * mu[0] + mu[1] * mu[1] - 1.);
      // std::cout << mu[0] <<" ";//<<mu[1]<<" ";


//       if (counter == 0) mu[0] = 0.8;
//
//       for (unsigned K = 0; K < DIM; K++) {
//         if (counter > 0 && norm2Xz.value() > 0.) {
//           mu[K] += (1. / norm2Xz.value()) * XzBarXz_Bar[K].value();
//         }
//         //if(counter % 2 == 0) mu[K]*=1.01;
//         //else mu[K]/=1.01;
//       }

      //std::cout << mu[0] <<" "<< mu[1]<<" ";


      adept::adouble M[DIM][dim];

      M[0][0] = muMinus + muPlus * (normal[1] * normal[1] + normal[2] * normal[2]) * solNx_uv[0][0]
                - muPlus * normal[0] * (solNx_uv[1][0] * normal[1] + solNx_uv[2][0] * normal[2])
                - muSqr * (solNx_uv[2][1] * normal[1] - solNx_uv[1][1] * normal[2]);

      M[1][0] = muMinus + muPlus * (normal[2] * normal[2] + normal[0] * normal[0]) * solNx_uv[1][0]
                - muPlus * normal[1] * (solNx_uv[2][0] * normal[2] + solNx_uv[0][0] * normal[0])
                - muSqr * (solNx_uv[0][1] * normal[2] - solNx_uv[2][1] * normal[0]);

      M[2][0] = muMinus + muPlus * (normal[0] * normal[0] + normal[1] * normal[1]) * solNx_uv[2][0]
                - muPlus * normal[2] * (solNx_uv[0][0] * normal[0] + solNx_uv[1][0] * normal[1])
                - muSqr * (solNx_uv[1][1] * normal[0] - solNx_uv[0][1] * normal[1]);

      M[0][1] = muMinus + muPlus * (normal[1] * normal[1] + normal[2] * normal[2]) * solNx_uv[0][1]
                - muPlus * normal[0] * (solNx_uv[2][1] * normal[2] + solNx_uv[1][1] * normal[1])
                + muSqr * (solNx_uv[2][0] * normal[1] - solNx_uv[1][0] * normal[2]);

      M[1][1] = muMinus + muPlus * (normal[2] * normal[2] + normal[0] * normal[0]) * solNx_uv[1][1]
                - muPlus * normal[1] * (solNx_uv[0][1] * normal[0] + solNx_uv[2][1] * normal[2])
                + muSqr * (solNx_uv[0][0] * normal[2] - solNx_uv[2][0] * normal[0]);

      M[2][1] = muMinus + muPlus * (normal[0] * normal[0] + normal[1] * normal[1]) * solNx_uv[2][1]
                - muPlus * normal[2] * (solNx_uv[1][1] * normal[1] + solNx_uv[0][1] * normal[0])
                + muSqr * (solNx_uv[1][0] * normal[0] - solNx_uv[0][0] * normal[1]);

      adept::adouble DnXmDxdotN = 0.;
      for(unsigned K = 0; K < DIM; K++) {
        DnXmDxdotN += (solDxg[K]) * normalMSqrtDetg[K];
      }

      // Lagrange multiplier equation (with trick).
      for(unsigned i = 0; i < nLDofs; i++) {
        aResL[i] += phiL[i] * (DnXmDxdotN + eps * solL[i]) * Area2; // no2
      }

      // Implement the Conformal Minimization equations.
      for(unsigned K = 0; K < DIM; K++) {
        for(unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;
          for(unsigned j = 0; j < dim; j++) {
            term1 += M[K][j] * phix_uv[j][i];
          }
          // Conformal energy equation (with trick).
          aResDx[K][i] += (term1 + solLg * phix[i] * normalMSqrtDetg[K]) * Area2;
        }
      }
    } // end GAUSS POINT LOOP

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    //copy the value of the adept::adoube aRes in double Res and store

    for(int K = 0; K < DIM; K++) {
      for(int i = 0; i < nxDofs; i++) {
        Res[ K * nxDofs + i] = - aResDx[K][i].value(); // for X
      }
    }

    for(int i = 0; i < nLDofs; i++) {
      Res[DIM * nxDofs + i] = - aResL[i].value(); // for Lambda
    }

    RES->add_vector_blocked(Res, SYSDOF);

    // Resize Jacobian.
    Jac.resize((DIM * nxDofs + nLDofs) * (DIM * nxDofs + nLDofs));

    // Define the dependent variables.
    for(int K = 0; K < DIM; K++) {
      s.dependent(&aResDx[K][0], nxDofs);
    }
    s.dependent(&aResL[0], nLDofs);

    // Define the independent variables.
    for(int K = 0; K < DIM; K++) {
      s.independent(&solDx[K][0], nxDofs);
    }
    s.independent(&solL[0], nLDofs);

    // Get the jacobian matrix (ordered by row).
    s.jacobian(&Jac[0], true);

    KK->add_matrix_blocked(Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } //end ELEMENT LOOP for each process.

  RES->close();
  KK->close();

  counter++;

} // end AssembleConformalMinimization.




void AssembleShearMinimization(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< LinearImplicitSystem> ("conformal");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh *msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem *el = msh->GetMeshElements();  // pointer to the elem object in msh (level)

  MultiLevelSolution *mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  SparseMatrix *KK = pdeSys->_KK;  // pointer to the global stiffness matrix object in pdeSys (level)
  NumericVector *RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension();

  std::vector < double > phi;  // local test function for velocity
  std::vector <adept::adouble> phi_x; // local test function first order partial derivatives
  adept::adouble weight; // gauss point weight

  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  std::vector < unsigned > solDxIndex(dim);
  solDxIndex[0] = mlSol->GetIndex("Dx1");  // get the position of "DX" in the ml_sol object
  solDxIndex[1] = mlSol->GetIndex("Dx2");  // get the position of "DY" in the ml_sol object
  if(dim == 3) solDxIndex[2] = mlSol->GetIndex("Dx3");   // get the position of "DY" in the ml_sol object

  unsigned solType;
  solType = mlSol->GetSolutionType(solDxIndex[0]);   // get the finite element type for "U"

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector < unsigned > solDxPdeIndex(dim);
  solDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex("Dx1");    // get the position of "Dx1" in the pdeSys object
  solDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex("Dx2");    // get the position of "Dx2" in the pdeSys object
  if(dim == 3) solDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex("Dx3");    // get the position of "Dx3" in the pdeSys object

  std::vector < std::vector < adept::adouble > > solDx(dim);  // local Y solution
  std::vector < std::vector < adept::adouble > > x(dim);

  std::vector< int > SYSDOF; // local to global pdeSys dofs

  vector< double > Res; // local redidual vector
  std::vector< adept::adouble > aResDx[dim]; // local redidual vector

  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)

  KK->zero();  // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual

  // element loop: each process loops only on the elements that owns
  for(int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nxDofs  = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs

    for(unsigned k = 0; k < dim; k++) {
      solDx[k].resize(nxDofs);
      x[k].resize(nxDofs);
    }

    // resize local arrays
    SYSDOF.resize(dim * nxDofs);
    Res.resize(dim * nxDofs);        //resize

    for(unsigned k = 0; k < dim; k++) {
      aResDx[k].assign(nxDofs, 0.);   //resize and zet to zero
    }


    // local storage of global mapping and solution
    for(unsigned i = 0; i < nxDofs; i++) {
      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof(i, iel, solType);
      for(unsigned k = 0; k < dim; k++) {
        solDx[k][i] = (*sol->_Sol[solDxIndex[k]])(iDDof);
        // Global-to-global mapping between NDx solution node and pdeSys dof.
        SYSDOF[ k * nxDofs + i] = pdeSys->GetSystemDof(solDxIndex[k], solDxPdeIndex[k], i, iel);
      }
    }

    // start a new recording of all the operations involving adept variables.
    s.new_recording();
    for(unsigned i = 0; i < nxDofs; i++) {
      unsigned iXDof  = msh->GetSolutionDof(i, iel, xType);
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->GetTopology()->_Sol[k])(iXDof);
      }
    }

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielGeom][solType]->Jacobian(x, ig, weight, phi, phi_x);

      std::vector < std::vector < adept::adouble > > gradSolDx(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolDx[k].assign(dim, 0.);
      }

      for(unsigned i = 0; i < nxDofs; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolDx[k][j] += (x[k][i] + solDx[k][i]) * phi_x[i * dim + j];
          }
        }
      }

      for(unsigned i = 0; i < nxDofs; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          adept::adouble term = 0.;
          term  +=  phi_x[i * dim + k] * (gradSolDx[k][k]);
          aResDx[k][i] += term * weight;
        }
      }
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store


    for(int k = 0; k < dim; k++) {
      for(int i = 0; i < nxDofs; i++) {
        Res[ k * nxDofs + i] = -aResDx[k][i].value();
      }
    }

    RES->add_vector_blocked(Res, SYSDOF);

    Jac.resize((dim * nxDofs) * (dim * nxDofs));

    // define the dependent variables

    for(int k = 0; k < dim; k++) {
      s.dependent(&aResDx[k][0], nxDofs);
    }

    // define the dependent variables

    for(int k = 0; k < dim; k++) {
      s.independent(&solDx[k][0], nxDofs);
    }

    // get the jacobian matrix (ordered by row)
    s.jacobian(&Jac[0], true);

    KK->add_matrix_blocked(Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();

  counter++;

  // ***************** END ASSEMBLY *******************
}

// Building the Conformal Minimization system.
void AssembleConformalO1Minimization(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  Extract pointers to the several objects that we are going to use.
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem> ("conformal");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Pointers to the mesh (level) object and elem object in mesh (level).
  Mesh *msh = ml_prob._ml_msh->GetLevel(level);
  elem *el = msh->GetMeshElements();

  // Pointers to the multilevel solution, solution (level) and equation (level).
  MultiLevelSolution *mlSol = ml_prob._ml_sol;

  if(counter > 0 && !stopIterate) {
    UpdateMu(*mlSol);
  }

  Solution *sol = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level];

  // Pointers to global stiffness matrix and residual vector in pdeSys (level).
  SparseMatrix *KK = pdeSys->_KK;
  NumericVector *RES = pdeSys->_RES;

  // Convenience variables to keep track of the dimension.
  const unsigned  dim = 2;
  const unsigned  DIM = 3;

  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // Setting the reference elements to be equilateral triangles.
  std::vector < std::vector < double > > xT(2);
  xT[0].resize(3);
  xT[0][0] = -0.5;
  xT[0][1] = 0.5;
  xT[0][2] = 0.;

  xT[1].resize(3);
  xT[1][0] = 0.;
  xT[1][1] = 0.;
  xT[1][2] = sqrt(3.) / 2.;

  std::vector<double> phi_uv0;
  std::vector<double> phi_uv1;

  std::vector< double > stdVectorPhi;
  std::vector< double > stdVectorPhi_uv;

  // Extract positions of Dx in ml_sol object.
//   unsigned solDxIndex[DIM];
//   solDxIndex[0] = mlSol->GetIndex ("Dx1");
//   solDxIndex[1] = mlSol->GetIndex ("Dx2");
//   solDxIndex[2] = mlSol->GetIndex ("Dx3");

  // Extract finite element type for the solution.


  // Local solution vectors for X, Dx, Xhat, XC.
  std::vector < double > solx[DIM];
  //std::vector < double > solDx[DIM];
  std::vector < double > xhat[DIM];
  std::vector < double > xc[DIM];

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get the poitions of Y in the ml_sol object.
  unsigned solNDxIndex[DIM];
  solNDxIndex[0] = mlSol->GetIndex("Dx1");
  solNDxIndex[1] = mlSol->GetIndex("Dx2");
  solNDxIndex[2] = mlSol->GetIndex("Dx3");

  unsigned solType;
  solType = mlSol->GetSolutionType(solNDxIndex[0]);

  // Get the positions of Y in the pdeSys object.
  unsigned solNDxPdeIndex[DIM];
  solNDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex("Dx1");
  solNDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex("Dx2");
  solNDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex("Dx3");

  // Local solution vectors for Nx and NDx.
  std::vector < adept::adouble > solNDx[DIM];
  std::vector < adept::adouble > solNx[DIM];

  // Get the position of "Lambda1" in the ml_sol object.
  unsigned solLIndex;
  solLIndex = mlSol->GetIndex("Lambda1");

  // Get the finite element type for "Lambda1".
  unsigned solLType;
  solLType = mlSol->GetSolutionType(solLIndex);

  // Get the position of "Lambda1" in the pdeSys object.
  unsigned solLPdeIndex;
  solLPdeIndex = mlPdeSys->GetSolPdeIndex("Lambda1");

  // Local Lambda1 solution.
  std::vector < adept::adouble > solL;


  std::vector < unsigned > solMuIndex(dim);
  solMuIndex[0] = mlSol->GetIndex("mu1");
  solMuIndex[1] = mlSol->GetIndex("mu2");
  unsigned solType1 = mlSol->GetSolutionType(solMuIndex[0]);

  std::vector < std::vector < double > > solMu(dim);

  // Local-to-global pdeSys dofs.
  std::vector < int > SYSDOF;

  // Local residual vectors.
  vector< double > Res;
  std::vector< adept::adouble > aResNDx[3];
  std::vector< adept::adouble > aResL;

  // Local Jacobian matrix (ordered by column).
  vector < double > Jac;

  KK->zero();  // Zero all the entries of the Global Matrix
  RES->zero(); // Zero all the entries of the Global Residual

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for(int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    // Numer of solution element dofs.
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nxDofs  = msh->GetElementDofNumber(iel, solType);
    unsigned nLDofs  = msh->GetElementDofNumber(iel, solLType);
    unsigned nDofs1  = msh->GetElementDofNumber(iel, solType1);

    // Resize local arrays.
    for(unsigned K = 0; K < DIM; K++) {

      xhat[K].resize(nxDofs);

      //solDx[K].resize (nxDofs);
      solx[K].resize(nxDofs);

      solNDx[K].resize(nxDofs);
      solNx[K].resize(nxDofs);

    }
    solL.resize(nLDofs);

    for(unsigned k = 0; k < dim; k++) {
      solMu[k].resize(nDofs1);
    }

    // Resize local arrays
    SYSDOF.resize(DIM * nxDofs + nLDofs);
    Res.resize(DIM * nxDofs + nLDofs);

    for(unsigned K = 0; K < DIM; K++) {
      aResNDx[K].assign(nxDofs, 0.);
    }
    aResL.assign(nLDofs, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof(i, iel, solType);
      unsigned iXDof  = msh->GetSolutionDof(i, iel, xType);

      for(unsigned K = 0; K < DIM; K++) {
        xhat[K][i] = (*msh->GetTopology()->_Sol[K])(iXDof);
        //solDx[K][i] = (*sol->_Sol[solDxIndex[K]]) (iDDof);
        solx[K][i] = xhat[K][i];// + solDx[K][i];
        solNDx[K][i] = (*sol->_Sol[solNDxIndex[K]])(iDDof);

        // Global-to-global mapping between NDx solution node and pdeSys dof.
        SYSDOF[ K * nxDofs + i] =
          pdeSys->GetSystemDof(solNDxIndex[K], solNDxPdeIndex[K], i, iel);
      }
    }

    // Local storage of global mapping and solution.
    for(unsigned i = 0; i < nLDofs; i++) {

      // Global-to-local mapping between Lambda solution node and solution dof.
      unsigned iLDof = msh->GetSolutionDof(i, iel, solLType);
      solL[i] = (*sol->_Sol[solLIndex])(iLDof);

      // Global-to-global mapping between Lambda solution node and pdeSys dof.
      SYSDOF[DIM * nxDofs + i] =
        pdeSys->GetSystemDof(solLIndex, solLPdeIndex, i, iel);
    }

    for(unsigned i = 0; i < nDofs1; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solType1);
      for(unsigned k = 0; k < dim; k++) {
        solMu[k][i] = (*sol->_Sol[solMuIndex[k]])(iDof);
      }
    }

    // start a new recording of all the operations involving adept variables.
    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function

      const double *phix_uv[dim]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // Get Gauss point weight, test function, and first order derivatives.
      if(ielGeom == QUAD) {
        phix = msh->_finiteElement[ielGeom][solType]->GetPhi(ig);


        phix_uv[0] = msh->_finiteElement[ielGeom][solType]->GetDPhiDXi(ig);
        phix_uv[1] = msh->_finiteElement[ielGeom][solType]->GetDPhiDEta(ig);

        weight = msh->_finiteElement[ielGeom][solType]->GetGaussWeight(ig);
      }

      // Special adjustments for triangles.
      else {
        msh->_finiteElement[ielGeom][solType]->Jacobian(xT, ig, weight, stdVectorPhi, stdVectorPhi_uv);

        phix = &stdVectorPhi[0];
///////// WHAT ABOUT PHIL? ////////////////////

        phi_uv0.resize(nxDofs);
        phi_uv1.resize(nxDofs);


        for(unsigned i = 0; i < nxDofs; i++) {
          phi_uv0[i] = stdVectorPhi_uv[i * dim];
          phi_uv1[i] = stdVectorPhi_uv[i * dim + 1];
        }

        phix_uv[0] = &phi_uv0[0];
        phix_uv[1] = &phi_uv1[0];

      }

      const double *phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi(ig);  // local test function
      const double *phi1 = msh->_finiteElement[ielGeom][solType1]->GetPhi(ig);  // local test function

      // Initialize and compute values of x, Dx, NDx, x_uv at the Gauss points.
      //double solDxg[3] = {0., 0., 0.};
      adept::adouble solNDxg[3] = {0., 0., 0.};

      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solNx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for(unsigned K = 0; K < DIM; K++) {
        for(unsigned i = 0; i < nxDofs; i++) {
          //solDxg[K] += phix[i] * solDx[K][i];
          solNDxg[K] += phix[i] * solNDx[K][i];
        }
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * solx[K][i];
            solNx_uv[K][j]   += phix_uv[j][i] * (xhat[K][i] + solNDx[K][i]);
          }
        }
      }

      ///////// ADDED THIS /////////
      adept::adouble solLg = 0.;
      for(unsigned i = 0; i < nLDofs; i++) {
        solLg += phiL[i] * solL[i];
      }

      // Compute the metric, metric determinant, and area element.
      double g[dim][dim] = {{0., 0.}, {0., 0.}};
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      double Area = weight * sqrt(detg);
      double Area2 = weight; // Trick to give equal weight to each element.

      // Compute the metric inverse.
      double gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;

      // Compute components of the unit normal N.
      double normal[DIM];
      normal[0] = 0.;//(solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) / sqrt(detg);
      normal[1] = 0.;//(solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt(detg);
      normal[2] = 1.;//(solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt(detg);

      //
      // // Discretize the equation \delta CD = 0 on the basis d/du, d/dv.
      // adept::adouble V[DIM];
      // V[0] = solNx_uv[0][1] - normal[1] * solNx_uv[2][0] + normal[2] * solNx_uv[1][0];
      // V[1] = solNx_uv[1][1] - normal[2] * solNx_uv[0][0] + normal[0] * solNx_uv[2][0];
      // V[2] = solNx_uv[2][1] - normal[0] * solNx_uv[1][0] + normal[1] * solNx_uv[0][0];
      //
      // adept::adouble W[DIM];
      // W[0] = solNx_uv[0][0] + normal[1] * solNx_uv[2][1] - normal[2] * solNx_uv[1][1];
      // W[1] = solNx_uv[1][0] + normal[2] * solNx_uv[0][1] - normal[0] * solNx_uv[2][1];
      // W[2] = solNx_uv[2][0] + normal[0] * solNx_uv[1][1] - normal[1] * solNx_uv[0][1];
      //
      // adept::adouble M[DIM][dim];
      // M[0][0] = W[0] - normal[2] * V[1] + normal[1] * V[2];
      // M[1][0] = W[1] - normal[0] * V[2] + normal[2] * V[0];
      // M[2][0] = W[2] - normal[1] * V[0] + normal[0] * V[1];
      //
      // M[0][1] = V[0] + normal[2] * W[1] - normal[1] * W[2];
      // M[1][1] = V[1] + normal[0] * W[2] - normal[2] * W[0];
      // M[2][1] = V[2] + normal[1] * W[0] - normal[0] * W[1];


      double mu[2] = {0., 0.};

      for(unsigned i = 0; i < nDofs1; i++) {
        for(unsigned k = 0; k < 2; k++) {
          mu[k] += phi1[i] * solMu[k][i];
        }
      }
      double muPlus = pow(1 + mu[0], 2) + mu[1] * mu[1];
      double muMinus = pow(1 - mu[0], 2) + mu[1] * mu[1];
      double mu1Sqr = mu[1] * mu[1];
      double sqrMin = ((1 - mu[0]) * (1 - mu[0]));
      double sqrPlus = ((1 + mu[0]) * (1 + mu[0]));
      double xMin = mu[1] * (1 - mu[0]);
      double xPlus = mu[1] * (1 + mu[0]);
      double muSqr = (mu[0] * mu[0] + mu[1] * mu[1] - 1.);
      // std::cout << mu[0] <<" ";//<<mu[1]<<" ";


//       if (counter == 0) mu[0] = 0.8;
//
//       for (unsigned K = 0; K < DIM; K++) {
//         if (counter > 0 && norm2Xz.value() > 0.) {
//           mu[K] += (1. / norm2Xz.value()) * XzBarXz_Bar[K].value();
//         }
//         //if(counter % 2 == 0) mu[K]*=1.01;
//         //else mu[K]/=1.01;
//       }

      //std::cout << mu[0] <<" "<< mu[1]<<" ";

      adept::adouble M2[DIM][dim];
      double muSqrPlus = 1 + mu[0] * mu[0] + mu[1] * mu[1];
      double muSqrMinus = mu[0] * mu[0] + mu[1] * mu[1] - 1;

      M2[0][0] = muSqrPlus * solNx_uv[0][0] + muSqrMinus * solNx_uv[1][1];
      M2[1][0] = muSqrPlus * solNx_uv[1][0] - muSqrMinus * solNx_uv[0][1];
      M2[0][1] = muSqrPlus * solNx_uv[0][1] - muSqrMinus * solNx_uv[1][0];
      M2[1][1] = muSqrPlus * solNx_uv[1][1] + muSqrMinus * solNx_uv[0][0];


      adept::adouble M1[DIM][dim];

      M1[0][0] = (muMinus + muPlus * (normal[1] * normal[1] + normal[2] * normal[2])) * solNx_uv[0][0]
                 - muPlus * normal[0] * (solNx_uv[1][0] * normal[1] + solNx_uv[2][0] * normal[2])
                 - 2 * muSqr * (solNx_uv[2][1] * normal[1] - solNx_uv[1][1] * normal[2]);

      M1[1][0] = (muMinus + muPlus * (normal[2] * normal[2] + normal[0] * normal[0])) * solNx_uv[1][0]
                 - muPlus * normal[1] * (solNx_uv[2][0] * normal[2] + solNx_uv[0][0] * normal[0])
                 - 2 * muSqr * (solNx_uv[0][1] * normal[2] - solNx_uv[2][1] * normal[0]);

      M1[2][0] = (muMinus + muPlus * (normal[0] * normal[0] + normal[1] * normal[1])) * solNx_uv[2][0]
                 - muPlus * normal[2] * (solNx_uv[0][0] * normal[0] + solNx_uv[1][0] * normal[1])
                 - 2 * muSqr * (solNx_uv[1][1] * normal[0] - solNx_uv[0][1] * normal[1]);

      M1[0][1] = (muMinus + muPlus * (normal[1] * normal[1] + normal[2] * normal[2])) * solNx_uv[0][1]
                 - muPlus * normal[0] * (solNx_uv[2][1] * normal[2] + solNx_uv[1][1] * normal[1])
                 + 2 * muSqr * (solNx_uv[2][0] * normal[1] - solNx_uv[1][0] * normal[2]);

      M1[1][1] = (muMinus + muPlus * (normal[2] * normal[2] + normal[0] * normal[0])) * solNx_uv[1][1]
                 - muPlus * normal[1] * (solNx_uv[0][1] * normal[0] + solNx_uv[2][1] * normal[2])
                 + 2 * muSqr * (solNx_uv[0][0] * normal[2] - solNx_uv[2][0] * normal[0]);

      M1[2][1] = (muMinus + muPlus * (normal[0] * normal[0] + normal[1] * normal[1])) * solNx_uv[2][1]
                 - muPlus * normal[2] * (solNx_uv[1][1] * normal[1] + solNx_uv[0][1] * normal[0])
                 + 2 * muSqr * (solNx_uv[1][0] * normal[0] - solNx_uv[0][0] * normal[1]);


      adept::adouble N[DIM][dim];

      N[0][0] = (+ (sqrMin + mu1Sqr  * (normal[1] * normal[1] + normal[2] * normal[2])) * solNx_uv[0][0]
                 - (xMin   + xPlus   * (normal[1] * normal[1] + normal[2] * normal[2])) * solNx_uv[0][1])
                + (- (solNx_uv[1][0] * normal[1] + solNx_uv[2][0] * normal[2]) * mu1Sqr * normal[0]
                   + (solNx_uv[1][1] * normal[1] + solNx_uv[2][1] * normal[2]) * xPlus  * normal[0]
                   + (solNx_uv[1][1] * normal[2] - solNx_uv[2][1] * normal[1]) * muSqr);

      N[1][0] = (+ (sqrMin + mu1Sqr  * (normal[2] * normal[2] + normal[0] * normal[0])) * solNx_uv[1][0]
                 - (xMin   + xPlus   * (normal[2] * normal[2] + normal[0] * normal[0])) * solNx_uv[1][1])
                + (- (solNx_uv[2][0] * normal[2] + solNx_uv[0][0] * normal[0]) * mu1Sqr * normal[1]
                   + (solNx_uv[2][1] * normal[2] + solNx_uv[0][1] * normal[0]) * xPlus  * normal[1]
                   + (solNx_uv[2][1] * normal[0] - solNx_uv[0][1] * normal[2]) * muSqr);

      N[2][0] = (+ (sqrMin + mu1Sqr  * (normal[0] * normal[0] + normal[1] * normal[1])) * solNx_uv[2][0]
                 - (xMin   + xPlus   * (normal[0] * normal[0] + normal[1] * normal[1])) * solNx_uv[2][1])
                + (- (solNx_uv[0][0] * normal[0] + solNx_uv[1][0] * normal[1]) * mu1Sqr * normal[2]
                   + (solNx_uv[0][1] * normal[0] + solNx_uv[1][1] * normal[1]) * xPlus  * normal[2]
                   + (solNx_uv[0][1] * normal[1] - solNx_uv[1][1] * normal[0]) * muSqr);


      N[0][1] = (+ (mu1Sqr + sqrPlus * (normal[1] * normal[1] + normal[2] * normal[2])) * solNx_uv[0][1]
                 - (xMin   + xPlus   * (normal[1] * normal[1] + normal[2] * normal[2])) * solNx_uv[0][0])
                + (- (solNx_uv[1][1] * normal[1] + solNx_uv[2][1] * normal[2]) * sqrPlus * normal[0]
                   + (solNx_uv[1][0] * normal[1] + solNx_uv[2][0] * normal[2]) * xPlus   * normal[0]
                   - (solNx_uv[1][0] * normal[2] - solNx_uv[2][0] * normal[1]) * muSqr);

      N[1][1] = (+ (mu1Sqr + sqrPlus * (normal[2] * normal[2] + normal[0] * normal[0])) * solNx_uv[1][1]
                 - (xMin   + xPlus   * (normal[2] * normal[2] + normal[0] * normal[0])) * solNx_uv[1][0])
                + (- (solNx_uv[2][1] * normal[2] + solNx_uv[0][1] * normal[0]) * sqrPlus * normal[1]
                   + (solNx_uv[2][0] * normal[2] + solNx_uv[0][0] * normal[0]) * xPlus   * normal[1]
                   - (solNx_uv[2][0] * normal[0] - solNx_uv[0][0] * normal[2]) * muSqr);

      N[2][1] = (+ (mu1Sqr + sqrPlus * (normal[0] * normal[0] + normal[1] * normal[1])) * solNx_uv[2][1]
                 - (xMin   + xPlus   * (normal[0] * normal[0] + normal[1] * normal[1])) * solNx_uv[2][0])
                + (- (solNx_uv[0][1] * normal[0] + solNx_uv[1][1] * normal[1]) * sqrPlus * normal[2]
                   + (solNx_uv[0][0] * normal[0] + solNx_uv[1][0] * normal[1]) * xPlus   * normal[2]
                   - (solNx_uv[0][0] * normal[1] - solNx_uv[1][0] * normal[0]) * muSqr);


      adept::adouble V[DIM];
      V[0] = (1 - mu[0]) * solNx_uv[0][0] - (1 + mu[0]) * solNx_uv[1][1] + mu[1] * (solNx_uv[1][0] - solNx_uv[0][1]);
      V[1] = (1 - mu[0]) * solNx_uv[1][0] + (1 + mu[0]) * solNx_uv[0][1] - mu[1] * (solNx_uv[0][0] + solNx_uv[1][1]);
      V[2] = 0;


      adept::adouble M[DIM][dim];

      M[0][0] = (1 - mu[0]) * V[0] - mu[1] * V[1];
      M[1][0] = (1 - mu[0]) * V[1] + mu[1] * V[0];
      M[2][0] = 0;
      //M[0][0] = (1 - mu1) * V[0] - mu2 * V[1];
      //M[1][0] = (1 - mu1) * V[1] + mu2 * V[0];

      M[0][1] = (1 + mu[0]) * V[1] - mu[1] * V[0];
      M[1][1] = - (1 + mu[0]) * V[0] - mu[1] * V[1];
      M[2][1] = 0;

      // std::cout << 00 << " "<< M1[0][0]<<std::endl;

//       if(iel == 4 && ig == 1) {
//
//         std::cout << "derivatives " << std::endl;
//         std::cout << 1 << " " << solNx_uv[0][0] << std::endl;
//         std::cout << 2 << " " << solNx_uv[1][1] << std::endl;
//         std::cout << 3 << " " << solNx_uv[1][0] << std::endl;
//         std::cout << 4 << " " << solNx_uv[0][1] << std::endl;
//
//         std::cout << "symmetric " << std::endl;
//         std::cout << "00" << " " << M1[0][0] << " " << M2[0][0] << std::endl;
//         std::cout << "01" << " " << M1[1][0] << " " << M2[1][0] << std::endl;
//         std::cout << "01" << " " << M1[0][1] << " " << M2[0][1] << std::endl;
//         std::cout << "11" << " " << M1[1][1] << " " << M2[1][1] << std::endl;
//
//         std::cout << "asymmetric " << std::endl;
//         std::cout << "00" << " " << N[0][0] << " " << M[0][0] << std::endl;
//         std::cout << "01" << " " << N[1][0] << " " << M[1][0] << std::endl;
//         std::cout << "01" << " " << N[0][1] << " " << M[0][1] << std::endl;
//         std::cout << "11" << " " << N[1][1] << " " << M[1][1] << std::endl;

//         std::cout << "00" << " " << N[0][0] << " " << M[0][0] << std::endl;
//         std::cout << "10" << " " << N[1][0] << " " << M[1][0] << std::endl;
//         //std::cout << "20" << " "<< N[2][0]<<" "<< M[2][0]<<std::endl;
//         std::cout << "01" << " " << N[0][1] << " " << M[0][1] << std::endl;
//         std::cout << "11" << " " << N[1][1] << " " << M[1][1] << std::endl;
//         //std::cout << "21" << " "<< N[2][1]<<" "<< M[2][1]<<std::endl;

//     }



//       adept::adouble Q[DIM][dim];
//       Q[0][0] = (gi[1][1] * W[0]
//                  + gi[0][0] * (normal[1] * V[2] - normal[2] * V[1])
//                  + gi[0][1] * (normal[2] * W[1] - normal[1] * W[2] - V[0]));
//
//       Q[1][0] = (gi[1][1] * W[1]
//                  + gi[0][0] * (normal[2] * V[0] - normal[0] * V[2])
//                  + gi[0][1] * (normal[0] * W[2] - normal[2] * W[0] - V[1]));
//
//       Q[2][0] = (gi[1][1] * W[2]
//                  + gi[0][0] * (normal[0] * V[1] - normal[1] * V[0])
//                  + gi[0][1] * (normal[1] * W[0] - normal[0] * W[1] - V[2]));
//
//       Q[0][1] = (gi[0][0] * V[0]
//                  + gi[1][1] * (normal[2] * W[1] - normal[1] * W[2])
//                  + gi[0][1] * (normal[1] * V[2] - normal[2] * V[1] - W[0]));
//
//       Q[1][1] = (gi[0][0] * V[1]
//                  + gi[1][1] * (normal[0] * W[2] - normal[2] * W[0])
//                  + gi[0][1] * (normal[2] * V[0] - normal[0] * V[2] - W[1]));
//
//       Q[2][1] = (gi[0][0] * V[2]
//                  + gi[1][1] * (normal[1] * W[0] - normal[0] * W[1])
//                  + gi[0][1] * (normal[0] * V[1] - normal[1] * V[0] - W[2]));


      // Compute new X minus old X dot N, for "reparametrization".
      adept::adouble DnXmDxdotN = 0.;
      for(unsigned K = 0; K < DIM; K++) {
        DnXmDxdotN += (/*solDxg[K]*/ - solNDxg[K]) * normal[K];
      }


      // Implement the Conformal Minimization equations.
      for(unsigned K = 0; K < DIM; K++) {
        for(unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;

          for(unsigned j = 0; j < dim; j++) {
            term1 += N[K][j] * phix_uv[j][i]; //asymmetric
            //term1 += M1[K][j] * phix_uv[j][i]; //symmetric
            //term1 += Q[K][j] * phix_uv[j][i];
          }

          // Conformal energy equation (with trick).
          aResNDx[K][i] += term1 * Area2
                           // + timederiv * (solNDxg[K] - solDxg[K]) * phix[i] * Area2
                           + solL[0] * phix[i] * normal[K] * Area; //2 occasionally better???
        }
      }

      // Lagrange multiplier equation (with trick).
      for(unsigned i = 0; i < nLDofs; i++) {
        aResL[i] += phiL[i] * (DnXmDxdotN + eps * solL[i]) * Area; // no2
      }
      //aResL[0] += (DnXmDxdotN + eps * solL[0]) * Area;

//       if(iel == 4 && ig == 1) {
//         for(unsigned i = 0; i < nxDofs; i++) {
//           std::cout <<  mu[0] << " " << mu[1] << " " << aResNDx[0][i] << " " << aResNDx[1][i] << "\n";
//         }
//       }

    } // end GAUSS POINT LOOP

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    //copy the value of the adept::adoube aRes in double Res and store

    for(int K = 0; K < DIM; K++) {
      for(int i = 0; i < nxDofs; i++) {
        Res[ K * nxDofs + i] = -aResNDx[K][i].value();
      }
    }

    for(int i = 0; i < nLDofs; i++) {
      Res[DIM * nxDofs + i] = - aResL[i].value();
    }
    RES->add_vector_blocked(Res, SYSDOF);

    // Resize Jacobian.
    Jac.resize((DIM * nxDofs + nLDofs) * (DIM * nxDofs + nLDofs));

    // Define the dependent variables.
    for(int K = 0; K < DIM; K++) {
      s.dependent(&aResNDx[K][0], nxDofs);
    }
    s.dependent(&aResL[0], nLDofs);


    // Define the independent variables.
    for(int K = 0; K < DIM; K++) {
      s.independent(&solNDx[K][0], nxDofs);
    }
    s.independent(&solL[0], nLDofs);

    // Get the jacobian matrix (ordered by row).
    s.jacobian(&Jac[0], true);

    KK->add_matrix_blocked(Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } //end ELEMENT LOOP for each process.

  RES->close();
  KK->close();

  counter++;

} // end AssembleConformalMinimization.
