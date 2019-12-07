
// Building the Conformal Minimization system.
void AssembleConformalMinimization(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  //  Extract pointers to the several objects that we are going to use.
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem> ("nProj");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Pointers to the mesh (level) object and elem object in mesh (level).
  Mesh *msh = ml_prob._ml_msh->GetLevel(level);
  elem *el = msh->el;

  // Pointers to the multilevel solution, solution (level) and equation (level).
  MultiLevelSolution *mlSol = ml_prob._ml_sol;
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

  unsigned solENVNIndex = mlSol->GetIndex("ENVN");
  unsigned solENVNType = mlSol->GetSolutionType(solENVNIndex);

  // Extract positions of Dx in ml_sol object.
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol->GetIndex("Dx1");
  solDxIndex[1] = mlSol->GetIndex("Dx2");
  solDxIndex[2] = mlSol->GetIndex("Dx3");

  // Extract finite element type for the solution.
  unsigned solType;
  solType = mlSol->GetSolutionType(solDxIndex[0]);

  // Local solution vectors for X, Dx, Xhat, XC.
  std::vector < double > solx[DIM];
  std::vector < double > solNx[DIM];
  std::vector < double > solDx[DIM];
  std::vector < double > xhat[DIM];
  std::vector < double > xc[DIM];

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get the poitions of Y in the ml_sol object.
  unsigned solNDxIndex[DIM];
  solNDxIndex[0] = mlSol->GetIndex("nDx1");
  solNDxIndex[1] = mlSol->GetIndex("nDx2");
  solNDxIndex[2] = mlSol->GetIndex("nDx3");

  // Get the positions of Y in the pdeSys object.
  unsigned solNDxPdeIndex[DIM];
  solNDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex("nDx1");
  solNDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex("nDx2");
  solNDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex("nDx3");

  // Local solution vectors for Nx and NDx.
  std::vector < double > solNDx[DIM];

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
  std::vector < double > solL;

  // Local-to-global pdeSys dofs.
  std::vector < int > SYSDOF;
  // Local residual vectors.
  vector< double > Res;
  // Local Jacobian matrix (ordered by column).
  vector < double > Jac;

  KK->zero();  // Zero all the entries of the Global Matrix
  RES->zero(); // Zero all the entries of the Global Residual

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Numer of solution element dofs.
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nxDofs  = msh->GetElementDofNumber(iel, solType);
    unsigned nLDofs  = msh->GetElementDofNumber(iel, solLType);

    // Resize local arrays.
    for(unsigned K = 0; K < DIM; K++) {

      xhat[K].resize(nxDofs);

      solDx[K].resize(nxDofs);
      solx[K].resize(nxDofs);

      solNx[K].resize(nxDofs);

      solNDx[K].resize(nxDofs);
      solNx[K].resize(nxDofs);

      solL.resize(nLDofs);
    }

    // Resize local arrays
    unsigned sizeAll = DIM * nxDofs + nLDofs;

    SYSDOF.resize(sizeAll);
    Res.assign(sizeAll, 0.);
    Jac.assign(sizeAll * sizeAll, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof(i, iel, solType);
      unsigned iXDof  = msh->GetSolutionDof(i, iel, xType);

      for(unsigned K = 0; K < DIM; K++) {
        xhat[K][i] = (*msh->_topology->_Sol[K])(iXDof);
        solDx[K][i] = (*sol->_Sol[solDxIndex[K]])(iDDof);
        solx[K][i] = xhat[K][i] + solDx[K][i];
        solNDx[K][i] = (*sol->_Sol[solNDxIndex[K]])(iDDof);
        solNx[K][i] = xhat[K][i] + solNDx[K][i];

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

    if(ielGeom == TRI) {

      xT[0][1] = 0.5;
      std::vector < unsigned > ENVN(3);
      std::vector < double > angle(3);

      for(unsigned j = 0; j < 3; j++) {
        unsigned jnode  = msh->GetSolutionDof(j, iel, solENVNType);
        ENVN[j] = (*sol->_Sol[solENVNIndex])(jnode);
        angle[j] = 2 * M_PI / ENVN[j];
      }


      if(conformalTriangleType == 1) {  //this works with moo two levels
        ChangeTriangleConfiguration1(ENVN, angle);
      }
      else if(conformalTriangleType == 2) {  //this works with mao
        ChangeTriangleConfiguration2(ENVN, angle);
      }
      else { //no change
        angle.assign(3, M_PI / 3.);
      }

      double l = xT[0][1] - xT[0][0];
      double d = l * sin(angle[0]) * sin(angle[1]) / sin(angle[0] + angle[1]);
      double scale = sqrt((sqrt(3.) / 2.) / (l * d));
      l = l * scale;
      d = d * scale;
      xT[0][1] = xT[0][0] + l;
      xT[0][2] = xT[0][0] + d / tan(angle[0]);
      xT[1][2] = d;

      //std::cout << l << " " << d<<" "<< angle[0] << " " << angle[1] <<" "<< angle[2] << " " << l * d <<" "<< xT[0][2]<< " " << xT[1][2]<<  std::endl;
    }


    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
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

      // Initialize and compute values of x, Dx, NDx, x_uv at the Gauss points.
      double solDxg[3] = {0., 0., 0.};
      double solNDxg[3] = {0., 0., 0.};

      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solMx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solNx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for(unsigned K = 0; K < DIM; K++) {
        for(unsigned i = 0; i < nxDofs; i++) {
          solDxg[K] += phix[i] * solDx[K][i];
          solNDxg[K] += phix[i] * solNDx[K][i];
        }
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * solx[K][i];
            solMx_uv[K][j]   += phix_uv[j][i] * (xhat[K][i] + 0.5 * ((1. + !O2conformal) * solDx[K][i] + O2conformal * solNDx[K][i]));
            solNx_uv[K][j]   += phix_uv[j][i] * (xhat[K][i] + solNDx[K][i]);
          }
        }
      }

      ///////// ADDED THIS /////////
      double solLg = 0.;
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
      normal[0] = (solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) / sqrt(detg);
      normal[1] = (solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt(detg);
      normal[2] = (solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt(detg);

      double normalMSqrtDetg[DIM];
      normalMSqrtDetg[0] = (solMx_uv[1][0] * solMx_uv[2][1] - solMx_uv[2][0] * solMx_uv[1][1]);
      normalMSqrtDetg[1] = (solMx_uv[2][0] * solMx_uv[0][1] - solMx_uv[0][0] * solMx_uv[2][1]);
      normalMSqrtDetg[2] = (solMx_uv[0][0] * solMx_uv[1][1] - solMx_uv[1][0] * solMx_uv[0][1]);

      // Computing the "reduced Jacobian" g^{ij}X_j .
      double Jir[dim][DIM] = {{0., 0., 0.}, {0., 0., 0.}};
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned J = 0; J < DIM; J++) {
          for(unsigned k = 0; k < dim; k++) {
            Jir[i][J] += gi[i][k] * solx_uv[J][k];
          }
        }
      }

      // Initializing tangential gradients of X and W (new, middle, old).
      double solx_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      double solNx_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      for(unsigned I = 0; I < DIM; I++) {
        for(unsigned J = 0; J < DIM; J++) {
          for(unsigned k = 0; k < dim; k++) {
            solx_Xtan[I][J] += solx_uv[I][k] * Jir[k][J];
            solNx_Xtan[I][J] += solNx_uv[I][k] * Jir[k][J];
          }
        }
      }

      // Define and compute gradients of test functions for X and W.
      std::vector < double > phix_Xtan[DIM];

      for(unsigned J = 0; J < DIM; J++) {
        phix_Xtan[J].assign(nxDofs, 0.);

        for(unsigned inode  = 0; inode < nxDofs; inode++) {
          for(unsigned k = 0; k < dim; k++) {
            phix_Xtan[J][inode] += phix_uv[k][inode] * Jir[k][J];
          }
        }
      }

      // Discretize the equation \delta CD = 0 on the basis d/du, d/dv.
      double V[DIM];
      V[0] = solNx_uv[0][1] - normal[1] * solNx_uv[2][0] + normal[2] * solNx_uv[1][0];
      V[1] = solNx_uv[1][1] - normal[2] * solNx_uv[0][0] + normal[0] * solNx_uv[2][0];
      V[2] = solNx_uv[2][1] - normal[0] * solNx_uv[1][0] + normal[1] * solNx_uv[0][0];

      double W[DIM];
      W[0] = solNx_uv[0][0] + normal[1] * solNx_uv[2][1] - normal[2] * solNx_uv[1][1];
      W[1] = solNx_uv[1][0] + normal[2] * solNx_uv[0][1] - normal[0] * solNx_uv[2][1];
      W[2] = solNx_uv[2][0] + normal[0] * solNx_uv[1][1] - normal[1] * solNx_uv[0][1];

      double Q[DIM][dim];
      Q[0][0] = (+ gi[1][1] * W[0]
                 + gi[0][0] * (normal[1] * V[2] - normal[2] * V[1])
                 + gi[0][1] * (normal[2] * W[1] - normal[1] * W[2] - V[0]));

      Q[1][0] = (+ gi[1][1] * W[1]
                 + gi[0][0] * (normal[2] * V[0] - normal[0] * V[2])
                 + gi[0][1] * (normal[0] * W[2] - normal[2] * W[0] - V[1]));

      Q[2][0] = (+ gi[1][1] * W[2]
                 + gi[0][0] * (normal[0] * V[1] - normal[1] * V[0])
                 + gi[0][1] * (normal[1] * W[0] - normal[0] * W[1] - V[2]));

      Q[0][1] = (+ gi[0][0] * V[0]
                 + gi[1][1] * (normal[2] * W[1] - normal[1] * W[2])
                 + gi[0][1] * (normal[1] * V[2] - normal[2] * V[1] - W[0]));

      Q[1][1] = (+ gi[0][0] * V[1]
                 + gi[1][1] * (normal[0] * W[2] - normal[2] * W[0])
                 + gi[0][1] * (normal[2] * V[0] - normal[0] * V[2] - W[1]));

      Q[2][1] = (+ gi[0][0] * V[2]
                 + gi[1][1] * (normal[1] * W[0] - normal[0] * W[1])
                 + gi[0][1] * (normal[0] * V[1] - normal[1] * V[0] - W[2]));

      // Compute new X minus old X dot N, for "reparametrization".
      double DnXmDxdotNSqrtDetg = 0.;
      for(unsigned K = 0; K < DIM; K++) {
        DnXmDxdotNSqrtDetg += (solDxg[K] - solNDxg[K]) * normalMSqrtDetg[K];
      }

      // Implement the Conformal Minimization equations.
      for(unsigned K = 0; K < DIM; K++) {
        for(unsigned i = 0; i < nxDofs; i++) {
          double M1 = 0.;

          for(unsigned k = 0; k < dim; k++) {
            M1 +=  Q[K][k] * phix_uv[k][i];
          }

          //Residual (Conformal Minimization + Lagrange Multiplier)
          double term = M1 * Area + solLg * phix[i] * normalMSqrtDetg[K] * Area2;
          Res[K * nxDofs + i] -= term;

          unsigned irow = K * nxDofs + i;
          unsigned istart = irow * sizeAll;
          // Jacobian (part1 Lagrange Multiplier: delta NDX * normal) 
          for(unsigned j = 0; j < nLDofs; j++) {
            Jac[istart + DIM * nxDofs + j] += phiL[j] * phix[i] * normalMSqrtDetg[K] * Area2;
          }

        }
      }

      for(unsigned j = 0; j < nxDofs; j++) {
        double DQ[DIM][DIM][dim];

        DQ[0][0][0] = ((+gi[1][1] + gi[0][0] * (normal[1] * normal[1] + normal[2] * normal[2]))          * phix_uv[0][j] +
                       (-gi[0][1] * (normal[0] * normal[0]))                                             * phix_uv[1][j]);
        DQ[0][1][0] = ((-gi[0][0] * normal[0] * normal[1])                                               * phix_uv[0][j] +
                       (-gi[0][1] * normal[0] * normal[1] - gi[1][1] * normal[2] - gi[0][0] * normal[2]) * phix_uv[1][j]);
        DQ[0][2][0] = ((-gi[0][0] * normal[0] * normal[2])                                               * phix_uv[0][j] +
                       (-gi[0][1] * normal[0] * normal[2] + gi[1][1] * normal[1] + gi[0][0] * normal[1]) * phix_uv[1][j]);

        DQ[0][0][1] = ((+gi[0][0] + gi[1][1] * (normal[1] * normal[1] + normal[2] * normal[2]))          * phix_uv[1][j] +
                       (-gi[0][1] * (normal[0] * normal[0]))                                             * phix_uv[0][j]);
        DQ[0][1][1] = ((-gi[1][1] * normal[0] * normal[1])                                               * phix_uv[1][j] +
                       (-gi[0][1] * normal[0] * normal[1] + gi[1][1] * normal[2] + gi[0][0] * normal[2]) * phix_uv[0][j]);
        DQ[0][2][1] = ((-gi[1][1] * normal[0] * normal[2])                                               * phix_uv[1][j] +
                       (-gi[0][1] * normal[0] * normal[2] - gi[1][1] * normal[1] - gi[0][0] * normal[1]) * phix_uv[0][j]);

        DQ[1][0][0] = ((-gi[0][0] * normal[1] * normal[0])                                               * phix_uv[0][j] +
                       (-gi[0][1] * normal[1] * normal[0] + gi[1][1] * normal[2] + gi[0][0] * normal[2]) * phix_uv[1][j]);
        DQ[1][1][0] = ((+gi[1][1] + gi[0][0] * (normal[2] * normal[2] + normal[0] * normal[0]))          * phix_uv[0][j] +
                       (-gi[0][1] * (normal[1] * normal[1]))                                             * phix_uv[1][j]);
        DQ[1][2][0] = ((-gi[0][0] * normal[1] * normal[2])                                               * phix_uv[0][j] +
                       (-gi[0][1] * normal[1] * normal[2] - gi[1][1] * normal[0] - gi[0][0] * normal[0]) * phix_uv[1][j]);

        DQ[1][0][1] = ((-gi[1][1] * normal[1] * normal[0])                                               * phix_uv[1][j] +
                       (-gi[0][1] * normal[1] * normal[0] - gi[1][1] * normal[2] - gi[0][0] * normal[2]) * phix_uv[0][j]);
        DQ[1][1][1] = ((+gi[0][0] + gi[1][1] * (normal[2] * normal[2] + normal[0] * normal[0]))          * phix_uv[1][j] +
                       (-gi[0][1] * (normal[1] * normal[1]))                                             * phix_uv[0][j]);
        DQ[1][2][1] = ((-gi[1][1] * normal[1] * normal[2])                                               * phix_uv[1][j] +
                       (-gi[0][1] * normal[1] * normal[2] + gi[1][1] * normal[0] + gi[0][0] * normal[0]) * phix_uv[0][j]);

        DQ[2][0][0] = ((-gi[0][0] * normal[2] * normal[0])                                               * phix_uv[0][j] +
                       (-gi[0][1] * normal[2] * normal[0] - gi[1][1] * normal[1] - gi[0][0] * normal[1]) * phix_uv[1][j]);
        DQ[2][1][0] = ((-gi[0][0] * normal[2] * normal[1])                                               * phix_uv[0][j] +
                       (-gi[0][1] * normal[2] * normal[1] + gi[1][1] * normal[0] + gi[0][0] * normal[0]) * phix_uv[1][j]);
        DQ[2][2][0] = ((+gi[1][1] + gi[0][0] * (normal[0] * normal[0] + normal[1] * normal[1]))          * phix_uv[0][j] +
                       (-gi[0][1] * (normal[2] * normal[2]))                                             * phix_uv[1][j]);

        DQ[2][0][1] = ((-gi[1][1] * normal[2] * normal[0])                                               * phix_uv[1][j] +
                       (-gi[0][1] * normal[2] * normal[0] + gi[1][1] * normal[1] + gi[0][0] * normal[1]) * phix_uv[0][j]);
        DQ[2][1][1] = ((-gi[1][1] * normal[2] * normal[1])                                               * phix_uv[1][j] +
                       (-gi[0][1] * normal[2] * normal[1] - gi[1][1] * normal[0] - gi[0][0] * normal[0]) * phix_uv[0][j]);
        DQ[2][2][1] = ((+gi[0][0] + gi[1][1] * (normal[0] * normal[0] + normal[1] * normal[1]))          * phix_uv[1][j] +
                       (-gi[0][1] * (normal[2] * normal[2]))                                             * phix_uv[0][j]);

        double DnormalMSqrtDetg[DIM][DIM];

        DnormalMSqrtDetg[0][0] =  0.;
        DnormalMSqrtDetg[0][1] =  0.5 * (solMx_uv[2][1] * phix_uv[0][j] - solMx_uv[2][0] * phix_uv[1][j]);
        DnormalMSqrtDetg[0][2] = -0.5 * (solMx_uv[1][1] * phix_uv[0][j] - solMx_uv[1][0] * phix_uv[1][j]);

        DnormalMSqrtDetg[1][0] = -0.5 * (solMx_uv[2][1] * phix_uv[0][j] - solMx_uv[2][0] * phix_uv[1][j]);
        DnormalMSqrtDetg[1][1] =  0.;
        DnormalMSqrtDetg[1][2] =  0.5 * (solMx_uv[0][1] * phix_uv[0][j] - solMx_uv[0][0] * phix_uv[1][j]);

        DnormalMSqrtDetg[2][0] =  0.5 * (solMx_uv[1][1] * phix_uv[0][j] - solMx_uv[1][0] * phix_uv[1][j]);
        DnormalMSqrtDetg[2][1] = -0.5 * (solMx_uv[0][1] * phix_uv[0][j] - solMx_uv[0][0] * phix_uv[1][j]);
        DnormalMSqrtDetg[2][2] =  0.;

        // Jacobian (Conformal Minimization + part2 Lagrange Multiplier: NDX * delta normal) 
        for(unsigned i = 0; i < nxDofs; i++) {
          for(unsigned K = 0; K < DIM; K++) {
            unsigned irow = K * nxDofs + i;
            unsigned istart = irow * sizeAll;
            for(unsigned J = 0; J < DIM; J++) {
              double term = 0.;
              for(unsigned k = 0; k < dim; k++) {
                term += DQ[K][J][k] * phix_uv[k][i];
              }
              Jac[istart + J * nxDofs + j] += term * Area + solLg * phix[i] * DnormalMSqrtDetg[K][J] * Area2;
            }
          }
        }
         
        // Jacobian (constraint: -\delta NDx * normal + (DX - NDx) * delta normal) 
        for(unsigned i = 0; i < nLDofs; i++) {
          unsigned irow = DIM * nxDofs + i;
          unsigned istart = irow * sizeAll;
          for(unsigned K = 0; K < DIM; K++) {
            double term = 0.;
            for(unsigned J = 0; J < DIM; J++) {
              term += (solDxg[J] - solNDxg[J]) * DnormalMSqrtDetg[J][K];
            }
            Jac[istart + K * nxDofs + j] += phiL[i] * (term - phix[j] * normalMSqrtDetg[K]) * Area2;
          }
        }
      }

      // Residual (constraint + eps) and Jacobian (eps)
      for(unsigned i = 0; i < nLDofs; i++) {
        unsigned irow = DIM * nxDofs + i;
        double term = phiL[i] * (DnXmDxdotNSqrtDetg * Area2 + eps * solLg * Area);
        Res[irow] -= term;
        unsigned istart = irow * sizeAll;
        for(unsigned j = 0; j < nLDofs; j++) {
          Jac[istart + DIM * nxDofs + j] += eps * phiL[i] * phiL[j] *  Area;
        }
      }

    } // end GAUSS POINT LOOP

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    //copy the value of the adept::adoube aRes in double Res and store

    RES->add_vector_blocked(Res, SYSDOF);
    KK->add_matrix_blocked(Jac, SYSDOF, SYSDOF);


  } //end ELEMENT LOOP for each process.

  RES->close();
  KK->close();

//     //VecView ( (static_cast<PetscVector*> (RES))->vec(),  PETSC_VIEWER_STDOUT_SELF);
//     MatView ( (static_cast<PetscMatrix*> (KK))->mat(), PETSC_VIEWER_STDOUT_SELF);
//
//     PetscViewer    viewer;
//     PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//     PetscObjectSetName ( (PetscObject) viewer, "PWilmore matrix");
//     PetscViewerPushFormat (viewer, PETSC_VIEWER_DRAW_LG);
//     MatView ( (static_cast<PetscMatrix*> (KK))->mat(), viewer);
//     double a;
//     std::cin >> a;

} // end AssembleO2ConformalMinimization.
