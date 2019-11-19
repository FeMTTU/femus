
using namespace femus;

// Building the Conformal Minimization system.
void AssembleConformalMinimization (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  Extract pointers to the several objects that we are going to use.
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem> ("nProj");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Pointers to the mesh (level) object and elem object in mesh (level).
  Mesh *msh = ml_prob._ml_msh->GetLevel (level);
  elem *el = msh->el;

  // Pointers to the multilevel solution, solution (level) and equation (level).
  MultiLevelSolution *mlSol = ml_prob._ml_sol;
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel (level);
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
  std::vector < std::vector < double > > xT (2);
  xT[0].resize (3);
  xT[0][0] = -0.5;
  xT[0][1] = 0.5;
  xT[0][2] = 0.;

  xT[1].resize (3);
  xT[1][0] = 0.;
  xT[1][1] = 0.;
  xT[1][2] = sqrt (3.) / 2.;

  std::vector<double> phi_uv0;
  std::vector<double> phi_uv1;

  std::vector< double > stdVectorPhi;
  std::vector< double > stdVectorPhi_uv;

  // Extract positions of Dx in ml_sol object.
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol->GetIndex ("Dx1");
  solDxIndex[1] = mlSol->GetIndex ("Dx2");
  solDxIndex[2] = mlSol->GetIndex ("Dx3");

  // Extract finite element type for the solution.
  unsigned solType;
  solType = mlSol->GetSolutionType (solDxIndex[0]);

  // Local solution vectors for X, Dx, Xhat, XC.
  std::vector < double > solx[DIM];
  std::vector < double > solDx[DIM];
  std::vector < double > xhat[DIM];
  std::vector < double > xc[DIM];

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get the poitions of Y in the ml_sol object.
  unsigned solNDxIndex[DIM];
  solNDxIndex[0] = mlSol->GetIndex ("nDx1");
  solNDxIndex[1] = mlSol->GetIndex ("nDx2");
  solNDxIndex[2] = mlSol->GetIndex ("nDx3");

  // Get the positions of Y in the pdeSys object.
  unsigned solNDxPdeIndex[DIM];
  solNDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("nDx1");
  solNDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("nDx2");
  solNDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("nDx3");

  // Local solution vectors for Nx and NDx.
  std::vector < adept::adouble > solNDx[DIM];
  std::vector < adept::adouble > solNx[DIM];

  // Get the position of "Lambda1" in the ml_sol object.
  unsigned solLIndex;
  solLIndex = mlSol->GetIndex ("Lambda1");

  // Get the finite element type for "Lambda1".
  unsigned solLType;
  solLType = mlSol->GetSolutionType (solLIndex);

  // Get the position of "Lambda1" in the pdeSys object.
  unsigned solLPdeIndex;
  solLPdeIndex = mlPdeSys->GetSolPdeIndex ("Lambda1");

  // Local Lambda1 solution.
  std::vector < adept::adouble > solL;

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
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Numer of solution element dofs.
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solType);
    unsigned nLDofs  = msh->GetElementDofNumber (iel, solLType);

    // Resize local arrays.
    for (unsigned K = 0; K < DIM; K++) {

      xhat[K].resize (nxDofs);

      solDx[K].resize (nxDofs);
      solx[K].resize (nxDofs);

      solNDx[K].resize (nxDofs);
      solNx[K].resize (nxDofs);

      solL.resize (nLDofs);
    }

    // Resize local arrays
    SYSDOF.resize (DIM * nxDofs + nLDofs);
    Res.resize (DIM * nxDofs + nLDofs);

    for (unsigned K = 0; K < DIM; K++) {
      aResNDx[K].assign (nxDofs, 0.);
    }
    aResL.assign (nLDofs, 0.);

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof (i, iel, solType);
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);

      for (unsigned K = 0; K < DIM; K++) {
        xhat[K][i] = (*msh->_topology->_Sol[K]) (iXDof);
        solDx[K][i] = (*sol->_Sol[solDxIndex[K]]) (iDDof);
        solx[K][i] = xhat[K][i] + solDx[K][i];
        solNDx[K][i] = (*sol->_Sol[solNDxIndex[K]]) (iDDof);

        // Global-to-global mapping between NDx solution node and pdeSys dof.
        SYSDOF[ K * nxDofs + i] =
          pdeSys->GetSystemDof (solNDxIndex[K], solNDxPdeIndex[K], i, iel);
      }
    }

    // Local storage of global mapping and solution.
    for (unsigned i = 0; i < nLDofs; i++) {

      // Global-to-local mapping between Lambda solution node and solution dof.
      unsigned iLDof = msh->GetSolutionDof (i, iel, solLType);
      solL[i] = (*sol->_Sol[solLIndex]) (iLDof);

      // Global-to-global mapping between Lambda solution node and pdeSys dof.
      SYSDOF[DIM * nxDofs + i] =
        pdeSys->GetSystemDof (solLIndex, solLPdeIndex, i, iel);
    }

    // start a new recording of all the operations involving adept variables.
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phiL;  // local test function
      const double *phix_uv[dim]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // Get Gauss point weight, test function, and first order derivatives.
      if (ielGeom == QUAD) {
        phix = msh->_finiteElement[ielGeom][solType]->GetPhi (ig);
        phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi (ig);

        phix_uv[0] = msh->_finiteElement[ielGeom][solType]->GetDPhiDXi (ig);
        phix_uv[1] = msh->_finiteElement[ielGeom][solType]->GetDPhiDEta (ig);

        weight = msh->_finiteElement[ielGeom][solType]->GetGaussWeight (ig);
      }

      // Special adjustments for triangles.
      else {
        msh->_finiteElement[ielGeom][solType]->Jacobian (xT, ig, weight, stdVectorPhi, stdVectorPhi_uv);

        phix = &stdVectorPhi[0];
///////// WHAT ABOUT PHIL? ////////////////////

        phi_uv0.resize (nxDofs);
        phi_uv1.resize (nxDofs);


        for (unsigned i = 0; i < nxDofs; i++) {
          phi_uv0[i] = stdVectorPhi_uv[i * dim];
          phi_uv1[i] = stdVectorPhi_uv[i * dim + 1];
        }

        phix_uv[0] = &phi_uv0[0];
        phix_uv[1] = &phi_uv1[0];

        phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi (ig);

      }

      // Initialize and compute values of x, Dx, NDx, x_uv at the Gauss points.
      double solDxg[3] = {0., 0., 0.};
      adept::adouble solNDxg[3] = {0., 0., 0.};

      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solNx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          solDxg[K] += phix[i] * solDx[K][i];
          solNDxg[K] += phix[i] * solNDx[K][i];
        }
        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * solx[K][i];
            solNx_uv[K][j]   += phix_uv[j][i] * (xhat[K][i] + solNDx[K][i]);
          }
        }
      }

      ///////// ADDED THIS /////////
      adept::adouble solLg = 0.;
      for (unsigned i = 0; i < nLDofs; i++) {
        solLg += phiL[i] * solL[i];
      }

      // Compute the metric, metric determinant, and area element.
      double g[dim][dim] = {{0., 0.}, {0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      double Area = weight * sqrt (detg);
      double Area2 = weight; // Trick to give equal weight to each element.

      // Compute the metric inverse.
      double gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;

      // Compute components of the unit normal N.
      double normal[DIM];
      normal[0] = (solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) / sqrt (detg);
      normal[1] = (solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);
      normal[2] = (solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);


      // Discretize the equation \delta CD = 0 on the basis d/du, d/dv.
      adept::adouble V[DIM];
      V[0] = solNx_uv[0][1] - normal[1] * solNx_uv[2][0] + normal[2] * solNx_uv[1][0];
      V[1] = solNx_uv[1][1] - normal[2] * solNx_uv[0][0] + normal[0] * solNx_uv[2][0];
      V[2] = solNx_uv[2][1] - normal[0] * solNx_uv[1][0] + normal[1] * solNx_uv[0][0];

      adept::adouble W[DIM];
      W[0] = solNx_uv[0][0] + normal[1] * solNx_uv[2][1] - normal[2] * solNx_uv[1][1];
      W[1] = solNx_uv[1][0] + normal[2] * solNx_uv[0][1] - normal[0] * solNx_uv[2][1];
      W[2] = solNx_uv[2][0] + normal[0] * solNx_uv[1][1] - normal[1] * solNx_uv[0][1];

      // adept::adouble M[DIM][dim];
      // M[0][0] = W[0] - normal[2] * V[1] + normal[1] * V[2];
      // M[1][0] = W[1] - normal[0] * V[2] + normal[2] * V[0];
      // M[2][0] = W[2] - normal[1] * V[0] + normal[0] * V[1];
      //
      // M[0][1] = V[0] + normal[2] * W[1] - normal[1] * W[2];
      // M[1][1] = V[1] + normal[0] * W[2] - normal[2] * W[0];
      // M[2][1] = V[2] + normal[1] * W[0] - normal[0] * W[1];

      adept::adouble Q[DIM][dim];
      Q[0][0] = (gi[1][1] * W[0]
                 + gi[0][0] * (normal[1] * V[2] - normal[2] * V[1])
                 + gi[0][1] * (normal[2] * W[1] - normal[1] * W[2] - V[0]));

      Q[1][0] = (gi[1][1] * W[1]
                 + gi[0][0] * (normal[2] * V[0] - normal[0] * V[2])
                 + gi[0][1] * (normal[0] * W[2] - normal[2] * W[0] - V[1]));

      Q[2][0] = (gi[1][1] * W[2]
                 + gi[0][0] * (normal[0] * V[1] - normal[1] * V[0])
                 + gi[0][1] * (normal[1] * W[0] - normal[0] * W[1] - V[2]));

      Q[0][1] = (gi[0][0] * V[0]
                 + gi[1][1] * (normal[2] * W[1] - normal[1] * W[2])
                 + gi[0][1] * (normal[1] * V[2] - normal[2] * V[1] - W[0]));

      Q[1][1] = (gi[0][0] * V[1]
                 + gi[1][1] * (normal[0] * W[2] - normal[2] * W[0])
                 + gi[0][1] * (normal[2] * V[0] - normal[0] * V[2] - W[1]));

      Q[2][1] = (gi[0][0] * V[2]
                 + gi[1][1] * (normal[1] * W[0] - normal[0] * W[1])
                 + gi[0][1] * (normal[0] * V[1] - normal[1] * V[0] - W[2]));


      // Compute new X minus old X dot N, for "reparametrization".
      adept::adouble DnXmDxdotN = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        DnXmDxdotN += (solDxg[K] - solNDxg[K]) * normal[K];
      }


      // Implement the Conformal Minimization equations.
      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;

          for (unsigned j = 0; j < dim; j++) {
            //term1 +=  M[K][j] * phix_uv[j][i];
            term1 += Q[K][j] * phix_uv[j][i];
          }

          // Conformal energy equation (with trick).
          aResNDx[K][i] += term1 * Area
                           // + timederiv * (solNDxg[K] - solDxg[K]) * phix[i] * Area2
                           + solL[0] * phix[i] * normal[K] * Area; //2 occasionally better???
        }
      }

      // Lagrange multiplier equation (with trick).
      for (unsigned i = 0; i < nLDofs; i++) {
        aResL[i] += phiL[i] * (DnXmDxdotN + eps * solL[i]) * Area; // no2
      }
      //aResL[0] += (DnXmDxdotN + eps * solL[0]) * Area;

    } // end GAUSS POINT LOOP

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    //copy the value of the adept::adoube aRes in double Res and store

    for (int K = 0; K < DIM; K++) {
      for (int i = 0; i < nxDofs; i++) {
        Res[ K * nxDofs + i] = -aResNDx[K][i].value();
      }
    }

    for (int i = 0; i < nLDofs; i++) {
      Res[DIM * nxDofs + i] = - aResL[i].value();
    }
    RES->add_vector_blocked (Res, SYSDOF);

    // Resize Jacobian.
    Jac.resize ( (DIM * nxDofs + nLDofs) * (DIM * nxDofs + nLDofs));

    // Define the dependent variables.
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResNDx[K][0], nxDofs);
    }
    s.dependent (&aResL[0], nLDofs);


    // Define the independent variables.
    for (int K = 0; K < DIM; K++) {
      s.independent (&solNDx[K][0], nxDofs);
    }
    s.independent (&solL[0], nLDofs);

    // Get the jacobian matrix (ordered by row).
    s.jacobian (&Jac[0], true);

    KK->add_matrix_blocked (Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } //end ELEMENT LOOP for each process.

  RES->close();
  KK->close();

} // end AssembleConformalMinimization.


void AssembleShearMinimization (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem> ("nProj");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh *msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem *el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution *mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object

  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix *KK = pdeSys->_KK;  // pointer to the global stiffness matrix object in pdeSys (level)
  NumericVector *RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = 2;
  const unsigned  DIM = 3;
  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol->GetIndex ("Dx1"); // get the position of "DX" in the ml_sol object
  solDxIndex[1] = mlSol->GetIndex ("Dx2"); // get the position of "DY" in the ml_sol object
  solDxIndex[2] = mlSol->GetIndex ("Dx3"); // get the position of "DZ" in the ml_sol object
  unsigned solType;
  solType = mlSol->GetSolutionType (solDxIndex[0]);  // get the finite element type for "U"
  std::vector < double > solx[DIM];  // surface coordinates
  std::vector < double > solDx[DIM];  // surface coordinates
  std::vector < double > solOldDx[DIM];  // surface coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  unsigned solNDxIndex[DIM];
  solNDxIndex[0] = mlSol->GetIndex ("nDx1");   // get the position of "Y1" in the ml_sol object
  solNDxIndex[1] = mlSol->GetIndex ("nDx2");   // get the position of "Y2" in the ml_sol object
  solNDxIndex[2] = mlSol->GetIndex ("nDx3");   // get the position of "Y3" in the ml_sol object
  unsigned solNDxPdeIndex[DIM];
  solNDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("nDx1");   // get the position of "Y1" in the pdeSys object
  solNDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("nDx2");   // get the position of "Y2" in the pdeSys object
  solNDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("nDx3");   // get the position of "Y3" in the pdeSys object
  std::vector < adept::adouble > solNDx[DIM]; // local Y solution


  unsigned solLIndex;
  solLIndex = mlSol->GetIndex ("Lambda1");   // get the position of "lambda" in the ml_sol object
  unsigned solLType;
  solLType = mlSol->GetSolutionType (solLIndex);  // get the finite element type for "lambda"
  unsigned solLPdeIndex;
  solLPdeIndex = mlPdeSys->GetSolPdeIndex ("Lambda1");   // get the position of "lambda" in the pdeSys object
  std::vector < adept::adouble > solL; // local lambda solution

  std::vector< int > SYSDOF; // local to global pdeSys dofs

  vector< double > Res; // local redidual vector
  std::vector< adept::adouble > aResNDx[3]; // local redidual vector
  std::vector< adept::adouble > aResL; // local redidual vector


  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)

  KK->zero();  // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solType);   // number of solution element dofs
    unsigned nLDofs  = msh->GetElementDofNumber (iel, solLType);   // number of solution element dofs


    for (unsigned K = 0; K < DIM; K++) {
      solDx[K].resize (nxDofs);
      solOldDx[K].resize (nxDofs);
      solx[K].resize (nxDofs);
      solNDx[K].resize (nxDofs);
      solL.resize (nLDofs);
    }

    // resize local arrays
    SYSDOF.resize (DIM * nxDofs + nLDofs);
    Res.resize (DIM * nxDofs + nLDofs);       //resize

    for (unsigned K = 0; K < DIM; K++) {
      aResNDx[K].assign (nxDofs, 0.);  //resize and zet to zero
    }
    aResL.assign (nLDofs, 0.);


    // local storage of global mapping and solution
    for (unsigned i = 0; i < nxDofs; i++) {
      unsigned iDDof = msh->GetSolutionDof (i, iel, solType); // global to local mapping between solution node and solution dof
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);
      for (unsigned K = 0; K < DIM; K++) {
        solx[K][i] = (*msh->_topology->_Sol[K]) (iXDof) + (*sol->_Sol[solDxIndex[K]]) (iDDof);
        solDx[K][i] = (*sol->_Sol[solDxIndex[K]]) (iDDof);
        solOldDx[K][i] = (*sol->_SolOld[solDxIndex[K]]) (iDDof);
        solNDx[K][i] = (*sol->_Sol[solNDxIndex[K]]) (iDDof);
        SYSDOF[ K * nxDofs + i] = pdeSys->GetSystemDof (solNDxIndex[K], solNDxPdeIndex[K], i, iel); // global to global mapping between solution node and pdeSys dof
      }
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nLDofs; i++) {
      unsigned iLDof = msh->GetSolutionDof (i, iel, solLType); // global to local mapping between solution node and solution dof
      solL[i] = (*sol->_Sol[solLIndex]) (iLDof); // global to local solution
      SYSDOF[DIM * nxDofs + i] = pdeSys->GetSystemDof (solLIndex, solLPdeIndex, i, iel);  // global to global mapping between solution node and pdeSys dof
    }


    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phix_uv[dim]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // *** get gauss point weight, test function and test function partial derivatives ***
      phix = msh->_finiteElement[ielGeom][solType]->GetPhi (ig);
      phix_uv[0] = msh->_finiteElement[ielGeom][solType]->GetDPhiDXi (ig); //derivative in u
      phix_uv[1] = msh->_finiteElement[ielGeom][solType]->GetDPhiDEta (ig); //derivative in v

      weight = msh->_finiteElement[ielGeom][solType]->GetGaussWeight (ig);

      double solOldDxg[3] = {0., 0., 0.};
      double solDxg[3] = {0., 0., 0.};
      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      adept::adouble solNDx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solDxOld_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      adept::adouble solNDxg[3] = {0., 0., 0.};


      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          solOldDxg[K] += phix[i] * solOldDx[K][i];
          solDxg[K] += phix[i] * solDx[K][i];
          solNDxg[K] += phix[i] * solNDx[K][i];
        }
        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]     += phix_uv[j][i] * solx[K][i];
            solDxOld_uv[K][j] += phix_uv[j][i] * solOldDx[K][i];
            solNDx_uv[K][j]   += phix_uv[j][i] * solNDx[K][i];
          }
        }
      }


      double g[dim][dim] = {{0., 0.}, {0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];

      double normal[DIM];
      normal[0] = (solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) / sqrt (detg);
      normal[1] = (solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);
      normal[2] = (solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);

      adept::adouble DnXmDxdotN = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        DnXmDxdotN += (solDxg[K] - solNDxg[K]) * normal[K];
      }

      double Area = weight * sqrt (detg);

      double gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;


      double Jir[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            Jir[i][J] += gi[i][k] * solx_uv[J][k];
          }
        }
      }

      adept::adouble solNDx_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solDxOld_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            solNDx_Xtan[I][J] += solNDx_uv[I][k] * Jir[k][J];
            solDxOld_Xtan[I][J] += solDxOld_uv[I][k] * Jir[k][J];
          }
        }
      }

      std::vector < adept::adouble > phix_Xtan[DIM];

      for (unsigned J = 0; J < DIM; J++) {
        phix_Xtan[J].assign (nxDofs, 0.);
        for (unsigned inode  = 0; inode < nxDofs; inode++) {
          for (unsigned k = 0; k < dim; k++) {
            phix_Xtan[J][inode] += phix_uv[k][inode] * Jir[k][J];
          }
        }
      }

      for (unsigned K = 0; K < DIM; K++) {
        //Energy equation
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;
          for (unsigned J = 0; J < DIM; J++) {
            if (J != K) {
              term1 += (solNDx_Xtan[K][J] + solNDx_Xtan[J][K]) * phix_Xtan[J][i];
            }
            else {
              term1 += 0.5 * (solNDx_Xtan[K][J] + solNDx_Xtan[J][K]) * phix_Xtan[J][i];
            }
          }
          aResNDx[K][i] += (term1 + solL[0] * phix[i] * normal[K]) * Area;
        }
      }

      aResL[0] += DnXmDxdotN * Area;




    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store


    for (int K = 0; K < DIM; K++) {
      for (int i = 0; i < nxDofs; i++) {
        Res[ K * nxDofs + i] = -aResNDx[K][i].value();
      }
    }

    for (int i = 0; i < nLDofs; i++) {
      Res[DIM * nxDofs + i] = - aResL[i].value();
    }

    RES->add_vector_blocked (Res, SYSDOF);



    Jac.resize ( (DIM * nxDofs + nLDofs) * (DIM * nxDofs + nLDofs));

    // define the dependent variables

    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResNDx[K][0], nxDofs);
    }
    s.dependent (&aResL[0], nLDofs);


    // define the dependent variables

    for (int K = 0; K < DIM; K++) {
      s.independent (&solNDx[K][0], nxDofs);
    }
    s.independent (&solL[0], nLDofs);

    // get the jacobian matrix (ordered by row)
    s.jacobian (&Jac[0], true);

    KK->add_matrix_blocked (Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();

  // ***************** END ASSEMBLY *******************
}



void AssembleSphereConformalMinimization (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  Extract pointers to the several objects that we are going to use.
  NonLinearImplicitSystem* mlPdeSys =
    &ml_prob.get_system< NonLinearImplicitSystem> ("nProj");

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Pointers to the mesh (level) object and elem object in mesh (level).
  Mesh *msh = ml_prob._ml_msh->GetLevel (level);
  elem *el = msh->el;

  // Pointers to the multilevel solution, solution (level) and equation (level).
  MultiLevelSolution *mlSol = ml_prob._ml_sol;
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel (level);
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
  std::vector < std::vector < double > > xT (2);

  xT[0].resize (3);
  xT[0][0] = -0.5;
  xT[0][1] = 0.5;
  xT[0][2] = 0.;

  xT[1].resize (3);
  xT[1][0] = 0.;
  xT[1][1] = 0.;
  xT[1][2] = sqrt (3.) / 2.;

  std::vector<double> phi_uv0;
  std::vector<double> phi_uv1;

  std::vector< double > stdVectorPhi;
  std::vector< double > stdVectorPhi_uv;

  // Extract positions of Dx in ml_sol object.
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol->GetIndex ("Dx1");
  solDxIndex[1] = mlSol->GetIndex ("Dx2");
  solDxIndex[2] = mlSol->GetIndex ("Dx3");


  unsigned solYIndex[DIM];
  solYIndex[0] = mlSol->GetIndex ("Y1");
  solYIndex[1] = mlSol->GetIndex ("Y2");
  solYIndex[2] = mlSol->GetIndex ("Y3");

  // Extract finite element type for the solution.
  unsigned solType;
  solType = mlSol->GetSolutionType (solDxIndex[0]);

  // Local solution vectors for X, Dx, Xhat, XC.
  std::vector < double > solY[DIM];
  std::vector < double > solx[DIM];
  std::vector < double > solDx[DIM];
  std::vector < double > xhat[DIM];


  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get the poitions of Y in the ml_sol object.
  unsigned solNDxIndex[DIM];
  solNDxIndex[0] = mlSol->GetIndex ("nDx1");
  solNDxIndex[1] = mlSol->GetIndex ("nDx2");
  solNDxIndex[2] = mlSol->GetIndex ("nDx3");

  // Get the positions of Y in the pdeSys object.
  unsigned solNDxPdeIndex[DIM];
  solNDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("nDx1");
  solNDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("nDx2");
  solNDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("nDx3");

  // Local solution vectors for Nx and NDx.
  std::vector < adept::adouble > solNDx[DIM];
  std::vector < adept::adouble > solNx[DIM];

  // Get the position of "Lambda1" in the ml_sol object.
  unsigned solLIndex;
  solLIndex = mlSol->GetIndex ("Lambda1");

  // Get the finite element type for "Lambda1".
  unsigned solLType;
  solLType = mlSol->GetSolutionType (solLIndex);

  // Get the position of "Lambda1" in the pdeSys object.
  unsigned solLPdeIndex;
  solLPdeIndex = mlPdeSys->GetSolPdeIndex ("Lambda1");

  // Local Lambda1 solution.
  std::vector < adept::adouble > solL;

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
  for (int iel = msh->_elementOffset[iproc];
       iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Numer of solution element dofs.
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solType);
    unsigned nLDofs  = msh->GetElementDofNumber (iel, solLType);

    // Resize local arrays.
    for (unsigned K = 0; K < DIM; K++) {

      xhat[K].resize (nxDofs);

      solDx[K].resize (nxDofs);
      solx[K].resize (nxDofs);

      solY[K].resize (nxDofs);

      solNDx[K].resize (nxDofs);
      solNx[K].resize (nxDofs);

      solL.resize (nLDofs);
    }

    SYSDOF.resize (DIM * nxDofs + nLDofs);
    Res.resize (DIM * nxDofs + nLDofs);

    for (unsigned K = 0; K < DIM; K++) {
      aResNDx[K].assign (nxDofs, 0.);
    }
    aResL.assign (nLDofs, 0.);


    // local storage of global mapping and solution
    for (unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof (i, iel, solType);
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);

      for (unsigned K = 0; K < DIM; K++) {
        xhat[K][i] = (*msh->_topology->_Sol[K]) (iXDof);
        solDx[K][i] = (*sol->_Sol[solDxIndex[K]]) (iDDof);

        solY[K][i] = (*sol->_Sol[solYIndex[K]]) (iDDof);

        solx[K][i] = xhat[K][i] + solDx[K][i];
        solNDx[K][i] = (*sol->_Sol[solNDxIndex[K]]) (iDDof);

        // Global-to-global mapping between NDx solution node and pdeSys dof.
        SYSDOF[ K * nxDofs + i] =
          pdeSys->GetSystemDof (solNDxIndex[K], solNDxPdeIndex[K], i, iel);
      }
    }

    // Local storage of global mapping and solution.
    for (unsigned i = 0; i < nLDofs; i++) {

      // Global-to-local mapping between Lambda solution node and solution dof.
      unsigned iLDof = msh->GetSolutionDof (i, iel, solLType);
      solL[i] = (*sol->_Sol[solLIndex]) (iLDof);

      // Global-to-global mapping between Lambda solution node and pdeSys dof.
      SYSDOF[DIM * nxDofs + i] =
        pdeSys->GetSystemDof (solLIndex, solLPdeIndex, i, iel);
    }

    // start a new recording of all the operations involving adept variables.
    s.new_recording();

    // begin GAUSS POINT LOOP
    for (unsigned ig = 0; ig < msh->
         _finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phiL;  // local test function
      const double *phix_uv[dim]; // first order derivatives in (u,v)

      double weight; // gauss point weight

      // Get Gauss point weight, test function, and first order derivatives.
      if (ielGeom == QUAD) {
        phix = msh->_finiteElement[ielGeom][solType]->GetPhi (ig);
        phix_uv[0] = msh->_finiteElement[ielGeom][solType]->GetDPhiDXi (ig);
        phix_uv[1] = msh->_finiteElement[ielGeom][solType]->GetDPhiDEta (ig);
        weight = msh->_finiteElement[ielGeom][solType]->GetGaussWeight (ig);
        phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi (ig);
      }

      // Special adjustments for triangles.
      else {
        msh->_finiteElement[ielGeom][solType]->Jacobian (xT, ig,
                                                         weight, stdVectorPhi, stdVectorPhi_uv);
        phix = &stdVectorPhi[0];

        phi_uv0.resize (nxDofs);
        phi_uv1.resize (nxDofs);

        for (unsigned i = 0; i < nxDofs; i++) {
          phi_uv0[i] = stdVectorPhi_uv[i * dim];
          phi_uv1[i] = stdVectorPhi_uv[i * dim + 1];
        }
        phix_uv[0] = &phi_uv0[0];
        phix_uv[1] = &phi_uv1[0];

        phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi (ig);

      }

      // Initialize and compute values of x, Dx, NDx, x_uv at the Gauss points.
      double solYg[3] = {0., 0., 0.};
      double solxg[3] = {0., 0., 0.};

      double solDxg[3] = {0., 0., 0.};
      adept::adouble solNDxg[3] = {0., 0., 0.};

      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble X_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};



      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          solxg[K] += phix[i] * solx[K][i];
          solYg[K] += phix[i] * solY[K][i];
          solDxg[K] += phix[i] * solDx[K][i];
          solNDxg[K] += phix[i] * solNDx[K][i];
        }

        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * solx[K][i];
            X_uv[K][j]   += phix_uv[j][i] * (xhat[K][i] + solNDx[K][i]);
          }
        }
      }

      adept::adouble solLg = 0.;
      for (unsigned i = 0; i < nLDofs; i++) {
        solLg += phiL[i] * solL[i];
      }

      // Compute the metric, metric determinant, and area element.
      double g[dim][dim] = {{0., 0.}, {0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      double Area = weight * sqrt (detg);
      double Area2 = weight; // Trick to give equal weight to each element.

      // Compute components of the unit normal N.
      double normal[DIM];
      normal[0] = (solx_uv[1][0] * solx_uv[2][1]
                   - solx_uv[2][0] * solx_uv[1][1]) / sqrt (detg);
      normal[1] = (solx_uv[2][0] * solx_uv[0][1]
                   - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);
      normal[2] = (solx_uv[0][0] * solx_uv[1][1]
                   - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);

      double a = 0;
      for (unsigned I = 0; I < DIM; I++) {
        a -= solYg[I] * normal[I];
      }
      double signa = (a > 0) ? 1 : -1;
      a = 2. / (a + 0 * signa * eps);

      //       std::cout << solxg[0] << " " << solxg[1] <<" "<< solxg[2]<<std::endl;
      //       std::cout << normal[0] << " " << normal[1] <<" "<< normal[2]<<std::endl;
      //       std::cout << a<<std::endl;

      // Discretize the equation \delta CD = 0 on the basis d/du, d/dv.
      adept::adouble V[DIM];
      V[0] = X_uv[0][1] - normal[1] * X_uv[2][0] + normal[2] * X_uv[1][0];
      V[1] = X_uv[1][1] - normal[2] * X_uv[0][0] + normal[0] * X_uv[2][0];
      V[2] = X_uv[2][1] - normal[0] * X_uv[1][0] + normal[1] * X_uv[0][0];

      adept::adouble W[DIM];
      W[0] = X_uv[0][0] + normal[1] * X_uv[2][1] - normal[2] * X_uv[1][1];
      W[1] = X_uv[1][0] + normal[2] * X_uv[0][1] - normal[0] * X_uv[2][1];
      W[2] = X_uv[2][0] + normal[0] * X_uv[1][1] - normal[1] * X_uv[0][1];

      adept::adouble M[DIM][dim];
      M[0][0] = W[0] - normal[2] * V[1] + normal[1] * V[2];
      M[1][0] = W[1] - normal[0] * V[2] + normal[2] * V[0];
      M[2][0] = W[2] - normal[1] * V[0] + normal[0] * V[1];

      M[0][1] = V[0] + normal[2] * W[1] - normal[1] * W[2];
      M[1][1] = V[1] + normal[0] * W[2] - normal[2] * W[0];
      M[2][1] = V[2] + normal[1] * W[0] - normal[0] * W[1];

      // Compute new X minus old X dot N, for "reparametrization".
      adept::adouble DnXmDxdotN = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        DnXmDxdotN += (solDxg[K] - solNDxg[K]) * normal[K];
      }

      // Compute the metric inverse.
      double gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;

      // Compute the "reduced Jacobian" g^{ij}X_j .
      double Jir[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            Jir[i][J] += gi[i][k] * solx_uv[J][k];
          }
        }
      }

      // Initializing tangential gradients of X and W (new, middle, old).
      double solx_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble NX_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble XminusXold_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      // adept::adouble solW_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

      // Computing tangential gradients defined above.
      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            NX_Xtan[I][J] += X_uv[I][k] * Jir[k][J];
            solx_Xtan[I][J] += solx_uv[I][k] * Jir[k][J];
            // solW_Xtan[I][J] += solW_uv[I][k] * Jir[k][J];
          }
        }
      }

      // Define and compute gradients of test functions for X and W.
      std::vector < adept::adouble > phiY_Xtan[DIM];
      std::vector < adept::adouble > phix_Xtan[DIM];

      for (unsigned J = 0; J < DIM; J++) {
        phix_Xtan[J].assign (nxDofs, 0.);
        phiY_Xtan[J].assign (nxDofs, 0.);

        for (unsigned inode  = 0; inode < nxDofs; inode++) {
          for (unsigned k = 0; k < dim; k++) {
            phix_Xtan[J][inode] += phix_uv[k][inode] * Jir[k][J];
          }
        }

        for (unsigned inode  = 0; inode < nxDofs; inode++) {
          for (unsigned k = 0; k < dim; k++) {
            phiY_Xtan[J][inode] += phix_uv[k][inode] * Jir[k][J];
          }
        }
      }

      // Implement the Conformal Minimization equations.
      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;

          for (unsigned j = 0; j < dim; j++) {
            term1 +=  M[K][j] * phix_uv[j][i];
            // term3 += solx_uv[K][j] * phix_uv[j][i];
          }


          adept::adouble term2 = 0.;
          //           for (unsigned J = 0; J < DIM; J++) {
          //             term2 += (a * normal[J] - solDxg[J] + solNDxg [J]) * phix[i];
          //           }

          term2 += (a * normal[K] - solDxg[K] + solNDxg [K]) * phix[i];

          // Conformal energy equation (with trick).
          aResNDx[K][i] += term1 * Area2
                           + solLg * term2 * Area;  //no2
        }

        // Lagrange multiplier equation (with trick).
        for (unsigned i = 0; i < nLDofs; i++) {
          adept::adouble term1 = 0.;
          adept::adouble term2 = 0.;
          for (unsigned J = 0; J < DIM; J++) {
            term1 += (a * normal[J] - solDxg[J] + solNDxg [J]) * (a * normal[J] - solDxg[J] + solNDxg [J]) ;
            term2 += (a * normal[J]) * (a * normal[J]) ;
          }
          aResL[i] += (phiL[i] * (term1 - term2) + eps * 0 * solLg * phiL[i]) * Area;    // no2
        }

        //         for(unsigned i = 0; i< nLDofs; i++){
        //           aResL[i] += phiL[i] * (DnXmDxdotN + eps * solL[i]) * Area2; // no2
        //         }

      }

    } // end GAUSS POINT LOOP

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    //copy the value of the adept::adoube aRes in double Res and store

    for (int K = 0; K < DIM; K++) {
      for (int i = 0; i < nxDofs; i++) {
        Res[ K * nxDofs + i] = -aResNDx[K][i].value(); // for X
      }
    }

    for (int i = 0; i < nLDofs; i++) {
      Res[DIM * nxDofs + i] = - aResL[i].value(); // for Lambda
    }

    RES->add_vector_blocked (Res, SYSDOF);

    // Resize Jacobian.
    Jac.resize ( (DIM * nxDofs + nLDofs) * (DIM * nxDofs + nLDofs));

    // Define the dependent variables.
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResNDx[K][0], nxDofs);
    }
    s.dependent (&aResL[0], nLDofs);

    // Define the independent variables.
    for (int K = 0; K < DIM; K++) {
      s.independent (&solNDx[K][0], nxDofs);
    }
    s.independent (&solL[0], nLDofs);

    // Get the jacobian matrix (ordered by row).
    s.jacobian (&Jac[0], true);

    KK->add_matrix_blocked (Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } //end ELEMENT LOOP for each process.

  RES->close();
  KK->close();

} // end AssembleSphereConformalMinimization.


// Building the Conformal Minimization system.
void AssembleO2ConformalMinimizationEA (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  Extract pointers to the several objects that we are going to use.
  NonLinearImplicitSystem* mlPdeSys =
    &ml_prob.get_system< NonLinearImplicitSystem> ("nProj");

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Pointers to the mesh (level) object and elem object in mesh (level).
  Mesh *msh = ml_prob._ml_msh->GetLevel (level);
  elem *el = msh->el;

  // Pointers to the multilevel solution, solution (level) and equation (level).
  MultiLevelSolution *mlSol = ml_prob._ml_sol;
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel (level);
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
  std::vector < std::vector < double > > xT (2);

  xT[0].resize (3);
  xT[0][0] = -0.5;
  xT[0][1] = 0.5;
  xT[0][2] = 0.;

  xT[1].resize (3);
  xT[1][0] = 0.;
  xT[1][1] = 0.;
  xT[1][2] = sqrt (3.) / 2.;

  std::vector<double> phi_uv0;
  std::vector<double> phi_uv1;

  std::vector< double > stdVectorPhi;
  std::vector< double > stdVectorPhi_uv;

  // Extract positions of Dx in ml_sol object.
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol->GetIndex ("Dx1");
  solDxIndex[1] = mlSol->GetIndex ("Dx2");
  solDxIndex[2] = mlSol->GetIndex ("Dx3");

  // Extract finite element type for the solution.
  unsigned solType;
  solType = mlSol->GetSolutionType (solDxIndex[0]);

  // Local solution vectors for X, Dx, Xhat, XC.
  std::vector < double > solx[DIM];
  std::vector < double > solDx[DIM];
  std::vector < double > xhat[DIM];
  std::vector < double > xc[DIM];

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get the poitions of Y in the ml_sol object.
  unsigned solNDxIndex[DIM];
  solNDxIndex[0] = mlSol->GetIndex ("nDx1");
  solNDxIndex[1] = mlSol->GetIndex ("nDx2");
  solNDxIndex[2] = mlSol->GetIndex ("nDx3");

  // Get the positions of Y in the pdeSys object.
  unsigned solNDxPdeIndex[DIM];
  solNDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("nDx1");
  solNDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("nDx2");
  solNDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("nDx3");

  // Local solution vectors for Nx and NDx.
  std::vector < adept::adouble > solNDx[DIM];
  std::vector < adept::adouble > solNx[DIM];

  // Get the position of "Lambda1" in the ml_sol object.
  unsigned solLIndex;
  solLIndex = mlSol->GetIndex ("Lambda1");

  // Get the finite element type for "Lambda1".
  unsigned solLType;
  solLType = mlSol->GetSolutionType (solLIndex);

  // Get the position of "Lambda1" in the pdeSys object.
  unsigned solLPdeIndex;
  solLPdeIndex = mlPdeSys->GetSolPdeIndex ("Lambda1");

  // Local Lambda1 solution.
  std::vector < adept::adouble > solL;

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
  for (int iel = msh->_elementOffset[iproc];
       iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Numer of solution element dofs.
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solType);
    unsigned nLDofs  = msh->GetElementDofNumber (iel, solLType);

    // Resize local arrays.
    for (unsigned K = 0; K < DIM; K++) {

      xhat[K].resize (nxDofs);
      xc[K].assign (nxDofs, 0.);

      solDx[K].resize (nxDofs);
      solx[K].resize (nxDofs);

      solNDx[K].resize (nxDofs);
      solNx[K].resize (nxDofs);

      solL.resize (nLDofs);
    }

    SYSDOF.resize (DIM * nxDofs + nLDofs);
    Res.resize (DIM * nxDofs + nLDofs);

    for (unsigned K = 0; K < DIM; K++) {
      aResNDx[K].assign (nxDofs, 0.);
    }
    aResL.assign (nLDofs, 0.);


    // local storage of global mapping and solution
    for (unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof (i, iel, solType);
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);

      for (unsigned K = 0; K < DIM; K++) {
        xhat[K][i] = (*msh->_topology->_Sol[K]) (iXDof);
        solDx[K][i] = (*sol->_Sol[solDxIndex[K]]) (iDDof);
        solx[K][i] = xhat[K][i] + solDx[K][i];
        solNDx[K][i] = (*sol->_Sol[solNDxIndex[K]]) (iDDof);

        // Global-to-global mapping between NDx solution node and pdeSys dof.
        SYSDOF[ K * nxDofs + i] =
          pdeSys->GetSystemDof (solNDxIndex[K], solNDxPdeIndex[K], i, iel);
      }
    }

    // Local storage of global mapping and solution.
    for (unsigned i = 0; i < nLDofs; i++) {

      // Global-to-local mapping between Lambda solution node and solution dof.
      unsigned iLDof = msh->GetSolutionDof (i, iel, solLType);
      solL[i] = (*sol->_Sol[solLIndex]) (iLDof);

      // Global-to-global mapping between Lambda solution node and pdeSys dof.
      SYSDOF[DIM * nxDofs + i] =
        pdeSys->GetSystemDof (solLIndex, solLPdeIndex, i, iel);
    }

    // start a new recording of all the operations involving adept variables.
    s.new_recording();

    // begin GAUSS POINT LOOP
    for (unsigned ig = 0; ig < msh->
         _finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phiL;  // local test function
      const double *phix_uv[dim]; // first order derivatives in (u,v)

      double weight; // gauss point weight

      // Get Gauss point weight, test function, and first order derivatives.
      if (ielGeom == QUAD) {
        phix = msh->_finiteElement[ielGeom][solType]->GetPhi (ig);
        phix_uv[0] = msh->_finiteElement[ielGeom][solType]->GetDPhiDXi (ig);
        phix_uv[1] = msh->_finiteElement[ielGeom][solType]->GetDPhiDEta (ig);
        weight = msh->_finiteElement[ielGeom][solType]->GetGaussWeight (ig);
        phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi (ig);
      }

      // Special adjustments for triangles.
      else {
        msh->_finiteElement[ielGeom][solType]->Jacobian (xT, ig,
                                                         weight, stdVectorPhi, stdVectorPhi_uv);
        phix = &stdVectorPhi[0];

        phi_uv0.resize (nxDofs);
        phi_uv1.resize (nxDofs);

        for (unsigned i = 0; i < nxDofs; i++) {
          phi_uv0[i] = stdVectorPhi_uv[i * dim];
          phi_uv1[i] = stdVectorPhi_uv[i * dim + 1];
        }
        phix_uv[0] = &phi_uv0[0];
        phix_uv[1] = &phi_uv1[0];

        phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi (ig);

      }

      // Initialize and compute values of x, Dx, NDx, x_uv at the Gauss points.
      double solDxg[3] = {0., 0., 0.};
      adept::adouble solNDxg[3] = {0., 0., 0.};

      adept::adouble solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solNx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};



      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          solDxg[K] += phix[i] * solDx[K][i];
          solNDxg[K] += phix[i] * solNDx[K][i];
        }

        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * (xhat[K][i] + 0.5 * (solDx[K][i] + solNDx[K][i]));
            solNx_uv[K][j]   += phix_uv[j][i] * (xhat[K][i] + (solNDx[K][i]));
          }
        }
      }

      adept::adouble solLg = 0.;
      for (unsigned i = 0; i < nLDofs; i++) {
        solLg += phiL[i] * solL[i];
      }

      // Compute the metric, metric determinant, and area element.
      adept::adouble g[dim][dim] = {{0., 0.}, {0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      adept::adouble detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      adept::adouble Area = weight * sqrt (detg);
      double Area2 = weight; // Trick to give equal weight to each element.

      // Computing the metric inverse
      adept::adouble gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;

      // Computing the "reduced Jacobian" g^{ij}X_j .
      adept::adouble Jir[dim][DIM] = {{0., 0., 0.}, {0., 0., 0.}};
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            Jir[i][J] += gi[i][k] * solx_uv[J][k];
          }
        }
      }

      // Initializing tangential gradients of X and W (new, middle, old).
      adept::adouble solNx_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble solx_Xtan[DIM][DIM] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

      // Computing tangential gradients defined above.
      for (unsigned I = 0; I < DIM; I++) {
        for (unsigned J = 0; J < DIM; J++) {
          for (unsigned k = 0; k < dim; k++) {
            solNx_Xtan[I][J] += solNx_uv[I][k] * Jir[k][J];
            solx_Xtan[I][J] += solx_uv[I][k] * Jir[k][J];
          }
        }
      }

      // Define and compute gradients of test functions for X.
      std::vector < adept::adouble > phix_Xtan[DIM];

      for (unsigned J = 0; J < DIM; J++) {
        phix_Xtan[J].assign (nxDofs, 0.);

        for (unsigned inode  = 0; inode < nxDofs; inode++) {
          for (unsigned k = 0; k < dim; k++) {
            phix_Xtan[J][inode] += phix_uv[k][inode] * Jir[k][J];
          }
        }
      }

      // Compute components of the unit normal N.
      adept::adouble normal[DIM];
      normal[0] = (solx_uv[1][0] * solx_uv[2][1]
                   - solx_uv[2][0] * solx_uv[1][1]) / sqrt (detg);
      normal[1] = (solx_uv[2][0] * solx_uv[0][1]
                   - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);
      normal[2] = (solx_uv[0][0] * solx_uv[1][1]
                   - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);

      // Discretize the equation \delta CD = 0 on the basis d/du, d/dv.
      adept::adouble V[DIM];
      V[0] = solNx_uv[0][1] - normal[1] * solNx_uv[2][0] + normal[2] * solNx_uv[1][0];
      V[1] = solNx_uv[1][1] - normal[2] * solNx_uv[0][0] + normal[0] * solNx_uv[2][0];
      V[2] = solNx_uv[2][1] - normal[0] * solNx_uv[1][0] + normal[1] * solNx_uv[0][0];

      adept::adouble W[DIM];
      W[0] = solNx_uv[0][0] + normal[1] * solNx_uv[2][1] - normal[2] * solNx_uv[1][1];
      W[1] = solNx_uv[1][0] + normal[2] * solNx_uv[0][1] - normal[0] * solNx_uv[2][1];
      W[2] = solNx_uv[2][0] + normal[0] * solNx_uv[1][1] - normal[1] * solNx_uv[0][1];

      adept::adouble M[DIM][dim];
      M[0][0] = W[0] - normal[2] * V[1] + normal[1] * V[2] + /*(1. / sqrt (detg)) * */ (
                  (solNx_uv[1][0] * solNx_uv[1][1] + solNx_uv[2][0] * solNx_uv[2][1]) * V[0]
                  - (solNx_uv[1][1] * solNx_uv[1][1] + solNx_uv[2][1] * solNx_uv[2][1]) * W[0]
                  + solNx_uv[0][1] * (solNx_uv[1][1] * W[1] + solNx_uv[2][1] * W[2])
                  - solNx_uv[0][0] * (solNx_uv[1][1] * V[1] + solNx_uv[2][1] * V[2])
                );
      M[1][0] = W[1] - normal[0] * V[2] + normal[2] * V[0] + /*(1. / sqrt (detg)) * */ (
                  (solNx_uv[2][0] * solNx_uv[2][1] + solNx_uv[0][0] * solNx_uv[0][1]) * V[1]
                  - (solNx_uv[2][1] * solNx_uv[2][1] + solNx_uv[0][1] * solNx_uv[0][1]) * W[1]
                  + solNx_uv[1][1] * (solNx_uv[2][1] * W[2] + solNx_uv[0][1] * W[0])
                  - solNx_uv[1][0] * (solNx_uv[2][1] * V[2] + solNx_uv[0][1] * V[0])
                );
      M[2][0] = W[2] - normal[1] * V[0] + normal[0] * V[1] + /*(1. / sqrt (detg)) **/ (
                  (solNx_uv[0][0] * solNx_uv[0][1] + solNx_uv[1][0] * solNx_uv[1][1]) * V[2]
                  - (solNx_uv[0][1] * solNx_uv[0][1] + solNx_uv[1][1] * solNx_uv[1][1]) * W[2]
                  + solNx_uv[2][1] * (solNx_uv[0][1] * W[0] + solNx_uv[1][1] * W[1])
                  - solNx_uv[2][0] * (solNx_uv[0][1] * V[0] + solNx_uv[1][1] * V[1])
                );
      M[0][1] = V[0] + normal[2] * W[1] - normal[1] * W[2] + /*(1. / sqrt (detg)) **/ (
                  (solNx_uv[1][1] * solNx_uv[1][0] + solNx_uv[2][1] * solNx_uv[2][0]) * W[0]
                  - (solNx_uv[1][0] * solNx_uv[1][0] + solNx_uv[2][0] * solNx_uv[2][0]) * V[0]
                  + solNx_uv[0][0] * (solNx_uv[1][0] * V[1] + solNx_uv[2][0] * V[2])
                  - solNx_uv[0][1] * (solNx_uv[1][0] * W[1] + solNx_uv[2][0] * W[2])
                );
      M[1][1] = V[1] + normal[0] * W[2] - normal[2] * W[0] + /*(1. / sqrt (detg)) **/ (
                  (solNx_uv[2][1] * solNx_uv[2][0] + solNx_uv[0][1] * solNx_uv[0][0]) * W[1]
                  - (solNx_uv[2][0] * solNx_uv[2][0] + solNx_uv[0][0] * solNx_uv[0][0]) * V[1]
                  + solNx_uv[1][0] * (solNx_uv[2][0] * V[2] + solNx_uv[0][0] * V[0])
                  - solNx_uv[1][1] * (solNx_uv[2][0] * W[2] + solNx_uv[0][0] * W[0])
                );
      M[2][1] = V[2] + normal[1] * W[0] - normal[0] * W[1] + /*(1. / sqrt (detg)) **/ (
                  (solNx_uv[0][1] * solNx_uv[0][0] + solNx_uv[1][1] * solNx_uv[1][0]) * W[2]
                  - (solNx_uv[0][0] * solNx_uv[0][0] + solNx_uv[1][0] * solNx_uv[1][0]) * V[2]
                  + solNx_uv[2][0] * (solNx_uv[0][0] * V[0] + solNx_uv[1][0] * V[1])
                  - solNx_uv[2][1] * (solNx_uv[0][0] * W[0] + solNx_uv[1][0] * W[1])
                );

      adept::adouble G = 0;
      G = normal[0] * (solNx_uv[2][1] * W[1] - solNx_uv[2][0] * V[1] + solNx_uv[1][0] * V[2] - solNx_uv[1][1] * W[2])
          + normal[1] * (solNx_uv[0][1] * W[2] - solNx_uv[0][0] * V[2] + solNx_uv[2][0] * V[0] - solNx_uv[2][1] * W[0])
          + normal[2] * (solNx_uv[1][1] * W[0] - solNx_uv[1][0] * V[0] + solNx_uv[0][0] * V[1] - solNx_uv[0][1] * W[1]);

      // Compute new X minus old X dot N, for "reparametrization".
      adept::adouble DnXmDxdotN = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        DnXmDxdotN += (solDxg[K] - solNDxg[K]) * normal[K];
      }

      // // Implement the curvature equation Y = \Delta X .
      // for (unsigned K = 0; K < DIM; K++) {
      //   for (unsigned i = 0; i < nxDofs; i++) {
      //     adept::adouble term2 = 0.;
      //
      //     for (unsigned J = 0; J < DIM; J++) {
      //       term2 +=  solNx_Xtan[K][J] * phix_Xtan[J][i];
      //     }
      //     aResx[K][i] += (solYg[K] * phix[i] + term1) * Area;
      //   }

      // Implement the Conformal Minimization equations.
      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;
          adept::adouble term1new = 0.;
          adept::adouble gxgp = 0.;

          for (unsigned j = 0; j < dim; j++) {
            for (unsigned l = 0; l < dim; l++) {
              term1new = gi[l][j] * M[K][l] * phix_uv[j][i];
            }
            term1 +=  M[K][j] * phix_uv[j][i];
          }

          for (unsigned J = 0; J < DIM; J++) {
            gxgp +=  solx_Xtan[K][J] * phix_Xtan[J][i];
          }

          // Conformal energy equation (with trick).
          aResNDx[K][i] += term1 * Area2
                           + G * gxgp * Area
                           + solLg * phix[i] * normal[K] * Area;  //no2
        }
      }

      // Lagrange multiplier equation (with trick).
      for (unsigned i = 0; i < nLDofs; i++) {
        aResL[i] += phiL[i] * (DnXmDxdotN + eps * solL[i]) * Area; // no2
      }

    } // end GAUSS POINT LOOP

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    //copy the value of the adept::adoube aRes in double Res and store

    for (int K = 0; K < DIM; K++) {
      for (int i = 0; i < nxDofs; i++) {
        Res[ K * nxDofs + i] = -aResNDx[K][i].value(); // for X
      }
    }

    for (int i = 0; i < nLDofs; i++) {
      Res[DIM * nxDofs + i] = - aResL[i].value(); // for Lambda
    }

    RES->add_vector_blocked (Res, SYSDOF);

    // Resize Jacobian.
    Jac.resize ( (DIM * nxDofs + nLDofs) * (DIM * nxDofs + nLDofs));

    // Define the dependent variables.
    for (int K = 0; K < DIM; K++) {
      s.dependent (&aResNDx[K][0], nxDofs);
    }
    s.dependent (&aResL[0], nLDofs);

    // Define the independent variables.
    for (int K = 0; K < DIM; K++) {
      s.independent (&solNDx[K][0], nxDofs);
    }
    s.independent (&solL[0], nLDofs);

    // Get the jacobian matrix (ordered by row).
    s.jacobian (&Jac[0], true);

    KK->add_matrix_blocked (Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } //end ELEMENT LOOP for each process.

  RES->close();
  KK->close();

} // end AssembleO2ConformalMinimization.
