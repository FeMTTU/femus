 void GetElementNearVertexNumber(MultiLevelSolution &mlSol){

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* sol = mlSol.GetSolutionLevel (level);
  Mesh* msh = mlSol._mlMesh->GetLevel (level);

  unsigned DIM = 3u;
  unsigned solIndex = mlSol.GetIndex ("ENVN");
  
  unsigned solType = mlSol.GetSolutionType (solIndex);

  sol->_Sol[solIndex]->zero();
  
  unsigned iproc = msh->processor_id();
  
  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Number of solution element dofs.
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofs  = msh->GetElementDofNumber (iel, solType);
    
    for (unsigned i = 0; i < nDofs; i++) {

      unsigned iDof = msh->GetSolutionDof (i, iel, solType);

      sol->_Sol[solIndex]->add(iDof, 1);
      
    }
  }
   
   sol->_Sol[solIndex]->close();
  
}

// Eugenio's standard FEMuS function.
void CopyDisplacement (MultiLevelSolution &mlSol,  const bool &forward) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* solution  = mlSol.GetSolutionLevel (level);
  Mesh* msh = mlSol._mlMesh->GetLevel (level);

  unsigned DIM = 3u;
  vector < unsigned > solDxIndex (DIM);
  solDxIndex[0] = mlSol.GetIndex ("Dx1");
  solDxIndex[1] = mlSol.GetIndex ("Dx2");
  solDxIndex[2] = mlSol.GetIndex ("Dx3");

  vector < unsigned > solNDxIndex (DIM);
  solNDxIndex[0] = mlSol.GetIndex ("nDx1");
  solNDxIndex[1] = mlSol.GetIndex ("nDx2");
  solNDxIndex[2] = mlSol.GetIndex ("nDx3");

  if (forward) {
    for (unsigned i = 0; i < DIM; i++) {
      * (solution->_Sol[solNDxIndex[i]]) = * (solution->_Sol[solDxIndex[i]]);
    }
  }
  else {
    for (unsigned i = 0; i < DIM; i++) {
      * (solution->_Sol[solDxIndex[i]]) = * (solution->_Sol[solNDxIndex[i]]);
    }
  }
}

double max (const double &a , const double &b) {
  return (a > b) ? a : b;
}


void ChangeTriangleConfiguration1 (const std::vector<unsigned> & ENVN, std::vector <double> &angle) {
  double scale;
  if (ENVN[0] < ENVN[1] && ENVN[0] < ENVN[2]) {
    scale = (M_PI - angle[0]) / (angle[1] + angle [2]);
    angle[1] *= scale;
    angle[2] *= scale;
  }
  else if (ENVN[0] < ENVN[1] && ENVN[0] == ENVN[2]) {
    angle[1] = M_PI - 2. * angle[0];
  }
  else if (ENVN[0] <= ENVN[1]  && ENVN[0] > ENVN[2]) {
    scale = (M_PI - angle[2]) / (angle[1] + angle [0]);
    angle[1] *= scale;
    angle[0] *= scale;
  }
  else if (ENVN[0] == ENVN[1] && ENVN[0] < ENVN[2]) {
    angle[2] = M_PI - 2. * angle[0];
  }
  else if (ENVN[0] == ENVN[1] && ENVN[0] == ENVN[2]) {
    angle[0] = angle[1] = angle[2] =  M_PI / 3.;
  }
  else if (ENVN[0] > ENVN[1] && ENVN[0] <= ENVN[2]) {
    scale = (M_PI - angle[1]) / (angle[0] + angle [2]);
    angle[0] *= scale;
    angle[2] *= scale;
  }
  else if (ENVN[0] > ENVN[1] && ENVN[0] > ENVN[2]) {
    if (ENVN[1] < ENVN[2]) {
      scale = (M_PI - angle[1]) / (angle[0] + angle [2]);
      angle[0] *= scale;
      angle[2] *= scale;
    }
    else if (ENVN[1] == ENVN[2]) {
      angle[0] = M_PI - 2. * angle[1];
    }
    else if (ENVN[1] > ENVN[2]) {
      scale = (M_PI - angle[2]) / (angle[0] + angle [1]);
      angle[0] *= scale;
      angle[1] *= scale;
    }
  }
}


void ChangeTriangleConfiguration2 (const std::vector<unsigned> & ENVN, std::vector <double> &angle) {
  unsigned type = 3; // there are 2 or 3 leading angles
  if (ENVN[0] < ENVN[1]) { // 0 leads on 1
    if (ENVN[0] < ENVN[2]) type = 0; // 0 is leading angle
    else if (ENVN[0] > ENVN[2]) type = 2; // 2 is leading angle
  }
  else if (ENVN[0] > ENVN[1]) { // 1 leads on 0
    if (ENVN[1] < ENVN[2]) type = 1; // 1 is leading angle
    else if (ENVN[1] > ENVN[2]) type = 2; // 2 is leading angle
  }
  else { // 0 equals 1
    if (ENVN[0] > ENVN[2]) type = 2; // 2 is leading angle
  }

  double scale;
  if (type == 0) {
    scale = (M_PI - angle[0]) / (angle[1] + angle [2]);
    angle[1] *= scale;
    angle[2] *= scale;
  }
  else if (type == 1) {
    scale = (M_PI - angle[1]) / (angle[0] + angle [2]);
    angle[0] *= scale;
    angle[2] *= scale;
  }
  else if (type == 2) {
    scale = (M_PI - angle[2]) / (angle[1] + angle [0]);
    angle[1] *= scale;
    angle[0] *= scale;
  }
  else {
    scale = M_PI / (angle[0] + angle[1] + angle[2]);
    angle[0] *= scale;
    angle[1] *= scale;
    angle[2] *= scale;
  }
}

double GetPWillmoreEnergy (MultiLevelSolution &mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* sol  = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);
  elem* el =  msh->el;

  // Convenience variables to encode the dimension.
  const unsigned dim = 2;
  const unsigned DIM = 3;

  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // Extract the solution vector; get solDx positions in the ml_sol object.
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol.GetIndex ("Dx1");
  solDxIndex[1] = mlSol.GetIndex ("Dx2");
  solDxIndex[2] = mlSol.GetIndex ("Dx3");

  // Extract the finite element type for solx.
  unsigned solxType;
  solxType = mlSol.GetSolutionType (solDxIndex[0]);

  // Define solx and solxOld.
  std::vector < double > solx[DIM];
  std::vector < double > solxOld[DIM];

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get positions of Y in the ml_sol object.
  unsigned solYIndex[DIM];
  solYIndex[0] = mlSol.GetIndex ("Y1");
  solYIndex[1] = mlSol.GetIndex ("Y2");
  solYIndex[2] = mlSol.GetIndex ("Y3");

  // Extract the finite element type for Y.
  unsigned solYType;
  solYType = mlSol.GetSolutionType (solYIndex[0]);


  // Define solY and solYOld.
  std::vector < double > solY[DIM];
  std::vector < double > solYOld[DIM];

  // Initialize area, volume, P-Willmore energy.

  double energy = 0.;

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Number of solution element dofs.
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nxDofs  = msh->GetElementDofNumber (iel, solxType);
    unsigned nYDofs  = msh->GetElementDofNumber (iel, solYType);

    // Resize solution vectors.
    for (unsigned K = 0; K < DIM; K++) {
      solx[K].resize (nxDofs);
      solxOld[K].resize (nxDofs);
      solY[K].resize (nYDofs);
      solYOld[K].resize (nYDofs);
    }


    // Loop which handles local storage of global mapping and solution X.
    for (unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof (i, iel, solxType);
      unsigned iXDof  = msh->GetSolutionDof (i, iel, xType);

      for (unsigned K = 0; K < DIM; K++) {
        solxOld[K][i] = (*msh->_topology->_Sol[K]) (iXDof) + (*sol->_SolOld[solDxIndex[K]]) (iDDof);
        solx[K][i] = (*msh->_topology->_Sol[K]) (iXDof) + (*sol->_Sol[solDxIndex[K]]) (iDDof);
      }
    }

    // Loop which handles local storage of global mapping and solution Y.
    for (unsigned i = 0; i < nYDofs; i++) {

      // Global-to-local mapping between solution node and solution dof.
      unsigned iYDof = msh->GetSolutionDof (i, iel, solYType);
      for (unsigned K = 0; K < DIM; K++) {

        // Global-to-local solutions.
        solYOld[K][i] = (*sol->_SolOld[solYIndex[K]]) (iYDof);
        solY[K][i] = (*sol->_Sol[solYIndex[K]]) (iYDof);

      }
    }

    // begin GAUSS POINT LOOP
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solxType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phix_uv[dim]; // first order derivatives in (u,v)

      const double *phiY;  // local test function


      double weight; // gauss point weight

      //Extract Gauss point weight, test functions, and their partial derivatives.
      // "0" is derivative in u, "1" is derivative in v.
      weight = msh->_finiteElement[ielGeom][solxType]->GetGaussWeight (ig);

      phix = msh->_finiteElement[ielGeom][solxType]->GetPhi (ig);
      phix_uv[0] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDXi (ig);
      phix_uv[1] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDEta (ig);

      phiY = msh->_finiteElement[ielGeom][solYType]->GetPhi (ig);


      double solYg[3] = {0., 0., 0.};
      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for (unsigned K = 0; K < DIM; K++) {
        for (unsigned i = 0; i < nYDofs; i++) {
          solYg[K] += phiY[i] * 0.5 * (solYOld[K][i] + solY[K][i]);
        }
        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * 0.5 * (solx[K][i] + solxOld[K][i]);
          }
        }
      }

      // Computing the metric, metric determinant, and area element.
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

      // Computing the unit normal vector N.
      double normal[DIM];
      normal[0] = normalSign * (solx_uv[1][0] * solx_uv[2][1]
      - solx_uv[2][0] * solx_uv[1][1]) / sqrt (detg);
      normal[1] = normalSign * (solx_uv[2][0] * solx_uv[0][1]
      - solx_uv[0][0] * solx_uv[2][1]) / sqrt (detg);
      normal[2] = normalSign * (solx_uv[0][0] * solx_uv[1][1]
      - solx_uv[1][0] * solx_uv[0][1]) / sqrt (detg);

      // Computing Y.N and |Y|^2, which are essentially 2H and 4H^2.
      double YdotN = 0.;
      double YdotY = 0.;
      for (unsigned K = 0; K < DIM; K++) {
        YdotN += solYg[K] * normal[K];
        YdotY += solYg[K] * solYg[K];
      }
      // double signYdotN = (YdotN.value() >= 0.) ? 1. : -1.;
      double signYdotN = 1.;

      // Some necessary quantities when working with polynomials.
      double sumP3 = 0.;
      for (unsigned p = 0; p < 3; p++) {
        double signP = (P[p] % 2u == 0) ? 1. : signYdotN;
        sumP3 += signP * ap[p] * pow (YdotY, P[p] / 2.);
      }

      energy += sumP3 * Area;

    } // end GAUSS POINT LOOP.


  } // End ELEMENT LOOP for each process.


  double energyAll;
  MPI_Reduce (&energy, &energyAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  return energyAll;

}
