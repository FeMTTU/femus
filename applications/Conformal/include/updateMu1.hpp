


void UpdateMu(MultiLevelSolution& mlSol) {

  //MultiLevelSolution*  mlSol = ml_prob._ml_sol;
  unsigned level = mlSol.GetMLMesh()->GetNumberOfLevels() - 1u;

  Solution* sol = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol.GetMLMesh()->GetLevel(level);
  elem* el = msh->GetMeshElements();

  //unsigned  dim = msh->GetDimension();
  unsigned dim = 2;
  unsigned DIM = 3;

  std::vector < unsigned > indexDx(DIM);
  indexDx[0] = mlSol.GetIndex("Dx1");
  indexDx[1] = mlSol.GetIndex("Dx2");
  indexDx[2] = mlSol.GetIndex("Dx3");
  unsigned solTypeDx = mlSol.GetSolutionType(indexDx[0]);

  std::vector < unsigned > indexMu(dim);
  indexMu[0] = mlSol.GetIndex("mu1");
  indexMu[1] = mlSol.GetIndex("mu2");


  unsigned indexW1 = mlSol.GetIndex("weight1");
  unsigned solType1 = mlSol.GetSolutionType(indexMu[0]);


  std::vector< double > dof1;

  std::vector < std::vector < double > > solx(DIM);
  std::vector < std::vector < double > > solNx(DIM);

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexMu[k]]->zero();
  }
  sol->_Sol[indexW1]->zero();

  double weight; // gauss point weight

  unsigned iproc = msh->processor_id();
  unsigned nprocs = msh->n_processors();

// Setting the reference elements to be equilateral triangles.
  std::vector < std::vector < double > > xT(2);
  xT[0].resize(7);
  xT[0][0] = -0.5;
  xT[0][1] = 0.5;
  xT[0][2] = 0.;
  xT[0][3] = 0.5 * (xT[0][0] + xT[0][1]);
  xT[0][4] = 0.5 * (xT[0][1] + xT[0][2]);
  xT[0][5] = 0.5 * (xT[0][2] + xT[0][0]);
  xT[0][6] = 0.33333333333333333333 * (xT[0][0] + xT[0][1] + xT[0][2]);

  xT[1].resize(7);
  xT[1][0] = 0.;
  xT[1][1] = 0.;
  xT[1][2] = sqrt(3.) / 2.;
  xT[1][3] = 0.5 * (xT[1][0] + xT[1][1]);
  xT[1][4] = 0.5 * (xT[1][1] + xT[1][2]);
  xT[1][5] = 0.5 * (xT[1][2] + xT[1][0]);
  xT[1][6] = 0.33333333333333333333 * (xT[1][0] + xT[1][1] + xT[1][2]);

  double angles[2][4] = {
    {0., 0.5 * M_PI, M_PI, 1.5 * M_PI}, // for square
    {0., 2. / 3. * M_PI, 4. / 3 * M_PI} // for equilateral triangle
  };

  unsigned solENVNIndex = mlSol.GetIndex("ENVN");
  unsigned solENVNType = mlSol.GetSolutionType(solENVNIndex);

  std::vector<double> phi_uv0;
  std::vector<double> phi_uv1;

  std::vector< double > stdVectorPhi;
  std::vector< double > stdVectorPhi_uv;

  for(int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs1  = msh->GetElementDofNumber(iel, solType1);
    unsigned nDofsDx  = msh->GetElementDofNumber(iel, solTypeDx);

    dof1.resize(nDofs1);

    for(int K = 0; K < DIM; K++) {
      solx[K].resize(nDofsDx);
      solNx[K].resize(nDofsDx);
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs1; i++) {
      dof1[i] = msh->GetSolutionDof(i, iel, solType1);
    }
    // local storage of coordinates
    for(unsigned i = 0; i < nDofsDx; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solTypeDx);
      unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
      for(unsigned K = 0; K < DIM; K++) {
        solx[K][i] = (*msh->GetTopology()->_Sol[K])(xDof);
        solNx[K][i] = solx[K][i]  + (*sol->_Sol[indexDx[K]])(idof);
      }
    }


//     if(ielGeom == TRI) {
// 
//       xT[0][1] = 0.5;
//       std::vector < unsigned > ENVN(3);
//       std::vector < double > angle(3);
// 
//       for(unsigned j = 0; j < 3; j++) {
//         unsigned jnode  = msh->GetSolutionDof(j, iel, solENVNType);
//         ENVN[j] = (*sol->_Sol[solENVNIndex])(jnode);
//         angle[j] = 2 * M_PI / ENVN[j];
//       }
// 
// 
//       if(conformalTriangleType == 1) {  //this works with moo two levels
//         ChangeTriangleConfiguration1(ENVN, angle);
//       }
//       else if(conformalTriangleType == 2) {  //this works with mao
//         ChangeTriangleConfiguration2(ENVN, angle);
//       }
//       else { //no change
//         angle.assign(3, M_PI / 3.);
//       }
// 
//       double l = xT[0][1] - xT[0][0];
//       double d = l * sin(angle[0]) * sin(angle[1]) / sin(angle[0] + angle[1]);
//       double scale = sqrt((sqrt(3.) / 2.) / (l * d));
//       l = l * scale;
//       d = d * scale;
//       xT[0][1] = xT[0][0] + l;
//       xT[0][2] = xT[0][0] + d / tan(angle[0]);
//       xT[1][2] = d;
// 
//       xT[0][3] = 0.5 * (xT[0][0] + xT[0][1]);
//       xT[0][4] = 0.5 * (xT[0][1] + xT[0][2]);
//       xT[0][5] = 0.5 * (xT[0][2] + xT[0][0]);
//       xT[0][6] = 0.33333333333333333333 * (xT[0][0] + xT[0][1] + xT[0][2]);
// 
//       xT[1][3] = 0.5 * (xT[1][0] + xT[1][1]);
//       xT[1][4] = 0.5 * (xT[1][1] + xT[1][2]);
//       xT[1][5] = 0.5 * (xT[1][2] + xT[1][0]);
//       xT[1][6] = 0.33333333333333333333 * (xT[1][0] + xT[1][1] + xT[1][2]);
// 
//       //std::cout << l << " " << d<<" "<< angle[0] << " " << angle[1] <<" "<< angle[2] << " " << l * d <<" "<< xT[0][2]<< " " << xT[1][2]<<  std::endl;
//     }


// *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solTypeDx]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function

      const double *phix_uv[dim]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // Get Gauss point weight, test function, and first order derivatives.
      if(ielGeom == QUAD) {
        phix = msh->_finiteElement[ielGeom][solTypeDx]->GetPhi(ig);


        phix_uv[0] = msh->_finiteElement[ielGeom][solTypeDx]->GetDPhiDXi(ig);
        phix_uv[1] = msh->_finiteElement[ielGeom][solTypeDx]->GetDPhiDEta(ig);

        weight = msh->_finiteElement[ielGeom][solTypeDx]->GetGaussWeight(ig);
      }

      // Special adjustments for triangles.
      else {
        msh->_finiteElement[ielGeom][solTypeDx]->Jacobian(xT, ig, weight, stdVectorPhi, stdVectorPhi_uv);

        phix = &stdVectorPhi[0];
        phi_uv0.resize(nDofsDx);
        phi_uv1.resize(nDofsDx);


        for(unsigned i = 0; i < nDofsDx; i++) {
          phi_uv0[i] = stdVectorPhi_uv[i * dim];
          phi_uv1[i] = stdVectorPhi_uv[i * dim + 1];
        }

        phix_uv[0] = &phi_uv0[0];
        phix_uv[1] = &phi_uv1[0];

      }
      const double *phi1 = msh->_finiteElement[ielGeom][solType1]->GetPhi(ig);  // local test function

      // Initialize and compute values of x, Dx, NDx, x_uv at the Gauss points.
      double solxg[3] = {0., 0., 0.};
      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      double solNx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for(unsigned K = 0; K < DIM; K++) {
        for(unsigned i = 0; i < nDofsDx; i++) {
          solxg[K] += phix[i] * solx[K][i];
        }
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nDofsDx; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * solx[K][i];
            solNx_uv[K][j]    += phix_uv[j][i] * solNx[K][i];
          }
        }
      }

      // Compute the metric, metric determinant, and area element.
      double g[2][2] = {{0., 0.}, {0., 0.}};
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];

      double normal[DIM];
      normal[0] = (solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) / sqrt(detg);
      normal[1] = (solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt(detg);
      normal[2] = (solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt(detg);


//       //Analytic for the cylinder
//       normal[0] = 0;
//       normal[1] = solxg[1] / sqrt(solxg[1] * solxg[1] + solxg[2] * solxg[2]);
//       normal[2] = solxg[2] / sqrt(solxg[1] * solxg[1] + solxg[2] * solxg[2]);

      double dxPlus[DIM];
      dxPlus[0] = solNx_uv[0][0] + solNx_uv[1][1] * normal[2] - solNx_uv[2][1] * normal[1];
      dxPlus[1] = solNx_uv[1][0] + solNx_uv[2][1] * normal[0] - solNx_uv[0][1] * normal[2];
      dxPlus[2] = solNx_uv[2][0] + solNx_uv[0][1] * normal[1] - solNx_uv[1][1] * normal[0];

      double sdxPlus[DIM];
      sdxPlus[0] = solNx_uv[0][1] - solNx_uv[1][0] * normal[2] + solNx_uv[2][0] * normal[1];
      sdxPlus[1] = solNx_uv[1][1] - solNx_uv[2][0] * normal[0] + solNx_uv[0][0] * normal[2];
      sdxPlus[2] = solNx_uv[2][1] - solNx_uv[0][0] * normal[1] + solNx_uv[1][0] * normal[0];

      double dxMinus[DIM];
      dxMinus[0] = solNx_uv[0][0] - solNx_uv[1][1] * normal[2] + solNx_uv[2][1] * normal[1];
      dxMinus[1] = solNx_uv[1][0] - solNx_uv[2][1] * normal[0] + solNx_uv[0][1] * normal[2];
      dxMinus[2] = solNx_uv[2][0] - solNx_uv[0][1] * normal[1] + solNx_uv[1][1] * normal[0];

      double norm2dxPlus = 0;
      double norm2sdxPlus = 0;
      double rhsmu1 = 0;
      double rhsmu2 = 0;

      //double dxsdxp = 0.;
      for(unsigned K = 0; K < DIM; K++) {
        norm2dxPlus += dxPlus[K] * dxPlus[K];
        norm2sdxPlus += sdxPlus[K] * sdxPlus[K];
        rhsmu1 += dxPlus[K] * dxMinus[K];
        rhsmu2 += sdxPlus[K] * dxMinus[K];
      }

      // Comment out for working code
      double mu[2] = {0., 0.};
      mu[0] = rhsmu1 / norm2dxPlus;
      mu[1] = rhsmu2 / norm2sdxPlus;

      for(unsigned i = 0; i < nDofs1; i++) {
        sol->_Sol[indexW1]->add(dof1[i], phi1[i] * weight);
        for(unsigned k = 0; k < dim; k++) {
          sol->_Sol[indexMu[k]]->add(dof1[i], mu[k] * phi1[i] * weight);
        }
      } // end phi_i loop
    } // end gauss point loop
  } //end element loop for each process*/

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexMu[k]]->close();
  }
  sol->_Sol[indexW1]->close();

  for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {

    double weight = (*sol->_Sol[indexW1])(i);

    double mu[2];
    for(unsigned k = 0; k < dim; k++) {
      mu[k] = (*sol->_Sol[indexMu[k]])(i) / weight;
      sol->_Sol[indexMu[k]]->set(i, mu[k]);
    }
  }
  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexMu[k]]->close();
  }




  double MuNormAverageBefore;
  {
    //Norm before the smoothing
    double MuNormLocalSum = 0.;
    for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {
      MuNormLocalSum += sqrt(pow((*sol->_Sol[indexMu[0]])(i), 2) + pow((*sol->_Sol[indexMu[1]])(i), 2));
    }
    MPI_Allreduce(&MuNormLocalSum, &MuNormAverageBefore, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MuNormAverageBefore /= msh->_dofOffset[solType1][nprocs];

    std::cout << " before " << MuNormAverageBefore << std::endl;
  }

  std::vector < unsigned > indexMuEdge(dim);
  indexMuEdge[0] = mlSol.GetIndex("mu1Edge");
  indexMuEdge[1] = mlSol.GetIndex("mu2Edge");
  unsigned indexCntEdge = mlSol.GetIndex("cntEdge");
  unsigned solType2 = mlSol.GetSolutionType(indexMuEdge[0]);

  unsigned smoothMax = 6;//2 + counter/10;
  //if(smoothMax > 10) smoothMax = 10;
  std::cout << "Max Number of Smoothing = " << smoothMax << std::endl;

  for(unsigned ismooth = 0; ismooth < smoothMax; ismooth++) {

    for(unsigned k = 0; k < dim; k++) {
      sol->_Sol[indexMuEdge[k]]->zero();
    }
    sol->_Sol[indexCntEdge]->zero();

    for(int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned idx = (ielGeom == QUAD) ? 0 : 1;

      double mu[2];
      for(unsigned k = 0; k < 2; k++) {
        mu[k] = (*sol->_Sol[indexMu[k]])(iel);
      }

//       if(ielGeom == TRI) {
// 
//         xT[0][1] = 0.5;
//         std::vector < unsigned > ENVN(3);
//         std::vector < double > angle(3);
// 
//         for(unsigned j = 0; j < 3; j++) {
//           unsigned jnode  = msh->GetSolutionDof(j, iel, solENVNType);
//           ENVN[j] = (*sol->_Sol[solENVNIndex])(jnode);
//           angle[j] = 2 * M_PI / ENVN[j];
//         }
// 
// 
//         if(conformalTriangleType == 1) {  //this works with moo two levels
//           ChangeTriangleConfiguration1(ENVN, angle);
//         }
//         else if(conformalTriangleType == 2) {  //this works with mao
//           ChangeTriangleConfiguration2(ENVN, angle);
//         }
//         else { //no change
//           angle.assign(3, M_PI / 3.);
//         }
// 
//         double l = xT[0][1] - xT[0][0];
//         double d = l * sin(angle[0]) * sin(angle[1]) / sin(angle[0] + angle[1]);
//         double scale = sqrt((sqrt(3.) / 2.) / (l * d));
//         l = l * scale;
//         d = d * scale;
//         xT[0][1] = xT[0][0] + l;
//         xT[0][2] = xT[0][0] + d / tan(angle[0]);
//         xT[1][2] = d;
// 
//         angles[idx][1] = atan2(xT[1][2], xT[0][2] - xT[0][1]);
//         angles[idx][2] = atan2(-xT[1][2], -0.5 + xT[0][2]);
// 
//       }

      unsigned nDofs0  = msh->GetElementDofNumber(iel, 0);
      for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {

        unsigned idof = msh->GetSolutionDof(nDofs0 + iface, iel, solType2);

        //double weight = (iface % 2) ? -1. : 1.;
        //sol->_Sol[indexMuEdge[0]]->add(idof, weight * mu[0]);
        //sol->_Sol[indexMuEdge[1]]->add(idof, weight * mu[1]);

        double a = cos(angles[idx][iface]);
        double b = sin(angles[idx][iface]);

        double mu0s = (a * a - b * b) * mu[0] + 2. * a * b * mu[1];
        double mu1s = (a * a - b * b) * mu[1] - 2. * a * b * mu[0];

        sol->_Sol[indexMuEdge[0]]->add(idof, mu0s);
        sol->_Sol[indexMuEdge[1]]->add(idof, mu1s);

        sol->_Sol[indexCntEdge]->add(idof, 1);
      }
    }
    for(unsigned k = 0; k < 2; k++) {
      sol->_Sol[indexMuEdge[k]]->close();
    }
    sol->_Sol[indexCntEdge]->close();

    for(int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned idx = (ielGeom == QUAD) ? 0 : 1;

//       if(ielGeom == TRI) {
// 
//         xT[0][1] = 0.5;
//         std::vector < unsigned > ENVN(3);
//         std::vector < double > angle(3);
// 
//         for(unsigned j = 0; j < 3; j++) {
//           unsigned jnode  = msh->GetSolutionDof(j, iel, solENVNType);
//           ENVN[j] = (*sol->_Sol[solENVNIndex])(jnode);
//           angle[j] = 2 * M_PI / ENVN[j];
//         }
// 
// 
//         if(conformalTriangleType == 1) {  //this works with moo two levels
//           ChangeTriangleConfiguration1(ENVN, angle);
//         }
//         else if(conformalTriangleType == 2) {  //this works with mao
//           ChangeTriangleConfiguration2(ENVN, angle);
//         }
//         else { //no change
//           angle.assign(3, M_PI / 3.);
//         }
// 
//         double l = xT[0][1] - xT[0][0];
//         double d = l * sin(angle[0]) * sin(angle[1]) / sin(angle[0] + angle[1]);
//         double scale = sqrt((sqrt(3.) / 2.) / (l * d));
//         l = l * scale;
//         d = d * scale;
//         xT[0][1] = xT[0][0] + l;
//         xT[0][2] = xT[0][0] + d / tan(angle[0]);
//         xT[1][2] = d;
// 
//         angles[idx][1] = atan2(xT[1][2], xT[0][2] - xT[0][1]);
//         angles[idx][2] = atan2(-xT[1][2], -0.5 + xT[0][2]);
// 
//       }

      double mu[2] = {0., 0.};
      double cnt = 0.;

      unsigned nDofs0  = msh->GetElementDofNumber(iel, 0);
      for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {

        unsigned idof = msh->GetSolutionDof(nDofs0 + iface, iel, solType2);
        //double sign = (iface % 2) ? -1. : 1.;
//         mu[0] += sign * (*sol->_Sol[indexMuEdge[0]])(idof);
//         mu[1] += sign * (*sol->_Sol[indexMuEdge[1]])(idof);


        double mu0s = (*sol->_Sol[indexMuEdge[0]])(idof);
        double mu1s = (*sol->_Sol[indexMuEdge[1]])(idof);


        double a = cos(angles[idx][iface]);
        double b = sin(angles[idx][iface]);

        mu[0] += (a * a - b * b) * mu0s - 2. * a * b * mu1s;
        mu[1] += (a * a - b * b) * mu1s + 2. * a * b * mu0s;

        cnt += (*sol->_Sol[indexCntEdge])(idof);
      }

      for(unsigned k = 0; k < 2; k++) {
        sol->_Sol[indexMu[k]]->set(iel, mu[k] / cnt);
      }

    }
    for(unsigned k = 0; k < 2; k++) {
      sol->_Sol[indexMu[k]]->close();
    }
  }


  double MuNormAverageAfter;
  {
    //BEGIN mu update
    double MuNormLocalSum = 0.;
    for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {
      MuNormLocalSum += sqrt(pow((*sol->_Sol[indexMu[0]])(i), 2) + pow((*sol->_Sol[indexMu[1]])(i), 2));
    }
    MPI_Allreduce(&MuNormLocalSum, &MuNormAverageAfter, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MuNormAverageAfter /= msh->_dofOffset[solType1][nprocs];
  }
  std::cout << " after " << MuNormAverageAfter << std::endl;

  for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {

    double mu[2];
    for(unsigned k = 0; k < 2; k++) {
      mu[k] = (*sol->_Sol[indexMu[k]])(i);
    }

    double norm = sqrt(mu[0] * mu[0] + mu[1] * mu[1]);
    double cosTheta = mu[0] / norm;
    double sinTheta = mu[1] / norm;

    sol->_Sol[indexMu[0]]->set(i, MuNormAverageBefore * cosTheta);
    sol->_Sol[indexMu[1]]->set(i, MuNormAverageBefore * sinTheta);

  }
  for(unsigned k = 0; k < 2; k++) {
    sol->_Sol[indexMu[k]]->close();
  }





}
