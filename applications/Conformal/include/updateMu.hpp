 
void UpdateMu(MultiLevelSolution& mlSol) {

  //MultiLevelSolution*  mlSol = ml_prob._ml_sol;
  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1u;

  Solution* sol = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);
  elem* el = msh->el;

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

  unsigned indexMuN1 = mlSol.GetIndex("muN1"); //piecewice linear discontinuous

  unsigned indexW1 = mlSol.GetIndex("weight1");
  unsigned solType1 = mlSol.GetSolutionType(indexMu[0]);

  unsigned indexTheta1 = mlSol.GetIndex("theta1");
  unsigned indexTheta2 = mlSol.GetIndex("theta2");

  unsigned indexPhi1 = mlSol.GetIndex("phi1");
  unsigned indexPhi2 = mlSol.GetIndex("phi2");

  std::vector< double > dof1;

  std::vector < std::vector < double > > solx(DIM);
  std::vector < std::vector < double > > xHat(DIM);

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexMu[k]]->zero();
  }
  sol->_Sol[indexW1]->zero();

  std::vector < double > phi;  // local test function for velocity
  std::vector < double > phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  unsigned iproc = msh->processor_id();
  unsigned nprocs = msh->n_processors();

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

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs1  = msh->GetElementDofNumber(iel, solType1);
    unsigned nDofsDx  = msh->GetElementDofNumber(iel, solTypeDx);

    dof1.resize(nDofs1);

    for(int K = 0; K < DIM; K++) {
      xHat[K].resize(nDofsDx);
      solx[K].resize(nDofsDx);
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
        xHat[K][i] = (*msh->_topology->_Sol[K])(xDof);
        solx[K][i] = xHat[K][i] + (*sol->_Sol[indexDx[K]])(idof);
      }
    }

//####### DIDN'T CHANGE THIS
    // for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solTypeDx]->GetGaussPointNumber(); ig++) {

    // msh->_finiteElement[ielGeom][solTypeDx]->Jacobian(xHat, ig, weight, phi, phi_x);
    //
    // // std::vector < std::vector < double > > gradSolx(dim);
    // //
    // // for(unsigned  k = 0; k < dim; k++) {
    // //   gradSolx[k].assign(dim, 0.);
    // // }
    //
    // for(unsigned i = 0; i < nDofsDx; i++) {
    //   for(unsigned j = 0; j < dim; j++) {
    //     for(unsigned  k = 0; k < dim; k++) {
    //       gradSolx[k][j] += solx[k][i] * phi_x[i * dim + j];
    //     }
    //   }
    // }
// ##############

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
      //double solDxg[3] = {0., 0., 0.};
      double solxg[3] = {0., 0., 0.};
      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for(unsigned K = 0; K < DIM; K++) {
        for(unsigned i = 0; i < nDofsDx; i++) {
          solxg[K] += phix[i] * solx[K][i];
        }
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nDofsDx; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * solx[K][i];
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
      //double Area = weight * sqrt(detg);
      //double Area2 = weight; // Trick to give equal weight to each element.

      // std::cout << detg << " ";

      double normal[DIM];
      normal[0] = (solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) / sqrt(detg);
      normal[1] = (solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt(detg);
      normal[2] = (solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt(detg);


//       //Analytic for the cylinder
//       normal[0] = 0;
//       normal[1] = solxg[1] / sqrt(solxg[1] * solxg[1] + solxg[2] * solxg[2]);
//       normal[2] = solxg[2] / sqrt(solxg[1] * solxg[1] + solxg[2] * solxg[2]);

      double dxPlus[DIM];
      dxPlus[0] = solx_uv[0][0] + solx_uv[1][1] * normal[2] - solx_uv[2][1] * normal[1];
      dxPlus[1] = solx_uv[1][0] + solx_uv[2][1] * normal[0] - solx_uv[0][1] * normal[2];
      dxPlus[2] = solx_uv[2][0] + solx_uv[0][1] * normal[1] - solx_uv[1][1] * normal[0];

      double sdxPlus[DIM];
      sdxPlus[0] = solx_uv[0][1] - solx_uv[1][0] * normal[2] + solx_uv[2][0] * normal[1];
      sdxPlus[1] = solx_uv[1][1] - solx_uv[2][0] * normal[0] + solx_uv[0][0] * normal[2];
      sdxPlus[2] = solx_uv[2][1] - solx_uv[0][0] * normal[1] + solx_uv[1][0] * normal[0];

      double dxMinus[DIM];
      dxMinus[0] = solx_uv[0][0] - solx_uv[1][1] * normal[2] + solx_uv[2][1] * normal[1];
      dxMinus[1] = solx_uv[1][0] - solx_uv[2][1] * normal[0] + solx_uv[0][1] * normal[2];
      dxMinus[2] = solx_uv[2][0] - solx_uv[0][1] * normal[1] + solx_uv[1][1] * normal[0];

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

        //dxsdxp += dxMinus[K] * dxMinus[K];
      }

      // double norm2Xz = (1. / 4.) * (pow((gradSolx[0][0] + gradSolx[1][1]), 2) + pow((gradSolx[1][0] - gradSolx[0][1]), 2));
      // double XzBarXz_Bar[2];
      //
      // XzBarXz_Bar[0] = (1. / 4.) * (pow(gradSolx[0][0], 2) + pow(gradSolx[1][0], 2) - pow(gradSolx[0][1], 2) - pow(gradSolx[1][1], 2));
      // XzBarXz_Bar[1] = (1. / 2.) * (gradSolx[0][0] * gradSolx[0][1] + gradSolx[1][0] * gradSolx[1][1]);

      // Comment out for working code

      double mu[2] = {0., 0.};
      mu[0] = rhsmu1 / norm2dxPlus;
      mu[1] = rhsmu2 / norm2sdxPlus;

//       if(iel == 4) {
//         std::cout << mu[0] << " " << mu[1] << " " << norm2dxPlus << " " << norm2sdxPlus << "\n";
//       }

      // for(unsigned k = 0; k < 2; k++) {
      //   if(norm2Xz > 0.) {
      //     mu[k] += (1. / norm2Xz) * XzBarXz_Bar[k];
      //   }
      // }

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

  sol->_Sol[indexTheta1]->zero();
  sol->_Sol[indexPhi1]->zero();

  for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {

    double weight = (*sol->_Sol[indexW1])(i);

    //std::cout << weight << " ";

    double mu[2];
    for(unsigned k = 0; k < dim; k++) {
      mu[k] = (*sol->_Sol[indexMu[k]])(i) / weight;
      sol->_Sol[indexMu[k]]->set(i, mu[k]);
    }
    sol->_Sol[indexMuN1]->set(i, sqrt(mu[0] * mu[0] + mu[1] * mu[1]));

    if(muSmoothingType == 1) {
      sol->_Sol[indexTheta1]->set(i, atan2(mu[1], fabs(mu[0])));
      sol->_Sol[indexPhi1]->set(i, atan2(mu[0], fabs(mu[1])));
    }
    else if(muSmoothingType == 2) { //invariant on quad element orientation
      sol->_Sol[indexTheta1]->set(i, atan(mu[1] / mu[0]));
      sol->_Sol[indexPhi1]->set(i, atan(mu[0] / mu[1]));
    }
//     else if(muSmoothingType == 3) { //invariant on quad element orientation
//       sol->_Sol[indexTheta1]->set(i, sin(2. * atan2(mu[1], mu[0])));
//       sol->_Sol[indexPhi1]->set(i, cos(2.*atan2(mu[0], mu[1])));
//     }

  }
  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexMu[k]]->close();
  }
  sol->_Sol[indexMuN1]->close();
  sol->_Sol[indexTheta1]->close();
  sol->_Sol[indexPhi1]->close();

  double norm = sol->_Sol[indexMuN1]->linfty_norm();
  std::cout << norm << std::endl;

  //BEGIN Iterative smoothing element -> nodes -> element

  for(unsigned smooth = 0; smooth < 1; smooth++) {
    unsigned indexW2 = mlSol.GetIndex("weight2");
    unsigned indexMuN2 = mlSol.GetIndex("muN2");  //smooth ni norm

    unsigned solType2 = mlSol.GetSolutionType(indexMuN2);

    std::vector< double > dof2;
    std::vector< double > sol1;
    std::vector< double > solTheta1;
    std::vector< double > solPhi1;

    sol->_Sol[indexMuN2]->zero();
    sol->_Sol[indexW2]->zero();
    sol->_Sol[indexTheta2]->zero();
    sol->_Sol[indexPhi2]->zero();

    //std::vector < double > phi2;  // local test function for velocity
    //std::vector < double > phi2_x; // local test function first order partial derivatives
    //double weight2; // gauss point weight

    for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofs1  = msh->GetElementDofNumber(iel, solType1);
      unsigned nDofs2  = msh->GetElementDofNumber(iel, solType2);

      sol1.resize(nDofs1);
      solTheta1.resize(nDofs1);
      solPhi1.resize(nDofs1);

      dof2.resize(nDofs2);

      for(int K = 0; K < DIM; K++) {
        xHat[K].resize(nDofs2);
      }

      for(unsigned i = 0; i < nDofs1; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solType1);
        sol1[i] = (*sol->_Sol[indexMuN1])(idof);
        solTheta1[i] = (*sol->_Sol[indexTheta1])(idof);
        solPhi1[i] = (*sol->_Sol[indexPhi1])(idof);
      }

      // local storage of global mapping and solution
      for(unsigned i = 0; i < nDofs2; i++) {
        dof2[i] = msh->GetSolutionDof(i, iel, solType2);
      }
      // local storage of coordinates
      for(unsigned i = 0; i < nDofs2; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
        for(unsigned K = 0; K < DIM; K++) {
          xHat[K][i] = (*msh->_topology->_Sol[K])(xDof);
        }
      }

      for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType2]->GetGaussPointNumber(); ig++) {
//         msh->_finiteElement[ielGeom][solType2]->Jacobian(xHat, ig, weight2, phi2, phi2_x);

        double *phi2 = msh->_finiteElement[ielGeom][solType2]->GetPhi(ig);
        double weight2 = msh->_finiteElement[ielGeom][solType2]->GetGaussWeight(ig);
        double *phi1 = msh->_finiteElement[ielGeom][solType1]->GetPhi(ig);
        double sol1g = 0.;
        double solTheta1g = 0.;
        double solPhi1g = 0.;
        for(unsigned i = 0; i < nDofs1; i++) {
          sol1g += phi1[i] * sol1[i];
          solTheta1g += phi1[i] * solTheta1[i];
          solPhi1g += phi1[i] * solPhi1[i];
        }

        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofs2; i++) {
          sol->_Sol[indexW2]->add(dof2[i], phi2[i] * weight2);
          sol->_Sol[indexMuN2]->add(dof2[i], sol1g * phi2[i] * weight2);
          sol->_Sol[indexTheta2]->add(dof2[i], solTheta1g * phi2[i] * weight2);
          sol->_Sol[indexPhi2]->add(dof2[i], solPhi1g * phi2[i] * weight2);
        } // end phi_i loop
      } // end gauss point loop

    } //end element loop for each process*/
    sol->_Sol[indexW2]->close();
    sol->_Sol[indexMuN2]->close();
    sol->_Sol[indexTheta2]->close();
    sol->_Sol[indexPhi2]->close();

    for(unsigned i = msh->_dofOffset[solType2][iproc]; i < msh->_dofOffset[solType2][iproc + 1]; i++) {
      double weight = (*sol->_Sol[indexW2])(i);
      double value = (*sol->_Sol[indexMuN2])(i);
      sol->_Sol[indexMuN2]->set(i, value / weight);

      value = (*sol->_Sol[indexTheta2])(i);
      sol->_Sol[indexTheta2]->set(i, value / weight);

      value = (*sol->_Sol[indexPhi2])(i);
      sol->_Sol[indexPhi2]->set(i, value / weight);
    }
    sol->_Sol[indexMuN2]->close();
    sol->_Sol[indexTheta2]->close();
    sol->_Sol[indexPhi2]->close();


    sol->_Sol[indexMuN1]->zero();
    sol->_Sol[indexTheta1]->zero();
    sol->_Sol[indexPhi1]->zero();

    std::vector< double > sol2;
    std::vector< double > solTheta2;
    std::vector< double > solPhi2;

    for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofs1  = msh->GetElementDofNumber(iel, solType1);
      unsigned nDofs2  = msh->GetElementDofNumber(iel, solType2);

      dof1.resize(nDofs1);

      for(int k = 0; k < dim; k++) {
        //xHat[k].resize(nDofs2);
        sol2.resize(nDofs2);
        solTheta2.resize(nDofs2);
        solPhi2.resize(nDofs2);
      }
      for(int K = 0; K < DIM; K++) {
        xHat[K].resize(nDofs2);
      }

      // local storage of global mapping and solution
      for(unsigned i = 0; i < nDofs1; i++) {
        dof1[i] = msh->GetSolutionDof(i, iel, solType1);
      }
      // local storage of coordinates
      for(unsigned i = 0; i < nDofs2; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solType2);
        unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
        for(unsigned K = 0; K < DIM; K++) {
          xHat[K][i] = (*msh->_topology->_Sol[K])(xDof);
        }
        sol2[i] = (*sol->_Sol[indexMuN2])(idof);
        solTheta2[i] = (*sol->_Sol[indexTheta2])(idof);
        solPhi2[i] = (*sol->_Sol[indexPhi2])(idof);

      }

      for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solTypeDx]->GetGaussPointNumber(); ig++) {

        //msh->_finiteElement[ielGeom][solTypeDx]->Jacobian(xHat, ig, weight, phi, phi_x);

        weight = msh->_finiteElement[ielGeom][solTypeDx]->GetGaussWeight(ig);
        double *phi = msh->_finiteElement[ielGeom][solTypeDx]->GetPhi(ig);

        double sol2g = 0;
        double solTheta2g = 0;
        double solPhi2g = 0;

        for(unsigned i = 0; i < nDofs2; i++) {
          sol2g += sol2[i] * phi[i];
          solTheta2g += solTheta2[i] * phi[i];
          solPhi2g += solPhi2[i] * phi[i];
        }

        double *phi1 = msh->_finiteElement[ielGeom][solType1]->GetPhi(ig);

        for(unsigned i = 0; i < nDofs1; i++) {
          sol->_Sol[indexMuN1]->add(dof1[i], sol2g * phi1[i] * weight);
          sol->_Sol[indexTheta1]->add(dof1[i], solTheta2g * phi1[i] * weight);
          sol->_Sol[indexPhi1]->add(dof1[i], solPhi2g * phi1[i] * weight);
        } // end phi_i loop

      } // end gauss point loop

    } //end element loop for each process*/

    sol->_Sol[indexMuN1]->close();
    sol->_Sol[indexTheta1]->close();
    sol->_Sol[indexPhi1]->close();

    for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {
      double weight = (*sol->_Sol[indexW1])(i);
      double radius = (*sol->_Sol[indexMuN1])(i) / weight;
      sol->_Sol[indexMuN1]->set(i,  radius);

      double theta = (*sol->_Sol[indexTheta1])(i) / weight;
      sol->_Sol[indexTheta1]->set(i, theta);

      double phi = (*sol->_Sol[indexPhi1])(i) / weight;
      sol->_Sol[indexPhi1]->set(i, phi);


      double mu[2];
      for(unsigned k = 0; k < dim; k++) {
        mu[k] = (*sol->_Sol[indexMu[k]])(i);
      }
    }
    sol->_Sol[indexMuN1]->close();
    sol->_Sol[indexTheta1]->close();
    sol->_Sol[indexPhi1]->close();
  }
  //END Iterative smoothing element -> nodes -> element


  //BEGIN mu update
  double MuNormLocalSum = 0.;
  for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {
    MuNormLocalSum += (*sol->_Sol[indexMuN1])(i);
  }

  double MuNormAverage;
  MPI_Allreduce(&MuNormLocalSum, &MuNormAverage, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MuNormAverage /= msh->_dofOffset[solType1][nprocs];

  for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {

    double theta = (*sol->_Sol[indexTheta1])(i);
    double phi = (*sol->_Sol[indexPhi1])(i);

    double mu[2];
    for(unsigned k = 0; k < dim; k++) {
      mu[k] = (*sol->_Sol[indexMu[k]])(i);
    }

    if(muSmoothingType == 1) {
      if(fabs(theta) < fabs(phi)) {
        if(mu[0] < 0) {
          if(mu[1] < 0) theta = -theta - M_PI;
          else theta = M_PI - theta;
        }
      }
      else {
        if(mu[1] < 0) {
          theta = -M_PI / 2 + phi;
        }
        else {
          theta =  M_PI / 2 - phi;
        }
      }
    }
    else if(muSmoothingType == 2) {
      double thetaOriginal = atan2(mu[1], mu[0]);
      if(thetaOriginal < -0.75 * M_PI) {
        theta = theta - M_PI;
      }
      else if(thetaOriginal < -0.25 * M_PI) {
        theta = -phi - 0.5 * M_PI;
      }
      else if(thetaOriginal < 0.25 * M_PI) {
        theta = theta;
      }
      else if(thetaOriginal < 0.75 * M_PI) {
        theta = -phi + 0.5 * M_PI;
      }
      else {
        theta = theta + M_PI;
      }
    }
//     else if(muSmoothingType == 3) {
//       double thetaOriginal = atan2(mu[1], mu[0]);
//       if(thetaOriginal < -7. / 8. * M_PI) {
//         theta = 0.5 * asin(theta) - M_PI;
//       }
//       else if(thetaOriginal < -5. / 8. * M_PI) {
//         theta = -0.5 * acos(phi) - 0.5 * M_PI;
//       }
//       else if(thetaOriginal < -3. / 8. * M_PI) {
//         theta = -0.5 * asin(theta) - 0.5 * M_PI;
//       }
//       else if(thetaOriginal < -1. / 8. * M_PI) {
//         theta =  0.5 * acos(phi) - 0.5 * M_PI;
//       }
//       else if(thetaOriginal < 1. / 8. * M_PI) {
//         theta = 0.5 * asin(theta);
//       }
//       else if(thetaOriginal < 3. / 8. * M_PI) {
//         theta = -0.5 * acos(phi) + 0.5 * M_PI;
//       }
//       else if(thetaOriginal < 5. / 8. * M_PI) {
//         theta = -0.5 * asin(theta) + 0.5 * M_PI;
//       }
//       else if(thetaOriginal < 7. / 8. * M_PI) {
//         theta = 0.5 * acos(phi) + 0.5 * M_PI;
//       }
//       else {
//         theta = 0.5 * asin(theta) + M_PI;
//       }
//     }

    sol->_Sol[indexTheta1]->set(i, theta);
    sol->_Sol[indexMu[0]]->set(i, MuNormAverage * cos(theta));
    sol->_Sol[indexMu[1]]->set(i, MuNormAverage * sin(theta));

  }

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexMu[k]]->close();
    sol->_Sol[indexTheta1]->close();
  }
  //END mu update
}
