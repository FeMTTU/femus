
using namespace femus;


//BEGIN Stochastic Input Parameters

unsigned pIndex = 4;
unsigned qIndex = 5;

int numberOfEigPairs = 2; //dimension of the stochastic variable
double stdDeviationInput = 0.2;  //standard deviation of the normal distribution (it is the same as the standard deviation of the covariance function in GetEigenPair)
double meanInput = 0.;
double amin = 1. / 100.; // for the KL expansion
std::vector < std::pair<double, double> > eigenvalues (numberOfEigPairs);
//END Stochastic Input Parameters



void AssembleSysSG (MultiLevelProblem& ml_prob) {

  uq &myuq = FemusInit::_uq;

  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("SG");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh* msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix* KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + ! (dim - 1));       // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27

  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  const std::vector < std::vector < std::vector < double > > > & integralMatrix = myuq.GetIntegralMatrix (qIndex, pIndex);

  const std::vector < std::vector < std::vector < double > > >  &G = myuq.GetStochasticMassMatrix (qIndex, pIndex, numberOfEigPairs);

  const std::vector < std::vector <unsigned> > &Jq = myuq.GetIndexSet (qIndex, numberOfEigPairs);
  const std::vector < std::vector <unsigned> > &Jp = myuq.GetIndexSet (pIndex, numberOfEigPairs);

  unsigned maxPolyOrder = (qIndex > pIndex) ? qIndex : pIndex;



  unsigned n1 = 2 * pIndex + qIndex + 1;
  n1 = (n1 % 2 == 0) ? n1 / 2 : (n1 + 1) / 2;
  unsigned numberOfQuadraturePoints = (n1 <= 16) ? n1 : 16;
  if (n1 > 16) {
    std::cout <<
              "------------------------------- WARNING: less quadrature points than needed were employed in function AssembleSysSG -------------------------------"
              << std::endl;
    std::cout << " Needed : " << n1 << " , " << " Used : " << 16 << std::endl;
  }

  const std::vector < std::vector < double > >  &hermitePoly = myuq.GetHermitePolynomial (numberOfQuadraturePoints, maxPolyOrder);

  //solution Index
  std::vector <unsigned> soluIndex (Jp.size());
  for (unsigned i = 0; i < Jp.size(); i++) {
    char name[10];
    sprintf (name, "uSG%d", i);
    soluIndex[i] = mlSol->GetIndex (name);   // get the position of "u" in the ml_sol object
  }
  unsigned soluType = mlSol->GetSolutionType (soluIndex[0]);


  //solution PdeIndex
  std::vector <unsigned> soluPdeIndex (Jp.size());
  for (unsigned i = 0; i < Jp.size(); i++) {
    char name[10];
    sprintf (name, "uSG%d", i);
    soluPdeIndex[i] = mlPdeSys->GetSolPdeIndex (name);   // get the position of "u" in the pdeSys object
  }


  //eigenfunction Index
  std::vector <unsigned> eigfIndex (numberOfEigPairs);
  for (unsigned i = 0; i < numberOfEigPairs; i++) {
    char name[10];
    sprintf (name, "egnf%d", i);
    eigfIndex[i] = mlSol->GetIndex (name);   // get the position of "u" in the ml_sol object
  }

  vector < double > KLexpansion; // local solution

  vector < vector < double > >  solu (Jp.size()); // local solution

  vector < vector < double > > x (dim);   // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx;
  double weight; // gauss point weight

  vector< int > l2GMap; // local to global mapping
  vector< double > Res; // local redidual vector
  vector < double > Jac;

  KK->zero(); // Set to zero all the entries of the Global Matrix

  //BEGIN terms for the coefficient obtained projecting the KL (computation later)
//   std::vector < std::vector < double > > productTerms(Jq.size());
//   for(unsigned q1 = 0; q1 < Jq.size(); q1 ++) {
//
//     productTerms[q1].resize(numberOfEigPairs);
//
//     unsigned numberOfQuadraturePointsForProjection = ((qIndex + 2) % 2 == 0) ? ((qIndex + 2) / 2) : ((qIndex + 3) / 2) ;
//     std::vector < std::vector < double > >  HermitePolyProjection;
//     ComputeHermitePoly(HermitePolyProjection, numberOfQuadraturePointsForProjection, qIndex + 1);
//     for(unsigned i = 0; i < numberOfEigPairs; i++) {
//       double termWithY = 0.;
//       for(unsigned j = 0; j < numberOfQuadraturePointsForProjection; j++) {
//         termWithY +=  HermiteQuadrature[numberOfQuadraturePointsForProjection - 1][1][j]
//                       * HermitePolyProjection[Jq[q1][i]][j] * HermiteQuadrature[numberOfQuadraturePointsForProjection - 1][0][j];
//       }
//
//       productTerms[q1][i] = termWithY;
// //       std::cout << " termWithY = "  << termWithY << " ";
//
//       for(unsigned j = 0; j < numberOfEigPairs; j++) {
//         if(j != i) {
//           productTerms[q1][i] *= integralMatrix[Jq[q1][j]][0][0];
//         }
//       }
// //       std::cout << " productTerms[" << q1 << "][" << i << "] = " << productTerms[q1][i] << " ";
//     }
// //         std::cout << "------------------" <<  std::endl;
//   }
  //END terms for coefficient obtained projecting the KL (computation later)

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofu  = msh->GetElementDofNumber (iel, soluType);   // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber (iel, xType);   // number of coordinate element dofs

    // resize local arrays
    for (unsigned j = 0; j < Jp.size(); j++) {
      solu[j].resize (nDofu);
    }

    for (int i = 0; i < dim; i++) {
      x[i].resize (nDofx);
    }

    l2GMap.resize (nDofu * Jp.size());
    Res.assign (nDofu *  Jp.size(), 0.);
    Jac.assign (nDofu *  Jp.size() * nDofu * Jp.size(), 0.);
    KLexpansion.resize (nDofu);

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, soluType);   // global to global mapping between solution node and solution dof
      for (unsigned j = 0; j < Jp.size(); j++) {
        solu[j][i] = (*sol->_Sol[soluIndex[j]]) (solDof);     // global extraction and local storage for the solution
        l2GMap[j * nDofu + i] = pdeSys->GetSystemDof (soluIndex[j], soluPdeIndex[j], i, iel);   // global to global mapping between solution node and pdeSys dof
      }
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim]) (xDof);     // global extraction and local storage for the element coordinates
      }
    }

    std::vector <double> eigVectorGauss (numberOfEigPairs);

    //  *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielGeom][soluType]->Jacobian (x, ig, weight, phi, phi_x, phi_xx);

      for (unsigned i = 0; i < numberOfEigPairs; i++) {
        unsigned solDof = msh->GetSolutionDof (i, iel, soluType);
        eigVectorGauss[i] = 0.;
        for (unsigned j = 0; j < nDofu; j++) {
          eigVectorGauss[i] += (*sol->_Sol[eigfIndex[i]]) (solDof) * phi[i];
        }
      }


      vector< double > aStochasticGauss (Jq.size());

//BEGIN coefficient obtained projecting the exponential of the KL
      for (unsigned q1 = 0; q1 < Jq.size(); q1 ++) {
//         std::vector <double> aStochasticTerm1(numberOfEigPairs);
        std::vector <double> aStochasticTerm2 (numberOfEigPairs);

        unsigned numberOfQuadraturePointsForProjection = 16;

        const double *hermiteQuadraturePoints = myuq.GetHermiteQuadraturePoints (numberOfQuadraturePointsForProjection);
        const double *hermiteQuadratureWeights = myuq.GetHermiteQuadratureWeights (numberOfQuadraturePointsForProjection);

        const std::vector < std::vector < double > >  &hermitePolyProjection = myuq.GetHermitePolynomial (numberOfQuadraturePointsForProjection, qIndex);

        for (unsigned i = 0; i < numberOfEigPairs; i++) {
          aStochasticTerm2[i] = 0.;
          for (unsigned j = 0; j < numberOfQuadraturePointsForProjection; j++) {
//             aStochasticTerm1[i] += hermitePoly[Jq[q1][i]][j] * hermiteQuadrature[numberOfQuadraturePoints - 1][0][j];
            aStochasticTerm2[i] += exp (sqrt (eigenvalues[i].first) * eigVectorGauss[i] * hermiteQuadraturePoints[j])
                                   * hermitePolyProjection[Jq[q1][i]][j] * hermiteQuadratureWeights[j];
          }
        }

        double aS1 = 1.;
        double aS2 = 1.;
        for (unsigned i = 0; i < numberOfEigPairs; i++) {
//           std::cout << "------------------- " << Jq[q1][i] << " ";
          aS1 *= integralMatrix[Jq[q1][i]][0][0];
          aS2 *= aStochasticTerm2[i];
        }
//         std::cout << " stochastic term 1= " << aS1 << " " << "stochastic term 2= " << aS2 << std::endl;

        aStochasticGauss[q1] = amin * aS1 + aS2; //a_q(x_ig)
        if (fabs (aStochasticGauss[q1]) > 10.) std::cout << " coeff =  " << aStochasticGauss[q1] << std::endl;
      }
//END coefficient obtained projecting the exponential of the KL

//BEGIN coefficient obtained projecting the  KL
//       for(unsigned q1 = 0; q1 < Jq.size(); q1 ++) {
//         double aS1 = meanInput;
//         double aS2 = 0.;
//         for(unsigned i = 0; i < numberOfEigPairs; i++) {
//           aS1 *= integralMatrix[Jq[q1][i]][0][0];
//           aS2 += sqrt(eigenvalues[i].first) * eigVectorGauss[i] * productTerms[q1][i];
//         }
//
//         aStochasticGauss[q1] = aS1 + aS2;
// //         std::cout << " coeff =  " << aStochasticGauss[q1] << std::endl;
//
//       }
//END coefficient obtained projecting the  KL


      vector < vector < double > > laplace (nDofu);
      for (unsigned i = 0; i < nDofu; i++) {
        laplace[i].assign (nDofu, 0.);
        for (unsigned j = 0; j < nDofu; j++) {
          for (unsigned kdim = 0; kdim < dim; kdim++) {
            laplace[i][j] += (phi_x[i * dim + kdim] * phi_x[j * dim + kdim]) * weight;
          }
        }
      }


      for (unsigned p1 = 0; p1 < Jp.size(); p1++) {

        double srcTermStoch = 1.;
        for (unsigned i = 0; i < numberOfEigPairs; i++) {
          srcTermStoch *= integralMatrix[0][Jp[p1][i]][0];
        }

        for (unsigned i = 0; i < nDofu; i++) {
          double resU = 1. * phi[i] * srcTermStoch * weight;
          for (unsigned p2 = 0; p2 < Jp.size(); p2++) {
            for (unsigned j = 0; j < nDofu; j++) {
              double AG = 0;
              for (unsigned q1 = 0; q1 < Jq.size(); q1++) {
                AG += aStochasticGauss[q1] * laplace[i][j] * G[q1][p1][p2];
              }
              Jac[ (p1 * nDofu + i) * (Jp.size() * nDofu) +  p2 * nDofu + j] -= AG;
              resU +=  AG * solu[p2][j];
            }
          }
          Res[ p1 * nDofu + i ] += resU;
        }
      }
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    RES->add_vector_blocked (Res, l2GMap);

    //store K in the global matrix KK
    KK->add_matrix_blocked (Jac, l2GMap, l2GMap);

  } //end element loop for each process

  RES->close();

  KK->close();

//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//   PetscObjectSetName((PetscObject)viewer, "SG matrix");
//   PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(), viewer);
//   double a;
//   std::cin >> a;
//  abort();
// ***************** END ASSEMBLY *******************
}







