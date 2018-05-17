#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace femus;


unsigned factorial(unsigned n) {

  unsigned fac = 1.;

  if(n >= 0 && n < 13) {

    if(n > 0) {
      while(n > 1) {
        fac *= n;
        n--;
      }
    }
  }

  else {
    std::cout << " input of factorial function must be an integer that is non-negative or less than 13" << std::endl;
    abort();
  }

  return fac;

}


//BEGIN Stochastic Input Parameters
unsigned pIndex = 3;
unsigned qIndex = 2;

int numberOfEigPairs = 2; //dimension of the stochastic variable
double stdDeviationInput = 1.;  //standard deviation of the normal distribution (it is the same as the standard deviation of the covariance function in GetEigenPair)
double amin = 1. / 100; // for the KL expansion
std::vector < std::pair<double, double> > eigenvalues(numberOfEigPairs);
//END Stochastic Input Parameters


 const double HermiteQuadrature[10][2][10] = { //Number of quadrature points, first row: weights, second row: coordinates
  {{1.77245385091}, {0.}},
  { {0.8862269254527577, 0.8862269254527577},
    { -0.7071067811865476, 0.7071067811865476}
  },
  { {0.29540897515091946, 1.1816359006036772, 0.29540897515091946},
    { -1.224744871391589, 0., 1.224744871391589}
  },
  { {0.08131283544724523, 0.8049140900055128, 0.8049140900055128, 0.08131283544724523},
    { -1.6506801238857844, -0.5246476232752904, 0.5246476232752904, 1.6506801238857844}
  },
  { {0.019953242059045917, 0.39361932315224113, 0.9453087204829418, 0.39361932315224113, 0.019953242059045917},
    { -2.0201828704560856, -0.9585724646138185, 0., 0.9585724646138185, 2.0201828704560856}
  },
  { {0.004530009905508841, 0.15706732032285653, 0.7246295952243923, 0.7246295952243923, 0.15706732032285653, 0.004530009905508841},
    { -2.3506049736744923, -1.335849074013697, -0.4360774119276165, 0.4360774119276165, 1.335849074013697, 2.3506049736744923}
  },
  { {
      0.0009717812450995206, 0.05451558281912702, 0.4256072526101277, 0.8102646175568072, 0.4256072526101277, 0.05451558281912702,
      0.0009717812450995206
    },
    { -2.6519613568352334, -1.6735516287674714, -0.8162878828589647, 0., 0.8162878828589647, 1.6735516287674714, 2.6519613568352334}
  },
  { {
      0.00019960407221136756, 0.017077983007413464, 0.20780232581489194, 0.6611470125582413, 0.6611470125582413, 0.20780232581489194,
      0.017077983007413464, 0.00019960407221136756
    },
    {
      -2.930637420257244, -1.981656756695843, -1.1571937124467802, -0.3811869902073221, 0.3811869902073221, 1.1571937124467802,
      1.981656756695843, 2.930637420257244
    }
  },
  { {
      0.00003960697726326427, 0.004943624275536942, 0.08847452739437658, 0.4326515590025556, 0.7202352156060509, 0.4326515590025556,
      0.08847452739437658, 0.004943624275536942, 0.00003960697726326427
    },
    {
      -3.1909932017815277, -2.266580584531843, -1.468553289216668, -0.7235510187528376, 0., 0.7235510187528376, 1.468553289216668,
      2.266580584531843, 3.1909932017815277
    }
  },
  { {
      7.640432855232643e-6, 0.001343645746781235, 0.033874394455481044, 0.24013861108231477, 0.6108626337353258, 0.6108626337353258,
      0.24013861108231477, 0.033874394455481044, 0.001343645746781235, 7.640432855232643e-6
    },
    {
      -3.4361591188377374, -2.5327316742327897, -1.7566836492998819, -1.0366108297895136, -0.3429013272237046, 0.3429013272237046,
      1.0366108297895136, 1.7566836492998819, 2.5327316742327897, 3.4361591188377374
    }
  }
};

void MultidimensionalHermiteQuadrature(std::vector < std::vector <unsigned> > & MultidimHermitePoints, const unsigned & numberOfQuadraturePoints,
                                       const unsigned & numberOfEigPairs) {

  unsigned tensorProductDim = pow(numberOfQuadraturePoints + 1, numberOfEigPairs);
  
  MultidimHermitePoints.resize(tensorProductDim);
  for(unsigned i = 0; i < tensorProductDim; i++) {
    MultidimHermitePoints[i].assign(1 + numberOfEigPairs, 1.); //first entry is weight, second is coordinate
  }

  unsigned index = 0;
  unsigned counters[numberOfEigPairs + 1];
  memset(counters, 0, sizeof(counters));

  while(!counters[numberOfEigPairs]) {

    for(unsigned j = 0; j < numberOfEigPairs; j++) {
      unsigned indexx = counters[numberOfEigPairs - 1 - j];
      MultidimHermitePoints[index][0] *= HermiteQuadrature[numberOfQuadraturePoints][0][indexx];
      MultidimHermitePoints[index][j+1] =  HermiteQuadrature[numberOfQuadraturePoints][j+1][indexx];
    }
    index++;

     unsigned i;
    for(i = 0; counters[i] == numberOfQuadraturePoints; i++) { // inner loops that are at maxval restart at zero
      counters[i] = 0;
    }
    ++counters[i];  // the innermost loop that isn't yet at maxval, advances by 1
    
  }
};

void EvaluateHermitePoly(std::vector < std::vector < double > >  & HermitePoly, const unsigned & numberOfQuadraturePoints, const unsigned & maxPolyOrder) {

  if(numberOfQuadraturePoints < 1 || numberOfQuadraturePoints > 10) {
    std::cout << "The selected order of integraiton has not been implemented yet, choose an integer in [1,10]" << std::endl;
    abort();
  }

  else {
    HermitePoly.resize(maxPolyOrder + 1);
    for(unsigned i = 0; i < maxPolyOrder + 1; i++) {
      HermitePoly[i].resize(numberOfQuadraturePoints);
    }

    for(unsigned j = 0; j < numberOfQuadraturePoints; j++) {

      double x = HermiteQuadrature[numberOfQuadraturePoints - 1][1][j];

      HermitePoly[0][j] = 1. ;

      if(maxPolyOrder > 0) {
        HermitePoly[1][j] = x ;
        if(maxPolyOrder > 1) {
          HermitePoly[2][j] = (pow(x, 2) - 1.) / sqrt(2) ;
          if(maxPolyOrder > 2) {
            HermitePoly[3][j] = (pow(x, 3) - 3. * x) / sqrt(6) ;
            if(maxPolyOrder > 3) {
              HermitePoly[4][j] = (pow(x, 4) - 6. * x * x + 3.) / sqrt(24) ;
              if(maxPolyOrder > 4) {
                HermitePoly[5][j] = (pow(x, 5) - 10. * pow(x, 3) + 15. * x) / sqrt(120) ;
                if(maxPolyOrder > 5) {
                  HermitePoly[6][j] = (pow(x, 6) - 15. * pow(x, 4) + 45. * pow(x, 2) - 15.) / sqrt(720) ;
                  if(maxPolyOrder > 6) {
                    HermitePoly[7][j] = (pow(x, 7) - 21. * pow(x, 5) + 105. * pow(x, 3) -  105. * x) / sqrt(5040) ;
                    if(maxPolyOrder > 7) {
                      HermitePoly[8][j] = (pow(x, 8) - 28. * pow(x, 6) + 210. * pow(x, 4) - 420. * pow(x, 4) + 105.) / sqrt(40320) ;
                      if(maxPolyOrder > 8) {
                        HermitePoly[9][j] = (pow(x, 9) - 36. * pow(x, 7) + 378. * pow(x, 5) - 1260. * pow(x, 3) + 945. * x) / sqrt(362880);
                        if(maxPolyOrder > 9) {
                          HermitePoly[10][j] = (pow(x, 10) - 45. * pow(x, 8) + 630. * pow(x, 6) - 3150. * pow(x, 4) + 4725. * pow(x, 2) - 945.) / sqrt(3628800);
                          if(maxPolyOrder > 10) {
                            std::cout << "Polynomial order is too big. For now, it has to be not greater than 10." << std::endl;
                            abort();
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }

    }

  }

};

void ComputeIndexSetJp(std::vector < std::vector <unsigned> > & Jp, const unsigned & p, const unsigned & numberOfEigPairs) { //p is max poly degree

  unsigned dimJp = factorial(numberOfEigPairs + p) / (factorial(numberOfEigPairs) * factorial(p));

  Jp.resize(dimJp);
  for(unsigned i = 0; i < dimJp; i++) {
    Jp[i].resize(numberOfEigPairs);
  }

  unsigned index = 0;
  unsigned counters[numberOfEigPairs + 1];
  memset(counters, 0, sizeof(counters));

  while(!counters[numberOfEigPairs]) {

//     for(unsigned i = numberOfEigPairs; i-- > 0;) {
//       std::cout << counters[i] << " ";
//     }
//     std::cout << std::endl;

    unsigned entrySum = 0;
    for(unsigned j = 0; j < numberOfEigPairs; j++) {
      entrySum += counters[j];
    }

    if(entrySum <= p) {
      for(unsigned j = 0; j < numberOfEigPairs; j++) {
        Jp[index][j] = counters[numberOfEigPairs - 1 - j];
      }
      index++;
    }
    unsigned i;
    for(i = 0; counters[i] == p; i++) { // inner loops that are at maxval restart at zero
      counters[i] = 0;
    }
    ++counters[i];  // the innermost loop that isn't yet at maxval, advances by 1
  }
};

void EvaluateIntegralsMatrix(const unsigned & q0, const unsigned & p0, std::vector < std::vector < std::vector < double > > > &integralsMatrix) {

  unsigned maxPolyOrder = (q0 > p0) ? q0 : p0;

  std::vector < std::vector < double > >  HermitePoly;

  unsigned n1 = 2 * p0 + q0 + 1;
  n1 = (n1 % 2 == 0) ? n1 / 2 : (n1 + 1) / 2;
  unsigned numberOfQuadraturePoints = (n1 <= 10) ? n1 : 10;

  EvaluateHermitePoly(HermitePoly, numberOfQuadraturePoints, maxPolyOrder);

  unsigned q = q0 + 1;
  unsigned p = p0 + 1;

  integralsMatrix.resize(q);
  for(unsigned q1 = 0; q1 < q; q1++) {
    integralsMatrix[q1].resize(p);
    for(unsigned p1 = 0; p1 < p; p1++) {
      integralsMatrix[q1][p1].assign(p, 0.);
      for(unsigned p2 = 0; p2 < p; p2++) {
        integralsMatrix[q1][p1][p2]  = 0.;
        for(unsigned i = 0; i < numberOfQuadraturePoints; i++) {
          double w = HermiteQuadrature[numberOfQuadraturePoints - 1][0][i];
          integralsMatrix[q1][p1][p2]  +=  w * HermitePoly[q1][i] * HermitePoly[p1][i] * HermitePoly[p2][i];
        }
      }
    }
  }

};

void EvaluateStochasticMassMatrices(const unsigned & q0, const unsigned & p0, std::vector < std::vector < std::vector < double > > > & G,
                                    const unsigned & numberOfEigPairs, const std::vector < std::vector < std::vector < double > > > & integralsMatrix) {

  std::vector < std::vector <unsigned> > Jq;
  std::vector < std::vector <unsigned> > Jp;
  ComputeIndexSetJp(Jq, q0, numberOfEigPairs);
  ComputeIndexSetJp(Jp, p0, numberOfEigPairs);

  G.resize(Jq.size());
  for(unsigned q1 = 0; q1 < Jq.size(); q1++) {
    G[q1].resize(Jp.size());
    for(unsigned p1 = 0; p1 < Jp.size(); p1++) {
      G[q1][p1].assign(Jp.size(), 1.);
      for(unsigned p2 = 0; p2 < Jp.size(); p2++) {
        for(unsigned i = 0; i < numberOfEigPairs; i++) {
          G[q1][p1][p2] *= integralsMatrix[Jq[q1][i]][Jp[p1][i]][Jp[p2][i]];
        }
        G[q1][p1][p2] = (fabs(G[q1][p1][p2]) < 1.e-14) ? 0. : G[q1][p1][p2];
      }
    }
  }

};


void AssembleSysSG(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("SG");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix* KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  std::vector < std::vector < std::vector < double > > > integralsMatrix;
  EvaluateIntegralsMatrix(qIndex, pIndex, integralsMatrix);

  std::vector < std::vector < std::vector < double > > >  G; //vector with stochastic mass matrices
  EvaluateStochasticMassMatrices(qIndex, pIndex, G, numberOfEigPairs, integralsMatrix);

  std::vector < std::vector <unsigned> > Jq;
  ComputeIndexSetJp(Jq, qIndex, numberOfEigPairs);
  std::vector < std::vector <unsigned> > Jp;
  ComputeIndexSetJp(Jp, pIndex, numberOfEigPairs);


  unsigned maxPolyOrder = (qIndex > pIndex) ? qIndex : pIndex;

  std::vector < std::vector < double > >  HermitePoly;

  unsigned n1 = 2 * pIndex + qIndex + 1;
  n1 = (n1 % 2 == 0) ? n1 / 2 : (n1 + 1) / 2;
  unsigned numberOfQuadraturePoints = (n1 <= 10) ? n1 : 10;

  EvaluateHermitePoly(HermitePoly, numberOfQuadraturePoints, maxPolyOrder);

  //solution Index
  std::vector <unsigned> soluIndex(Jp.size());
  for(unsigned i = 0; i < Jp.size(); i++) {
    char name[10];
    sprintf(name, "uSG%d", i);
    soluIndex[i] = mlSol->GetIndex(name);    // get the position of "u" in the ml_sol object
  }
  unsigned soluType = mlSol->GetSolutionType(soluIndex[0]);


  //solution PdeIndex
  std::vector <unsigned> soluPdeIndex(Jp.size());
  for(unsigned i = 0; i < Jp.size(); i++) {
    char name[10];
    sprintf(name, "uSG%d", i);
    soluPdeIndex[i] = mlPdeSys->GetSolPdeIndex(name);    // get the position of "u" in the pdeSys object
  }


  //eigenfunction Index
  std::vector <unsigned> eigfIndex(numberOfEigPairs);
  for(unsigned i = 0; i < numberOfEigPairs; i++) {
    char name[10];
    sprintf(name, "egnf%d", i);
    eigfIndex[i] = mlSol->GetIndex(name);    // get the position of "u" in the ml_sol object
  }

  vector < double > KLexpansion; // local solution

  vector < vector < double > >  solu(Jp.size()); // local solution

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx;
  double weight; // gauss point weight

  vector< int > l2GMap; // local to global mapping
  vector< double > Res; // local redidual vector
  vector < double > Jac;

  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    for(unsigned j = 0; j < Jp.size(); j++) {
      solu[j].resize(nDofu);
    }

    for(int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    l2GMap.resize(nDofu * Jp.size());
    Res.assign(nDofu *  Jp.size(), 0.);
    Jac.assign(nDofu *  Jp.size() * nDofu * Jp.size(), 0.);
    KLexpansion.resize(nDofu);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      for(unsigned j = 0; j < Jp.size(); j++) {
        solu[j][i] = (*sol->_Sol[soluIndex[j]])(solDof);      // global extraction and local storage for the solution
        l2GMap[j * nDofu + i] = pdeSys->GetSystemDof(soluIndex[j], soluPdeIndex[j], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for(unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    std::vector <double> eigVectorGauss(numberOfEigPairs);

    //  *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

      for(unsigned i = 0; i < numberOfEigPairs; i++) {
        unsigned solDof = msh->GetSolutionDof(i, iel, soluType);
        eigVectorGauss[i] = 0.;
        for(unsigned j = 0; j < nDofu; j++) {
          eigVectorGauss[i] += (*sol->_Sol[eigfIndex[i]])(solDof) * phi[i];
        }
      }


      vector< double > aStochasticGauss(Jq.size());
      for(unsigned q1 = 0; q1 < Jq.size(); q1 ++) {
        std::vector <double> aStochasticTerm1(numberOfEigPairs);
        std::vector <double> aStochasticTerm2(numberOfEigPairs);

        for(unsigned i = 0; i < numberOfEigPairs; i++) {
          aStochasticTerm2[i] = 0.;
          for(unsigned j = 0; j < numberOfQuadraturePoints; j++) {
//             aStochasticTerm1[i] += HermitePoly[Jq[q1][i]][j] * HermiteQuadrature[numberOfQuadraturePoints - 1][0][j];
            aStochasticTerm2[i] += exp(sqrt(eigenvalues[i].first) * eigVectorGauss[i] * HermiteQuadrature[numberOfQuadraturePoints - 1][1][j])
                                   * HermitePoly[Jq[q1][i]][j] * HermiteQuadrature[numberOfQuadraturePoints - 1][0][j];
          }
        }


        double aS1 = 1.;
        double aS2 = 1.;
        for(unsigned i = 0; i < numberOfEigPairs; i++) {
          aS1 *= integralsMatrix[Jq[q1][i]][0][0];
          aS2 *= aStochasticTerm2[i];
        }

        aStochasticGauss[q1] = amin * aS1 + aS2; //a_q(x_ig)
        if(fabs(aStochasticGauss[q1]) > 10.) std::cout << aStochasticGauss[q1] << " ";
      }
      // evaluate the solution, the solution derivatives and the coordinates in the gauss point

      vector < vector < double > > laplace(nDofu);
      for(unsigned i = 0; i < nDofu; i++) {
        laplace[i].assign(nDofu, 0.);
        for(unsigned j = 0; j < nDofu; j++) {
          for(unsigned kdim = 0; kdim < dim; kdim++) {
            laplace[i][j] += (phi_x[i * dim + kdim] * phi_x[j * dim + kdim]) * weight;
          }
        }
      }


      for(unsigned p1 = 0; p1 < Jp.size(); p1++) {

        double srcTermStoch = 1.;
        for(unsigned i = 0; i < numberOfEigPairs; i++) {
          srcTermStoch *= integralsMatrix[0][Jp[p1][i]][0];
        }

        for(unsigned i = 0; i < nDofu; i++) {
          double resU = 1. * phi[i] * srcTermStoch * weight;
          for(unsigned p2 = 0; p2 < Jp.size(); p2++) {
            for(unsigned j = 0; j < nDofu; j++) {
              double AG = 0;
              for(unsigned q1 = 0; q1 < Jq.size(); q1++) {
                AG += aStochasticGauss[q1] * laplace[i][j] * G[q1][p1][p2];
              }
              Jac[(p1 * nDofu + i) * (Jp.size() * nDofu) +  p2 * nDofu + j] += AG;
              resU -=  AG * solu[p2][j];
            }
          }
          Res[ p1 * nDofu + i ] += resU;
        }
      }
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    RES->add_vector_blocked(Res, l2GMap);

    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac, l2GMap, l2GMap);

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





