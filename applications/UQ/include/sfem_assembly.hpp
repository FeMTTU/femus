#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace femus;

// double GetExactSolutionValue(const std::vector < double >& x) {
//   double pi = acos(-1.);
//   return cos(pi * x[0]) * cos(pi * x[1]);
// };
//
// void GetExactSolutionGradient(const std::vector < double >& x, vector < double >& solGrad) {
//   double pi = acos(-1.);
//   solGrad[0]  = -pi * sin(pi * x[0]) * cos(pi * x[1]);
//   solGrad[1] = -pi * cos(pi * x[0]) * sin(pi * x[1]);
// };


int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


const double HermiteQuadrature[15][2][15] = { //Number of quadrature points, first row: weights, second row: coordinates
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
  },
  {
    {
      1.4395603937142598e-6, 0.0003468194663233449, 0.011911395444911533, 0.11722787516770858, 0.42935975235612506, 0.6547592869145917, 0.42935975235612506,
      0.11722787516770858, 0.011911395444911533, 0.0003468194663233449, 1.4395603937142598e-6
    },
    {
      -3.6684708465595826, -2.783290099781652, -2.0259480158257555, -1.3265570844949328, -0.6568095668820998, 0, 0.6568095668820998, 1.3265570844949328, 2.0259480158257555,
      2.783290099781652, 3.6684708465595826
    }
  },
  {
    {
      2.658551684356303e-7, 0.00008573687043587891, 0.003905390584629065, 0.05160798561588393, 0.260492310264161, 0.5701352362624797, 0.5701352362624797,
      0.260492310264161, 0.05160798561588393, 0.003905390584629065, 0.00008573687043587891, 2.658551684356303e-7
    },
    {
      -3.889724897869782, -3.0206370251208896, -2.2795070805010598, -1.5976826351526048, -0.9477883912401638, -0.31424037625435913, 0.31424037625435913,
      0.9477883912401638, 1.5976826351526048, 2.2795070805010598, 3.0206370251208896, 3.889724897869782
    }
  },
  {
    {
      4.825731850073113e-8, 0.000020430360402707097, 0.0012074599927193834, 0.02086277529616993, 0.14032332068702352, 0.4216162968985431, 0.6043931879211616,
      0.4216162968985431, 0.14032332068702352, 0.02086277529616993, 0.0012074599927193834, 0.000020430360402707097, 4.825731850073113e-8
    },
    {
      -4.10133759617864, -3.24660897837241, -2.519735685678238, -1.8531076516015121, -1.2200550365907483, -0.6057638791710601, 0, 0.6057638791710601, 1.2200550365907483,
      1.8531076516015121, 2.519735685678238, 3.24660897837241, 4.10133759617864
    }
  },
  {
    {
      8.628591168125167e-9, 4.7164843550189165e-6, 0.00035509261355192383, 0.007850054726457932, 0.06850553422346523, 0.2731056090642466, 0.5364059097120905,
      0.5364059097120905, 0.2731056090642466, 0.06850553422346523, 0.007850054726457932, 0.00035509261355192383, 4.7164843550189165e-6, 8.628591168125167e-9
    },
    {
      -4.304448570473632, -3.4626569336022706, -2.7484707249854026, -2.095183258507717, -1.4766827311411408, -0.8787137873293994, -0.2917455106725621, 0.2917455106725621,
      0.8787137873293994, 1.4766827311411408, 2.095183258507717, 2.7484707249854026, 3.4626569336022706, 4.304448570473632
    }
  },
  {
    {
      1.5224758042535108e-9, 1.059115547711064e-6, 0.00010000444123249993, 0.0027780688429127837, 0.030780033872546103, 0.1584889157959358, 0.41202868749889865,
      0.5641003087264175, 0.41202868749889865, 0.1584889157959358, 0.030780033872546103, 0.0027780688429127837, 0.00010000444123249993, 1.059115547711064e-6,
      1.5224758042535108e-9
    },
    {
      -4.499990707309392, -3.6699503734044527, -2.967166927905603, -2.3257324861738575, -1.7199925751864888, -1.1361155852109206, -0.5650695832555758, 0,
      0.5650695832555758, 1.1361155852109206, 1.7199925751864888, 2.3257324861738575, 2.967166927905603, 3.6699503734044527, 4.499990707309392
    }
  }
};




int numberOfEigPairs = 1; //dimension of the stochastic variable
std::vector < std::pair<double, double> > eigenvalues(numberOfEigPairs);

double amin = 1. / 100;

double stdDeviationInput = 1./*0.2*/;  //standard deviation of the normal distribution (it is the same as the standard deviation of the covariance function in GetEigenPair)

boost::mt19937 rng; // I don't seed it on purpouse (it's not relevant)

boost::normal_distribution<> nd(0.0, stdDeviationInput);

boost::variate_generator < boost::mt19937&,

      boost::normal_distribution<> > var_nor(rng, nd);

void EvaluateHermitePoly(std::vector < std::vector < double > >  & HermitePoly, const unsigned & numberOfQuadraturePoints, const unsigned & maxPolyOrder)
{

  if(numberOfQuadraturePoints < 1 || numberOfQuadraturePoints > 15) {
    std::cout << "The selected order of integraiton has not been implemented yet, choose an integer in [1,15]" << std::endl;
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

double GetExactSolutionLaplace(const std::vector < double >& x)
{
  double pi = acos(-1.);
  return -pi * pi * cos(pi * x[0]) * cos(pi * x[1]) - pi * pi * cos(pi * x[0]) * cos(pi * x[1]);
};


void ComputeIndexSetJp(std::vector < std::vector <unsigned> > & Jp, const unsigned & p)   //p is max poly degree
{

  unsigned dimJp = factorial(numberOfEigPairs + p) / (factorial(numberOfEigPairs) * factorial(p));
//   std::cout << dimJp <<std::endl;
  Jp.resize(dimJp);
  for(unsigned i = 0; i < dimJp; i++) {
    Jp[i].resize(numberOfEigPairs);
  }

  unsigned index = 0;
  unsigned counters[numberOfEigPairs + 1];
  memset(counters, 0, sizeof(counters));

  while(!counters[numberOfEigPairs]) {

//     for(unsigned i = numberOfEigPairs; i-- > 0;) {
//       std::cout << counters[i] <<" ";
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

void EvaluateStochasticMassMatrices(const unsigned & q0, const unsigned & p0, std::vector < std::vector < std::vector < double > > > & G)
{
  
  unsigned maxPolyOrder = (q0 > p0) ? q0 : p0;
  
  std::vector < std::vector < double > >  HermitePoly;

  unsigned n1 = 2 * p0 + q0 + 1;
  n1 = (n1 % 2 == 0 ) ? n1 / 2 : (n1 + 1) / 2;
  unsigned numberOfQuadraturePoints = (n1 <= 15) ? n1 : 15;

  EvaluateHermitePoly(HermitePoly, numberOfQuadraturePoints, maxPolyOrder);

  unsigned q = q0 + 1;
  unsigned p = p0 + 1;
  
  std::vector < std::vector < std::vector < double > > > integralsMatrix;
  integralsMatrix.resize(q);
  for(unsigned q1 = 0; q1 < q; q1++) {
    integralsMatrix[q1].resize(p);
    for(unsigned p1 = 0; p1 < p; p1++) {
      integralsMatrix[q1][p1].assign(p,0.);
      for(unsigned p2 = 0; p2 < p; p2++) {
        integralsMatrix[q1][p1][p2]  = 0.;
        for(unsigned i = 0; i < numberOfQuadraturePoints; i++) {
          double w = HermiteQuadrature[numberOfQuadraturePoints - 1][0][i];
          integralsMatrix[q1][p1][p2]  +=  w * HermitePoly[q1][i] * HermitePoly[p1][i] * HermitePoly[p2][i];
        }
      }
    }
  }
  
  std::vector < std::vector <unsigned> > Jq;
  std::vector < std::vector <unsigned> > Jp;

  ComputeIndexSetJp(Jq, q0);
  ComputeIndexSetJp(Jp, p0);

  G.resize(Jq.size());
  for(unsigned q1 = 0; q1 < Jq.size(); q1++) {
    G[q1].resize(Jp.size());
    for(unsigned p1 = 0; p1 < Jp.size(); p1++) {
      G[q1][p1].assign(Jp.size(), 1.);
      for(unsigned p2 = 0; p2 < Jp.size(); p2++) {
        for(unsigned i = 0; i < numberOfEigPairs; i++) {
          G[q1][p1][p2] *= integralsMatrix[Jq[q1][i]][Jp[p1][i]][Jp[p2][i]];
        }
        G[q1][p1][p2] = ( fabs(G[q1][p1][p2]) < 1.e-14 ) ? 0. : G[q1][p1][p2];
      }
    }
  }

};

void AssembleUQSys(MultiLevelProblem& ml_prob)
{
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object


  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("UQ");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  char name[10];
  std::vector <unsigned> eigfIndex(numberOfEigPairs);

  for(unsigned i = 0; i < numberOfEigPairs; i++) {
    sprintf(name, "egnf%d", i);
    eigfIndex[i] = mlSol->GetIndex(name);    // get the position of "u" in the ml_sol object
  }


  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solu; // local solution
  solu.reserve(maxSize);


  vector < double > KLexpansion; // local solution
  KLexpansion.reserve(maxSize);


  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for(unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  phi_xx.reserve(maxSize * dim2);

  vector< adept::adouble > aRes; // local redidual vector
  aRes.reserve(maxSize);

  vector< int > l2GMap; // local to global mapping
  l2GMap.reserve(maxSize);
  vector< double > Res; // local redidual vector
  Res.reserve(maxSize);
  vector < double > Jac;
  Jac.reserve(maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  std::vector <double> yOmega(numberOfEigPairs, 0.);
  for(unsigned eig = 0; eig < numberOfEigPairs; eig++) {
    if(iproc == 0) {
      yOmega[eig] = var_nor();
    }
    MPI_Bcast(&yOmega[eig], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    std::cout << yOmega[eig] << " ";
  }
  std::cout << std::endl;






  // element loop: each process loops only on the elements that owns
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    l2GMap.resize(nDofu);
    solu.resize(nDofu);
    KLexpansion.resize(nDofu);

    for(int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    aRes.resize(nDofu);    //resize
    std::fill(aRes.begin(), aRes.end(), 0);    //set aRes to zero

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      KLexpansion[i] = 0.;
      for(unsigned j = 0; j < numberOfEigPairs; j++) {
        KLexpansion[i] += sqrt(eigenvalues[j].first) * (*sol->_Sol[eigfIndex[j]])(solDof) * yOmega[j];
      }
      l2GMap[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for(unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }


    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point

      double KLexpansion_gss = 0.;
      vector < adept::adouble > gradSolu_gss(dim, 0.);
      vector < double > x_gss(dim, 0.);

      for(unsigned i = 0; i < nDofu; i++) {

        KLexpansion_gss += phi[i] * KLexpansion[i];

        for(unsigned jdim = 0; jdim < dim; jdim++) {
          gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
          x_gss[jdim] += x[jdim][i] * phi[i];
        }
      }

      double aCoeff = amin + exp(KLexpansion_gss);
//       std::cout << "COEEEEEEEEEEEEEEEEEEEEF" << aCoeff << std::endl;

      // *** phi_i loop ***
      for(unsigned i = 0; i < nDofu; i++) {

        adept::adouble laplace = 0.;

        for(unsigned jdim = 0; jdim < dim; jdim++) {
          laplace   +=  aCoeff * phi_x[i * dim + jdim] * gradSolu_gss[jdim];
        }

        double srcTerm = 1./*- GetExactSolutionLaplace(x_gss)*/ ;
        aRes[i] += (srcTerm * phi[i] + laplace) * weight;

      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    Res.resize(nDofu);    //resize

    for(int i = 0; i < nDofu; i++) {
      Res[i] = - aRes[i].value();
    }

    RES->add_vector_blocked(Res, l2GMap);



    // define the dependent variables
    s.dependent(&aRes[0], nDofu);

    // define the independent variables
    s.independent(&solu[0], nDofu);

    // get the jacobian matrix (ordered by row major )
    Jac.resize(nDofu * nDofu);    //resize
    s.jacobian(&Jac[0], true);

    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac, l2GMap, l2GMap);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}



