
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/binomial.hpp>

using namespace femus;

class uq {

  public:

    /// Get Hermite quadrature point coordinates
    const double* GetHermiteQuadraturePoints (const unsigned &numberOfQuadraturePoints);

    /// Get Hermite quadrature weights
    const double* GetHermiteQuadratureWeights (const unsigned &numberOfQuadraturePoints);

    ////////////////////////////////////////////

    /// Compute the tensor product at the key < numberOfQuadraturePoints,  numberOfEigPairs >
    void ComputeTensorProductSet (std::vector < std::vector <unsigned> > & Tp,
                                  const unsigned & numberOfQuadraturePoints,
                                  const unsigned & numberOfEigPairs);

    /// Return the tensor product set at the key < numberOfQuadraturePoints , numberOfEigPairs >
    const std::vector < std::vector <unsigned> > & GetTensorProductSet (const unsigned & numberOfQuadraturePoints,
        const unsigned & numberOfEigPairs);

    /// Erase the tensor product set at the key < numberOfQuadraturePoints , numberOfEigPairs >, if stored
    void EraseTensorProductSet (const unsigned & numberOfQuadraturePoints, const unsigned & numberOfEigPairs);

    /// Clear all stored tensor product sets
    void ClearTensorProductSet();

    ////////////////////////////////////////////

    /// Compute the Hermite Polynomial at the key < numberOfQuadraturePoints, maxPolyOrder >
    void ComputeHermitePoly (std::vector < std::vector < double > >  & HermitePoly,
                             const unsigned &numberOfQuadraturePoints, const unsigned & maxPolyOrder);


    /// Return the Hermite polynomial at the key < numberOfQuadraturePoints , maxPolyOrder >
    const std::vector < std::vector <double> > & GetHermitePolynomial (const unsigned & numberOfQuadraturePoints,
                                                                       const unsigned & maxPolyOrder);

    /// Erase the Hermite polynomial at the key < numberOfQuadraturePoints , maxPolyOrder >, if stored
    void EraseHermitePolynomial (const unsigned & numberOfQuadraturePoints, const unsigned & maxPolyOrder);

    /// Clear all stored Hermite polynomials
    void ClearHermitePolynomial();

    /////////////////////////////////////////////////

    void ComputeHermitePolyHistogram (std::vector < std::vector < double > >  & HermitePoly, const unsigned & pIndex,
                                       std::vector<double> & samplePoints, const unsigned & numberOfEigPairs);

    /////////////////////////////////////////////////

    /// Compute the Index Set Jp at the key < p, numberOfEigPairs>
    void ComputeIndexSetJp (std::vector < std::vector <unsigned> > & Jp, const unsigned & p,
                            const unsigned & numberOfEigPairs);

    /// Return the Index Set at the key < p, numberOfEigPairs>
    const std::vector < std::vector <unsigned> > & GetIndexSet (const unsigned & p, const unsigned & numberOfEigPairs);

    /// Erase the Index Set at the key < p, numberOfEigPairs>, if stored
    void EraseIndexSet (const unsigned & numberOfQuadraturePoints, const unsigned & maxPolyOrder);

    /// Clear all Index Sets Jp
    void ClearIndexSet();

    /////////////////////////////////////////////////

    /// Compute the Integral Matrix at the key < q0, p0>
    void ComputeIntegralMatrix (std::vector < std::vector < std::vector < double > > > &integralMatrix,
                                  const unsigned & q0, const unsigned & p0);
    
    /// Get the Integral Matrix at the key < q0, p0>
    const std::vector < std::vector < std::vector < double > > > &GetIntegralMatrix(const unsigned & q0, const unsigned & p0);
    
    /// Erase the Integral Matrix at the key < q0, p0>
    void EraseIntegralMatrix (const unsigned & q0, const unsigned & p0);;
    
    /// Clear all Integral Matrices
    void ClearIntegralMatrix();
    
    /////////////////////////////////////////////////
    
    /// Compute the Integral Matrix at the key < q0, p0, numberOfEigPairs>
    void ComputeStochasticMassMatrix (std::vector < std::vector < std::vector < double > > > & G, 
                                      const unsigned & q0, const unsigned & p0, const unsigned & numberOfEigPairs);
    
    /// Return the Stochastic Mass Matrix at the key < q0, p0, numberOfEigPairs>
    std::vector < std::vector < std::vector < double > > > & GetStochasticMassMatrix (const unsigned & q0, const unsigned & p0, 
                                                                                      const unsigned & numberOfEigPairs);
                                                                        
    
    /// Erase the Stochastic Mass Matrix at the key < q0, p0, numberOfEigPairs>
    void EraseStochasticMassMatrix (const unsigned & q0, const unsigned & p0, 
                                    const unsigned & numberOfEigPairs);
    
    /// Clear all stored Stochastic Mass Matrices 
    void ClearStochasticMassMatrix();

    /////////////////////////////////////////////////
    
    void ComputeMultivariateHermitePoly (std::vector < std::vector < double > >  & multivariateHermitePoly,
                                          std::vector < double > & multivariateHermiteQuadratureWeights,
                                          const unsigned & numberOfQuadraturePoints, const unsigned & p,
                                          const unsigned & numberOfEigPairs);

  private:
    static const double _HermiteQuadrature[16][2][16];
    std::map<std::pair<unsigned, unsigned>, std::vector < std::vector <unsigned> > > _Jp;
    std::map<std::pair<unsigned, unsigned>, std::vector < std::vector <unsigned> > > _Tp;
    std::map<std::pair<unsigned, unsigned>, std::vector < std::vector <double> > > _HermitePoly;
    std::map<std::pair<unsigned, unsigned>, std::vector < std::vector < std::vector < double > > > > _integralMatrix;
    std::map<std::pair < std::pair<unsigned, unsigned>, unsigned >, std::vector < std::vector < std::vector < double > > > > _stochasticMassMatrix;
};


////////////////////////////////////////////////

const double uq::_HermiteQuadrature[16][2][16] = { //Number of quadrature points, first row: weights, second row: coordinates
  {{1.}, {0.}},
  { {0.5000000000000002, 0.5000000000000002},
    { -1., 1.}
  },
  { {0.1666666666666667, 0.6666666666666666, 0.1666666666666667},
    { -1.7320508075688772, 0, 1.7320508075688772}
  },
  { {0.04587585476806852, 0.4541241452319316, 0.4541241452319316, 0.04587585476806852},
    { -2.3344142183389773, -0.7419637843027258, 0.7419637843027258, 2.3344142183389773}
  },
  { {0.011257411327720693, 0.22207592200561266, 0.5333333333333333, 0.22207592200561266, 0.011257411327720693},
    { -2.8569700138728056, -1.355626179974266, 0, 1.355626179974266, 2.8569700138728056}
  },
  { {0.0025557844020562544, 0.08861574604191463, 0.40882846955602936, 0.40882846955602936, 0.08861574604191463, 0.0025557844020562544},
    { -3.324257433552119, -1.8891758777537107, -0.6167065901925941, 0.6167065901925941, 1.8891758777537107, 3.324257433552119}
  },
  { {
      0.0005482688559722186, 0.03075712396758655, 0.2401231786050128, 0.4571428571428571, 0.2401231786050128, 0.03075712396758655,
      0.0005482688559722186
    },
    { -3.7504397177257425, -2.366759410734541, -1.1544053947399682, 0, 1.1544053947399682, 2.366759410734541, 3.7504397177257425}
  },
  { {
      0.00011261453837536777, 0.009635220120788286, 0.11723990766175918, 0.3730122576790776, 0.3730122576790776, 0.11723990766175918,
      0.009635220120788286, 0.00011261453837536777
    },
    {
      -4.1445471861258945, -2.8024858612875416, -1.636519042435108, -0.5390798113513751, 0.5390798113513751, 1.636519042435108,
      2.8024858612875416, 4.1445471861258945
    }
  },
  { {
      0.00002234584400774667, 0.0027891413212317783, 0.04991640676521788, 0.2440975028949394, 0.40634920634920635, 0.2440975028949394,
      0.04991640676521788, 0.0027891413212317783, 0.00002234584400774667
    },
    {
      -4.512745863399783, -3.20542900285647, -2.07684797867783, -1.0232556637891326, 0, 1.0232556637891326, 2.07684797867783,
      3.20542900285647, 4.512745863399783
    }
  },
  { {
      4.3106526307183275e-6, 0.0007580709343122191, 0.019111580500770314, 0.13548370298026785, 0.34464233493201907, 0.34464233493201907,
      0.13548370298026785, 0.019111580500770314, 0.0007580709343122191, 4.3106526307183275e-6
    },
    {
      -4.859462828332312, -3.581823483551927, -2.4843258416389546, -1.4659890943911582, -0.48493570751549764, 0.48493570751549764,
      1.4659890943911582, 2.4843258416389546, 3.581823483551927, 4.859462828332312
    }
  },
  {
    {
      8.121849790214966e-7, 0.00019567193027122425, 0.006720285235537278, 0.06613874607105787, 0.24224029987396994, 0.3694083694083694,
      0.24224029987396994, 0.06613874607105787, 0.006720285235537278, 0.00019567193027122425, 8.121849790214966e-7
    },
    {
      -5.1880012243748705, -3.936166607129977, -2.865123160643645, -1.876035020154846, -0.928868997381064, 0, 0.928868997381064,
      1.876035020154846, 2.865123160643645, 3.936166607129977, 5.1880012243748705
    }
  },
  {
    {
      1.4999271676371695e-7, 0.00004837184922590667, 0.0022033806875332014, 0.02911668791236416, 0.14696704804533003, 0.3216643615128302,
      0.3216643615128302, 0.14696704804533003, 0.02911668791236416, 0.0022033806875332014, 0.00004837184922590667, 1.4999271676371695e-7
    },
    {
      -5.500901704467748, -4.2718258479322815, -3.2237098287700974, -2.2594644510007993, -1.3403751971516167, -0.44440300194413895,
      0.44440300194413895, 1.3403751971516167, 2.2594644510007993, 3.2237098287700974, 4.2718258479322815, 5.500901704467748
    }
  },
  {
    {
      2.7226276428059343e-8, 0.000011526596527333888, 0.0006812363504429289, 0.01177056050599653, 0.07916895586045004, 0.23787152296413616,
      0.34099234099234094, 0.23787152296413616, 0.07916895586045004, 0.01177056050599653, 0.0006812363504429289, 0.000011526596527333888,
      2.7226276428059343e-8
    },
    {
      -5.8001672523865, -4.591398448936521, -3.5634443802816342, -2.620689973432215, -1.7254183795882392, -0.85667949351945, 0,
      0.85667949351945, 1.7254183795882392, 2.620689973432215, 3.5634443802816342, 4.591398448936521, 5.8001672523865
    }
  },
  {
    {
      4.868161257748387e-9, 2.6609913440676444e-6, 0.0002003395537607452, 0.004428919106947417, 0.03865010882425344, 0.1540833398425136,
      0.30263462681301934, 0.30263462681301934, 0.1540833398425136, 0.03865010882425344, 0.004428919106947417, 0.0002003395537607452,
      2.6609913440676444e-6, 4.868161257748387e-9
    },
    {
      -6.087409546901291, -4.896936397345565, -3.886924575059769, -2.9630365798386675, -2.088344745701944, -1.242688955485464,
      -0.41259045795460186, 0.41259045795460186, 1.242688955485464, 2.088344745701944, 2.9630365798386675, 3.886924575059769,
      4.896936397345565, 6.087409546901291
    }
  },
  {
    {
      8.589649899633355e-10, 5.975419597920624e-7, 0.0000564214640518904, 0.001567357503549961, 0.01736577449213762, 0.08941779539984443,
      0.23246229360973228, 0.31825951825951826, 0.23246229360973228, 0.08941779539984443, 0.01736577449213762, 0.001567357503549961,
      0.0000564214640518904, 5.975419597920624e-7, 8.589649899633355e-10
    },
    {
      -6.363947888829839, -5.190093591304781, -4.1962077112690155, -3.2890824243987664, -2.432436827009758, -1.6067100690287297,
      -0.799129068324548, 0, 0.799129068324548, 1.6067100690287297, 2.432436827009758, 3.2890824243987664, 4.1962077112690155,
      5.190093591304781, 6.363947888829839
    }
  },
  {
    {
      1.4978147231618547e-10, 1.3094732162868277e-7, 0.000015300032162487367, 0.0005259849265739115, 0.0072669376011847515, 0.047284752354014054,
      0.15833837275094978, 0.28656852123801213, 0.28656852123801213, 0.15833837275094978, 0.047284752354014054, 0.0072669376011847515,
      0.0005259849265739115, 0.000015300032162487367, 1.3094732162868277e-7, 1.4978147231618547e-10
    },
    {
      -6.630878198393129, -5.472225705949343, -4.492955302520011, -3.600873624171548, -2.7602450476307014, -1.9519803457163334,
      -1.1638291005549648, -0.3867606045005573, 0.3867606045005573, 1.1638291005549648, 1.9519803457163334, 2.7602450476307014,
      3.600873624171548, 4.492955302520011, 5.472225705949343, 6.630878198393129
    }
  }
};

/// Get Hermite quadrature point coordinates
const double* uq::GetHermiteQuadraturePoints (const unsigned &numberOfQuadraturePoints) {
  if (numberOfQuadraturePoints < 1 || numberOfQuadraturePoints > 16) {
    std::cout << "Wrong Number of Quadature points (= " << numberOfQuadraturePoints << " ) in function uq::GetHermiteQuadaturePoints(...)" << std::endl;
    abort();
  }
  return _HermiteQuadrature[numberOfQuadraturePoints - 1][1];
}

/// Get Hermite quadrature weights
const double* uq::GetHermiteQuadratureWeights (const unsigned &numberOfQuadraturePoints) {
  if (numberOfQuadraturePoints < 1 || numberOfQuadraturePoints > 16) {
    std::cout << "Wrong Number of Quadature points (= " << numberOfQuadraturePoints << " ) in function uq::GetHermiteQuadatureWeights(...)" << std::endl;
    abort();
  }
  return _HermiteQuadrature[numberOfQuadraturePoints - 1][0];
}

////////////////////////////////////////////////

/// Compute the tensor product set corresponding to the pair < numberOfQuadraturePoints, numberOfEigPairs >
void uq::ComputeTensorProductSet (std::vector < std::vector <unsigned> > & Tp,
                                  const unsigned & numberOfQuadraturePoints,
                                  const unsigned & numberOfEigPairs) { //p is max poly degree

  unsigned tensorProductDim = pow (numberOfQuadraturePoints, numberOfEigPairs);

  std::cout << "tensorProductDim = " << tensorProductDim << std::endl;

  Tp.resize (tensorProductDim);
  for (unsigned i = 0; i < tensorProductDim; i++) {
    Tp[i].resize (numberOfEigPairs);
  }

  unsigned index = 0;
  unsigned counters[numberOfEigPairs + 1];
  memset (counters, 0, sizeof (counters));

  while (!counters[numberOfEigPairs]) {

    for (unsigned j = 0; j < numberOfEigPairs; j++) {
      Tp[index][j] = counters[numberOfEigPairs - 1 - j];
      std::cout << " Tp[" << index << "][" << j << "]= " << Tp[index][j] ;
    }
    std::cout << std::endl;
    index++;

    unsigned i;
    for (i = 0; counters[i] == numberOfQuadraturePoints - 1; i++) { // inner loops that are at maxval restart at zero
      counters[i] = 0;
    }
    ++counters[i];  // the innermost loop that isn't yet at maxval, advances by 1
  }
}

/// Return the tensor product set in the key < numberOfQuadraturePoints , numberOfEigPairs >
const std::vector < std::vector <unsigned> > & uq::GetTensorProductSet (const unsigned & numberOfQuadraturePoints,
                                                                        const unsigned & numberOfEigPairs) {

  std::pair<unsigned, unsigned> TpIndex = std::make_pair (numberOfQuadraturePoints, numberOfEigPairs);

  std::map<std::pair<unsigned, unsigned>, std::vector < std::vector <unsigned> > >::iterator it;

  it = _Tp.find (TpIndex);
  if (it == _Tp.end()) {
    ComputeTensorProductSet (_Tp[TpIndex], numberOfQuadraturePoints, numberOfEigPairs);
  }
  return _Tp[TpIndex];
};

/// Erase the tensor product set in the key < numberOfQuadraturePoints , numberOfEigPairs >
void uq::EraseTensorProductSet (const unsigned & numberOfQuadraturePoints, const unsigned & numberOfEigPairs) {
  _Tp.erase (std::make_pair (numberOfQuadraturePoints, numberOfEigPairs));
}

/// Clear all stored tensor product sets
void uq::ClearTensorProductSet() {
  _Tp.clear();
}

////////////////////////////////////////////////

/// Compute the Hermite Polynomial corresponding corresponding to the pair < numberOfQuadraturePoints, maxPolyOrder >
void uq::ComputeHermitePoly (std::vector < std::vector < double > >  & HermitePoly,
                             const unsigned &numberOfQuadraturePoints, const unsigned & maxPolyOrder) {

  if (numberOfQuadraturePoints < 1 || numberOfQuadraturePoints > 16) {
    std::cout << "The selected order of integraiton has not been implemented yet, choose an integer in [1,16]" << std::endl;
    abort();
  }

  else {
    HermitePoly.resize (maxPolyOrder + 1);
    for (unsigned i = 0; i < maxPolyOrder + 1; i++) {
      HermitePoly[i].resize (numberOfQuadraturePoints);
    }

    const double* HermiteQuadraturePoints = GetHermiteQuadraturePoints (numberOfQuadraturePoints);

    for (unsigned j = 0; j < numberOfQuadraturePoints; j++) {

      double x = HermiteQuadraturePoints[j];//_HermiteQuadrature[numberOfQuadraturePoints - 1][1][j];

      HermitePoly[0][j] = 1. ;

      if (maxPolyOrder > 0) {
        HermitePoly[1][j] = x ;
        if (maxPolyOrder > 1) {
          HermitePoly[2][j] = (pow (x, 2) - 1.) / sqrt (2) ;
          if (maxPolyOrder > 2) {
            HermitePoly[3][j] = (pow (x, 3) - 3. * x) / sqrt (6) ;
            if (maxPolyOrder > 3) {
              HermitePoly[4][j] = (pow (x, 4) - 6. * x * x + 3.) / sqrt (24) ;
              if (maxPolyOrder > 4) {
                HermitePoly[5][j] = (pow (x, 5) - 10. * pow (x, 3) + 15. * x) / sqrt (120) ;
                if (maxPolyOrder > 5) {
                  HermitePoly[6][j] = (pow (x, 6) - 15. * pow (x, 4) + 45. * pow (x, 2) - 15.) / sqrt (720) ;
                  if (maxPolyOrder > 6) {
                    HermitePoly[7][j] = (pow (x, 7) - 21. * pow (x, 5) + 105. * pow (x, 3) -  105. * x) / sqrt (5040) ;
                    if (maxPolyOrder > 7) {
                      HermitePoly[8][j] = (pow (x, 8) - 28. * pow (x, 6) + 210. * pow (x, 4) - 420. * pow (x, 2) + 105.) / sqrt (40320) ;
                      if (maxPolyOrder > 8) {
                        HermitePoly[9][j] = (pow (x, 9) - 36. * pow (x, 7) + 378. * pow (x, 5) - 1260. * pow (x, 3) + 945. * x) / sqrt (362880);
                        if (maxPolyOrder > 9) {
                          HermitePoly[10][j] = (pow (x, 10) - 45. * pow (x, 8) + 630. * pow (x, 6) - 3150. * pow (x, 4) + 4725. * pow (x, 2) - 945.) / sqrt (3628800);
                          if (maxPolyOrder > 10) {
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
}

/// Return the Hermite polynomial in the key < numberOfQuadraturePoints , maxPolyOrder >
const std::vector < std::vector <double> > & uq::GetHermitePolynomial (const unsigned & numberOfQuadraturePoints,
                                                                       const unsigned & maxPolyOrder) {

  std::pair<unsigned, unsigned> HermitePolyIndex = std::make_pair (numberOfQuadraturePoints, maxPolyOrder);

  std::map<std::pair<unsigned, unsigned>, std::vector < std::vector <double> > >::iterator it;

  it = _HermitePoly.find (HermitePolyIndex);
  if (it == _HermitePoly.end()) {
    ComputeHermitePoly (_HermitePoly[HermitePolyIndex], numberOfQuadraturePoints, maxPolyOrder);
  }
  return _HermitePoly[HermitePolyIndex];
};

/// Erase the Hermite polynomial in the key < numberOfQuadraturePoints , maxPolyOrder >
void uq::EraseHermitePolynomial (const unsigned & numberOfQuadraturePoints, const unsigned & maxPolyOrder) {
  _HermitePoly.erase (std::make_pair (numberOfQuadraturePoints, maxPolyOrder));
}

/// Clear all stored Hermite polynomials
void uq::ClearHermitePolynomial() {
  _HermitePoly.clear();
}

////////////////////////////////////////////////

void uq::ComputeHermitePolyHistogram (std::vector < std::vector < double > >  & HermitePoly, const unsigned & pIndex,
                                       std::vector<double> & samplePoints, const unsigned & numberOfEigPairs) {

  HermitePoly.resize (pIndex + 1);
  for (unsigned i = 0; i < pIndex + 1; i++) {
    HermitePoly[i].resize (samplePoints.size());
  }

  for (unsigned j = 0; j < numberOfEigPairs; j++) {

    double x = samplePoints[j];

    HermitePoly[0][j] = 1. ;

    if (pIndex > 0) {
      HermitePoly[1][j] = x ;
      if (pIndex > 1) {
        HermitePoly[2][j] = (pow (x, 2) - 1.) / sqrt (2) ;
        if (pIndex > 2) {
          HermitePoly[3][j] = (pow (x, 3) - 3. * x) / sqrt (6) ;
          if (pIndex > 3) {
            HermitePoly[4][j] = (pow (x, 4) - 6. * x * x + 3.) / sqrt (24) ;
            if (pIndex > 4) {
              HermitePoly[5][j] = (pow (x, 5) - 10. * pow (x, 3) + 15. * x) / sqrt (120) ;
              if (pIndex > 5) {
                HermitePoly[6][j] = (pow (x, 6) - 15. * pow (x, 4) + 45. * pow (x, 2) - 15.) / sqrt (720) ;
                if (pIndex > 6) {
                  HermitePoly[7][j] = (pow (x, 7) - 21. * pow (x, 5) + 105. * pow (x, 3) -  105. * x) / sqrt (5040) ;
                  if (pIndex > 7) {
                    HermitePoly[8][j] = (pow (x, 8) - 28. * pow (x, 6) + 210. * pow (x, 4) - 420. * pow (x, 4) + 105.) / sqrt (40320) ;
                    if (pIndex > 8) {
                      HermitePoly[9][j] = (pow (x, 9) - 36. * pow (x, 7) + 378. * pow (x, 5) - 1260. * pow (x, 3) + 945. * x) / sqrt (362880);
                      if (pIndex > 9) {
                        HermitePoly[10][j] = (pow (x, 10) - 45. * pow (x, 8) + 630. * pow (x, 6) - 3150. * pow (x, 4) + 4725. * pow (x, 2) - 945.) / sqrt (3628800);
                        if (pIndex > 10) {
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

////////////////////////////////////////////////

void uq::ComputeIndexSetJp (std::vector < std::vector <unsigned> > & Jp, 
                            const unsigned & p, const unsigned & numberOfEigPairs) { //p is max poly degree


  unsigned dimJp = static_cast <unsigned> (boost::math::binomial_coefficient<double> (numberOfEigPairs + p, p));

  //long unsigned dimJp = factorial(numberOfEigPairs + p) / (factorial(numberOfEigPairs) * factorial(p));

  Jp.resize (dimJp);
  for (unsigned i = 0; i < dimJp; i++) {
    Jp[i].resize (numberOfEigPairs);
  }

  unsigned index = 0;
  unsigned counters[numberOfEigPairs + 1];
  memset (counters, 0, sizeof (counters));

  while (!counters[numberOfEigPairs]) {

//     for(unsigned i = numberOfEigPairs; i-- > 0;) {
//       std::cout << counters[i] << " ";
//     }
//     std::cout << std::endl;

    unsigned entrySum = 0;
    for (unsigned j = 0; j < numberOfEigPairs; j++) {
      entrySum += counters[j];
    }

    if (entrySum <= p) {
      for (unsigned j = 0; j < numberOfEigPairs; j++) {
        Jp[index][j] = counters[numberOfEigPairs - 1 - j];
        std::cout << " Jp[" << index << "][" << j << "]= " << Jp[index][j] ;
      }
      std::cout << std::endl;
      index++;
    }
    unsigned i;
    for (i = 0; counters[i] == p; i++) { // inner loops that are at maxval restart at zero
      counters[i] = 0;
    }
    ++counters[i];  // the innermost loop that isn't yet at maxval, advances by 1
  }
}

/// Return the Index Set corresponding to the pair < p, numberOfEigPairs>
const std::vector < std::vector <unsigned> > & uq::GetIndexSet (const unsigned & p, const unsigned & numberOfEigPairs) {

  std::pair<unsigned, unsigned> JpIndex = std::make_pair (p, numberOfEigPairs);

  std::map<std::pair<unsigned, unsigned>, std::vector < std::vector <unsigned> > >::iterator it;

  it = _Jp.find (JpIndex);
  if (it == _Jp.end()) {
    ComputeIndexSetJp (_Jp[JpIndex], p, numberOfEigPairs);
  }
  return _Jp[JpIndex];
}

/// Erase the Index Set corresponding to the pair < p, numberOfEigPairs>
void uq::EraseIndexSet (const unsigned & p, const unsigned & numberOfEigPairs) {
  _Jp.erase (std::make_pair (p, numberOfEigPairs));
}

/// Clear all Index Sets Jp
void uq::ClearIndexSet() {
  _Jp.clear();
}

////////////////////////////////////////////////

/// Compute the Integral Matrix at the key < q0, p0>
void uq::ComputeIntegralMatrix (std::vector < std::vector < std::vector < double > > > &integralMatrix,
                                  const unsigned & q0, const unsigned & p0) {

  unsigned maxPolyOrder = (q0 > p0) ? q0 : p0;



  unsigned n1 = 2 * p0 + q0 + 1;
  n1 = (n1 % 2 == 0) ? n1 / 2 : (n1 + 1) / 2;
  unsigned numberOfQuadraturePoints = (n1 <= 16) ? n1 : 16;
  if (n1 > 16) {
    std::cout <<
              "------------------------------- WARNING: less quadrature points than needed were employed in function ComputeIntegralsMatrix -------------------------------"
              << std::endl;
    std::cout << " Needed : " << n1 << " , " << " Used : " << 16 << std::endl;
  }

  const std::vector < std::vector < double > >  &HermitePoly = GetHermitePolynomial (numberOfQuadraturePoints, maxPolyOrder);

  unsigned q = q0 + 1;
  unsigned p = p0 + 1;

  integralMatrix.resize (q);
  for (unsigned q1 = 0; q1 < q; q1++) {
    integralMatrix[q1].resize (p);
    for (unsigned p1 = 0; p1 < p; p1++) {
      integralMatrix[q1][p1].assign (p, 0.);
      for (unsigned p2 = 0; p2 < p; p2++) {
        integralMatrix[q1][p1][p2]  = 0.;
        for (unsigned i = 0; i < numberOfQuadraturePoints; i++) {
          double w = _HermiteQuadrature[numberOfQuadraturePoints - 1][0][i];
          integralMatrix[q1][p1][p2]  +=  w * HermitePoly[q1][i] * HermitePoly[p1][i] * HermitePoly[p2][i];
        }
      }
    }
  }

//   for(unsigned q1 = 0; q1 < q; q1++) {
//     for(unsigned p1 = 0; p1 < p; p1++) {
//       for(unsigned p2 = 0; p2 < p; p2++) {
//         std::cout << "integralMatrix[" << q1 << "][" << p1 << "][" << p2 << "]=" << integralMatrix[q1][p1][p2] << std::endl;
//       }
//     }
//   }

};

/// Get the Integral Matrix at the key < q0, p0>
const std::vector < std::vector < std::vector < double > > > & uq::GetIntegralMatrix(const unsigned & q0, 
                                                                                     const unsigned & p0){
  std::pair<unsigned, unsigned> integralMatrixIndex = std::make_pair (q0,p0);
  
  std::map<std::pair<unsigned, unsigned>, std::vector < std::vector < std::vector <double> > > >::iterator it;
  
  it = _integralMatrix.find (integralMatrixIndex);
  if (it == _integralMatrix.end()) {
    ComputeIntegralMatrix (_integralMatrix[integralMatrixIndex], q0, p0);
  }
  return _integralMatrix[integralMatrixIndex];
}
                                                                                    

/// Erase the Integral Matrix at the key < q0, p0>
void uq::EraseIntegralMatrix (const unsigned & q0, const unsigned & p0){
  _integralMatrix.erase (std::make_pair (q0,p0));
}

/// Clear all Integral Matrices
void uq::ClearIntegralMatrix(){
  _integralMatrix.clear ();
}

///////////////////////////////////////////

/// Compute the Integral Matrix at the key < q0, p0, numberOfEigPairs>
void uq::ComputeStochasticMassMatrix (std::vector < std::vector < std::vector < double > > > & G, 
                                        const unsigned & q0, const unsigned & p0, const unsigned & numberOfEigPairs) {
  
  const std::vector < std::vector < std::vector < double > > > & integralMatrix = GetIntegralMatrix(q0,p0);

  const std::vector < std::vector <unsigned> > &Jq = GetIndexSet (q0, numberOfEigPairs);
  const std::vector < std::vector <unsigned> > &Jp = GetIndexSet (p0, numberOfEigPairs);


  G.resize (Jq.size());
  for (unsigned q1 = 0; q1 < Jq.size(); q1++) {
    G[q1].resize (Jp.size());
    for (unsigned p1 = 0; p1 < Jp.size(); p1++) {
      G[q1][p1].assign (Jp.size(), 1.);
      for (unsigned p2 = 0; p2 < Jp.size(); p2++) {
        for (unsigned i = 0; i < numberOfEigPairs; i++) {
          G[q1][p1][p2] *= integralMatrix[Jq[q1][i]][Jp[p1][i]][Jp[p2][i]];
        }
        G[q1][p1][p2] = (fabs (G[q1][p1][p2]) < 1.e-14) ? 0. : G[q1][p1][p2];
      }
    }
  }

//   for(unsigned q1 = 0; q1 < Jq.size(); q1++) {
//     for(unsigned p1 = 0; p1 < Jp.size(); p1++) {
//       for(unsigned p2 = 0; p2 < Jp.size(); p2++) {
//         std::cout << "G[" << q1 << "][" << p1 << "][" << p2 << "]=" << G[q1][p1][p2] << std::endl;
//       }
//     }
//   }

}

/// Return the Stochastic Mass Matrix at the key < q0, p0, numberOfEigPairs>
std::vector < std::vector < std::vector < double > > > & uq::GetStochasticMassMatrix (const unsigned & q0, const unsigned & p0, 
                                                                                      const unsigned & numberOfEigPairs){
                             
  std::pair < std::pair<unsigned, unsigned>,unsigned> stochasticMassMatrixIndex = std::make_pair ( std::make_pair (q0, p0), numberOfEigPairs );
  
  std::map<std::pair < std::pair<unsigned, unsigned>,unsigned>, std::vector < std::vector < std::vector <double> > > >::iterator it;
  
  it = _stochasticMassMatrix.find (stochasticMassMatrixIndex);
  if (it == _stochasticMassMatrix.end()) {
    ComputeStochasticMassMatrix (_stochasticMassMatrix[stochasticMassMatrixIndex], q0, p0, numberOfEigPairs);
  }
  return _stochasticMassMatrix[stochasticMassMatrixIndex];
  
  
                                                                                      }


/// Erase the Stochastic Mass Matrix at the key < q0, p0, numberOfEigPairs>
void uq::EraseStochasticMassMatrix (const unsigned & q0, const unsigned & p0, 
                                const unsigned & numberOfEigPairs){
  _stochasticMassMatrix.erase( std::make_pair ( std::make_pair (q0, p0), numberOfEigPairs ) );
}

/// Clear all stored Stochastic Mass Matrices 
void uq::ClearStochasticMassMatrix(){
  _stochasticMassMatrix.clear();
}


///////////////////////////////////////////


void uq::ComputeMultivariateHermitePoly (
    std::vector < std::vector < double > >  & multivariateHermitePoly,
    std::vector < double > & multivariateHermiteQuadratureWeights,
    const unsigned & numberOfQuadraturePoints, const unsigned & p, const unsigned & numberOfEigPairs) {
    
  const std::vector < std::vector <unsigned> > & Jp = GetIndexSet (p,numberOfEigPairs);
  const std::vector < std::vector <unsigned> > & Tp = GetTensorProductSet (numberOfQuadraturePoints, numberOfEigPairs);

  multivariateHermiteQuadratureWeights.assign (Tp.size(), 1.);

  multivariateHermitePoly.resize (Jp.size());
  for (unsigned i = 0; i < Jp.size(); i++) {
    multivariateHermitePoly[i].assign (Tp.size(), 1.);
  }

  const std::vector < std::vector < double > >  &hermitePoly = GetHermitePolynomial (numberOfQuadraturePoints, p);

  for (unsigned j = 0; j < Tp.size(); j++) {
    for (unsigned k = 0; k < numberOfEigPairs; k++) {
      multivariateHermiteQuadratureWeights[j] *= _HermiteQuadrature[numberOfQuadraturePoints - 1][0][Tp[j][k]] ;
      for (unsigned i = 0; i < Jp.size(); i++) {
        multivariateHermitePoly[i][j] *= hermitePoly[Jp[i][k]][Tp[j][k]] ;
      }
    }
  }
  
};






//BEGIN Stochastic Input Parameters

unsigned pIndex = 4;
unsigned qIndex = 5;

int numberOfEigPairs = 2; //dimension of the stochastic variable
double stdDeviationInput = 0.2;  //standard deviation of the normal distribution (it is the same as the standard deviation of the covariance function in GetEigenPair)
double meanInput = 0.;
double amin = 1. / 100.; // for the KL expansion
std::vector < std::pair<double, double> > eigenvalues (numberOfEigPairs);
//END Stochastic Input Parameters

uq myuq;

void AssembleSysSG (MultiLevelProblem& ml_prob) {



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

  const std::vector < std::vector < double > >  &HermitePoly = myuq.GetHermitePolynomial (numberOfQuadraturePoints, maxPolyOrder);

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

        const double *HermiteQuadraturePoints = myuq.GetHermiteQuadraturePoints (numberOfQuadraturePointsForProjection);
        const double *HermiteQuadratureWeights = myuq.GetHermiteQuadratureWeights (numberOfQuadraturePointsForProjection);

        const std::vector < std::vector < double > >  &HermitePolyProjection = myuq.GetHermitePolynomial (numberOfQuadraturePointsForProjection, qIndex);

        for (unsigned i = 0; i < numberOfEigPairs; i++) {
          aStochasticTerm2[i] = 0.;
          for (unsigned j = 0; j < numberOfQuadraturePointsForProjection; j++) {
//             aStochasticTerm1[i] += HermitePoly[Jq[q1][i]][j] * HermiteQuadrature[numberOfQuadraturePoints - 1][0][j];
            aStochasticTerm2[i] += exp (sqrt (eigenvalues[i].first) * eigVectorGauss[i] * HermiteQuadraturePoints[j])
                                   * HermitePolyProjection[Jq[q1][i]][j] * HermiteQuadratureWeights[j];
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







