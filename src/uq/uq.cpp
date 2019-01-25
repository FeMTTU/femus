
#include "uq.hpp"

namespace femus {

  const double uq::_hermiteQuadrature[16][2][16] = { //Number of quadrature points, first row: weights, second row: coordinates
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

  const double uq::_legendreQuadrature[16][2][16] = { //Number of quadrature points, first row: weights, second row: coordinates
    {{1.}, {0.}},
    { {0.4999999999999996, 0.4999999999999996},
      {-0.5773502691896257, 0.5773502691896257 }
    },
    { {0.2777777777777779, 0.4444444444444444, 0.2777777777777779},
      {-0.7745966692414834, 0, 0.7745966692414834}
    },
    { {0.17392742256872698, 0.32607257743127316, 0.32607257743127316, 0.17392742256872698},
      {-0.8611363115940526, -0.33998104358485626, 0.33998104358485626, 0.8611363115940526}
    },
    { {0.11846344252809478, 0.2393143352496832, 0.28444444444444444, 0.2393143352496832, 0.11846344252809478},
      {-0.906179845938664, -0.5384693101056831, 0., 0.5384693101056831, 0.906179845938664}
    },
    { {0.08566224618958442, 0.1803807865240689, 0.23395696728634563, 0.23395696728634563, 0.1803807865240689, 0.08566224618958442},
      {-0.932469514203152, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645, 0.932469514203152}
    },
    { {
        0.06474248308443509, 0.1398526957446383, 0.19091502525255943, 0.2089795918367347, 0.19091502525255943, 0.1398526957446383, 0.06474248308443509
      },
      { -0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0., 0.4058451513773972, 0.7415311855993945, 0.9491079123427585}
    },
    { {
        0.05061426814518826, 0.11119051722668719, 0.15685332293894363, 0.18134189168918102, 0.18134189168918102, 0.15685332293894363, 0.11119051722668719, 0.05061426814518826
      },
      {
        -0.9602898564975363, -0.7966664774136267, -0.525532409916329, -0.1834346424956498, 0.1834346424956498, 0.525532409916329, 0.7966664774136267, 0.9602898564975363
        }
    },
    { {
        0.04063719418078653, 0.09032408034742867, 0.13030534820146764, 0.15617353852000135, 0.1651196775006299, 0.15617353852000135, 0.13030534820146764, 0.09032408034742867,
        0.04063719418078653
      },
      {
        -0.9681602395076261, -0.8360311073266358, -0.6133714327005904, -0.3242534234038089, 0., 0.3242534234038089, 0.6133714327005904, 0.8360311073266358, 0.9681602395076261
        }
    },
    { {
        0.033335672154343965, 0.07472567457529057, 0.10954318125799119, 0.13463335965499823, 0.14776211235737635, 0.14776211235737635, 0.13463335965499823, 0.10954318125799119,
        0.07472567457529057, 0.033335672154343965
      },
      {
        -0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.14887433898163122, 0.14887433898163122, 0.4333953941292472, 0.6794095682990244,
          0.8650633666889845, 0.9739065285171717
        }
    },
    {
      {
        0.02783428355808611, 0.0627901847324527, 0.09314510546386695, 0.11659688229599524, 0.13140227225512335, 0.1364625433889503, 0.13140227225512335, 0.11659688229599524,
        0.09314510546386695, 0.0627901847324527, 0.027834283558086117
      },
      {
        -0.978228658146057, -0.8870625997680953, -0.7301520055740494, -0.5190961292068118, -0.26954315595234496, 0., 0.26954315595234496, 0.5190961292068118, 0.7301520055740494,
          0.8870625997680953, 0.978228658146057
        }
    },
    {
      {
        0.02358766819325607, 0.053469662997659574, 0.08003916427167311, 0.10158371336153296, 0.11674626826917732, 0.12457352290670136, 0.12457352290670136, 0.11674626826917732,
        0.10158371336153296, 0.08003916427167311, 0.053469662997659574, 0.02358766819325607
      },
      {
        -0.9815606342467192, -0.9041172563704749, -0.7699026741943047, -0.5873179542866175, -0.3678314989981802, -0.1252334085114689, 0.1252334085114689, 0.3678314989981802,
          0.5873179542866175, 0.7699026741943047, 0.9041172563704749, 0.9815606342467192
        }
    },
    {
      {
        0.020242002382658098, 0.04606074991886419, 0.06943675510989344, 0.08907299038097279, 0.10390802376844416, 0.11314159013144862, 0.11627577661543695, 0.11314159013144862,
        0.10390802376844416, 0.08907299038097279, 0.06943675510989344, 0.04606074991886419, 0.020242002382658098
      },
      {
        -0.9841830547185881, -0.9175983992229779, -0.8015780907333099, -0.6423493394403402, -0.44849275103644687, -0.2304583159551348, 0, 0.2304583159551348, 0.44849275103644687, 0.6423493394403402, 0.8015780907333099, 0.9175983992229779, 0.9841830547185881
        }
    },
    {
      {
        0.017559730165875483, 0.04007904357987999, 0.06075928534395155, 0.07860158357909665, 0.092769198738969, 0.10259923186064775, 0.10763192673157891, 0.10763192673157891,
        0.10259923186064775, 0.092769198738969, 0.07860158357909665, 0.06075928534395155, 0.04007904357987999, 0.017559730165875483
      },
      {
        -0.9862838086968123, -0.9284348836635735, -0.827201315069765, -0.6872929048116855, -0.5152486363581541, -0.31911236892788974, -0.10805494870734367, 0.10805494870734367,
          0.31911236892788974, 0.5152486363581541, 0.6872929048116855, 0.827201315069765, 0.9284348836635735, 0.9862838086968123
        }
    },
    {
      {
        0.015376620998058374, 0.03518302374405419, 0.05357961023358571, 0.06978533896307716, 0.08313460290849681, 0.09308050000778098, 0.09921574266355576, 0.10128912096278064,
        0.09921574266355576, 0.09308050000778098, 0.08313460290849681, 0.06978533896307716, 0.05357961023358571, 0.03518302374405419, 0.015376620998058374
      },
      {
        -0.9879925180204854, -0.937273392400706, -0.8482065834104272, -0.7244177313601701, -0.5709721726085388, -0.3941513470775634, -0.20119409399743451, 0, 0.20119409399743451, 0.3941513470775634, 0.5709721726085388, 0.7244177313601701, 0.8482065834104272, 0.937273392400706, 0.9879925180204854
        }
    },
    {
      {
        0.01357622970587767, 0.031126761969323898, 0.04757925584124643, 0.06231448562776654, 0.07479799440828827, 0.08457825969750124, 0.09130170752246176, 0.09472530522753425,
        0.09472530522753425, 0.09130170752246176, 0.08457825969750124, 0.07479799440828827, 0.06231448562776654, 0.04757925584124643, 0.031126761969323898, 0.01357622970587767
      },
      {
        -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.755404408355003, -0.6178762444026438, -0.45801677765722737, -0.2816035507792589, -0.09501250983763744,
          0.09501250983763744, 0.2816035507792589, 0.45801677765722737, 0.6178762444026438, 0.755404408355003, 0.8656312023878318, 0.9445750230732326, 0.9894009349916499
        }
    }
  };

/// Get Hermite quadrature point coordinates
  const double* uq::GetHermiteQuadraturePoints (const unsigned &numberOfQuadraturePoints) {
    if (numberOfQuadraturePoints < 1 || numberOfQuadraturePoints > 16) {
      std::cout << "Wrong Number of Quadature points (= " << numberOfQuadraturePoints << " ) in function uq::GetHermiteQuadaturePoints(...)" << std::endl;
      abort();
    }
    return _hermiteQuadrature[numberOfQuadraturePoints - 1][1];
  }

  /// Get Legendre quadrature point coordinates
  const double* uq::GetLegendreQuadraturePoints (const unsigned &numberOfQuadraturePoints) {
    if (numberOfQuadraturePoints < 1 || numberOfQuadraturePoints > 16) {
      std::cout << "Wrong Number of Quadature points (= " << numberOfQuadraturePoints << " ) in function uq::GetLegendreQuadaturePoints(...)" << std::endl;
      abort();
    }
    return _legendreQuadrature[numberOfQuadraturePoints - 1][1];

  }

/// Get quadrature point coordinates
  const double* uq::GetQuadraturePoints (const unsigned &numberOfQuadraturePoints) {
    if (_quadratureType == UQ_HERMITE) {
      return GetHermiteQuadraturePoints (numberOfQuadraturePoints);
    }
    else {
      return GetLegendreQuadraturePoints (numberOfQuadraturePoints);
    }
  }


/// Get Hermite quadrature weights
  const double* uq::GetHermiteQuadratureWeights (const unsigned &numberOfQuadraturePoints) {
    if (numberOfQuadraturePoints < 1 || numberOfQuadraturePoints > 16) {
      std::cout << "Wrong Number of Quadature points (= " << numberOfQuadraturePoints << " ) in function uq::GetHermiteQuadatureWeights(...)" << std::endl;
      abort();
    }
    return _hermiteQuadrature[numberOfQuadraturePoints - 1][0];

  }

/// Get Legendre quadrature weights
  const double* uq::GetLegendreQuadratureWeights (const unsigned &numberOfQuadraturePoints) {
    if (numberOfQuadraturePoints < 1 || numberOfQuadraturePoints > 16) {
      std::cout << "Wrong Number of Quadature points (= " << numberOfQuadraturePoints << " ) in function uq::GetLegendreQuadatureWeights(...)" << std::endl;
      abort();
    }
    return _legendreQuadrature[numberOfQuadraturePoints - 1][0];

  }


/// Get quadrature weights
  const double* uq::GetQuadratureWeights (const unsigned &numberOfQuadraturePoints) {
    if (_quadratureType == UQ_HERMITE) {
      return GetHermiteQuadratureWeights (numberOfQuadraturePoints);
    }
    else {
      return GetLegendreQuadratureWeights (numberOfQuadraturePoints);
    }
  }

////////////////////////////////////////////////

/// Compute the tensor product set corresponding to the pair < numberOfQuadraturePoints, numberOfEigPairs >
  void uq::ComputeTensorProductSet (std::vector < std::vector <unsigned> > & Tp,
                                    const unsigned & numberOfQuadraturePoints,
                                    const unsigned & numberOfEigPairs) { //p is max poly degree

    unsigned tensorProductDim = pow (numberOfQuadraturePoints, numberOfEigPairs);

    if (_output) {
      std::cout << "tensorProductDim = " << tensorProductDim << std::endl;
    }

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
        if (_output) {
          std::cout << " Tp[" << index << "][" << j << "]= " << Tp[index][j] ;
        }
      }
      if (_output) {
        std::cout << std::endl;
      }
      index++;

      unsigned i;
      for (i = 0; counters[i] == numberOfQuadraturePoints - 1; i++) {   // inner loops that are at maxval restart at zero
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

////////////////////////////////////////////////  POLY STUFF

/// Compute the Hermite Polynomial corresponding corresponding to the pair < numberOfQuadraturePoints, maxPolyOrder >
  void uq::ComputeHermitePolynomial (std::vector < std::vector < double > >  & hermitePoly,
                                     const unsigned &numberOfQuadraturePoints, const unsigned & maxPolyOrder) {

    if (numberOfQuadraturePoints < 1 || numberOfQuadraturePoints > 16) {
      std::cout << "The selected order of integration has not been implemented yet, choose an integer in [1,16]" << std::endl;
      abort();
    }

    else {
      hermitePoly.resize (maxPolyOrder + 1);
      for (unsigned i = 0; i < maxPolyOrder + 1; i++) {
        hermitePoly[i].resize (numberOfQuadraturePoints);
      }

      const double* hermiteQuadraturePoints = GetHermiteQuadraturePoints (numberOfQuadraturePoints);

      for (unsigned j = 0; j < numberOfQuadraturePoints; j++) {

        double x = hermiteQuadraturePoints[j];

        hermitePoly[0][j] = 1. ;

        if (maxPolyOrder > 0) {
          hermitePoly[1][j] = x ;
          if (maxPolyOrder > 1) {
            hermitePoly[2][j] = (pow (x, 2) - 1.) / sqrt (2) ;
            if (maxPolyOrder > 2) {
              hermitePoly[3][j] = (pow (x, 3) - 3. * x) / sqrt (6) ;
              if (maxPolyOrder > 3) {
                hermitePoly[4][j] = (pow (x, 4) - 6. * x * x + 3.) / sqrt (24) ;
                if (maxPolyOrder > 4) {
                  hermitePoly[5][j] = (pow (x, 5) - 10. * pow (x, 3) + 15. * x) / sqrt (120) ;
                  if (maxPolyOrder > 5) {
                    hermitePoly[6][j] = (pow (x, 6) - 15. * pow (x, 4) + 45. * pow (x, 2) - 15.) / sqrt (720) ;
                    if (maxPolyOrder > 6) {
                      hermitePoly[7][j] = (pow (x, 7) - 21. * pow (x, 5) + 105. * pow (x, 3) -  105. * x) / sqrt (5040) ;
                      if (maxPolyOrder > 7) {
                        hermitePoly[8][j] = (pow (x, 8) - 28. * pow (x, 6) + 210. * pow (x, 4) - 420. * pow (x, 2) + 105.) / sqrt (40320) ;
                        if (maxPolyOrder > 8) {
                          hermitePoly[9][j] = (pow (x, 9) - 36. * pow (x, 7) + 378. * pow (x, 5) - 1260. * pow (x, 3) + 945. * x) / sqrt (362880);
                          if (maxPolyOrder > 9) {
                            hermitePoly[10][j] = (pow (x, 10) - 45. * pow (x, 8) + 630. * pow (x, 6) - 3150. * pow (x, 4) + 4725. * pow (x, 2) - 945.) / sqrt (3628800);
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

/// Compute the Legendre Polynomial corresponding corresponding to the pair < numberOfQuadraturePoints, maxPolyOrder >
  void uq::ComputeLegendrePolynomial (std::vector < std::vector < double > >  & legendrePoly,
                                      const unsigned &numberOfQuadraturePoints, const unsigned & maxPolyOrder) {

    if (numberOfQuadraturePoints < 1 || numberOfQuadraturePoints > 16) {
      std::cout << "The selected order of integration has not been implemented yet, choose an integer in [1,16]" << std::endl;
      abort();
    }

    else {
      legendrePoly.resize (maxPolyOrder + 1);
      for (unsigned i = 0; i < maxPolyOrder + 1; i++) {
        legendrePoly[i].resize (numberOfQuadraturePoints);
      }

      const double* legendreQuadraturePoints = GetLegendreQuadraturePoints (numberOfQuadraturePoints);

      for (unsigned j = 0; j < numberOfQuadraturePoints; j++) {

        double x = legendreQuadraturePoints[j];

        legendrePoly[0][j] = 1. ;

        if (maxPolyOrder > 0) {
          legendrePoly[1][j] = x * sqrt (3) ;
          if (maxPolyOrder > 1) {
            legendrePoly[2][j] = (- 1. + 3. * x * x) * 0.5 * sqrt (5) ;
            if (maxPolyOrder > 2) {
              legendrePoly[3][j] = (- 3. * x + 5. * x * x * x) * 0.5 * sqrt (7) ;
              if (maxPolyOrder > 3) {
                legendrePoly[4][j] = (sqrt (9) / 8.) * (3. - 30. * x * x + 35. * pow (x, 4)) ;
                if (maxPolyOrder > 4) {
                  legendrePoly[5][j] = (15. * x - 70. * pow (x, 3) + 63. * pow (x, 5)) * sqrt (11) / 8. ;
                  if (maxPolyOrder > 5) {
                    legendrePoly[6][j] = (- 5. + 105. * x * x - 315. * pow (x, 4) + 231. * pow (x, 6)) * sqrt (13) / 16. ;
                    if (maxPolyOrder > 6) {
                      legendrePoly[7][j] = (- 35. * x + 315. * pow (x, 3) - 693. * pow (x, 5) + 429. * pow (x, 7)) * sqrt (15) / 16.;
                      if (maxPolyOrder > 7) {
                        legendrePoly[8][j] = (35. - 1260. * x * x + 6930. * pow (x, 4) - 12012. * pow (x, 6) + 6435. * pow (x, 8)) *                sqrt (17)  / 128.;
                        if (maxPolyOrder > 8) {
                          legendrePoly[9][j] = (315. * x - 4620. * pow (x, 3) + 18018. * pow (x, 5) - 25740. * pow (x, 7) + 12155. * pow (x, 9)) * sqrt (19) / 128.;
                          if (maxPolyOrder > 9) {
                            legendrePoly[10][j] = (- 63. + 3465. * x * x - 30030. * pow (x, 4) + 90090. * pow (x, 6) - 109395. * pow (x, 8) + 46189 * pow (x, 10)) * sqrt (21) / 256. ;
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

/// Return the polynomial in the key < numberOfQuadraturePoints , maxPolyOrder >
  const std::vector < std::vector <double> > & uq::GetPolynomial (const unsigned & numberOfQuadraturePoints,
                                                                  const unsigned & maxPolyOrder) {

    std::pair<unsigned, unsigned> PolyIndex = std::make_pair (numberOfQuadraturePoints, maxPolyOrder);

    std::map<std::pair<unsigned, unsigned>, std::vector < std::vector <double> > >::iterator it;

    it = _polynomial.find (PolyIndex);
    if (it == _polynomial.end()) {
      if (_quadratureType == UQ_HERMITE) {
        ComputeHermitePolynomial (_polynomial[PolyIndex], numberOfQuadraturePoints, maxPolyOrder);
      }
      else {
        ComputeLegendrePolynomial (_polynomial[PolyIndex], numberOfQuadraturePoints, maxPolyOrder);
      }
    }
    return _polynomial[PolyIndex];
  };

/// Erase the polynomial in the key < numberOfQuadraturePoints , maxPolyOrder >
  void uq::ErasePolynomial (const unsigned & numberOfQuadraturePoints, const unsigned & maxPolyOrder) {
    _polynomial.erase (std::make_pair (numberOfQuadraturePoints, maxPolyOrder));
  }

/// Clear all stored polynomials
  void uq::ClearPolynomial() {
    _polynomial.clear();
  }



/// Compute the Hermite Polynomial Histogram at the key < pIndex, samplePoints, numberOfEigPairs >
  const std::vector < std::vector < double > >  & uq::GetHermitePolynomialHistogram (
    const unsigned & pIndex, const std::vector<double> & samplePoints, const unsigned & numberOfEigPairs) {

    _polynomialHistogram.resize (pIndex + 1);
    for (unsigned i = 0; i < pIndex + 1; i++) {
      _polynomialHistogram[i].resize (samplePoints.size());
    }

    for (unsigned j = 0; j < numberOfEigPairs; j++) {

      double x = samplePoints[j];

      _polynomialHistogram[0][j] = 1. ;

      if (pIndex > 0) {
        _polynomialHistogram[1][j] = x ;
        if (pIndex > 1) {
          _polynomialHistogram[2][j] = (pow (x, 2) - 1.) / sqrt (2) ;
          if (pIndex > 2) {
            _polynomialHistogram[3][j] = (pow (x, 3) - 3. * x) / sqrt (6) ;
            if (pIndex > 3) {
              _polynomialHistogram[4][j] = (pow (x, 4) - 6. * x * x + 3.) / sqrt (24) ;
              if (pIndex > 4) {
                _polynomialHistogram[5][j] = (pow (x, 5) - 10. * pow (x, 3) + 15. * x) / sqrt (120) ;
                if (pIndex > 5) {
                  _polynomialHistogram[6][j] = (pow (x, 6) - 15. * pow (x, 4) + 45. * pow (x, 2) - 15.) / sqrt (720) ;
                  if (pIndex > 6) {
                    _polynomialHistogram[7][j] = (pow (x, 7) - 21. * pow (x, 5) + 105. * pow (x, 3) -  105. * x) / sqrt (5040) ;
                    if (pIndex > 7) {
                      _polynomialHistogram[8][j] = (pow (x, 8) - 28. * pow (x, 6) + 210. * pow (x, 4) - 420. * pow (x, 4) + 105.) / sqrt (40320) ;
                      if (pIndex > 8) {
                        _polynomialHistogram[9][j] = (pow (x, 9) - 36. * pow (x, 7) + 378. * pow (x, 5) - 1260. * pow (x, 3) + 945. * x) / sqrt (362880);
                        if (pIndex > 9) {
                          _polynomialHistogram[10][j] = (pow (x, 10) - 45. * pow (x, 8) + 630. * pow (x, 6) - 3150. * pow (x, 4) + 4725. * pow (x, 2) - 945.) / sqrt (3628800);
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
    return _polynomialHistogram;
  }

/// Compute the Legendre Polynomial Histogram at the key < pIndex, samplePoints, numberOfEigPairs >
  const std::vector < std::vector < double > >  & uq::GetLegendrePolynomialHistogram (
    const unsigned & pIndex, const std::vector<double> & samplePoints, const unsigned & numberOfEigPairs) {

    _polynomialHistogram.resize (pIndex + 1);
    for (unsigned i = 0; i < pIndex + 1; i++) {
      _polynomialHistogram[i].resize (samplePoints.size());
    }

    for (unsigned j = 0; j < numberOfEigPairs; j++) {

      double x = samplePoints[j];

      _polynomialHistogram[0][j] = 1. ;

      if (pIndex > 0) {
        _polynomialHistogram[1][j] = x / sqrt (3) ;
        if (pIndex > 1) {
          _polynomialHistogram[2][j] = (- 1. + 3. * x * x) / (2 * sqrt (5)) ;
          if (pIndex > 2) {
            _polynomialHistogram[3][j] = (- 3. * x + 5. * x * x * x) / (2. * sqrt (7)) ;
            if (pIndex > 3) {
              _polynomialHistogram[4][j] = (1. / 24.) * (3. - 30. * x * x + 35. * pow (x, 4)) ;
              if (pIndex > 4) {
                _polynomialHistogram[5][j] = (15. * x - 70. * pow (x, 3) + 63. * pow (x, 5)) / (8. * sqrt (11)) ;
                if (pIndex > 5) {
                  _polynomialHistogram[6][j] = (- 5. + 105. * x * x - 315. * pow (x, 4) + 231. * pow (x, 6)) / (16. * sqrt (13)) ;
                  if (pIndex > 6) {
                    _polynomialHistogram[7][j] = (- 35. * x + 315. * pow (x, 3) - 693. * pow (x, 5) + 429. * pow (x, 7)) / (16.* sqrt (15));
                    if (pIndex > 7) {
                      _polynomialHistogram[8][j] = (35. - 1260. * x * x + 6930. * pow (x, 4) - 12012. * pow (x, 6) + 6435. * pow (x, 8)) / (128. *                sqrt (17));
                      if (pIndex > 8) {
                        _polynomialHistogram[9][j] = (315. * x - 4620. * pow (x, 3) + 18018. * pow (x, 5) - 25740. * pow (x, 7) + 12155. * pow (x, 9)) / (128. * sqrt (19));
                        if (pIndex > 9) {
                          _polynomialHistogram[10][j] = (- 63. + 3465. * x * x - 30030. * pow (x, 4) + 90090. * pow (x, 6) - 109395. * pow (x, 8) + 46189 * pow (x, 10)) / (256. * sqrt (21));
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
    return _polynomialHistogram;
  }

  const std::vector < std::vector < double > >  & uq::GetPolynomialHistogram (
    const unsigned & pIndex, const std::vector<double> & samplePoints, const unsigned & numberOfEigPairs) {

    if (_quadratureType == UQ_HERMITE) {
      return GetHermitePolynomialHistogram (pIndex, samplePoints, numberOfEigPairs);
    }
    else {
      return GetLegendrePolynomialHistogram (pIndex, samplePoints, numberOfEigPairs);
    }
  }

////////////////////////////////////////////////

  void uq::ComputeIndexSet (std::vector < std::vector <unsigned> > & Jp,
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
          if (_output) {
            std::cout << " Jp[" << index << "][" << j << "]= " << Jp[index][j] ;
          }
        }
        if (_output) {
          std::cout << std::endl;
        }
        index++;
      }
      unsigned i;
      for (i = 0; counters[i] == p; i++) {   // inner loops that are at maxval restart at zero
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
      ComputeIndexSet (_Jp[JpIndex], p, numberOfEigPairs);
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

    const std::vector < std::vector < double > >  &poly = GetPolynomial (numberOfQuadraturePoints, maxPolyOrder);

    unsigned q = q0 + 1;
    unsigned p = p0 + 1;

    const double* quadratureWeights = GetQuadratureWeights (numberOfQuadraturePoints);

    integralMatrix.resize (q);
    for (unsigned q1 = 0; q1 < q; q1++) {
      integralMatrix[q1].resize (p);
      for (unsigned p1 = 0; p1 < p; p1++) {
        integralMatrix[q1][p1].assign (p, 0.);
        for (unsigned p2 = 0; p2 < p; p2++) {
          integralMatrix[q1][p1][p2]  = 0.;
          for (unsigned i = 0; i < numberOfQuadraturePoints; i++) {
            double w = quadratureWeights[i];
            integralMatrix[q1][p1][p2]  +=  w * poly[q1][i] * poly[p1][i] * poly[p2][i];
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
  const std::vector < std::vector < std::vector < double > > > & uq::GetIntegralMatrix (const unsigned & q0, const unsigned & p0) {
    std::pair<unsigned, unsigned> integralMatrixIndex = std::make_pair (q0, p0);

    std::map<std::pair<unsigned, unsigned>, std::vector < std::vector < std::vector <double> > > >::iterator it;

    it = _integralMatrix.find (integralMatrixIndex);
    if (it == _integralMatrix.end()) {
      ComputeIntegralMatrix (_integralMatrix[integralMatrixIndex], q0, p0);
    }
    return _integralMatrix[integralMatrixIndex];
  }


/// Erase the Integral Matrix at the key < q0, p0>
  void uq::EraseIntegralMatrix (const unsigned & q0, const unsigned & p0) {
    _integralMatrix.erase (std::make_pair (q0, p0));
  }

/// Clear all Integral Matrices
  void uq::ClearIntegralMatrix() {
    _integralMatrix.clear ();
  }

///////////////////////////////////////////

/// Compute the Integral Matrix at the key < q0, p0, numberOfEigPairs>
  void uq::ComputeStochasticMassMatrix (std::vector < std::vector < std::vector < double > > > & G,
                                        const unsigned & q0, const unsigned & p0, const unsigned & numberOfEigPairs) {

    const std::vector < std::vector < std::vector < double > > > & integralMatrix = GetIntegralMatrix (q0, p0);

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
                                                                                        const unsigned & numberOfEigPairs) {

    std::pair < std::pair<unsigned, unsigned>, unsigned> stochasticMassMatrixIndex = std::make_pair (std::make_pair (q0, p0), numberOfEigPairs);

    std::map<std::pair < std::pair<unsigned, unsigned>, unsigned>, std::vector < std::vector < std::vector <double> > > >::iterator it;

    it = _stochasticMassMatrix.find (stochasticMassMatrixIndex);
    if (it == _stochasticMassMatrix.end()) {
      ComputeStochasticMassMatrix (_stochasticMassMatrix[stochasticMassMatrixIndex], q0, p0, numberOfEigPairs);
    }
    return _stochasticMassMatrix[stochasticMassMatrixIndex];


  }


/// Erase the Stochastic Mass Matrix at the key < q0, p0, numberOfEigPairs>
  void uq::EraseStochasticMassMatrix (const unsigned & q0, const unsigned & p0,
                                      const unsigned & numberOfEigPairs) {
    _stochasticMassMatrix.erase (std::make_pair (std::make_pair (q0, p0), numberOfEigPairs));
  }

/// Clear all stored Stochastic Mass Matrices
  void uq::ClearStochasticMassMatrix() {
    _stochasticMassMatrix.clear();
  }


/////////////////////////////////////////// MULTIVARIATE HERMITE STUFF


/// Compute the Multivariate  polynomials and weights at the key < numberOfQuadraturePoints, p,numberOfEigPairs>
  void uq::ComputeMultivariate (
    std::vector < std::vector < double > >  & multivariatePoly,
    std::vector < double > & multivariateQuadratureWeights,
    const unsigned & numberOfQuadraturePoints, const unsigned & p, const unsigned & numberOfEigPairs) {

    const std::vector < std::vector <unsigned> > & Jp = GetIndexSet (p, numberOfEigPairs);
    const std::vector < std::vector <unsigned> > & Tp = GetTensorProductSet (numberOfQuadraturePoints, numberOfEigPairs);

    multivariateQuadratureWeights.assign (Tp.size(), 1.);

    multivariatePoly.resize (Jp.size());
    for (unsigned i = 0; i < Jp.size(); i++) {
      multivariatePoly[i].assign (Tp.size(), 1.);
    }

    const std::vector < std::vector < double > >  &Polynomial = GetPolynomial (numberOfQuadraturePoints, p);

    const double* QuadratureWeights = GetQuadratureWeights (numberOfQuadraturePoints);

    for (unsigned j = 0; j < Tp.size(); j++) {
      for (unsigned k = 0; k < numberOfEigPairs; k++) {
        multivariateQuadratureWeights[j] *= QuadratureWeights[Tp[j][k]] ;
        for (unsigned i = 0; i < Jp.size(); i++) {
          multivariatePoly[i][j] *= Polynomial[Jp[i][k]][Tp[j][k]] ;
        }
      }
    }
  };

/// Return the Multivariate  polynomials at the key < numberOfQuadraturePoints, p,numberOfEigPairs>
  const std::vector < std::vector < double > >  & uq::GetMultivariatePolynomial (
    const unsigned & numberOfQuadraturePoints, const unsigned & p, const unsigned & numberOfEigPairs) {

    std::pair < std::pair<unsigned, unsigned>, unsigned> multivariateIndex = std::make_pair (std::make_pair (numberOfQuadraturePoints, p), numberOfEigPairs);

    std::map<std::pair < std::pair<unsigned, unsigned>, unsigned>, std::vector < std::vector <double> > >::iterator it;

    it = _multivariatePolynomial.find (multivariateIndex);
    if (it == _multivariatePolynomial.end()) {
      ComputeMultivariate (_multivariatePolynomial[multivariateIndex], _multivariateWeight[multivariateIndex], numberOfQuadraturePoints, p, numberOfEigPairs);
    }
    return _multivariatePolynomial[multivariateIndex];

  }

/// Return the Multivariate  weight at the key < numberOfQuadraturePoints, p,numberOfEigPairs>
  const std::vector < double > & uq::GetMultivariateWeights (
    const unsigned & numberOfQuadraturePoints, const unsigned & p, const unsigned & numberOfEigPairs) {

    std::pair < std::pair<unsigned, unsigned>, unsigned> multivariateIndex = std::make_pair (std::make_pair (numberOfQuadraturePoints, p), numberOfEigPairs);

    std::map<std::pair < std::pair<unsigned, unsigned>, unsigned>, std::vector < std::vector <double> > >::iterator it;

    it = _multivariatePolynomial.find (multivariateIndex);
    if (it == _multivariatePolynomial.end()) {
      ComputeMultivariate (_multivariatePolynomial[multivariateIndex], _multivariateWeight[multivariateIndex], numberOfQuadraturePoints, p, numberOfEigPairs);
    }
    return _multivariateWeight[multivariateIndex];

  }

/// Erase the Multivariate  polynomials and weights at the key < numberOfQuadraturePoints, p,numberOfEigPairs>
  void uq::EraseMultivariate (const unsigned & numberOfQuadraturePoints, const unsigned & p, const unsigned & numberOfEigPairs) {
    std::pair < std::pair<unsigned, unsigned>, unsigned> multivariateIndex = std::make_pair (std::make_pair (numberOfQuadraturePoints, p), numberOfEigPairs);
    _multivariatePolynomial.erase (multivariateIndex);
    _multivariateWeight.erase (multivariateIndex);
  }

/// Clear all Multivariate  polynomials and weights
  void uq::ClearMultivariate() {
    _multivariatePolynomial.clear();
    _multivariateWeight.clear();
  }

///////////////////////////////////////////

/// Set the standard output for some of the objects;
  void uq::SetOutput (const bool & output) {
    _output = output;
  };


}











