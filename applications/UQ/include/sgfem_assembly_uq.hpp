#include "UqQuadratureTypeEnum.hpp"
#include "MultiLevelSolution.hpp"

using namespace femus;

//THIS IS THE MOST UPDATED ASSEMBLY FOR SGM SIMULATIONS OF POISSON's EQUATION with HERMITE or LEGENDRE POLYNOMIALS

//BEGIN Stochastic Input Parameters

uq &myuq = FemusInit::_uqHermite;
// uq &myuq = FemusInit::_uqLegendre;

unsigned pIndex = 4;
unsigned qIndex = 5;

int  numberOfEigPairs = 2; //dimension of the stochastic variable
double stdDeviationInput = 0.08;  //standard deviation of the normal distribution (it is the same as the standard deviation of the covariance function in GetEigenPair)
double meanInput = 0.;
double amin = 1. / 100.; // for the KL expansion
std::vector < std::pair<double, double> > eigenvalues (numberOfEigPairs);

//END Stochastic Input Parameters

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

  const std::vector < std::vector < double > >  &poly = myuq.GetPolynomial (numberOfQuadraturePoints, maxPolyOrder);

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

  vector < vector < double > >  solu (Jp.size());   // local solution

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
        eigVectorGauss[i] = 0.;

        for (unsigned j = 0; j < nDofu; j++) {
          unsigned solDof = msh->GetSolutionDof (j, iel, soluType);
          eigVectorGauss[i] += (*sol->_Sol[eigfIndex[i]]) (solDof) * phi[j];
        }
      }


      vector< double > aStochastic (Jq.size());

//BEGIN coefficient obtained projecting the exponential of the KL
      for (unsigned q1 = 0; q1 < Jq.size(); q1 ++) {
        std::vector <double> aStochasticTerm2 (numberOfEigPairs);

        unsigned numberOfQuadraturePointsForProjection = 16;

        const double *quadraturePoints = myuq.GetQuadraturePoints (numberOfQuadraturePointsForProjection);
        const double *quadratureWeights = myuq.GetQuadratureWeights (numberOfQuadraturePointsForProjection);

        const std::vector < std::vector < double > >  &polyProjection = myuq.GetPolynomial (numberOfQuadraturePointsForProjection, qIndex);

        for (unsigned i = 0; i < numberOfEigPairs; i++) {
          aStochasticTerm2[i] = 0.;

          for (unsigned j = 0; j < numberOfQuadraturePointsForProjection; j++) {
            aStochasticTerm2[i] += exp (sqrt (eigenvalues[i].first) * eigVectorGauss[i] * quadraturePoints[j])
                                   * polyProjection[Jq[q1][i]][j] * quadratureWeights[j];
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

        aStochastic[q1] = amin * aS1 + aS2; //a_q(x_ig)

        if (fabs (aStochastic[q1]) > 10.) {
          std::cout << " coeff =  " << aStochastic[q1] << std::endl;
        }
      }
//END coefficient obtained projecting the exponential of the KL


//BEGIN random coefficient (not a random field), a(y1,y2,y3) = 1 + y1^2 + y2^2 + y3^2
//             for ( unsigned q1 = 0; q1 < Jq.size(); q1 ++ ) {
// 
//                 unsigned numberOfQuadraturePointsForProjection = 16;
// 
//                 const double *quadraturePoints = myuq.GetQuadraturePoints ( numberOfQuadraturePointsForProjection );
//                 const double *quadratureWeights = myuq.GetQuadratureWeights ( numberOfQuadraturePointsForProjection );
// 
//                 const std::vector < std::vector < double > >  &polyProjection = myuq.GetPolynomial ( numberOfQuadraturePointsForProjection, qIndex );
// 
//                 std::vector <double> aStochasticTerms ( numberOfEigPairs, 0. );
// 
//                 for ( unsigned i = 0; i < numberOfEigPairs; i++ ) {
//                     for ( unsigned j = 0; j < numberOfQuadraturePointsForProjection; j++ ) {
//                         aStochasticTerms[i] += quadraturePoints[j] * quadraturePoints[j] * polyProjection[Jq[q1][i]][j] * quadratureWeights[j];
//                     }
// 
//                     for ( unsigned k = 0; k < numberOfEigPairs; k++ ) {
//                         if ( k != i ) aStochasticTerms[i] *= integralMatrix[Jq[q1][k]][0][0];
//                     }
//                 }
// 
//                 double aS1 = 1.;
// 
//                 for ( unsigned i = 0; i < numberOfEigPairs; i++ ) {
//                     aS1 *= integralMatrix[Jq[q1][i]][0][0];
//                 }
// 
//                 aStochastic[q1] = aS1;
// 
//                 for ( unsigned i = 0; i < numberOfEigPairs; i++ ) {
//                     aStochastic[q1] += aStochasticTerms[i];
//                 }
// 
// 
//                 if ( fabs ( aStochastic[q1] ) > 10. ) {
//                     std::cout << " coeff =  " << aStochastic[q1] << std::endl;
//                 }
//             }
//END coefficient without KL, not a field, just dependent on y

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
                AG += aStochastic[q1] * laplace[i][j] * G[q1][p1][p2];
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







