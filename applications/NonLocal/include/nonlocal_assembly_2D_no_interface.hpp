#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>


//THIS IS THE 2D ASSEMBLY FOR THE NONLOCAL DIFFUSION PROBLEM (NO INTERFACE)

using namespace femus;

bool nonLocalAssembly = true;

double delta = 0.25;
double kappa = 1.;

void GetBoundaryFunctionValue (double &value, const std::vector < double >& x) {

  value = 0.;
//     value = x[0];
//   value = x[0] * x[0];
//     value = ( x[0] < 0. ) ? x[0] * x[0] * x[0] : 3 * x[0] * x[0] * x[0];
//     value = x[0] * x[0] * x[0] + x[1] * x[1] * x[1];
//     value = x[0] * x[0] * x[0] * x[0] + 0.1 * x[0] * x[0];
//     value = x[0] * x[0] * x[0] * x[0];
//     value =  2 * x[0] + x[0] * x[0] * x[0] * x[0] * x[0]; //this is 2x + x^5

}

void ReorderElement (std::vector < int > &dofs, std::vector < double > & sol, std::vector < std::vector < double > > & x);

void RectangleAndBallRelation (bool &theyIntersect, const std::vector<double> &ballCenter, const double &ballRadius, const std::vector < std::vector < double> > &elementCoordinates,  std::vector < std::vector < double> > &newCoordinates);

void RectangleAndBallRelation2 (bool & theyIntersect, const std::vector<double> &ballCenter, const double & ballRadius, const std::vector < std::vector < double> > &elementCoordinates, std::vector < std::vector < double> > &newCoordinates);

const elem_type *fem = new const elem_type_2D ("quad", "linear", "second");   //to use a different quadrature rule in the inner integral

const elem_type *femQuadrature = new const elem_type_2D ("quad", "linear", "eighth");   //to use a different quadrature rule in the inner integral

void AssembleNonLocalSys (MultiLevelProblem& ml_prob) {
  adept::Stack& s = FemusInit::_adeptStack;

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("NonLocal");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel (level);
  elem*                     el = msh->el;

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel (level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             KK = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + ! (dim - 1));
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the noumber of processes (for parallel computation)

  unsigned soluIndex = mlSol->GetIndex ("u");   // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType (soluIndex);   // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex ("u");   // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solu; // local solution for the local assembly (it uses adept)
  solu.reserve (maxSize);

  vector < double >  solu1; // local solution for the nonlocal assembly
  vector < double >  solu2; // local solution for the nonlocal assembly
  solu1.reserve (maxSize);
  solu2.reserve (maxSize);

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector < vector < double > > x1 (dim);
  vector < vector < double > > x2 (dim);

  for (unsigned k = 0; k < dim; k++) {
    x1[k].reserve (maxSize);
    x2[k].reserve (maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  phi.reserve (maxSize);
  phi_x.reserve (maxSize * dim);
  phi_xx.reserve (maxSize * dim2);

  vector< adept::adouble > aRes; // local redidual vector
  aRes.reserve (maxSize);

  vector< int > l2GMap1; // local to global mapping
  vector< int > l2GMap2; // local to global mapping
  l2GMap1.reserve (maxSize);
  l2GMap2.reserve (maxSize);

  vector< double > Res1; // local redidual vector
  Res1.reserve (maxSize);
  vector< double > Res2; // local redidual vector
  Res2.reserve (maxSize);

  vector < double > Jac11;
  Jac11.reserve (maxSize * maxSize);
  vector < double > Jac12;
  Jac12.reserve (maxSize * maxSize);

  vector < double > Jac21;
  Jac21.reserve (maxSize * maxSize);
  vector < double > Jac22;
  Jac22.reserve (maxSize * maxSize);


  KK->zero(); // Set to zero all the entries of the Global Matrix

  //BEGIN nonlocal assembly

  for (int kproc = 0; kproc < nprocs; kproc++) {
    for (int jel = msh->_elementOffset[kproc]; jel < msh->_elementOffset[kproc + 1]; jel++) {

      short unsigned jelGeom;
      short unsigned jelGroup;
      unsigned nDof2;
      //unsigned nDofx2;

      if (iproc == kproc) {
        jelGeom = msh->GetElementType (jel);
        jelGroup = msh->GetElementGroup (jel);
        nDof2  = msh->GetElementDofNumber (jel, soluType);
      }

      MPI_Bcast (&jelGeom, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD);
      MPI_Bcast (&jelGroup, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD);
      MPI_Bcast (&nDof2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);

      l2GMap2.resize (nDof2);
      solu2.resize (nDof2);

      for (int k = 0; k < dim; k++) {
        x2[k].resize (nDof2);
      }

      if (iproc == kproc) {
        for (unsigned j = 0; j < nDof2; j++) {
          l2GMap2[j] = pdeSys->GetSystemDof (soluIndex, soluPdeIndex, j, jel);
          unsigned solDof = msh->GetSolutionDof (j, jel, soluType);
          solu2[j] = (*sol->_Sol[soluIndex]) (solDof);
          unsigned xDof  = msh->GetSolutionDof (j, jel, xType);

          for (unsigned k = 0; k < dim; k++) {
            x2[k][j] = (*msh->_topology->_Sol[k]) (xDof);
          }
        }

        ReorderElement (l2GMap2, solu2, x2);
      }

      MPI_Bcast (&l2GMap2[0], nDof2, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      MPI_Bcast (&solu2[0], nDof2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);

      for (unsigned k = 0; k < dim; k++) {
        MPI_Bcast (& x2[k][0], nDof2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }

      for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

        bool midpointQuadrature = false;

        short unsigned ielGeom = msh->GetElementType (iel);
        short unsigned ielGroup = msh->GetElementGroup (iel);
        unsigned nDof1  = msh->GetElementDofNumber (iel, soluType);

        l2GMap1.resize (nDof1);
        solu1.resize (nDof1);

        Jac11.assign (nDof1 * nDof1, 0.);
        Jac12.assign (nDof1 * nDof2, 0.);
        Jac21.assign (nDof2 * nDof1, 0.);
        Jac22.assign (nDof2 * nDof2, 0.);
        Res1.assign (nDof1, 0.);
        Res2.assign (nDof2, 0.);

        for (int k = 0; k < dim; k++) {
          x1[k].resize (nDof1);
        }

        for (unsigned i = 0; i < nDof1; i++) {
          l2GMap1[i] = pdeSys->GetSystemDof (soluIndex, soluPdeIndex, i, iel);
          unsigned solDof = msh->GetSolutionDof (i, iel, soluType);
          solu1[i] = (*sol->_Sol[soluIndex]) (solDof);
          unsigned xDof  = msh->GetSolutionDof (i, iel, xType);

          for (unsigned k = 0; k < dim; k++) {
            x1[k][i] = (*msh->_topology->_Sol[k]) (xDof);
          }
        }

        ReorderElement (l2GMap1, solu1, x1);

        double sideLength = fabs (x1[0][0] - x1[0][1]);

        double leftBoundInterface = - sideLength;
        double rightBoundInterface = sideLength;

        unsigned igNumber = (midpointQuadrature) ? 4 : msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber();
        vector < vector < double > > xg1 (igNumber);
        vector <double> weight1 (igNumber);
        vector < vector <double> > phi1x (igNumber);

        if (midpointQuadrature) {

          for (unsigned ig = 0; ig < igNumber; ig++) {

            std::vector <double> xg1Local (dim);

            weight1[ig] = 0.25 * sideLength * sideLength;

            xg1[ig].assign (dim, 0.);

            unsigned midpointDof = ig + 4;
            unsigned xDof  = msh->GetSolutionDof (midpointDof, iel, xType);

            for (unsigned k = 0; k < dim; k++) {
              xg1[ig][k] = (*msh->_topology->_Sol[k]) (xDof);
//                                 std::cout<< xg1[ig][k] << std::endl;
            }

            for (unsigned k = 0; k < dim; k++) {
              xg1Local[k] = - 1. + 2. * (xg1[ig][k] - x1[k][k]) / (x1[k][k + 1] - x1[k][k]);
            }

            double weightTemp;
            msh->_finiteElement[ielGeom][soluType]->Jacobian (x1, xg1Local, weightTemp, phi1x[ig], phi_x);
          }
        }

        else {

          for (unsigned ig = 0; ig < igNumber; ig++) {
            msh->_finiteElement[ielGeom][soluType]->Jacobian (x1, ig, weight1[ig], phi1x[ig], phi_x);

            xg1[ig].assign (dim, 0.);

            for (unsigned i = 0; i < nDof1; i++) {
              for (unsigned k = 0; k < dim; k++) {
                xg1[ig][k] += x1[k][i] * phi1x[ig][i];
              }
            }
          }

        }

        double radius = delta;

        bool coarseIntersectionTest = true;

        for (unsigned k = 0; k < dim; k++) {
          double min = 1.0e10;
          min = (min < fabs (x1[k][k]   - x2[k][k])) ?    min :  fabs (x1[k][k]   - x2[k][k]);
          min = (min < fabs (x1[k][k]   - x2[k][k + 1])) ?  min :  fabs (x1[k][k]   - x2[k][k + 1]);
          min = (min < fabs (x1[k][k + 1] - x2[k][k])) ?    min :  fabs (x1[k][k + 1] - x2[k][k]);
          min = (min < fabs (x1[k][k + 1] - x2[k][k + 1])) ?  min :  fabs (x1[k][k + 1] - x2[k][k + 1]);

          if (min >= radius - 1.0e-10) {
            coarseIntersectionTest = false;
            break;
          }
        }

        if (coarseIntersectionTest) {

          bool ifAnyIntersection = false;

          for (unsigned ig = 0; ig < igNumber; ig++) {

            if (iel == jel) {
              for (unsigned i = 0; i < nDof1; i++) {
//                                 Res1[i] -= 0. * weight[ig] * phi1x[ig][i]; //Ax - f (so f = 0)
                Res1[i] -=  1. * weight1[ig]  * phi1x[ig][i]; //Ax - f (so f = 1)
//                 Res1[i] -=  - 2. * weight1[ig]  * phi1x[ig][i]; //Ax - f (so f = - 2)
//                                 Res1[i] -=  - 6. * xg1[ig][0] * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = - 6 x)
                // Res1[i] -=  - 6. * ( xg1[ig][0] + xg1[ig][1] ) * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = - 6 (x + y))
//                                 Res1[i] -= ( - 12. * xg1[ig][0] * xg1[ig][0] - 6. / 5. * radius * radius - 2. * radius ) * weight1[ig] * phi1x[ig][i];  //Ax - f (so f = - 12x^2 - 6/5 * delta^2 - 2 delta)
//                                      Res1[i] -=  - 20. * ( xg1[ig][0] * xg1[ig][0] * xg1[ig][0] ) * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = - 20 x^3 )
//                                 Res1[i] -=  - 12. * ( xg1[ig][0] * xg1[ig][0] ) * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = - 12 x^2 )
              }
            }

            std::vector< std::vector < double > > x2New;
            bool theyIntersect;
            RectangleAndBallRelation2 (theyIntersect, xg1[ig], radius, x2, x2New);

            if (theyIntersect) {

              ifAnyIntersection = true;

              unsigned jgNumber = msh->_finiteElement[jelGeom][soluType]->GetGaussPointNumber();
//                             unsigned jgNumber = fem->GetGaussPointNumber();

              for (unsigned jg = 0; jg < jgNumber; jg++) {

                vector <double>  phi2y;
                double weight2;

                msh->_finiteElement[jelGeom][soluType]->Jacobian (x2New, jg, weight2, phi2y, phi_x);
//                                 fem->Jacobian ( x2New, jg, weight2, phi2y, phi_x );

                std::vector< double > xg2 (dim, 0.);

                for (unsigned j = 0; j < nDof2; j++) {
                  for (unsigned k = 0; k < dim; k++) {
                    xg2[k] += x2New[k][j] * phi2y[j];
                  }
                }

                std::vector <double> xg2Local (dim);

                for (unsigned k = 0; k < dim; k++) {
                  xg2Local[k] = - 1. + 2. * (xg2[k] - x2[k][k]) / (x2[k][k + 1] - x2[k][k]);
                }

                double weightTemp;
                msh->_finiteElement[jelGeom][soluType]->Jacobian (x2, xg2Local, weightTemp, phi2y, phi_x);
//                                 fem->Jacobian ( x2, xg2Local, weightTemp, phi2y, phi_x );

                double kernel = 0.75 * kappa / (delta * delta * delta * delta) ;

                for (unsigned i = 0; i < nDof1; i++) {
                  for (unsigned j = 0; j < nDof1; j++) {
                    double jacValue11 = weight1[ig] * weight2 * kernel * (phi1x[ig][i]) * phi1x[ig][j];
                    Jac11[i * nDof1 + j] -= jacValue11;
                    Res1[i] +=  jacValue11 * solu1[j];
                  }

                  for (unsigned j = 0; j < nDof2; j++) {
                    double jacValue12 = - weight1[ig] * weight2 * kernel * (phi1x[ig][i]) * phi2y[j];
                    Jac12[i * nDof2 + j] -= jacValue12;
                    Res1[i] +=  jacValue12 * solu2[j];
                  }//endl j loop
                }

                for (unsigned i = 0; i < nDof2; i++) {
                  for (unsigned j = 0; j < nDof1; j++) {
                    double jacValue21 = weight1[ig] * weight2 * kernel * (- phi2y[i]) * phi1x[ig][j];
                    Jac21[i * nDof1 + j] -= jacValue21;
                    Res2[i] +=  jacValue21 * solu1[j];
                  }

                  for (unsigned j = 0; j < nDof2; j++) {
                    double jacValue22 = - weight1[ig] * weight2 * kernel * (- phi2y[i]) * phi2y[j];
                    Jac22[i * nDof2 + j] -= jacValue22;
                    Res2[i] +=  jacValue22 * solu2[j];
                  }//endl j loop
                } //endl i loop
              }//end jg loop
            }
          }//end ig loop

          if (ifAnyIntersection) {
            KK->add_matrix_blocked (Jac11, l2GMap1, l2GMap1);
            KK->add_matrix_blocked (Jac12, l2GMap1, l2GMap2);
            RES->add_vector_blocked (Res1, l2GMap1);

            KK->add_matrix_blocked (Jac21, l2GMap2, l2GMap1);
            KK->add_matrix_blocked (Jac22, l2GMap2, l2GMap2);
            RES->add_vector_blocked (Res2, l2GMap2);
          }
        }
      } //end iel loop
    } // end jel loop
  } //end kproc loop

  RES->close();

  KK->close();

//     Mat A = ( static_cast<PetscMatrix*> ( KK ) )->mat();
//     MatAssemblyBegin ( A, MAT_FINAL_ASSEMBLY );
//     MatAssemblyEnd ( A, MAT_FINAL_ASSEMBLY );
//     PetscViewer viewer;
//     MatView ( A, viewer );

//     Vec v = ( static_cast< PetscVector* > ( RES ) )->vec();
//     VecView(v,PETSC_VIEWER_STDOUT_WORLD);

  // ***************** END ASSEMBLY *******************
}




void AssembleLocalSys (MultiLevelProblem& ml_prob) {
  adept::Stack& s = FemusInit::_adeptStack;

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Local");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel (level);
  elem*                     el = msh->el;

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel (level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             KK = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + ! (dim - 1));
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the noumber of processes (for parallel computation)

  unsigned soluIndex = mlSol->GetIndex ("u_local");   // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType (soluIndex);   // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex ("u_local");   // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solu; // local solution for the local assembly (it uses adept)
  solu.reserve (maxSize);

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector < vector < double > > x1 (dim);

  for (unsigned k = 0; k < dim; k++) {
    x1[k].reserve (maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  phi.reserve (maxSize);
  phi_x.reserve (maxSize * dim);
  phi_xx.reserve (maxSize * dim2);

  vector< adept::adouble > aRes; // local redidual vector
  aRes.reserve (maxSize);

  vector< int > l2GMap1; // local to global mapping
  l2GMap1.reserve (maxSize);

  vector< double > Res1; // local redidual vector
  Res1.reserve (maxSize);

  vector < double > Jac11;
  Jac11.reserve (maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  //BEGIN local assembly

//   for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
//
//     short unsigned ielGroup = msh->GetElementGroup (iel);
//
//     if (ielGroup == 5 || ielGroup == 6) {   //5 and 6 are the boundary surfaces
//
//       unsigned nDofu  = msh->GetElementDofNumber (iel, soluType);
//       std::vector <double> dofCoordinates (dim);
//
//       for (unsigned i = 0; i < nDofu; i++) {
//         unsigned solDof = msh->GetSolutionDof (i, iel, soluType);
//         unsigned xDof = msh->GetSolutionDof (i, iel, xType);
//         sol->_Bdc[soluIndex]->set (solDof, 0.);
//
//         for (unsigned jdim = 0; jdim < dim; jdim++) {
//           dofCoordinates[jdim] = (*msh->_topology->_Sol[jdim]) (xDof);
//         }
//
//         double bdFunctionValue;
//         GetBoundaryFunctionValue (bdFunctionValue, dofCoordinates);
//         sol->_Sol[soluIndex]->set (solDof, bdFunctionValue);
//
//       }
//
//     }
//
//   }
//
//   sol->_Bdc[soluIndex]->close();
//   sol->_Sol[soluIndex]->close();


  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofu  = msh->GetElementDofNumber (iel, soluType);   // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber (iel, xType);   // number of coordinate element dofs

    // resize local arrays
    l2GMap1.resize (nDofu);
    solu.resize (nDofu);

    for (int i = 0; i < dim; i++) {
      x1[i].resize (nDofx);
    }

    aRes.resize (nDofu);   //resize
    std::fill (aRes.begin(), aRes.end(), 0);   //set aRes to zero

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, soluType);   // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex]) (solDof);     // global extraction and local storage for the solution
      l2GMap1[i] = pdeSys->GetSystemDof (soluIndex, soluPdeIndex, i, iel);   // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x1[jdim][i] = (*msh->_topology->_Sol[jdim]) (xDof);     // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian (x1, ig, weight, phi, phi_x, boost::none);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point

      vector < adept::adouble > gradSolu_gss (dim, 0.);
      vector < double > x_gss (dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        for (unsigned jdim = 0; jdim < dim; jdim++) {
          gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
          x_gss[jdim] += x1[jdim][i] * phi[i];
        }
      }


//       double aCoeff = 1.;
      double aCoeff = kappa;

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofu; i++) {

        adept::adouble laplace = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          laplace   +=  aCoeff * phi_x[i * dim + jdim] * gradSolu_gss[jdim];
        }

//                 double srcTerm =  12. * x_gss[0] * x_gss[0] ; // so f = - 12 x^2
//         double srcTerm =  2. ; // so f = - 2
        double srcTerm =  - 1. ; // so f = 1
        //double srcTerm =  0./*- GetExactSolutionLaplace(x_gss)*/ ;
        aRes[i] += (srcTerm * phi[i] + laplace) * weight;

      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    Res1.resize (nDofu);   //resize

    for (int i = 0; i < nDofu; i++) {
      Res1[i] = - aRes[i].value();
    }

    RES->add_vector_blocked (Res1, l2GMap1);

    // define the dependent variables
    s.dependent (&aRes[0], nDofu);

    // define the independent variables
    s.independent (&solu[0], nDofu);

    // get the jacobian matrix (ordered by row major )
    Jac11.resize (nDofu * nDofu);   //resize
    s.jacobian (&Jac11[0], true);

    //store K in the global matrix KK
    KK->add_matrix_blocked (Jac11, l2GMap1, l2GMap1);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  //END local assembly

  RES->close();

  KK->close();

//     Mat A = ( static_cast<PetscMatrix*> ( KK ) )->mat();
//     MatAssemblyBegin ( A, MAT_FINAL_ASSEMBLY );
//     MatAssemblyEnd ( A, MAT_FINAL_ASSEMBLY );
//     PetscViewer viewer;
//     MatView ( A, viewer );

//     Vec v = ( static_cast< PetscVector* > ( RES ) )->vec();
//     VecView(v,PETSC_VIEWER_STDOUT_WORLD);

  // ***************** END ASSEMBLY *******************
}



//BEGIN fine nonlocal assembly for h convergence computations


void AssembleNonLocalSysFine (MultiLevelProblem& ml_prob) {
  adept::Stack& s = FemusInit::_adeptStack;

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("NonLocalFine");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel (level);
  elem*                     el = msh->el;

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel (level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             KK = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + ! (dim - 1));
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the noumber of processes (for parallel computation)

  unsigned soluIndex = mlSol->GetIndex ("u_fine");   // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType (soluIndex);   // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex ("u_fine");   // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solu; // local solution for the local assembly (it uses adept)
  solu.reserve (maxSize);

  vector < double >  solu1; // local solution for the nonlocal assembly
  vector < double >  solu2; // local solution for the nonlocal assembly
  solu1.reserve (maxSize);
  solu2.reserve (maxSize);

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector < vector < double > > x1 (dim);
  vector < vector < double > > x2 (dim);

  for (unsigned k = 0; k < dim; k++) {
    x1[k].reserve (maxSize);
    x2[k].reserve (maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  phi.reserve (maxSize);
  phi_x.reserve (maxSize * dim);
  phi_xx.reserve (maxSize * dim2);

  vector< adept::adouble > aRes; // local redidual vector
  aRes.reserve (maxSize);

  vector< int > l2GMap1; // local to global mapping
  vector< int > l2GMap2; // local to global mapping
  l2GMap1.reserve (maxSize);
  l2GMap2.reserve (maxSize);

  vector< double > Res1; // local redidual vector
  Res1.reserve (maxSize);
  vector< double > Res2; // local redidual vector
  Res1.reserve (maxSize);

  vector < double > Jac11;
  Jac11.reserve (maxSize * maxSize);
  vector < double > Jac12;
  Jac12.reserve (maxSize * maxSize);

  vector < double > Jac21;
  Jac21.reserve (maxSize * maxSize);
  vector < double > Jac22;
  Jac22.reserve (maxSize * maxSize);


  KK->zero(); // Set to zero all the entries of the Global Matrix

  //BEGIN nonlocal assembly

  for (int kproc = 0; kproc < nprocs; kproc++) {
    for (int jel = msh->_elementOffset[kproc]; jel < msh->_elementOffset[kproc + 1]; jel++) {

      short unsigned jelGeom;
      short unsigned jelGroup;
      unsigned nDof2;
      //unsigned nDofx2;

      if (iproc == kproc) {
        jelGeom = msh->GetElementType (jel);
        jelGroup = msh->GetElementGroup (jel);
        nDof2  = msh->GetElementDofNumber (jel, soluType);
      }

      MPI_Bcast (&jelGeom, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD);
      MPI_Bcast (&jelGroup, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD);
      MPI_Bcast (&nDof2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);

      l2GMap2.resize (nDof2);
      solu2.resize (nDof2);

      for (int k = 0; k < dim; k++) {
        x2[k].resize (nDof2);
      }

      if (iproc == kproc) {
        for (unsigned j = 0; j < nDof2; j++) {
          l2GMap2[j] = pdeSys->GetSystemDof (soluIndex, soluPdeIndex, j, jel);
          unsigned solDof = msh->GetSolutionDof (j, jel, soluType);
          solu2[j] = (*sol->_Sol[soluIndex]) (solDof);
          unsigned xDof  = msh->GetSolutionDof (j, jel, xType);

          for (unsigned k = 0; k < dim; k++) {
            x2[k][j] = (*msh->_topology->_Sol[k]) (xDof);
          }
        }

        ReorderElement (l2GMap2, solu2, x2);
      }

      MPI_Bcast (&l2GMap2[0], nDof2, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      MPI_Bcast (&solu2[0], nDof2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);

      for (unsigned k = 0; k < dim; k++) {
        MPI_Bcast (& x2[k][0], nDof2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }

      for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

        bool midpointQuadrature = false;

        short unsigned ielGeom = msh->GetElementType (iel);
        short unsigned ielGroup = msh->GetElementGroup (iel);
        unsigned nDof1  = msh->GetElementDofNumber (iel, soluType);

        l2GMap1.resize (nDof1);
        solu1.resize (nDof1);

        Jac11.assign (nDof1 * nDof1, 0.);
        Jac12.assign (nDof1 * nDof2, 0.);
        Jac21.assign (nDof2 * nDof1, 0.);
        Jac22.assign (nDof2 * nDof2, 0.);
        Res1.assign (nDof1, 0.);
        Res2.assign (nDof2, 0.);

        for (int k = 0; k < dim; k++) {
          x1[k].resize (nDof1);
        }

        for (unsigned i = 0; i < nDof1; i++) {
          l2GMap1[i] = pdeSys->GetSystemDof (soluIndex, soluPdeIndex, i, iel);
          unsigned solDof = msh->GetSolutionDof (i, iel, soluType);
          solu1[i] = (*sol->_Sol[soluIndex]) (solDof);
          unsigned xDof  = msh->GetSolutionDof (i, iel, xType);

          for (unsigned k = 0; k < dim; k++) {
            x1[k][i] = (*msh->_topology->_Sol[k]) (xDof);
          }
        }

        ReorderElement (l2GMap1, solu1, x1);

        double sideLength = fabs (x1[0][0] - x1[0][1]);

        double leftBoundInterface = - sideLength;
        double rightBoundInterface = sideLength;

        unsigned igNumber = (midpointQuadrature) ? 4 : msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber();
        vector < vector < double > > xg1 (igNumber);
        vector <double> weight1 (igNumber);
        vector < vector <double> > phi1x (igNumber);

        if (midpointQuadrature) {

          for (unsigned ig = 0; ig < igNumber; ig++) {

            std::vector <double> xg1Local (dim);

            weight1[ig] = 0.25 * sideLength * sideLength;

            xg1[ig].assign (dim, 0.);

            unsigned midpointDof = ig + 4;
            unsigned xDof  = msh->GetSolutionDof (midpointDof, iel, xType);

            for (unsigned k = 0; k < dim; k++) {
              xg1[ig][k] = (*msh->_topology->_Sol[k]) (xDof);
//                                 std::cout<< xg1[ig][k] << std::endl;
            }

            for (unsigned k = 0; k < dim; k++) {
              xg1Local[k] = - 1. + 2. * (xg1[ig][k] - x1[k][k]) / (x1[k][k + 1] - x1[k][k]);
            }

            double weightTemp;
            msh->_finiteElement[ielGeom][soluType]->Jacobian (x1, xg1Local, weightTemp, phi1x[ig], phi_x);
          }
        }

        else {

          for (unsigned ig = 0; ig < igNumber; ig++) {
            msh->_finiteElement[ielGeom][soluType]->Jacobian (x1, ig, weight1[ig], phi1x[ig], phi_x);

            xg1[ig].assign (dim, 0.);

            for (unsigned i = 0; i < nDof1; i++) {
              for (unsigned k = 0; k < dim; k++) {
                xg1[ig][k] += x1[k][i] * phi1x[ig][i];
              }
            }
          }

        }

        double radius = delta;

        bool coarseIntersectionTest = true;

        for (unsigned k = 0; k < dim; k++) {
          double min = 1.0e10;
          min = (min < fabs (x1[k][k]   - x2[k][k])) ?    min :  fabs (x1[k][k]   - x2[k][k]);
          min = (min < fabs (x1[k][k]   - x2[k][k + 1])) ?  min :  fabs (x1[k][k]   - x2[k][k + 1]);
          min = (min < fabs (x1[k][k + 1] - x2[k][k])) ?    min :  fabs (x1[k][k + 1] - x2[k][k]);
          min = (min < fabs (x1[k][k + 1] - x2[k][k + 1])) ?  min :  fabs (x1[k][k + 1] - x2[k][k + 1]);

          if (min >= radius - 1.0e-10) {
            coarseIntersectionTest = false;
            break;
          }
        }

        if (coarseIntersectionTest) {

          bool ifAnyIntersection = false;

          for (unsigned ig = 0; ig < igNumber; ig++) {

            if (iel == jel) {
              for (unsigned i = 0; i < nDof1; i++) {
//                                 Res1[i] -= 0. * weight[ig] * phi1x[ig][i]; //Ax - f (so f = 0)
                Res1[i] -=  1. * weight1[ig]  * phi1x[ig][i]; //Ax - f (so f = 1)
//                 Res1[i] -=  - 2. * weight1[ig]  * phi1x[ig][i]; //Ax - f (so f = - 2 )
//                                 Res1[i] -=  - 6. * xg1[ig][0] * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = - 6 x)
                // Res1[i] -=  - 6. * ( xg1[ig][0] + xg1[ig][1] ) * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = - 6 (x + y))
//                                 Res1[i] -= ( - 12. * xg1[ig][0] * xg1[ig][0] - 6. / 5. * radius * radius - 2. * radius ) * weight1[ig] * phi1x[ig][i];  //Ax - f (so f = - 12x^2 - 6/5 * delta^2 - 2 delta)
//                                      Res1[i] -=  - 20. * ( xg1[ig][0] * xg1[ig][0] * xg1[ig][0] ) * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = - 20 x^3 )
//                                 Res1[i] -=  - 12. * ( xg1[ig][0] * xg1[ig][0] ) * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = - 12 x^2 )
              }
            }

            std::vector< std::vector < double > > x2New;
            bool theyIntersect;
            RectangleAndBallRelation2 (theyIntersect, xg1[ig], radius, x2, x2New);

            if (theyIntersect) {

              ifAnyIntersection = true;

              unsigned jgNumber = msh->_finiteElement[jelGeom][soluType]->GetGaussPointNumber();
//                             unsigned jgNumber = fem->GetGaussPointNumber();

              for (unsigned jg = 0; jg < jgNumber; jg++) {

                vector <double>  phi2y;
                double weight2;

                msh->_finiteElement[jelGeom][soluType]->Jacobian (x2New, jg, weight2, phi2y, phi_x);
//                                 fem->Jacobian ( x2New, jg, weight2, phi2y, phi_x );

                std::vector< double > xg2 (dim, 0.);

                for (unsigned j = 0; j < nDof2; j++) {
                  for (unsigned k = 0; k < dim; k++) {
                    xg2[k] += x2New[k][j] * phi2y[j];
                  }
                }

                std::vector <double> xg2Local (dim);

                for (unsigned k = 0; k < dim; k++) {
                  xg2Local[k] = - 1. + 2. * (xg2[k] - x2[k][k]) / (x2[k][k + 1] - x2[k][k]);
                }

                double weightTemp;
                msh->_finiteElement[jelGeom][soluType]->Jacobian (x2, xg2Local, weightTemp, phi2y, phi_x);
//                                 fem->Jacobian ( x2, xg2Local, weightTemp, phi2y, phi_x );

                double kernel = 0.75 * kappa / (delta * delta * delta * delta) ;


                for (unsigned i = 0; i < nDof1; i++) {
                  for (unsigned j = 0; j < nDof1; j++) {
                    double jacValue11 = weight1[ig] * weight2 * kernel * (phi1x[ig][i]) * phi1x[ig][j];
                    Jac11[i * nDof1 + j] -= jacValue11;
                    Res1[i] +=  jacValue11 * solu1[j];
                  }

                  for (unsigned j = 0; j < nDof2; j++) {
                    double jacValue12 = - weight1[ig] * weight2 * kernel * (phi1x[ig][i]) * phi2y[j];
                    Jac12[i * nDof2 + j] -= jacValue12;
                    Res1[i] +=  jacValue12 * solu2[j];
                  }//endl j loop
                }

                for (unsigned i = 0; i < nDof2; i++) {
                  for (unsigned j = 0; j < nDof1; j++) {
                    double jacValue21 = weight1[ig] * weight2 * kernel * (- phi2y[i]) * phi1x[ig][j];
                    Jac21[i * nDof1 + j] -= jacValue21;
                    Res2[i] +=  jacValue21 * solu1[j];
                  }

                  for (unsigned j = 0; j < nDof2; j++) {
                    double jacValue22 = - weight1[ig] * weight2 * kernel * (- phi2y[i]) * phi2y[j];
                    Jac22[i * nDof2 + j] -= jacValue22;
                    Res2[i] +=  jacValue22 * solu2[j];
                  }//endl j loop
                } //endl i loop
              }//end jg loop
            }
          }//end ig loop

          if (ifAnyIntersection) {
            KK->add_matrix_blocked (Jac11, l2GMap1, l2GMap1);
            KK->add_matrix_blocked (Jac12, l2GMap1, l2GMap2);
            RES->add_vector_blocked (Res1, l2GMap1);

            KK->add_matrix_blocked (Jac21, l2GMap2, l2GMap1);
            KK->add_matrix_blocked (Jac22, l2GMap2, l2GMap2);
            RES->add_vector_blocked (Res2, l2GMap2);
          }
        }
      } //end iel loop
    } // end jel loop
  } //end kproc loop

  RES->close();

  KK->close();

//     Mat A = ( static_cast<PetscMatrix*> ( KK ) )->mat();
//     MatAssemblyBegin ( A, MAT_FINAL_ASSEMBLY );
//     MatAssemblyEnd ( A, MAT_FINAL_ASSEMBLY );
//     PetscViewer viewer;
//     MatView ( A, viewer );

//     Vec v = ( static_cast< PetscVector* > ( RES ) )->vec();
//     VecView(v,PETSC_VIEWER_STDOUT_WORLD);

  // ***************** END ASSEMBLY *******************
}

//END



void RectangleAndBallRelation (bool & theyIntersect, const std::vector<double> &ballCenter, const double & ballRadius, const std::vector < std::vector < double> > &elementCoordinates, std::vector < std::vector < double> > &newCoordinates) {

  //theyIntersect = true : element and ball intersect
  //theyIntersect = false : element and ball are disjoint

  //elementCoordinates are the coordinates of the vertices of the element

  theyIntersect = false; //by default we assume the two sets are disjoint

  unsigned dim = 2;
  unsigned nDofs = elementCoordinates[0].size();

  std::vector< std::vector < double > > ballVerticesCoordinates (dim);
  newCoordinates.resize (dim);


  for (unsigned n = 0; n < dim; n++) {
    newCoordinates[n].resize (nDofs);
    ballVerticesCoordinates[n].resize (4);

    for (unsigned i = 0; i < nDofs; i++) {
      newCoordinates[n][i] = elementCoordinates[n][i]; //this is just an initalization, it will be overwritten
    }
  }

  double xMinElem = elementCoordinates[0][0];
  double yMinElem = elementCoordinates[1][0];
  double xMaxElem = elementCoordinates[0][2];
  double yMaxElem = elementCoordinates[1][2];


  for (unsigned i = 0; i < 4; i++) {
    if (elementCoordinates[0][i] < xMinElem) xMinElem = elementCoordinates[0][i];

    if (elementCoordinates[0][i] > xMaxElem) xMaxElem = elementCoordinates[0][i];

    if (elementCoordinates[1][i] < yMinElem) yMinElem = elementCoordinates[1][i];

    if (elementCoordinates[1][i] > yMaxElem) yMaxElem = elementCoordinates[1][i];
  }

  //bottom left corner of ball (south west)
  ballVerticesCoordinates[0][0] =  ballCenter[0] - ballRadius;
  ballVerticesCoordinates[1][0] =  ballCenter[1] - ballRadius;

  //top right corner of ball (north east)
  ballVerticesCoordinates[0][2] = ballCenter[0] + ballRadius;
  ballVerticesCoordinates[1][2] = ballCenter[1] + ballRadius;

  newCoordinates[0][0] = (ballVerticesCoordinates[0][0] >= xMinElem) ? ballVerticesCoordinates[0][0] : xMinElem;
  newCoordinates[1][0] = (ballVerticesCoordinates[1][0] >= yMinElem) ? ballVerticesCoordinates[1][0] : yMinElem;

  newCoordinates[0][2] = (ballVerticesCoordinates[0][2] >= xMaxElem) ? xMaxElem : ballVerticesCoordinates[0][2];
  newCoordinates[1][2] = (ballVerticesCoordinates[1][2] >= yMaxElem) ? yMaxElem : ballVerticesCoordinates[1][2];

  if (newCoordinates[0][0] < newCoordinates[0][2] && newCoordinates[1][0] < newCoordinates[1][2]) {   //ball and rectangle intersect

    theyIntersect = true;

    newCoordinates[0][1] = newCoordinates[0][2];
    newCoordinates[1][1] = newCoordinates[1][0];

    newCoordinates[0][3] = newCoordinates[0][0];
    newCoordinates[1][3] = newCoordinates[1][2];

    if (nDofs > 4) {   //TODO the quadratic case has not yet been debugged

      newCoordinates[0][4] = 0.5 * (newCoordinates[0][0] + newCoordinates[0][1]);
      newCoordinates[1][4] = newCoordinates[1][0];

      newCoordinates[0][5] = newCoordinates[0][1];
      newCoordinates[1][5] = 0.5 * (newCoordinates[1][1] + newCoordinates[1][2]);

      newCoordinates[0][6] = newCoordinates[0][4];
      newCoordinates[1][6] = newCoordinates[1][2];

      newCoordinates[0][7] = newCoordinates[0][0];
      newCoordinates[1][7] = newCoordinates[1][5];

      if (nDofs > 8) {

        newCoordinates[0][8] = newCoordinates[0][4];
        newCoordinates[1][8] = newCoordinates[1][5];

      }

    }

  }

}

const unsigned swap[4][9] = {
  {0, 1, 2, 3, 4, 5, 6, 7, 8},
  {3, 0, 1, 2, 7, 4, 5, 6, 8},
  {2, 3, 0, 1, 6, 7, 4, 5, 8},
  {1, 2, 3, 0, 5, 6, 7, 4, 8}
};

void ReorderElement (std::vector < int > &dofs, std::vector < double > & sol, std::vector < std::vector < double > > & x) {

  unsigned type = 0;

  if (fabs (x[0][0] - x[0][1]) > 1.e-10) {
    if (x[0][0] - x[0][1] > 0) {
      type = 2;
    }
  }

  else {
    type = 1;

    if (x[1][0] - x[1][1] > 0) {
      type = 3;
    }
  }

  if (type != 0) {
    std::vector < int > dofsCopy = dofs;
    std::vector < double > solCopy = sol;
    std::vector < std::vector < double > > xCopy = x;

    for (unsigned i = 0; i < dofs.size(); i++) {
      dofs[i] = dofsCopy[swap[type][i]];
      sol[i] = solCopy[swap[type][i]];

      for (unsigned k = 0; k < x.size(); k++) {
        x[k][i] = xCopy[k][swap[type][i]];
      }
    }
  }
}


void RectangleAndBallRelation2 (bool & theyIntersect, const std::vector<double> &ballCenter, const double & ballRadius, const std::vector < std::vector < double> > &elementCoordinates, std::vector < std::vector < double> > &newCoordinates) {

  theyIntersect = false; //by default we assume the two sets are disjoint

  unsigned dim = 2;
  unsigned nDofs = elementCoordinates[0].size();

  newCoordinates.resize (dim);

  for (unsigned i = 0; i < dim; i++) {
    newCoordinates[i].resize (nDofs);
  }

  double xMin = elementCoordinates[0][0];
  double xMax = elementCoordinates[0][2];
  double yMin = elementCoordinates[1][0];
  double yMax = elementCoordinates[1][2];

  if (xMin > xMax || yMin > yMax) {
    std::cout << "error" << std::endl;

    for (unsigned i = 0; i < nDofs; i++) {
      std::cout <<  elementCoordinates[0][i] << " " << elementCoordinates[1][i] << std::endl;
    }

    exit (0);
  }


  double xMinBall = ballCenter[0] - ballRadius;
  double xMaxBall = ballCenter[0] + ballRadius;
  double yMinBall = ballCenter[1] - ballRadius;
  double yMaxBall = ballCenter[1] + ballRadius;


  xMin = (xMin > xMinBall) ? xMin : xMinBall;
  xMax = (xMax < xMaxBall) ? xMax : xMaxBall;
  yMin = (yMin > yMinBall) ? yMin : yMinBall;
  yMax = (yMax < yMaxBall) ? yMax : yMaxBall;

  if (xMin < xMax && yMin < yMax) {   //ball and rectangle intersect

    theyIntersect = true;

    //std::cout<< xMin <<" "<<xMax<<" "<<yMin<<" "<<yMax<<std::endl;

    newCoordinates[0][0] = xMin;
    newCoordinates[0][1] = xMax;
    newCoordinates[0][2] = xMax;
    newCoordinates[0][3] = xMin;

    newCoordinates[1][0] = yMin;
    newCoordinates[1][1] = yMin;
    newCoordinates[1][2] = yMax;
    newCoordinates[1][3] = yMax;

    if (nDofs > 4) {   //TODO the quadratic case has not yet been debugged

      double xMid = 0.5 * (xMin + xMax);
      double yMid = 0.5 * (yMin + yMax);

      newCoordinates[0][4] = xMid;
      newCoordinates[0][5] = xMax;
      newCoordinates[0][6] = xMid;
      newCoordinates[0][7] = xMin;

      newCoordinates[1][4] = yMin;
      newCoordinates[1][5] = yMid;
      newCoordinates[1][6] = yMax;
      newCoordinates[1][7] = yMid;

      if (nDofs > 8) {

        newCoordinates[0][8] = xMid;
        newCoordinates[1][8] = yMid;

      }

    }

  }

}


