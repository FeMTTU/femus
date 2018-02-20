


using namespace femus;

void AssembleUQSys(MultiLevelProblem& ml_prob) {

  // ml_prob is the global object from/to where get/set all the data
  // level is the level of the PDE system to be assembled
  // levelMax is the Maximum level of the MultiLevelProblem
  // assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  //pointers and references

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FEM");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = ml_sol->GetSolutionLevel(level);     // pointer to the solution (level) object
  LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];  // pointer to the equation (level) object

  Mesh* mymsh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* myel = mymsh->el;   // pointer to the elem object in msh (level)
  SparseMatrix* myKK = myLinEqSolver->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* myRES =  myLinEqSolver->_RES;  // pointer to the global residual vector object in pdeSys (level)
  bool assembleMatrix = my_nnlin_impl_sys.GetAssembleMatrix();

// call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;
  if(assembleMatrix) s.continue_recording();
  else s.pause_recording();

  const unsigned dim = mymsh->GetDimension();

  // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  // data
  unsigned iproc  = mymsh->processor_id();

  vector < double > phi;
  vector < double > phi_hat;
  vector < adept::adouble> gradphi;
  vector < double > gradphi_hat;

  phi.reserve(max_size);
  phi_hat.reserve(max_size);

  gradphi.reserve(max_size * dim);
  gradphi_hat.reserve(max_size * dim);

  vector <vector < adept::adouble> > vx(dim); //vx is coordX in assembly of ex30
  vector <vector < double> > vx_hat(dim);

  vector< vector< adept::adouble > > SolUd(dim);      // local solution (displacement)

  vector< vector< int > > dofsVAR(dim);

  vector< vector< double > > Rhs(dim);     // local redidual vector
  vector< vector< adept::adouble > > aRhs(dim);     // local redidual vector

  vector < int > dofsAll;

  vector < double > Jac;

  adept::adouble weight;
  double weight_hat = 0.;

  vector < adept::adouble >* nullAdoublePointer = NULL;
  vector < double >* nullDoublePointer = NULL;

  //variable-name handling
  const char varname[3][5] = {"UX", "UY", "UW"};
  vector <unsigned> indexSolU(dim);
  vector <unsigned> indexPdeU(dim);
  unsigned solType = ml_sol->GetSolutionType(&varname[0][0]);


  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolU[ivar] = ml_sol->GetIndex(&varname[ivar][0]);
    indexPdeU[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);
  }

  start_time = clock();

  if(assembleMatrix) myKK->zero();


  for(int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = mymsh->GetElementType(iel);

    int material = mymsh->GetElementMaterial(iel);

    unsigned nDofsU = mymsh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    unsigned nDofs = dim * nDofsU ;//+ nDofsP;
    // resize local arrays
    std::vector <int> sysDof(nDofs);


    for(unsigned  k = 0; k < dim; k++) {
      SolUd[k].resize(nDofsU);
      vx[k].resize(nDofsU);
      vx_hat[k].resize(nDofsU);
    }


    for(unsigned  k = 0; k < dim; k++) {
      aRhs[k].resize(nDofsU);    //resize
      std::fill(aRhs[k].begin(), aRhs[k].end(), 0);    //set aRes to zero
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsU; i++) {
      unsigned idof = mymsh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
      unsigned idofX = mymsh->GetSolutionDof(i, iel, 2);    // global to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        SolUd[k][i] = (*mysolution->_Sol[indexSolU[k]])(idof);      // global extraction and local storage for the solution
        sysDof[i + k * nDofsU] = myLinEqSolver->GetSystemDof(indexSolU[k], indexPdeU[k], i, iel);    // global to global mapping between solution node and pdeSys dof
        vx_hat[k][i] = (*mymsh->_topology->_Sol[k])(idofX);
        vx[k][i] = vx_hat[k][i] + SolUd[k][i];
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    if(assembleMatrix) s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < mymsh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {


      mymsh->_finiteElement[ielt][solType]->Jacobian(vx, ig, weight, phi, gradphi, *nullAdoublePointer);
      mymsh->_finiteElement[ielt][solType]->Jacobian(vx_hat, ig, weight_hat, phi_hat, gradphi_hat, *nullDoublePointer);

      vector < adept::adouble > SolUgss(dim, 0);

      vector < vector < adept::adouble > > GradSolUgss(dim);
      vector < vector < adept::adouble > > GradSolUgssHat(dim);

      for(unsigned  k = 0; k < dim; k++) {
        GradSolUgssHat[k].resize(dim);
        std::fill(GradSolUgssHat[k].begin(), GradSolUgssHat[k].end(), 0);
      }

      for(unsigned i = 0; i < nDofsU; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            GradSolUgssHat[k][j] += gradphi_hat[i * dim + j] * SolUd[k][i];
          }
        }
      }


      for(unsigned  k = 0; k < dim; k++) {
        GradSolUgss[k].resize(dim);
        std::fill(GradSolUgss[k].begin(), GradSolUgss[k].end(), 0);
      }

      for(unsigned i = 0; i < nDofsU; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          SolUgss[k] += phi[i] * SolUd[k][i];
        }

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            GradSolUgss[k][j] += gradphi[i * dim + j] * SolUd[k][i];
          }
        }
      }

      adept::adouble F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      adept::adouble B[3][3];
      adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
      adept::adouble Cauchy[3][3];

      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          F[i][j] += GradSolUgssHat[i][j];
        }
      }

      adept::adouble J_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                              - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          B[i][j] = 0.;

          for(int k = 0; k < 3; k++) {
            //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
            B[i][j] += F[i][k] * F[j][k];
          }
        }
      }

      adept::adouble I1_B = B[0][0] + B[1][1] + B[2][2];

      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          //  Cauchy[i][j] = mu * (B[i][j] - I1_B * Id2th[i][j] / 3.) / pow(J_hat, 5. / 3.)
          //                 + K * (J_hat - 1.) * Id2th[i][j];  //Generalized Neo-Hookean solid, in Allan-Bower's book, for rubbers with very limited compressibility and K >> mu

          Cauchy[i][j] = 1.e5 * log(J_hat) / J_hat * Id2th[i][j] + 1.e5 / J_hat * (B[i][j] - Id2th[i][j]); //alternative formulation

        }
      }

      for(unsigned i = 0; i < nDofsU; i++) {
        adept::adouble CauchyDIR[3] = {0., 0., 0.};

        for(int idim = 0.; idim < dim; idim++) {
          for(int jdim = 0.; jdim < dim; jdim++) {
            CauchyDIR[idim] += gradphi[i * dim + jdim] * Cauchy[idim][jdim];
          }
        }

        for(int idim = 0; idim < dim; idim++) {
          aRhs[indexPdeU[idim]][i] += (phi[i] * 1000. / J_hat  - CauchyDIR[idim]
                                       - phi[i] * (-1. / (0.5) * SolUgss[idim]
                                           - (1. - 0.5) / 2. * SolUgss[idim])
                                      ) * weight;
        }
      }
    } // end gauss point loop


    //copy the value of the adept::adoube aRes in double Res and store them in RES
    std::vector<double> Rhs(nDofs);  //resize

    for(int i = 0; i < nDofsU; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        Rhs[ i +  k * nDofsU ] = -aRhs[k][i].value();
      }
    }

    myRES->add_vector_blocked(Rhs, sysDof);

    if(assembleMatrix) {
      Jac.resize(nDofs * nDofs);
      // define the dependent variables

      for(unsigned  k = 0; k < dim; k++) {
        s.dependent(&aRhs[k][0], nDofsU);
      }

      // define the independent variables
      for(unsigned  k = 0; k < dim; k++) {
        s.independent(&SolUd[k][0], nDofsU);
      }

      // get the and store jacobian matrix (row-major)
      s.jacobian(&Jac[0] , true);
      myKK->add_matrix_blocked(Jac, sysDof, sysDof);

      s.clear_independents();
      s.clear_dependents();
    }
  }
  //END building "soft" stiffness matrix


  myRES->close();

  if(assembleMatrix) {
    myKK->close();
  }

  // *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
  // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}
