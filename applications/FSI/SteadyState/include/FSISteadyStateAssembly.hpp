#ifndef __femus_include_FSISteadyStateAssembly_hpp__
#define __femus_include_FSISteadyStateAssembly_hpp__
#endif

#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "MultiLevelSolution.hpp"
#include "adept.h"


namespace femus {

  void FSISteadyStateAssembly(MultiLevelProblem& ml_prob) {

    clock_t AssemblyTime = 0;
    clock_t start_time, end_time;

    //pointers and references
    MonolithicFSINonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system < MonolithicFSINonLinearImplicitSystem>("Fluid-Structure-Interaction");
    const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();

    MultiLevelSolution*	        ml_sol	        = ml_prob._ml_sol;
    Solution*	                mysolution      = ml_sol->GetSolutionLevel(level);
    LinearEquationSolver*       myLinEqSolver   = my_nnlin_impl_sys._LinSolver[level];
    Mesh*		        mymsh		= ml_prob._ml_msh->GetLevel(level);
    elem*		        myel		= mymsh->el;
    SparseMatrix*	        myKK		= myLinEqSolver->_KK;
    NumericVector*	        myRES		= myLinEqSolver->_RES;

    bool assembleMatrix = my_nnlin_impl_sys.GetAssembleMatrix();
    adept::Stack& s = FemusInit::_adeptStack;
    if( assembleMatrix ) s.continue_recording();
    else s.pause_recording();

    const unsigned dim = mymsh->GetDimension();
    const unsigned max_size = static_cast < unsigned >(ceil(pow(3, dim)));

    // local objects
    vector < adept::adouble> SolVAR(2 * dim + 1);
    vector < vector < adept::adouble> > GradSolVAR(2 * dim);
    vector < vector < adept::adouble> > GradSolhatVAR(2 * dim);

    vector < vector < adept::adouble> > NablaSolVAR(2 * dim);
    vector < vector < adept::adouble> > NablaSolhatVAR(2 * dim);

    for (int i = 0; i < 2 * dim; i++) {
      GradSolVAR[i].resize(dim);
      GradSolhatVAR[i].resize(dim);

      NablaSolVAR[i].resize(3 * (dim - 1));
      NablaSolhatVAR[i].resize(3 * (dim - 1));
    }

    vector  < bool> solidmark;
    vector  < double > phi;
    vector  < double > phi_hat;
    vector  < adept::adouble> gradphi;
    vector  < double> gradphi_hat;
    vector  < adept::adouble> nablaphi;
    vector  < double> nablaphi_hat;

    phi.reserve(max_size);
    solidmark.reserve(max_size);
    phi_hat.reserve(max_size);
    gradphi.reserve(max_size * dim);
    gradphi_hat.reserve(max_size * dim);
    nablaphi.reserve(max_size * 3 * (dim - 1));
    nablaphi_hat.reserve(max_size * 3 * (dim - 1));

    const double* phi1;

    adept::adouble jacobian = 0.;
    double jacobianOverArea = 0.;
    double jacobian_hat = 0.;

    vector  < vector  <  adept::adouble> > vx(dim);
    vector  < vector  <  adept::adouble> > vx_face(dim);
    vector  < vector  <  double> > vx_hat(dim);

    for (int i = 0; i < dim; i++) {
      vx[i].reserve(max_size);
      vx_face[i].resize(9);
      vx_hat[i].reserve(max_size);
    }

    vector <  vector <  adept::adouble > > Soli(2 * dim + 1);
    vector <  vector <  int > > dofsVAR(2 * dim + 1);

    for (int i = 0; i < 2 * dim + 1; i++) {
      Soli[i].reserve(max_size);
      dofsVAR[i].reserve(max_size);
    }

    vector <  vector <  double > > Rhs(2 * dim + 1);
    vector <  vector <  adept::adouble > > aRhs(2 * dim + 1);

    for (int i = 0; i < 2 * dim + 1; i++) {
      aRhs[i].reserve(max_size);
      Rhs[i].reserve(max_size);
    }

    vector  <  int > dofsAll;
    dofsAll.reserve(max_size * (2 + dim + 1));

    vector  <  double > Jac;
    Jac.reserve(dim * max_size * (2 * dim + 1)*dim * max_size * (2 * dim + 1));

    // ------------------------------------------------------------------------
    // Physical parameters
    double rhof	 	= ml_prob.parameters.get < Fluid>("Fluid").get_density();
    double mu_lame 	= ml_prob.parameters.get < Solid>("Solid").get_lame_shear_modulus();
    double lambda_lame 	= ml_prob.parameters.get < Solid>("Solid").get_lame_lambda();
    double mus		= mu_lame / rhof;
    double IRe 		= ml_prob.parameters.get < Fluid>("Fluid").get_IReynolds_number();
    double lambda	= lambda_lame / rhof;
    double betans	= 1.;
    int    solid_model	= ml_prob.parameters.get < Solid>("Solid").get_physical_model();

    bool incompressible = (0.5  ==  ml_prob.parameters.get < Solid>("Solid").get_poisson_coeff()) ? 1 : 0;
    const bool penalty = ml_prob.parameters.get < Solid>("Solid").get_if_penalty();

    // gravity
    double _gravity[3] = {0., 0., 0.};

    // -----------------------------------------------------------------
    // space discretization parameters
    unsigned SolType2 = ml_sol->GetSolutionType(ml_sol->GetIndex("U"));
    unsigned SolType1 = ml_sol->GetSolutionType(ml_sol->GetIndex("P"));

    // mesh and procs
    unsigned nel    = mymsh->GetNumberOfElements();
    unsigned igrid  = mymsh->GetLevel();
    unsigned iproc  = mymsh->processor_id();

    //----------------------------------------------------------------------------------
    //variable-name handling
    const char varname[7][3] = {"DX", "DY", "DZ", "U", "V", "W", "P"};
    vector  < unsigned> indexVAR(2 * dim + 1);
    vector  < unsigned> indVAR(2 * dim + 1);
    vector  < unsigned> SolType(2 * dim + 1);

    for (unsigned ivar = 0; ivar < dim; ivar++) {
      indVAR[ivar] = ml_sol->GetIndex(&varname[ivar][0]);
      indVAR[ivar + dim] = ml_sol->GetIndex(&varname[ivar + 3][0]);
      SolType[ivar] = ml_sol->GetSolutionType(&varname[ivar][0]);
      SolType[ivar + dim] = ml_sol->GetSolutionType(&varname[ivar + 3][0]);
      indexVAR[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);
      indexVAR[ivar + dim] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar + 3][0]);
    }

    indexVAR[2 * dim] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[6][0]);
    indVAR[2 * dim] = ml_sol->GetIndex(&varname[6][0]);
    SolType[2 * dim] = ml_sol->GetSolutionType(&varname[6][0]);
    //----------------------------------------------------------------------------------

    int nprocs = mymsh->n_processors();
    NumericVector* area_elem_first;
    area_elem_first = NumericVector::build().release();

    if (nprocs == 1) {
      area_elem_first->init(nprocs, 1, false, SERIAL);
    }
    else {
      area_elem_first->init(nprocs, 1, false, PARALLEL);
    }

    area_elem_first->zero();
    double rapresentative_area = 1.;

    start_time = clock();

    if( assembleMatrix ) myKK->zero();

    // *** element loop ***
    for (int iel = mymsh->_elementOffset[iproc]; iel  <  mymsh->_elementOffset[iproc + 1]; iel++) {

      short unsigned ielt = mymsh->GetElementType(iel);
      unsigned nve        = mymsh->GetElementDofNumber(iel, SolType2);
      unsigned nve1       = mymsh->GetElementDofNumber(iel, SolType1);
      int flag_mat        = mymsh->GetElementMaterial(iel);

      // *******************************************************************************************************

      //initialization of everything in common between fluid and solid

      //Rhs
      for (int i = 0; i < 2 * dim; i++) {
        dofsVAR[i].resize(nve);
        Soli[indexVAR[i]].resize(nve);
        aRhs[indexVAR[i]].resize(nve);
      }

      dofsVAR[2 * dim].resize(nve1);
      Soli[indexVAR[2 * dim]].resize(nve1);
      aRhs[indexVAR[2 * dim]].resize(nve1);

      dofsAll.resize(0);

      // ----------------------------------------------------------------------------------------
      // coordinates, solutions, displacement and velocity dofs

      solidmark.resize(nve);

      for (int i = 0; i < dim; i++) {
        vx[i].resize(nve);
        vx_hat[i].resize(nve);
      }

      for (unsigned i = 0; i < nve; i++) {
        unsigned iDof = mymsh->GetSolutionDof(i, iel, SolType2);

        // flag nodes on the fluid-solid interface
        solidmark[i] = mymsh->GetSolidMark(iDof);

        for (int j = 0; j < dim; j++) {
          Soli[indexVAR[j]][i] = (*mysolution->_Sol[indVAR[j]])(iDof);
          Soli[indexVAR[j + dim]][i] = (*mysolution->_Sol[indVAR[j + dim]])(iDof);

          aRhs[indexVAR[j]][i] = 0.;
          aRhs[indexVAR[j + dim]][i] = 0.;

          //Fixed coordinates (Reference frame)
          vx_hat[j][i] = (*mymsh->_topology->_Sol[j])(iDof);
          // displacement dofs
          dofsVAR[j][i] = myLinEqSolver->GetSystemDof(indVAR[j], indexVAR[j], i, iel);
          // velocity dofs
          dofsVAR[j + dim][i] = myLinEqSolver->GetSystemDof(indVAR[j + dim], indexVAR[j + dim], i, iel);
        }
      }

      // pressure dofs
      for (unsigned i = 0; i < nve1; i++) {
        unsigned iDof = mymsh->GetSolutionDof(i, iel, SolType[2 * dim]);
        dofsVAR[2 * dim][i] = myLinEqSolver->GetSystemDof(indVAR[2 * dim], indexVAR[2 * dim], i, iel);
        Soli[indexVAR[2 * dim]][i] = (*mysolution->_Sol[indVAR[2 * dim]])(iDof);
        aRhs[indexVAR[2 * dim]][i] = 0.;
      }

      // compose the system dofs
      for (int idim = 0; idim < 2 * dim; idim++) {
        dofsAll.insert(dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end());
      }

      dofsAll.insert(dofsAll.end(), dofsVAR[2 * dim].begin(), dofsVAR[2 * dim].end());

      if( assembleMatrix ) s.new_recording();

      //Moving coordinates (Moving frame)
      for (unsigned idim = 0; idim < dim; idim++) {
        for (int j = 0; j < nve; j++) {
          vx[idim][j] = vx_hat[idim][j] + Soli[indexVAR[idim]][j];
        }
      }

      // Boundary integral
      {
        double tau = 0.;
        vector < adept::adouble> normal(dim, 0);

        // loop on faces
        for (unsigned jface = 0; jface < mymsh->GetElementFaceNumber(iel); jface++) {
          std::vector  <  double > xx(3, 0.);

          // look for boundary faces
          if (myel->GetFaceElementIndex(iel, jface) < 0) {
            unsigned int face = -(mymsh->el->GetFaceElementIndex(iel, jface) + 1);

            if ( !ml_sol->GetBdcFunction()(xx, "U", tau, face, 0.) && tau != 0.) {
                
              unsigned nve = mymsh->GetElementFaceDofNumber(iel, jface, SolType2);
              const unsigned felt = mymsh->GetElementFaceType(iel, jface);

              for (unsigned i = 0; i < nve; i++) {
                unsigned int ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);
                unsigned iDof = mymsh->GetSolutionDof(ilocal, iel, 2);

                for (unsigned idim = 0; idim < dim; idim++) {
                  vx_face[idim][i] = (*mymsh->_topology->_Sol[idim])(iDof) + Soli[indexVAR[idim]][ilocal];
                }
              }

              for (unsigned igs = 0; igs  <  mymsh->_finiteElement[felt][SolType2]->GetGaussPointNumber(); igs++) {
                mymsh->_finiteElement[felt][SolType2]->JacobianSur(vx_face, igs, jacobian, phi, gradphi, normal);

                // *** phi_i loop ***
                for (unsigned i = 0; i < nve; i++) {
                  adept::adouble value = -phi[i] * tau / rhof * jacobian;
                  unsigned int ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);

                  for (unsigned idim = 0; idim < dim; idim++) {
                    if ((!solidmark[ilocal])) { //if fluid node it goes to U, V, W
                      aRhs[indexVAR[dim + idim]][ilocal] +=  value * normal[idim];
                    }
                    else { //if solid node it goes to DX, DY, DZ
                      aRhs[indexVAR[idim]][ilocal] +=  value * normal[idim];
                    }
                  }
                }
              }
            }
          }
        }
      }

      // *** Gauss point loop ***
      double area = 1.;
      adept::adouble supg_tau;

      for (unsigned ig = 0; ig  <  mymsh->_finiteElement[ielt][SolType2]->GetGaussPointNumber(); ig++) {
        // *** get Jacobian and test function and test function derivatives in the moving frame***
        mymsh->_finiteElement[ielt][SolType2]->Jacobian(vx, ig, jacobian, phi, gradphi, nablaphi);
        mymsh->_finiteElement[ielt][SolType2]->Jacobian(vx_hat, ig, jacobian_hat, phi_hat, gradphi_hat, nablaphi_hat);
        phi1 = mymsh->_finiteElement[ielt][SolType1]->GetPhi(ig);

        if (flag_mat == 2 || iel  ==  mymsh->_elementOffset[iproc]) {
          if (ig  ==  0) {
            double GaussWeight = mymsh->_finiteElement[ielt][SolType2]->GetGaussWeight(ig);
            area = jacobian_hat / GaussWeight;

            if (iel == mymsh->_elementOffset[iproc]) {
              area_elem_first->add(mymsh->processor_id(), area);
              area_elem_first->close();
              rapresentative_area = area_elem_first->l1_norm() / nprocs;
            }
          }

          jacobianOverArea = jacobian_hat / area * rapresentative_area;
        }

        // ---------------------------------------------------------------------------
        // displacement and velocity
        for (int i = 0; i < 2 * dim; i++) {
          SolVAR[i] = 0.;

          for (int j = 0; j < dim; j++) {
            GradSolVAR[i][j] = 0.;
            GradSolhatVAR[i][j] = 0.;
          }

          for (int j = 0; j < 3 * (dim - 1); j++) {
            NablaSolVAR[i][j] = 0.;
            NablaSolhatVAR[i][j] = 0.;
          }

          for (unsigned inode = 0; inode < nve; inode++) {
            SolVAR[i] +=  phi[inode] * Soli[indexVAR[i]][inode];

            for (int j = 0; j < dim; j++) {
              GradSolVAR[i][j] +=  gradphi[inode * dim + j] * Soli[indexVAR[i]][inode];
              GradSolhatVAR[i][j] +=  gradphi_hat[inode * dim + j] * Soli[indexVAR[i]][inode];
            }

            for (int j = 0; j < 3 * (dim - 1); j++) {
              NablaSolVAR[i][j] +=  nablaphi[inode * 3 * (dim - 1) + j] * Soli[indexVAR[i]][inode];
              NablaSolhatVAR[i][j] +=  nablaphi_hat[inode * 3 * (dim - 1) + j] * Soli[indexVAR[i]][inode];
            }
          }
        }

        // pressure
        SolVAR[2 * dim] = 0.;

        for (unsigned inode = 0; inode < nve1; inode++) {
          adept::adouble soli = Soli[indexVAR[2 * dim]][inode];
          SolVAR[2 * dim] += phi1[inode] * soli;
        }

        // ---------------------------------------------------------------------------
        //BEGIN FLUID ASSEMBLY
        if (flag_mat == 2) {
          //BEGIN ALE + Momentum (Navier-Stokes)
          {
            for (unsigned i = 0; i < nve; i++) {

              //BEGIN redidual Laplacian ALE map in the reference domain
              adept::adouble LapmapVAR[3] = {0., 0., 0.};

              for (int idim = 0; idim < dim; idim++) {
                for (int jdim = 0; jdim < dim; jdim++) {
                  LapmapVAR[idim] += (GradSolhatVAR[idim][jdim] * gradphi_hat[i * dim + jdim]) ;
                }
              }

              for (int idim = 0; idim < dim; idim++) {
                aRhs[indexVAR[idim]][i] += (!solidmark[i]) * (-LapmapVAR[idim] * jacobianOverArea);
              }

              //END redidual Laplacian ALE map in reference domain

              //BEGIN redidual Navier-Stokes in moving domain
              adept::adouble LapvelVAR[3] = {0., 0., 0.};
              adept::adouble AdvaleVAR[3] = {0., 0., 0.};

              for (int idim = 0.; idim < dim; idim++) {
                for (int jdim = 0.; jdim < dim; jdim++) {
                  //LapvelVAR[idim] +=  GradSolVAR[dim+idim][jdim]*gradphi[i*dim+jdim];
                  LapvelVAR[idim] += (GradSolVAR[dim + idim][jdim] + GradSolVAR[dim + jdim][idim]) * gradphi[i * dim + jdim];
                  AdvaleVAR[idim] +=  SolVAR[dim + jdim] * GradSolVAR[dim + idim][jdim] * phi[i];
                }
              }

              for (int idim = 0; idim < dim; idim++) {
                adept::adouble value = (-AdvaleVAR[idim]      	           // advection term
                                        - IRe * LapvelVAR[idim]	   	 // viscous dissipation
                                        + SolVAR[2 * dim] * gradphi[i * dim + idim] // pressure gradient
                                       ) * jacobian;

                if ((!solidmark[i])) {
                  aRhs[indexVAR[dim + idim]][i] += value;
                }
                else {
                  aRhs[indexVAR[idim]][i] +=  value;
                }

                //END redidual Navier-Stokes in moving domain
              }
            }
          }
          //END ALE + Momentum (Navier-Stokes)

          //BEGIN continuity block
          {
            adept::adouble div_vel = 0.;

            for (int i = 0; i < dim; i++) {
              div_vel += GradSolVAR[dim + i][i];
            }

            for (unsigned i = 0; i < nve1; i++) {
              aRhs[indexVAR[2 * dim]][i] += -(-phi1[i] * div_vel) * jacobian;
            }
          }
          //END continuity block
        }
        //END FLUID ASSEMBLY

        //*******************************************************************************************************

        //BEGIN SOLID ASSEMBLY
        else {
            
            
          //BEGIN build Chauchy Stress in moving domain
          //physical quantity
          adept::adouble J_hat;
          adept::adouble I_e;
          adept::adouble Cauchy[3][3];
          adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};

          adept::adouble I1_B = 0.;
          adept::adouble I2_B = 0.;

          if (solid_model == 0) { // Saint-Venant
            adept::adouble e[3][3];

            //computation of the stress tensor
            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                e[i][j] = 0.5 * (GradSolhatVAR[i][j] + GradSolhatVAR[j][i]);
              }
            }

            I_e = 0;

            for (int i = 0; i < dim; i++) {
              I_e += e[i][i];
            }

            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                //incompressible
                Cauchy[i][j] = 2 * mus * e[i][j] - 2. * (incompressible) * mus * I_e * SolVAR[2 * dim] * Id2th[i][j];
                //+(penalty)*lambda*I_e*Id2th[i][j];
              }
            }
          }

          else { // hyperelastic non linear material
            adept::adouble F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
            adept::adouble B[3][3];

            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                F[i][j] += GradSolhatVAR[i][j];
              }
            }

            J_hat =   F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                      - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

            for (int I = 0; I < 3; ++I) {
              for (int J = 0; J < 3; ++J) {
                B[I][J] = 0.;

                for (int K = 0; K < 3; ++K) {
                  //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
                  B[I][J] += F[I][K] * F[J][K];
                }
              }
            }

            if (solid_model <= 4) { // Neo-Hookean
              I1_B = B[0][0] + B[1][1] + B[2][2];

              for (int I = 0; I < 3; ++I) {
                for (int J = 0; J < 3; ++J) {
                  if (1  ==  solid_model) Cauchy[I][J] = mus * B[I][J]
                        - (incompressible) * mus * I1_B * SolVAR[2 * dim] * Id2th[I][J]; 	//Wood-Bonet J_hat  =1;
                  else if (2  ==  solid_model) Cauchy[I][J] = mus / J_hat * B[I][J]
                        - (incompressible) * mus / J_hat * SolVAR[2 * dim] * Id2th[I][J]; //Wood-Bonet J_hat !=1;
                  else if (3  ==  solid_model) Cauchy[I][J] = mus * (B[I][J] - Id2th[I][J]) / J_hat
                        + lambda / J_hat * log(J_hat) * Id2th[I][J]; 	//Wood-Bonet penalty
                  else if (4  ==  solid_model) Cauchy[I][J] = mus * (B[I][J] - I1_B * Id2th[I][J] / 3.) / pow(J_hat, 5. / 3.)
                        + lambda * (J_hat - 1.) * Id2th[I][J];  	 //Allan-Bower

                }
              }
            }
            else if (5  ==  solid_model) {  //Mooney-Rivlin
              adept::adouble detB =   B[0][0] * (B[1][1] * B[2][2] - B[2][1] * B[1][2])
                                      - B[0][1] * (B[2][2] * B[1][0] - B[1][2] * B[2][0])
                                      + B[0][2] * (B[1][0] * B[2][1] - B[2][0] * B[1][1]);
              adept::adouble invdetB = 1. / detB;
              adept::adouble invB[3][3];

              invB[0][0] = (B[1][1] * B[2][2] - B[2][1] * B[1][2]) * invdetB;
              invB[1][0] = -(B[0][1] * B[2][2] - B[0][2] * B[2][1]) * invdetB;
              invB[2][0] = (B[0][1] * B[1][2] - B[0][2] * B[1][1]) * invdetB;
              invB[0][1] = -(B[1][0] * B[2][2] - B[1][2] * B[2][0]) * invdetB;
              invB[1][1] = (B[0][0] * B[2][2] - B[0][2] * B[2][0]) * invdetB;
              invB[2][1] = -(B[0][0] * B[1][2] - B[1][0] * B[0][2]) * invdetB;
              invB[0][2] = (B[1][0] * B[2][1] - B[2][0] * B[1][1]) * invdetB;
              invB[1][2] = -(B[0][0] * B[2][1] - B[2][0] * B[0][1]) * invdetB;
              invB[2][2] = (B[0][0] * B[1][1] - B[1][0] * B[0][1]) * invdetB;

              I1_B = B[0][0] + B[1][1] + B[2][2];
              I2_B = B[0][0] * B[1][1] + B[1][1] * B[2][2] + B[2][2] * B[0][0]
                     - B[0][1] * B[1][0] - B[1][2] * B[2][1] - B[2][0] * B[0][2];

              double C1 = mus / 3.;
              double C2 = C1 / 2.;

              for (int I = 0; I < 3; ++I) {
                for (int J = 0; J < 3; ++J) {
                  Cauchy[I][J] =  2.*(C1 * B[I][J] - C2 * invB[I][J])
                                  //- (2. / 3.) * (C1 * I1_B - C2 * I2_B) * SolVAR[2 * dim] * Id2th[I][J];
                                  - (incompressible) * SolVAR[2 * dim] * Id2th[I][J];
                }
              }

            }
          }

          //END build Chauchy Stress in moving domain

          //BEGIN v=0 + Momentum (Solid)
          {
            for (unsigned i = 0; i < nve; i++) {

              //BEGIN redidual v=0 in fixed domain
              for (int idim = 0; idim < dim; idim++) {
                aRhs[indexVAR[dim + idim]][i] += (-phi[i] * (-SolVAR[dim + idim])) * jacobian_hat;
              }

              //END redidual v=0 in fixed domain

              //BEGIN redidual Solid Momentum in moving domain
              adept::adouble CauchyDIR[3] = {0., 0., 0.};

              for (int idim = 0.; idim < dim; idim++) {
                for (int jdim = 0.; jdim < dim; jdim++) {
                  CauchyDIR[idim] += gradphi[i * dim + jdim] * Cauchy[idim][jdim];
                }
              }

              for (int idim = 0; idim < dim; idim++) {
                aRhs[indexVAR[idim]][i] += (phi[i] * _gravity[idim] - CauchyDIR[idim]) * jacobian;
              }

              //END redidual Solid Momentum in moving domain
            }
          }
          //END v=0 + Momentum (Solid)

          //BEGIN continuity block
          {
            for (unsigned i = 0; i < nve1; i++) {
                
              if (!penalty) {
                  
                if (0  ==  solid_model) {
                  aRhs[indexVAR[2 * dim]][i] += -(-phi1[i] * (I_e + (!incompressible) / lambda * SolVAR[2 * dim])) * jacobian_hat;
                }
                else if (1  ==  solid_model || 5  ==  solid_model) {
                  if (incompressible)  aRhs[indexVAR[2 * dim]][i] += phi1[i] * (J_hat - 1. + (!incompressible) / lambda * SolVAR[2 * dim]) * jacobian_hat;
//                   else                 aRhs[indexVAR[2 * dim]][i] +=  /*phi1[i]*/ /** 1.e10 */ SolVAR[2 * dim];
                  else                 aRhs[indexVAR[2 * dim]][i] =  Soli[indexVAR[2 * dim]][i];
                }
                else if (2  ==  solid_model) {
                  aRhs[indexVAR[2 * dim]][i] +=  -(-phi1[i] * (log(J_hat) / J_hat + (!incompressible) / lambda * SolVAR[2 * dim])) * jacobian_hat;
                }
                
              }
              else {
                if (3  ==  solid_model || 4  ==  solid_model) { // pressure = 0 in the solid
                   aRhs[indexVAR[2 * dim]][i] += - ( -phi1[i] * (SolVAR[2 * dim])) * jacobian_hat;
                }
                  
              }
            }
          }
          //END continuity block
        }

        //END SOLID ASSEMBLY
      }

      //BEGIN local to global assembly
      //copy adouble aRhs into double Rhs
      for (unsigned i = 0; i < 2 * dim; i++) {
        Rhs[indexVAR[i]].resize(nve);

        for (int j = 0; j < nve; j++) {
          Rhs[indexVAR[i]][j] = -aRhs[indexVAR[i]][j].value();
        }
      }

      Rhs[indexVAR[2 * dim]].resize(nve1);

      for (unsigned j = 0; j < nve1; j++) {
        Rhs[indexVAR[2 * dim]][j] = -aRhs[indexVAR[2 * dim]][j].value();
      }

      for (int i = 0; i < 2 * dim + 1; i++) {
        myRES->add_vector_blocked(Rhs[indexVAR[i]], dofsVAR[i]);
      }

      if( assembleMatrix ){
        //Store equations
        for (int i = 0; i < 2 * dim; i++) {
          s.dependent(&aRhs[indexVAR[i]][0], nve);
          s.independent(&Soli[indexVAR[i]][0], nve);
        }

        s.dependent(&aRhs[indexVAR[2 * dim]][0], nve1);
        s.independent(&Soli[indexVAR[2 * dim]][0], nve1);

        Jac.resize((2 * dim * nve + nve1) * (2 * dim * nve + nve1));

        s.jacobian(&Jac[0], true);

        myKK->add_matrix_blocked(Jac, dofsAll, dofsAll);
        s.clear_independents();
        s.clear_dependents();

      //END local to global assembly
      }
    } //end list of elements loop

    if( assembleMatrix ) myKK->close();
    myRES->close();

    delete area_elem_first;

    // *************************************
    end_time = clock();
    AssemblyTime += (end_time - start_time);
    // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

  }



  void FSISteadyStateAssemblyWithNoPivoting(MultiLevelProblem& ml_prob) {

    clock_t AssemblyTime = 0;
    clock_t start_time, end_time;

    //pointers and references
    MonolithicFSINonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system < MonolithicFSINonLinearImplicitSystem>("Fluid-Structure-Interaction");
    const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();

    MultiLevelSolution*         ml_sol          = ml_prob._ml_sol;
    Solution*                   mysolution      = ml_sol->GetSolutionLevel(level);
    LinearEquationSolver*       myLinEqSolver   = my_nnlin_impl_sys._LinSolver[level];
    Mesh*                       mymsh           = ml_prob._ml_msh->GetLevel(level);
    elem*                       myel            = mymsh->el;
    SparseMatrix*               myKK            = myLinEqSolver->_KK;
    NumericVector*              myRES           = myLinEqSolver->_RES;

    bool assembleMatrix = my_nnlin_impl_sys.GetAssembleMatrix();
    adept::Stack& s = FemusInit::_adeptStack;
    if( assembleMatrix ) s.continue_recording();
    else s.pause_recording();

    const unsigned dim = mymsh->GetDimension();
    const unsigned max_size = static_cast < unsigned >(ceil(pow(3, dim)));

    // local objects
    vector < adept::adouble> SolVAR(2 * dim + 1);
    vector < vector < adept::adouble> > GradSolVAR(2 * dim);
    vector < vector < adept::adouble> > GradSolhatVAR(2 * dim);

    vector < vector < adept::adouble> > NablaSolVAR(2 * dim);
    vector < vector < adept::adouble> > NablaSolhatVAR(2 * dim);

    for (int i = 0; i < 2 * dim; i++) {
      GradSolVAR[i].resize(dim);
      GradSolhatVAR[i].resize(dim);

      NablaSolVAR[i].resize(3 * (dim - 1));
      NablaSolhatVAR[i].resize(3 * (dim - 1));
    }

    vector  < bool> solidmark;
    vector  < double > phi;
    vector  < double > phi_hat;
    vector  < adept::adouble> gradphi;
    vector  < double> gradphi_hat;
    vector  < adept::adouble> nablaphi;
    vector  < double> nablaphi_hat;

    phi.reserve(max_size);
    solidmark.reserve(max_size);
    phi_hat.reserve(max_size);
    gradphi.reserve(max_size * dim);
    gradphi_hat.reserve(max_size * dim);
    nablaphi.reserve(max_size * 3 * (dim - 1));
    nablaphi_hat.reserve(max_size * 3 * (dim - 1));

    const double* phi1;

    adept::adouble jacobian = 0.;
    double jacobianOverArea = 0.;
    double jacobian_hat = 0.;

    vector  < vector  <  adept::adouble> > vx(dim);
    vector  < vector  <  adept::adouble> > vx_face(dim);
    vector  < vector  <  double> > vx_hat(dim);

    for (int i = 0; i < dim; i++) {
      vx[i].reserve(max_size);
      vx_face[i].resize(9);
      vx_hat[i].reserve(max_size);
    }

    vector <  vector <  adept::adouble > > Soli(2 * dim + 1);
    vector <  vector <  int > > dofsVAR(2 * dim + 1);

    for (int i = 0; i < 2 * dim + 1; i++) {
      Soli[i].reserve(max_size);
      dofsVAR[i].reserve(max_size);
    }

    vector <  vector <  double > > Rhs(2 * dim + 1);
    vector <  vector <  adept::adouble > > aRhs(2 * dim + 1);

    for (int i = 0; i < 2 * dim + 1; i++) {
      aRhs[i].reserve(max_size);
      Rhs[i].reserve(max_size);
    }

    vector  <  int > dofsAll;
    dofsAll.reserve(max_size * (2 + dim + 1));

    vector  <  double > Jac;
    Jac.reserve(dim * max_size * (2 * dim + 1)*dim * max_size * (2 * dim + 1));

    // ------------------------------------------------------------------------
    // Physical parameters
    double rhof         = ml_prob.parameters.get < Fluid>("Fluid").get_density();
    double mu_lame      = ml_prob.parameters.get < Solid>("Solid").get_lame_shear_modulus();
    double lambda_lame  = ml_prob.parameters.get < Solid>("Solid").get_lame_lambda();
    double mus          = mu_lame / rhof;
    double IRe          = ml_prob.parameters.get < Fluid>("Fluid").get_IReynolds_number();
    double lambda       = lambda_lame / rhof;
    double betans       = 1.;
    int    solid_model  = ml_prob.parameters.get < Solid>("Solid").get_physical_model();

    bool incompressible = (0.5  ==  ml_prob.parameters.get < Solid>("Solid").get_poisson_coeff()) ? 1 : 0;
    const bool penalty = ml_prob.parameters.get < Solid>("Solid").get_if_penalty();

    // gravity
    double _gravity[3] = {0., 0., 0.};

    // -----------------------------------------------------------------
    // space discretization parameters
    unsigned SolType2 = ml_sol->GetSolutionType(ml_sol->GetIndex("U"));
    unsigned SolType1 = ml_sol->GetSolutionType(ml_sol->GetIndex("P"));

    // mesh and procs
    unsigned nel    = mymsh->GetNumberOfElements();
    unsigned igrid  = mymsh->GetLevel();
    unsigned iproc  = mymsh->processor_id();

    //----------------------------------------------------------------------------------
    //variable-name handling
    const char varname[7][3] = {"DX", "DY", "DZ", "U", "V", "W", "P"};
    vector  < unsigned> indexVAR(2 * dim + 1);
    vector  < unsigned> indVAR(2 * dim + 1);
    vector  < unsigned> SolType(2 * dim + 1);

    for (unsigned ivar = 0; ivar < dim; ivar++) {
      indVAR[ivar] = ml_sol->GetIndex(&varname[ivar][0]);
      indVAR[ivar + dim] = ml_sol->GetIndex(&varname[ivar + 3][0]);
      SolType[ivar] = ml_sol->GetSolutionType(&varname[ivar][0]);
      SolType[ivar + dim] = ml_sol->GetSolutionType(&varname[ivar + 3][0]);
      indexVAR[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);
      indexVAR[ivar + dim] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar + 3][0]);
    }

    indexVAR[2 * dim] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[6][0]);
    indVAR[2 * dim] = ml_sol->GetIndex(&varname[6][0]);
    SolType[2 * dim] = ml_sol->GetSolutionType(&varname[6][0]);
    //----------------------------------------------------------------------------------

    int nprocs = mymsh->n_processors();
    NumericVector* area_elem_first;
    area_elem_first = NumericVector::build().release();

    if (nprocs == 1) {
      area_elem_first->init(nprocs, 1, false, SERIAL);
    }
    else {
      area_elem_first->init(nprocs, 1, false, PARALLEL);
    }

    area_elem_first->zero();
    double rapresentative_area = 1.;

    start_time = clock();

    if( assembleMatrix ) myKK->zero();

    // *** element loop ***
    for (int iel = mymsh->_elementOffset[iproc]; iel  <  mymsh->_elementOffset[iproc + 1]; iel++) {

      short unsigned ielt = mymsh->GetElementType(iel);
      unsigned nve        = mymsh->GetElementDofNumber(iel, SolType2);
      unsigned nve1       = mymsh->GetElementDofNumber(iel, SolType1);
      int flag_mat        = mymsh->GetElementMaterial(iel);

      // *******************************************************************************************************

      //initialization of everything in common between fluid and solid

      //Rhs
      for (int i = 0; i < 2 * dim; i++) {
        dofsVAR[i].resize(nve);
        Soli[indexVAR[i]].resize(nve);
        aRhs[indexVAR[i]].resize(nve);
      }

      dofsVAR[2 * dim].resize(nve1);
      Soli[indexVAR[2 * dim]].resize(nve1);
      aRhs[indexVAR[2 * dim]].resize(nve1);

      dofsAll.resize(0);

      // ----------------------------------------------------------------------------------------
      // coordinates, solutions, displacement and velocity dofs

      solidmark.resize(nve);

      for (int i = 0; i < dim; i++) {
        vx[i].resize(nve);
        vx_hat[i].resize(nve);
      }

      for (unsigned i = 0; i < nve; i++) {
        unsigned iDof = mymsh->GetSolutionDof(i, iel, SolType2);

        // flag nodes on the fluid-solid interface
        solidmark[i] = mymsh->GetSolidMark(iDof);

        for (int j = 0; j < dim; j++) {
          Soli[indexVAR[j]][i] = (*mysolution->_Sol[indVAR[j]])(iDof);
          Soli[indexVAR[j + dim]][i] = (*mysolution->_Sol[indVAR[j + dim]])(iDof);

          aRhs[indexVAR[j]][i] = 0.;
          aRhs[indexVAR[j + dim]][i] = 0.;

          //Fixed coordinates (Reference frame)
          vx_hat[j][i] = (*mymsh->_topology->_Sol[j])(iDof);
          // displacement dofs
          dofsVAR[j][i] = myLinEqSolver->GetSystemDof(indVAR[j], indexVAR[j], i, iel);
          // velocity dofs
          dofsVAR[j + dim][i] = myLinEqSolver->GetSystemDof(indVAR[j + dim], indexVAR[j + dim], i, iel);
        }
      }

      // pressure dofs
      for (unsigned i = 0; i < nve1; i++) {
        unsigned iDof = mymsh->GetSolutionDof(i, iel, SolType[2 * dim]);
        dofsVAR[2 * dim][i] = myLinEqSolver->GetSystemDof(indVAR[2 * dim], indexVAR[2 * dim], i, iel);
        Soli[indexVAR[2 * dim]][i] = (*mysolution->_Sol[indVAR[2 * dim]])(iDof);
        aRhs[indexVAR[2 * dim]][i] = 0.;
      }

      // compose the system dofs
      for (int idim = 0; idim < 2 * dim; idim++) {
        dofsAll.insert(dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end());
      }

      dofsAll.insert(dofsAll.end(), dofsVAR[2 * dim].begin(), dofsVAR[2 * dim].end());

      if( assembleMatrix ) s.new_recording();

      //Moving coordinates (Moving frame)
      for (unsigned idim = 0; idim < dim; idim++) {
        for (int j = 0; j < nve; j++) {
          vx[idim][j] = vx_hat[idim][j] + Soli[indexVAR[idim]][j];
        }
      }

      // Boundary integral
      {
        double tau = 0.;
        vector < adept::adouble> normal(dim, 0);

        // loop on faces
        for (unsigned jface = 0; jface < mymsh->GetElementFaceNumber(iel); jface++) {
          std::vector  <  double > xx(3, 0.);

          // look for boundary faces
          if (myel->GetFaceElementIndex(iel, jface) < 0) {
            unsigned int face = -(mymsh->el->GetFaceElementIndex(iel, jface) + 1);

            if ( !ml_sol->GetBdcFunction()(xx, "U", tau, face, 0.) && tau != 0.) {
              unsigned nve = mymsh->GetElementFaceDofNumber(iel, jface, SolType2);
              const unsigned felt = mymsh->GetElementFaceType(iel, jface);

              for (unsigned i = 0; i < nve; i++) {
                unsigned int ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);
                unsigned iDof = mymsh->GetSolutionDof(ilocal, iel, 2);

                for (unsigned idim = 0; idim < dim; idim++) {
                  vx_face[idim][i] = (*mymsh->_topology->_Sol[idim])(iDof) + Soli[indexVAR[idim]][ilocal];
                }
              }

              for (unsigned igs = 0; igs  <  mymsh->_finiteElement[felt][SolType2]->GetGaussPointNumber(); igs++) {
                mymsh->_finiteElement[felt][SolType2]->JacobianSur(vx_face, igs, jacobian, phi, gradphi, normal);

                // *** phi_i loop ***
                for (unsigned i = 0; i < nve; i++) {
                  adept::adouble value = -phi[i] * tau / rhof * jacobian;
                  unsigned int ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);

                  for (unsigned idim = 0; idim < dim; idim++) {
                    if ((!solidmark[ilocal])) { //if fluid node it goes to U, V, W
                      aRhs[indexVAR[dim + idim]][ilocal] +=  value * normal[idim];
                    }
                    else { //if solid node it goes to DX, DY, DZ
                      aRhs[indexVAR[dim + idim]][ilocal] +=  value * normal[idim];
                    }
                  }
                }
              }
            }
          }
        }
      }

      // *** Gauss point loop ***
      double area = 1.;
      adept::adouble supg_tau;

      for (unsigned ig = 0; ig  <  mymsh->_finiteElement[ielt][SolType2]->GetGaussPointNumber(); ig++) {
        // *** get Jacobian and test function and test function derivatives in the moving frame***
        mymsh->_finiteElement[ielt][SolType2]->Jacobian(vx, ig, jacobian, phi, gradphi, nablaphi);
        mymsh->_finiteElement[ielt][SolType2]->Jacobian(vx_hat, ig, jacobian_hat, phi_hat, gradphi_hat, nablaphi_hat);
        phi1 = mymsh->_finiteElement[ielt][SolType1]->GetPhi(ig);

        if (flag_mat == 2 || iel  ==  mymsh->_elementOffset[iproc]) {
          if (ig  ==  0) {
            double GaussWeight = mymsh->_finiteElement[ielt][SolType2]->GetGaussWeight(ig);
            area = jacobian_hat / GaussWeight;

            if (iel == mymsh->_elementOffset[iproc]) {
              area_elem_first->add(mymsh->processor_id(), area);
              area_elem_first->close();
              rapresentative_area = area_elem_first->l1_norm() / nprocs;
            }
          }

          jacobianOverArea = jacobian_hat / area * rapresentative_area;
        }

        // ---------------------------------------------------------------------------
        // displacement and velocity
        for (int i = 0; i < 2 * dim; i++) {
          SolVAR[i] = 0.;

          for (int j = 0; j < dim; j++) {
            GradSolVAR[i][j] = 0.;
            GradSolhatVAR[i][j] = 0.;
          }

          for (int j = 0; j < 3 * (dim - 1); j++) {
            NablaSolVAR[i][j] = 0.;
            NablaSolhatVAR[i][j] = 0.;
          }

          for (unsigned inode = 0; inode < nve; inode++) {
            SolVAR[i] +=  phi[inode] * Soli[indexVAR[i]][inode];

            for (int j = 0; j < dim; j++) {
              GradSolVAR[i][j] +=  gradphi[inode * dim + j] * Soli[indexVAR[i]][inode];
              GradSolhatVAR[i][j] +=  gradphi_hat[inode * dim + j] * Soli[indexVAR[i]][inode];
            }

            for (int j = 0; j < 3 * (dim - 1); j++) {
              NablaSolVAR[i][j] +=  nablaphi[inode * 3 * (dim - 1) + j] * Soli[indexVAR[i]][inode];
              NablaSolhatVAR[i][j] +=  nablaphi_hat[inode * 3 * (dim - 1) + j] * Soli[indexVAR[i]][inode];
            }
          }
        }

        // pressure
        SolVAR[2 * dim] = 0.;

        for (unsigned inode = 0; inode < nve1; inode++) {
          adept::adouble soli = Soli[indexVAR[2 * dim]][inode];
          SolVAR[2 * dim] += phi1[inode] * soli;
        }

        // ---------------------------------------------------------------------------
        //BEGIN FLUID ASSEMBLY
        if (flag_mat == 2) {
          //BEGIN ALE + Momentum (Navier-Stokes)
          {
            for (unsigned i = 0; i < nve; i++) {

              //BEGIN redidual Laplacian ALE map in the reference domain
              adept::adouble LapmapVAR[3] = {0., 0., 0.};

              for (int idim = 0; idim < dim; idim++) {
                for (int jdim = 0; jdim < dim; jdim++) {
                  LapmapVAR[idim] += (GradSolhatVAR[idim][jdim] * gradphi_hat[i * dim + jdim]) ;
                }
              }

              for (int idim = 0; idim < dim; idim++) {
                aRhs[indexVAR[idim]][i] += (!solidmark[i]) * (-LapmapVAR[idim] * jacobianOverArea);
              }

              //END redidual Laplacian ALE map in reference domain

              //BEGIN redidual Navier-Stokes in moving domain
              adept::adouble LapvelVAR[3] = {0., 0., 0.};
              adept::adouble AdvaleVAR[3] = {0., 0., 0.};

              for (int idim = 0.; idim < dim; idim++) {
                for (int jdim = 0.; jdim < dim; jdim++) {
                  //LapvelVAR[idim] +=  GradSolVAR[dim+idim][jdim]*gradphi[i*dim+jdim];
                  LapvelVAR[idim] += (GradSolVAR[dim + idim][jdim] + GradSolVAR[dim + jdim][idim]) * gradphi[i * dim + jdim];
                  AdvaleVAR[idim] +=  SolVAR[dim + jdim] * GradSolVAR[dim + idim][jdim] * phi[i];
                }
              }

              for (int idim = 0; idim < dim; idim++) {
                adept::adouble value = (-AdvaleVAR[idim]                   // advection term
                                        - IRe * LapvelVAR[idim]          // viscous dissipation
                                        + SolVAR[2 * dim] * gradphi[i * dim + idim] // pressure gradient
                                       ) * jacobian;

                if ((!solidmark[i])) {
                  aRhs[indexVAR[dim + idim]][i] += value;
                }
                else {
                  aRhs[indexVAR[dim+idim]][i] +=  value;
                }

                //END redidual Navier-Stokes in moving domain
              }
            }
          }
          //END ALE + Momentum (Navier-Stokes)

          //BEGIN continuity block
          {
            adept::adouble div_vel = 0.;

            for (int i = 0; i < dim; i++) {
              div_vel += GradSolVAR[dim + i][i];
            }

            for (unsigned i = 0; i < nve1; i++) {
              aRhs[indexVAR[2 * dim]][i] += -(-phi1[i] * div_vel) * jacobian;
            }
          }
          //END continuity block
        }
        //END FLUID ASSEMBLY

        //*******************************************************************************************************

        //BEGIN SOLID ASSEMBLY
        else {
          //BEGIN build Chauchy Stress in moving domain
          //physical quantity
          adept::adouble J_hat;
          adept::adouble I_e;
          adept::adouble Cauchy[3][3];
          adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};

          adept::adouble I1_B = 0.;
          adept::adouble I2_B = 0.;

          if (solid_model == 0) { // Saint-Venant
            adept::adouble e[3][3];

            //computation of the stress tensor
            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                e[i][j] = 0.5 * (GradSolhatVAR[i][j] + GradSolhatVAR[j][i]);
              }
            }

            I_e = 0;

            for (int i = 0; i < dim; i++) {
              I_e += e[i][i];
            }

            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                //incompressible
                Cauchy[i][j] = 2 * mus * e[i][j] - 2 * mus * I_e * SolVAR[2 * dim] * Id2th[i][j];
                //+(penalty)*lambda*I_e*Id2th[i][j];
              }
            }
          }

          else { // hyperelastic non linear material
            adept::adouble F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
            adept::adouble B[3][3];

            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                F[i][j] += GradSolhatVAR[i][j];
              }
            }

            J_hat =   F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                      - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

            for (int I = 0; I < 3; ++I) {
              for (int J = 0; J < 3; ++J) {
                B[I][J] = 0.;

                for (int K = 0; K < 3; ++K) {
                  //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
                  B[I][J] += F[I][K] * F[J][K];
                }
              }
            }

            if (solid_model <= 4) { // Neo-Hookean
              I1_B = B[0][0] + B[1][1] + B[2][2];

              for (int I = 0; I < 3; ++I) {
                for (int J = 0; J < 3; ++J) {
                  if (1  ==  solid_model) Cauchy[I][J] = mus * B[I][J]
                        - mus * I1_B * SolVAR[2 * dim] * Id2th[I][J];   //Wood-Bonet J_hat  =1;
                  else if (2  ==  solid_model) Cauchy[I][J] = mus / J_hat * B[I][J]
                        - mus / J_hat * SolVAR[2 * dim] * Id2th[I][J]; //Wood-Bonet J_hat !=1;
                  else if (3  ==  solid_model) Cauchy[I][J] = mus * (B[I][J] - Id2th[I][J]) / J_hat
                        + lambda / J_hat * log(J_hat) * Id2th[I][J];    //Wood-Bonet penalty
                  else if (4  ==  solid_model) Cauchy[I][J] = mus * (B[I][J] - I1_B * Id2th[I][J] / 3.) / pow(J_hat, 5. / 3.)
                        + lambda * (J_hat - 1.) * Id2th[I][J];           //Allan-Bower

                }
              }
            }
            else if (5  ==  solid_model) {  //Mooney-Rivlin
              adept::adouble detB =   B[0][0] * (B[1][1] * B[2][2] - B[2][1] * B[1][2])
                                      - B[0][1] * (B[2][2] * B[1][0] - B[1][2] * B[2][0])
                                      + B[0][2] * (B[1][0] * B[2][1] - B[2][0] * B[1][1]);
              adept::adouble invdetB = 1. / detB;
              adept::adouble invB[3][3];

              invB[0][0] = (B[1][1] * B[2][2] - B[2][1] * B[1][2]) * invdetB;
              invB[1][0] = -(B[0][1] * B[2][2] - B[0][2] * B[2][1]) * invdetB;
              invB[2][0] = (B[0][1] * B[1][2] - B[0][2] * B[1][1]) * invdetB;
              invB[0][1] = -(B[1][0] * B[2][2] - B[1][2] * B[2][0]) * invdetB;
              invB[1][1] = (B[0][0] * B[2][2] - B[0][2] * B[2][0]) * invdetB;
              invB[2][1] = -(B[0][0] * B[1][2] - B[1][0] * B[0][2]) * invdetB;
              invB[0][2] = (B[1][0] * B[2][1] - B[2][0] * B[1][1]) * invdetB;
              invB[1][2] = -(B[0][0] * B[2][1] - B[2][0] * B[0][1]) * invdetB;
              invB[2][2] = (B[0][0] * B[1][1] - B[1][0] * B[0][1]) * invdetB;

              I1_B = B[0][0] + B[1][1] + B[2][2];
              I2_B = B[0][0] * B[1][1] + B[1][1] * B[2][2] + B[2][2] * B[0][0]
                     - B[0][1] * B[1][0] - B[1][2] * B[2][1] - B[2][0] * B[0][2];

              double C1 = mus / 3.;
              double C2 = C1 / 2.;

              for (int I = 0; I < 3; ++I) {
                for (int J = 0; J < 3; ++J) {
                  Cauchy[I][J] =  2.*(C1 * B[I][J] - C2 * invB[I][J])
                                  //- (2. / 3.) * (C1 * I1_B - C2 * I2_B) * SolVAR[2 * dim] * Id2th[I][J];
                                  - SolVAR[2 * dim] * Id2th[I][J];
                }
              }

            }
          }

          //END build Chauchy Stress in moving domain

          //BEGIN v=0 + Momentum (Solid)
          {
            for (unsigned i = 0; i < nve; i++) {

              //BEGIN redidual v=0 in fixed domain
              for (int idim = 0; idim < dim; idim++) {
                aRhs[indexVAR[idim]][i] += (-phi[i] * (-SolVAR[dim + idim])) * jacobian_hat;
              }

              //END redidual v=0 in fixed domain

              //BEGIN redidual Solid Momentum in moving domain
              adept::adouble CauchyDIR[3] = {0., 0., 0.};

              for (int idim = 0.; idim < dim; idim++) {
                for (int jdim = 0.; jdim < dim; jdim++) {
                  CauchyDIR[idim] += gradphi[i * dim + jdim] * Cauchy[idim][jdim];
                }
              }

              for (int idim = 0; idim < dim; idim++) {
                aRhs[indexVAR[dim+idim]][i] += (phi[i] * _gravity[idim] - CauchyDIR[idim]) * jacobian;
              }

              //END redidual Solid Momentum in moving domain
            }
          }
          //END v=0 + Momentum (Solid)

          //BEGIN continuity block
          {
            for (unsigned i = 0; i < nve1; i++) {
              if (!penalty) {
                if (0  ==  solid_model) {
                  aRhs[indexVAR[2 * dim]][i] += -(-phi1[i] * (I_e + (!incompressible) / lambda * SolVAR[2 * dim])) * jacobian_hat;
                }
                else if (1  ==  solid_model || 5  ==  solid_model) {
                  aRhs[indexVAR[2 * dim]][i] += phi1[i] * (J_hat - 1. + (!incompressible) / lambda * SolVAR[2 * dim]) * jacobian_hat;
                }
                else if (2  ==  solid_model) {
                  aRhs[indexVAR[2 * dim]][i] +=  -(-phi1[i] * (log(J_hat) / J_hat + (!incompressible) / lambda * SolVAR[2 * dim])) * jacobian_hat;
                }
              }
              else if (3  ==  solid_model || 4  ==  solid_model) { // pressure = 0 in the solid
                aRhs[indexVAR[2 * dim]][i] += -(-phi1[i] * (SolVAR[2 * dim])) * jacobian_hat;
              }
            }
          }
          //END continuity block
        }

        //END SOLID ASSEMBLY
      }

      //BEGIN local to global assembly
      //copy adouble aRhs into double Rhs
      for (unsigned i = 0; i < 2 * dim; i++) {
        Rhs[indexVAR[i]].resize(nve);

        for (int j = 0; j < nve; j++) {
          Rhs[indexVAR[i]][j] = -aRhs[indexVAR[i]][j].value();
        }
      }

      Rhs[indexVAR[2 * dim]].resize(nve1);

      for (unsigned j = 0; j < nve1; j++) {
        Rhs[indexVAR[2 * dim]][j] = -aRhs[indexVAR[2 * dim]][j].value();
      }

      for (int i = 0; i < 2 * dim + 1; i++) {
        myRES->add_vector_blocked(Rhs[indexVAR[i]], dofsVAR[i]);
      }

      if( assembleMatrix ){
        //Store equations
        for (int i = 0; i < 2 * dim; i++) {
          s.dependent(&aRhs[indexVAR[i]][0], nve);
          s.independent(&Soli[indexVAR[i]][0], nve);
        }

        s.dependent(&aRhs[indexVAR[2 * dim]][0], nve1);
        s.independent(&Soli[indexVAR[2 * dim]][0], nve1);

        Jac.resize((2 * dim * nve + nve1) * (2 * dim * nve + nve1));

        s.jacobian(&Jac[0], true);

        myKK->add_matrix_blocked(Jac, dofsAll, dofsAll);
        s.clear_independents();
        s.clear_dependents();

      //END local to global assembly
      }
    } //end list of elements loop

    if( assembleMatrix ) myKK->close();
    myRES->close();

    delete area_elem_first;

    // *************************************
    end_time = clock();
    AssemblyTime += (end_time - start_time);
    // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

  }

  
  
  
 #define FE_DOMAIN  2

 #define POS_U   0

 #define WET_RIGID  5            //set GROUP_VOL_ELEMS 6
 #define WET_DEFORMABLE 6        //set GROUP_VOL_ELEMS 6
 #define DRY_RIGID_DEFORMABLE  7 //set GROUP_VOL_ELEMS 5

#define GROUP_VOL_ELEMS  6
 
 #define EPS_EDGE_LOCATION 1.e-5   //with e-5 it doesn't find the circle!!!


bool or_vector(const int current_face, const std::vector< int > all_face_flags) {
 
     bool is_face_there = false;
    
       for (unsigned l = 0; l < all_face_flags.size(); l++) {
           if (current_face == all_face_flags[l]) { is_face_there = true; break; }

       }
    
    return is_face_there;
}



 void allocate_element_faces(  std::vector< MyMatrix <int> > & _element_faces, /*const*/ MultiLevelMesh & ml_msh) {

  _element_faces[0].resize( ml_msh.GetLevel(0)->GetNumberOfElements(), NFC[0][1], -1);  /*NFC[0][1]: maximum possible number of faces*/

  
  for (unsigned lev = 1; lev < _element_faces.size(); lev++) {
      
    std::vector < double > coarseLocalizedAmrVector;
    ml_msh.GetLevel(lev - 1)->_topology->_Sol[ml_msh.GetLevel(lev - 1)->GetAmrIndex()]->localize_to_all(coarseLocalizedAmrVector);

//     ml_msh.GetLevel(lev - 1)->el->AllocateChildrenElement(ml_msh.GetLevel(lev)->GetRefIndex(), ml_msh.GetLevel(lev - 1) );  //I believe this was already done

    const unsigned n_elems = ml_msh.GetLevel(lev)->GetNumberOfElements();
       
       MyVector <unsigned> rowSizeElNearFace(n_elems);
       
    unsigned jel = 0;
    
    for (unsigned isdom = 0; isdom < ml_msh.GetLevel(lev)->n_processors(); isdom++) {
        
       ml_msh.GetLevel(lev)->GetElementArray()->GetElementTypeArray().broadcast(isdom);
       
      for (unsigned iel = ml_msh.GetLevel(lev)->GetElementArray()->GetElementTypeArray().begin(); 
                    iel < ml_msh.GetLevel(lev)->GetElementArray()->GetElementTypeArray().end(); iel++) {
          
        short unsigned elType = ml_msh.GetLevel(lev)->GetElementArray()->GetElementTypeArray()[iel];
      
        int increment = 1;
      
        if (static_cast < short unsigned >(coarseLocalizedAmrVector[iel] + 0.25) == 1) {
          increment = NRE[elType];
        }
        
        for (unsigned j = 0; j < increment; j++) {
          rowSizeElNearFace[jel + j] += NFC[elType][1];
        }
        
        jel += increment;
      }
      
      ml_msh.GetLevel(lev)->GetElementArray()->GetElementTypeArray().clearBroadcast();
    }
         _element_faces[lev] =   MyMatrix < int > (rowSizeElNearFace, -1); 

      }

     
}
// ==================================



 int find_faces_for_integration_based_on_face_center( const unsigned dim, const std::vector< double > x ) { 

     const double epsilon = EPS_EDGE_LOCATION;
     
     int face_flag = 0;      /*std::cout << "The face numbers start from 1" << std::endl;*/
     
     double x_offset = 0.;
     if (dim == 3) x_offset = 0.3;
     
     const double cyl_radius = 0.05;
     const double cyl_radius_squared = cyl_radius * cyl_radius;
     
     const double x2plusy2 = (x[0] - (x_offset + 0.2) ) * (x[0] - (x_offset + 0.2) ) + (x[1] - 0.2) * (x[1] - 0.2);
     
     const bool is_cylinder_surface = (x2plusy2 > (cyl_radius_squared - epsilon) &&
                                       x2plusy2 < (cyl_radius_squared + epsilon) );

     
     const bool is_part_of_cylinder_surface_where_flap_clamps = ( x[0] > (x_offset + 0.2) &&  x[1] > (0.19 - epsilon) &&  x[1] < (0.21 + epsilon) );

     const bool wet_deformable_is_above = ( x[1] > ( 0.21 - epsilon )  && x[1] < ( 0.21 + epsilon ) );
     
     const bool wet_deformable_is_below = ( x[1] > (0.19 - epsilon)  && x[1] < ( 0.19 + epsilon) );
     
     const bool wet_deformable_is_in_horizontal_range = ( x[0] > ( (x_offset + 0.248) - epsilon ) && x[0] < ( (x_offset + 0.6) + epsilon ) );

     const bool wet_deformable_is_in_vertical_range =   ( x[1] > ( 0.19 - epsilon)  && x[1] < (0.21 + epsilon) );   
     
     const bool wet_deformable_is_in_tip = ( x[0] > ( (x_offset + 0.6) - epsilon ) && x[0] < ( (x_offset + 0.6) + epsilon ) );
     
          if (dim == 2) {
         
     
              if (
                  is_cylinder_surface   &&   !(is_part_of_cylinder_surface_where_flap_clamps) 
                 ) { 
                     face_flag = WET_RIGID;
                  }
       
         else if ( ( wet_deformable_is_in_horizontal_range &&  wet_deformable_is_above ) 
            ||     ( wet_deformable_is_in_horizontal_range &&  wet_deformable_is_below )
            ||     ( wet_deformable_is_in_tip              &&  wet_deformable_is_in_vertical_range )
                 ) {
                     face_flag = WET_DEFORMABLE;
                  }
        
         else  if ( is_cylinder_surface &&  is_part_of_cylinder_surface_where_flap_clamps )  {  
                    face_flag = DRY_RIGID_DEFORMABLE; 
                  }

        
     }
     
     else if (dim == 3) { 
         
     const bool three_d_wet_deformable_is_in_z_range    = ( x[2] > ( 0.1 - epsilon)  && x[2] < (0.3 + epsilon) );
         
     const bool three_d_wet_deformable_is_in_side_left  = ( x[2] > ( 0.1 - epsilon)  && x[2] < (0.1 + epsilon) );

     const bool three_d_wet_deformable_is_in_side_right = ( x[2] > ( 0.3 - epsilon)  && x[2] < (0.3 + epsilon) );
     
     const bool three_d_is_part_of_cylinder_surface_where_flap_clamps = three_d_wet_deformable_is_in_z_range && ( x[0] > (x_offset + 0.2) &&  x[1] > (0.19 - epsilon) &&  x[1] < (0.21 + epsilon) );
     
         
              if (
                  is_cylinder_surface   &&   !(three_d_is_part_of_cylinder_surface_where_flap_clamps ) 
                 ) { 
                     face_flag = WET_RIGID;
                  }
       
      else if (        ( three_d_wet_deformable_is_in_z_range    && wet_deformable_is_in_horizontal_range &&  wet_deformable_is_above ) 
            ||     ( three_d_wet_deformable_is_in_z_range    && wet_deformable_is_in_horizontal_range &&  wet_deformable_is_below )
            ||     ( three_d_wet_deformable_is_in_z_range    && wet_deformable_is_in_tip              &&  wet_deformable_is_in_vertical_range )
            ||     ( three_d_wet_deformable_is_in_side_left  && wet_deformable_is_in_horizontal_range &&  wet_deformable_is_in_vertical_range )
            ||     ( three_d_wet_deformable_is_in_side_right && wet_deformable_is_in_horizontal_range &&  wet_deformable_is_in_vertical_range )
                 ) {
                     face_flag = WET_DEFORMABLE;
                  }
         
       else  if ( is_cylinder_surface &&  three_d_is_part_of_cylinder_surface_where_flap_clamps )  {  
                    face_flag = DRY_RIGID_DEFORMABLE; 
                  }
         
         
    }
     
     
     return face_flag;
     
 }

 
 
 
 // ==================================
 void fill_element_faces(  std::vector< MyMatrix <int> > & _element_faces, 
                           const MultiLevelMesh & ml_msh, 
                           const std::vector < int > & all_face_flags, 
                           const int group_outside_solid_to_which_all_faces_belong) {
     
   int dimension = ml_msh.GetLevel(0)->GetDimension();
  
   unsigned solType_coords = FE_DOMAIN;
   
   

        for (unsigned lev = 0; lev < _element_faces.size(); lev++) {
            
            unsigned face_count = 0;
                int count_volume_elements = 0;

        for (unsigned iel = 0; iel < _element_faces[lev].size(); iel++) {
            
       CurrentElem < double > geom_element(dimension, ml_msh.GetLevel(lev));
       
    geom_element.set_coords_at_dofs_and_geom_type(iel, solType_coords);
  
    int iel_group = ml_msh.GetLevel(lev)->GetElementGroup(iel);
    
    
      if (iel_group == group_outside_solid_to_which_all_faces_belong)  {
          
          count_volume_elements++;
    
        const unsigned int n_faces_iel = ml_msh.GetLevel(lev)->GetElementFaceNumber(iel);

          
        for (unsigned f = 0; f < n_faces_iel; f++) {
            
        const unsigned ielGeom_bdry = ml_msh.GetLevel(lev)->GetElementFaceType(iel, f);    
       

       geom_element.set_coords_at_dofs_bdry_3d(iel, f, solType_coords);
 
       geom_element.set_elem_center_bdry_3d(); //this is computed linearly, so you may offset!
       
       // here I use the Element Center of each face to locate the faces.
       // Alternatively, I could do a Node-based criterion and say that all Nodes of the face have to belong
       // Instead of doing it on all nodes, I just need to get the Dof of the center of the face
       //once I have that dof, I will get its coordinates
       
                   unsigned nv1 = ml_msh.GetLevel(lev)->GetElementFaceDofNumber(iel, f, solType_coords);  // only the face dofs

                  unsigned i_local_face_center = ml_msh.GetLevel(lev)->GetLocalFaceVertexIndex(iel, f, nv1 - 1/*iv*/);  //the center is always the last one
                  unsigned idof = ml_msh.GetLevel(lev)->GetSolutionDof(i_local_face_center, iel, solType_coords);
                  
       std::vector< double > face_center(3, 0.); 
       
//        face_center = geom_element.get_elem_center_bdry(); //old method
       
        for (unsigned d = 0; d < face_center.size(); d++) {
                        face_center[d] = (*ml_msh.GetLevel(lev)->_topology->_Sol[d])(idof);
        }
        
         const int face_flag_wet  = find_faces_for_integration_based_on_face_center (dimension, face_center);
         
         const int elem_near_face = ml_msh.GetLevel(lev)->GetElementArray()->GetFaceElementIndex(iel, f) - 1; //@todo have to subtract 1 because it was added before!

              bool face_already_found_from_near_elem = false;
              
         if (elem_near_face >= 0 ) {
              
             const unsigned int n_faces_near_elem = ml_msh.GetLevel(lev)->GetElementFaceNumber(elem_near_face);
               for (unsigned int v = 0; v < n_faces_near_elem; v++) {
                     if ( or_vector(_element_faces[lev][elem_near_face][v], all_face_flags) ) face_already_found_from_near_elem = true;
               }
         }
         

         if ( or_vector(face_flag_wet, all_face_flags) && !(face_already_found_from_near_elem) ) { //in this way the face is counted only once... but this already happens by reducing to the group
              face_count++;
           _element_faces[lev][iel][f] = face_flag_wet; 
         }
//             std::cout << _element_faces[lev][iel][f];

             }
             
           } //end group outside solid
            
            
        }  //end volume elems
        
            std::cout << "Volume elements in group " << group_outside_solid_to_which_all_faces_belong << ": " << count_volume_elements  << std::endl;
            
            std::cout << "Faces on which QoI integration is performed: " << face_count << std::endl;
            
        }


 }


// ==================================

 
  
  void Compute_normal_stress_interface(const MultiLevelProblem& ml_prob,
                       const unsigned level, 
                       const std::vector < int > all_face_flags,
                       const unsigned stress_component, 
                       /*const*/ std::vector< MyMatrix <int> > & element_faces, //@todo here the operator with const has to be added to MyMatrix
                       const MonolithicFSINonLinearImplicitSystem* mlPdeSys)    {
  
  
//   const MonolithicFSINonLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system< MonolithicFSINonLinearImplicitSystem > ("Fluid-Structure-Interaction");
//   const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->el;

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id();

  //=============== Geometry ========================================
   unsigned solType_coords = FE_DOMAIN;
 
  CurrentElem < double > geom_element(dim, msh);
    
  constexpr unsigned int space_dim = 3;
  
  std::vector< double > normal(space_dim, 0.);
 //***************************************************

  //=============== Integration ========================================
  double weight_vol = 0.;
  double weight_bdry = 0.;


  std::vector< std::vector< double > > deform_tensor_qp(dim);
   for (unsigned d = 0; d < deform_tensor_qp.size(); d++) {   deform_tensor_qp[d].resize(dim); }
  
 //*************** unknowns ***************************** 
 //*************************************************** 
   const unsigned int n_fluid_unknowns = dim + 1;
   const unsigned int pressure_index = n_fluid_unknowns - 1;
 
  vector < vector < double > > phi_u(n_fluid_unknowns);     
  vector < vector < double > > phi_u_x(n_fluid_unknowns);   
 
   for (unsigned d = 0; d < phi_u.size(); d++) {
       phi_u[d].reserve(max_size);
       phi_u_x[d].reserve(max_size * space_dim);
   }
   
   
   
 std::vector < std::string >  var_names(n_fluid_unknowns);
 var_names[0] = "U";
 var_names[1] = "V";
 if (dim == 3) var_names[2] =  "W";
 var_names[pressure_index] = "P";
 
  std::vector <  unsigned > solIndex_u(n_fluid_unknowns); 
  std::vector <  unsigned > solType_u(n_fluid_unknowns); 
  std::vector <  unsigned > nDof_u(n_fluid_unknowns); 
  
     for (unsigned d = 0; d < solIndex_u.size(); d++) {
         solIndex_u[d] = ml_sol->GetIndex(var_names[d].c_str());            // get the position of "state" in the ml_sol object
         solType_u[d]  = ml_sol->GetSolutionType(solIndex_u[d]);    // get the finite element type for "state"
    }

  vector < vector < double > > sol_u(n_fluid_unknowns); // local solution
   for (unsigned d = 0; d < sol_u.size(); d++) {  sol_u[d].reserve(max_size);   }
  
 //*************************************************** 
 //***************************************************
 //***************************************************
  vector < vector < double > > phi_u_bdry(n_fluid_unknowns);
  vector < vector < double > > phi_u_x_bdry(n_fluid_unknowns); 

   for (unsigned d = 0; d < phi_u_bdry.size(); d++) {   
       phi_u_bdry[d].reserve(max_size);
       phi_u_x_bdry[d].reserve(max_size * space_dim);
   }
   
  //volume shape functions at boundary
  vector < vector < double > > phi_u_vol_at_bdry(n_fluid_unknowns);
  vector < vector < double > > phi_u_x_vol_at_bdry(n_fluid_unknowns);
  for (unsigned d = 0; d < phi_u_vol_at_bdry.size(); d++) {   
       phi_u_vol_at_bdry[d].reserve(max_size);
     phi_u_x_vol_at_bdry[d].reserve(max_size * space_dim);
  }
   
  vector < vector < double > > sol_u_x_vol_at_bdry_gss(n_fluid_unknowns);
  for (unsigned d = 0; d < sol_u_x_vol_at_bdry_gss.size(); d++) {
      sol_u_x_vol_at_bdry_gss[d].resize(space_dim);
  }
 //***************************************************
  
  
  double integral_volume = 0.;
  std:vector < double > integral_norm_stress_component(dim, 0.);
  

 //*************************************************** 
     std::vector < std::vector < double > >  JacI_qp(space_dim);
     std::vector < std::vector < double > >  Jac_qp(dim);
    for (unsigned d = 0; d < Jac_qp.size(); d++) {   Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp.size(); d++) { JacI_qp[d].resize(dim); }
    
    double detJac_qp;

     std::vector < std::vector < double > >  JacI_qp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_qp_bdry(dim - 1);
    for (unsigned d = 0; d < Jac_qp_bdry.size(); d++) {   Jac_qp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp_bdry.size(); d++) { JacI_qp_bdry[d].resize(dim - 1); }
    
    double detJac_qp_bdry;
    
      //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
 //*************************************************** 
  
  
 //**** physical parameters ************************** 
    double mu_f = ml_prob.parameters.get<Fluid>("Fluid").get_viscosity();
 //*************************************************** 
 
    
  // element loop: each process loops only on its own elements
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    geom_element.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    const short unsigned ielGeom = geom_element.geom_type();

   
 //*********** state ********************************* 
   for (unsigned d = 0; d < sol_u.size(); d++) {   
    nDof_u[d]     = msh->GetElementDofNumber(iel, solType_u[d]);
    sol_u[d]    .resize(nDof_u[d]);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u[d].size(); i++) {
      unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u[d]);
      sol_u[d][i] = (*sol->_Sol[solIndex_u[d]])(solDof_u);
     }
   }
 //*********** state ********************************* 

 
 //********** ALL VARS ******************************* 
    int nDof_max    = msh->GetElementDofNumber(iel, solType_u[POS_U]);   //  TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
 //***************************************************

	       
	  // loop on faces of the current element

	  for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
          
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, jface);    
       

       geom_element.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);
 
       geom_element.set_elem_center_bdry_3d();

       const unsigned nve_bdry_u = msh->GetElementFaceDofNumber(iel,jface,solType_u[POS_U]);
       
// // // 	    // look for boundary faces
// // //             const int bdry_index = msh->el->GetFaceElementIndex(iel,jface);
// // //             
// // // 	    if( bdry_index < 0) {
// // // 	      unsigned int face = msh->el->GetBoundaryIndex(iel,jface);
// // // 	      
// // // 		
// // // // 	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"U",tau,face,0.) && tau!=0.){

       if ( or_vector(element_faces[level][iel][jface], all_face_flags) ) {
       
//          //========= check face center ================================================
//        std::vector< double > face_center(3, 0.); 
//        face_center = geom_element.get_elem_center_bdry(); //old method
//        
//        std::cout << "Face center ====" << std::endl;
//        
//        for (unsigned d = 0; d < dim; d++) {
//          std::cout << face_center[d] << " ";
//        }
//        std::cout << std::endl;
//          //========= check face center ================================================
	
		//============ initialize gauss quantities on the boundary ==========================================
                std::vector< double >                sol_u_bdry_gss(sol_u.size());
                std::vector< std::vector< double > > sol_u_x_bdry_gss(sol_u.size());
                       std::fill( sol_u_bdry_gss.begin(), sol_u_bdry_gss.end(), 0.); 
                   for (unsigned d = 0; d < sol_u_x_bdry_gss.size(); d++) {
                       sol_u_x_bdry_gss[d].resize(dim);
                       std::fill( sol_u_x_bdry_gss[d].begin(), sol_u_x_bdry_gss[d].end(), 0.); 
                }
		//============ initialize gauss quantities on the boundary ==========================================
		
		for(unsigned ig_bdry = 0; ig_bdry < ml_prob.GetQuadratureRule(ielGeom_bdry).GetGaussPointsNumber(); ig_bdry++) {
		  
    elem_all[ielGeom_bdry][solType_coords]->JacJacInv(geom_element.get_coords_at_dofs_bdry_3d(), ig_bdry, Jac_qp_bdry, JacI_qp_bdry, detJac_qp_bdry, space_dim);
	elem_all[ielGeom_bdry][solType_coords]->compute_normal(Jac_qp_bdry, normal);

    weight_bdry = detJac_qp_bdry * ml_prob.GetQuadratureRule(ielGeom_bdry).GetGaussWeightsPointer()[ig_bdry];
    
  for (unsigned d = 0; d < phi_u_bdry.size(); d++) {
      elem_all[ielGeom_bdry][solType_u[d]] ->shape_funcs_current_elem(ig_bdry, JacI_qp_bdry, phi_u_bdry[d], phi_u_x_bdry[d], boost::none, space_dim);
  }
    
    elem_all[ielGeom][solType_coords]->JacJacInv_vol_at_bdry_new(geom_element.get_coords_at_dofs_3d(), ig_bdry, jface, Jac_qp/*not_needed_here*/, JacI_qp, detJac_qp/*not_needed_here*/, space_dim);
    
  for (unsigned d = 0; d < phi_u_vol_at_bdry.size(); d++) {
    elem_all[ielGeom][solType_u[d]]->shape_funcs_vol_at_bdry_current_elem(ig_bdry, jface, JacI_qp, phi_u_vol_at_bdry[d], phi_u_x_vol_at_bdry[d], boost::none, space_dim);
  }
  
		  
		 //========== compute gauss quantities on the boundary ===============================================
          for (unsigned d = 0; d < sol_u_bdry_gss.size(); d++) {
		    sol_u_bdry_gss[d] = 0.;
                  std::fill(sol_u_x_bdry_gss[d].begin(), sol_u_x_bdry_gss[d].end(), 0.);
          }
          
          
	 for (int i_bdry = 0; i_bdry < nve_bdry_u; i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
			
         for (int e = 0; e < sol_u_bdry_gss.size(); e++) {
			sol_u_bdry_gss[e] +=  sol_u[e][i_vol] * phi_u_bdry[e][i_bdry];
                            for (int d = 0; d < dim; d++) {
			      sol_u_x_bdry_gss[e][d] += sol_u[e][i_vol] * phi_u_x_bdry[e][i_bdry * space_dim + d];
			    }
         }
	  }
         //========= compute gauss quantities on the boundary ================================================

//          //========= check normal computation ================================================
//        std::cout << "Normal " << std::endl;
//        for (unsigned d = 0; d < dim; d++) {
//          std::cout << normal[d] << " ";
//        }
//        std::cout << std::endl;
//          //========= check normal computation ================================================
       
       
//     compute gauss quantities on the boundary through VOLUME interpolation
          for (unsigned d = 0; d < sol_u_x_vol_at_bdry_gss.size(); d++) {
              std::fill(sol_u_x_vol_at_bdry_gss[d].begin(), sol_u_x_vol_at_bdry_gss[d].end(), 0.);
           }
           
         for (int e = 0; e < sol_u_x_vol_at_bdry_gss.size(); e++) {
		      for (int iv = 0; iv < nDof_u[e]; iv++)  {
			
            for (int d = 0; d < dim; d++) {
			      sol_u_x_vol_at_bdry_gss[e][d] += sol_u[e][iv] * phi_u_x_vol_at_bdry[e][iv * space_dim + d];
			    }
		       }
		       
           }
              
		      
       for (unsigned d = 0; d < deform_tensor_qp.size(); d++) {   std::fill( deform_tensor_qp[d].begin(), deform_tensor_qp[d].end(), 0.); }
       
       for (unsigned d = 0; d < deform_tensor_qp.size(); d++) {
           for (unsigned e = 0; e < deform_tensor_qp[d].size(); e++) {
               deform_tensor_qp[d][e] = 0.5 * (sol_u_x_vol_at_bdry_gss[d][e] + sol_u_x_vol_at_bdry_gss[e][d]);
           }
       }
       
    //--------       
       for (unsigned d = 0; d < deform_tensor_qp.size(); d++) {
           
          double normal_deform_tensor = 0.;
          
            for (unsigned e = 0; e < dim; e++) {
                  normal_deform_tensor += deform_tensor_qp[d][e] * normal[e]; 
               }
//     compute gauss quantities on the boundary through VOLUME interpolation


                 integral_norm_stress_component[d] +=  weight_bdry * mu_f * normal_deform_tensor - sol_u_bdry_gss[pressure_index] * normal[d] /*1.*/; /*normal[0] **/ /*normal[1] **/ /*normal[2] **/ 
                 
          }
    //--------       
            
           } //end gauss boundary
           
         } //end faces
// // // 	    } //end if boundary faces
	    
	  }  // loop over element faces   

//=====================================================================================================================  
//=====================================================================================================================  
//=====================================================================================================================  
  
  
   
      // *** Gauss point loop ***
      for (unsigned ig_vol = 0; ig_vol < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig_vol++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[ielGeom][solType_coords]->JacJacInv(geom_element.get_coords_at_dofs_3d(), ig_vol, Jac_qp, JacI_qp, detJac_qp, space_dim);
    weight_vol = detJac_qp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[ig_vol];

  for (unsigned d = 0; d < phi_u.size(); d++) {
    elem_all[ielGeom][solType_u[d]]                 ->shape_funcs_current_elem(ig_vol, JacI_qp, phi_u[d], phi_u_x[d], boost::none, space_dim);
  }
  
               integral_volume +=  weight_vol * 1.;
	  
      } // end gauss point loop
      
  } //end element loop

  
  ////////////////////////////////////////
        std::cout << "The normals on the fluid-solid interface are pointing from the fluid to the solid: " << std::endl;
 for (unsigned d = 0; d < dim; d++) {
      
       std::cout << "integral on processor: " << integral_norm_stress_component[d] << std::endl;

   std::vector<double> J(dim, 0.);
   
      MPI_Allreduce( &integral_norm_stress_component[d], &J[d], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  //THIS IS THE RIGHT ONE!!

    std::cout << "@@@@@@@@@@@@@@@@ functional value: " << J[d] << std::endl;
  }
  
//   std::cout << "The value of the integral_target is " << std::setw(11) << std::setprecision(10) << integral_target << std::endl;
 
return;
  
}

  
  

}
