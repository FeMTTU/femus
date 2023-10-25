#ifndef __femus_include_IncompressibleFSIAssembly_hpp__
#define __femus_include_IncompressibleFSIAssembly_hpp__

#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "MultiLevelSolution.hpp"
#include "adept.h"


namespace femus {

  void FSITimeDependentAssembly(MultiLevelProblem& ml_prob) {

    clock_t AssemblyTime = 0;
    clock_t start_time, end_time;

    adept::Stack& s = FemusInit::_adeptStack;

    //pointers and references

    //MonolithicFSINonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<MonolithicFSINonLinearImplicitSystem>("Fluid-Structure-Interaction");
    TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem>("Fluid-Structure-Interaction");
    const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();

    MultiLevelSolution*	 ml_sol		= ml_prob._ml_sol;
    Solution*	 mysolution		= ml_sol->GetSolutionLevel(level);

    LinearEquationSolver*  myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];
    Mesh*		mymsh		=  ml_prob._ml_msh->GetLevel(level);
    elem*		myel		=  mymsh->el;
    SparseMatrix*	myKK		=  myLinEqSolver->_KK;
    NumericVector*	 myRES		=  myLinEqSolver->_RES;

    const unsigned dim = mymsh->GetDimension();
    const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));
    
    double theta = 0.5;

    // local objects
    vector<adept::adouble> SolVAR(2 * dim + 1);
    vector<double> SolVAR_old(2 * dim + 1);

    vector<vector < adept::adouble > > GradSolVAR(2 * dim);
    vector<vector < double > > GradSolVAR_old(2 * dim);

    vector<vector < adept::adouble > > GradSolhatVAR(2 * dim);
    vector<vector < double > > GradSolhatVAR_old(2 * dim);

    for (int i = 0; i < 2 * dim; i++) {
      GradSolVAR[i].resize(dim);
      GradSolVAR_old[i].resize(dim);

      GradSolhatVAR[i].resize(dim);
      GradSolhatVAR_old[i].resize(dim);
    }

    vector < bool> solidmark;
    vector < double > phi;
    vector < double > phi_hat;
    vector < double > phi_old;
    vector < adept::adouble> gradphi;
    vector < double > gradphi_hat;
    vector < double > gradphi_old;
    vector < adept::adouble> nablaphi;
    vector < double > nablaphi_hat;
    vector < double > nablaphi_old;


    phi.reserve(max_size);
    solidmark.reserve(max_size);
    phi_hat.reserve(max_size);
    phi_old.reserve(max_size);

    gradphi.reserve(max_size * dim);
    gradphi_hat.reserve(max_size * dim);
    gradphi_old.reserve(max_size * dim);

    nablaphi.reserve(max_size * 3 * (dim - 1));
    nablaphi_hat.reserve(max_size * 3 * (dim - 1));
    nablaphi_old.reserve(max_size * 3 * (dim - 1));

    const double* phi1;

    adept::adouble Weight = 0.;
    adept::adouble Weight_nojac = 0.;
    double Weight_hat = 0.;
    double Weight_old = 0.;

    vector <vector < adept::adouble> > vx(dim);
    vector <vector < double> > vx_hat(dim);
    vector <vector < double> > vx_old(dim);

    vector <vector < adept::adouble > > vx_face(dim);
    vector <vector < double > > vx_face_old(dim);

    for (int i = 0; i < dim; i++) {
      vx[i].reserve(max_size);
      vx_hat[i].reserve(max_size);
      vx_old[i].reserve(max_size);

      vx_face[i].resize(9);
      vx_face_old[i].resize(9);
    }

    vector< vector< adept::adouble > > Soli(2 * dim + 1);
    vector< vector< double > > Soli_old(2 * dim + 1);
    vector< vector< int > > dofsVAR(2 * dim + 1);

    for (int i = 0; i < 2 * dim + 1; i++) {
      Soli[i].reserve(max_size);
      Soli_old[i].reserve(max_size);
      dofsVAR[i].reserve(max_size);
    }

    vector< vector< double > > Rhs(2 * dim + 1);
    vector< vector< adept::adouble > > aRhs(2 * dim + 1);

    for (int i = 0; i < 2 * dim + 1; i++) {
      aRhs[i].reserve(max_size);
      Rhs[i].reserve(max_size);
    }


    vector < int > dofsAll;
    dofsAll.reserve(max_size * (2 + dim + 1));

    vector < double > KKloc;
    KKloc.reserve(dim * max_size * (2 * dim + 1)*dim * max_size * (2 * dim + 1));

    vector < double > Jac;
    Jac.reserve(dim * max_size * (2 * dim + 1)*dim * max_size * (2 * dim + 1));

    // ------------------------------------------------------------------------
    // Physical parameters
    double rhof	 	= ml_prob.parameters.get<Fluid>("Fluid").get_density();
    double rhos	 	= ml_prob.parameters.get<Fluid>("Solid").get_density() / rhof;
    double mu_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_shear_modulus();
    double lambda_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_lambda();
    double mus		= mu_lame / rhof;
    double IRe 		= ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
    double lambda	= lambda_lame / rhof;
    double betans	= 1.;
    int    solid_model	= ml_prob.parameters.get<Solid>("Solid").get_physical_model();

    if (solid_model >= 2 && solid_model <= 4) {
      std::cout << "Error! Solid Model " << solid_model << "not implemented\n";
      abort();
    }

    bool incompressible = (0.5 == ml_prob.parameters.get<Solid>("Solid").get_poisson_coeff()) ? 1 : 0;
    const bool penalty = ml_prob.parameters.get<Solid>("Solid").get_if_penalty();

    // gravity
    double _gravity[3] = {0., 0., 0.};

    double dt =  my_nnlin_impl_sys.GetIntervalTime();
    double time =  my_nnlin_impl_sys.GetTime();
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
    vector <unsigned> indexVAR(2 * dim + 1);
    vector <unsigned> indVAR(2 * dim + 1);
    vector <unsigned> SolType(2 * dim + 1);

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

    myKK->zero();

    // *** element loop ***
    for (int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++) {

      short unsigned ielt = mymsh->GetElementType(iel);
      unsigned nve        = mymsh->GetElementDofNumber(iel, SolType2);
      unsigned nve1       = mymsh->GetElementDofNumber(iel, SolType1);
      int flag_mat        = mymsh->GetElementMaterial(iel);

      // *******************************************************************************************************

      //initialization of everything is in common fluid and solid

      //Rhs
      for (int i = 0; i < 2 * dim; i++) {
        dofsVAR[i].resize(nve);
        Soli[indexVAR[i]].resize(nve);
        Soli_old[indexVAR[i]].resize(nve);

        aRhs[indexVAR[i]].resize(nve);
        Rhs[indexVAR[i]].resize(nve);
      }

      dofsVAR[2 * dim].resize(nve1);
      Soli[indexVAR[2 * dim]].resize(nve1);
      Soli_old[indexVAR[2 * dim]].resize(nve1);
      aRhs[indexVAR[2 * dim]].resize(nve1);
      Rhs[indexVAR[2 * dim]].resize(nve1);

      dofsAll.resize(0);

      KKloc.resize((2 * dim * nve + nve1) * (2 * dim * nve + nve1));
      Jac.resize((2 * dim * nve + nve1) * (2 * dim * nve + nve1));

      // ----------------------------------------------------------------------------------------
      // coordinates, solutions, displacement, velocity dofs

      solidmark.resize(nve);

      for (int i = 0; i < dim; i++) {
        vx[i].resize(nve);
        vx_hat[i].resize(nve);
        vx_old[i].resize(nve);
      }

      for (unsigned i = 0; i < nve; i++) {
        unsigned idof = mymsh->GetSolutionDof(i, iel, SolType2);
        // flag to know if the node "idof" lays on the fluid-solid interface
        solidmark[i] = mymsh->GetSolidMark(idof); // to check

        for (int j = 0; j < dim; j++) {
          Soli[indexVAR[j]][i]     = (*mysolution->_Sol[indVAR[j]])(idof);
          Soli[indexVAR[j + dim]][i] = (*mysolution->_Sol[indVAR[j + dim]])(idof);

          Soli_old[indexVAR[j]][i]     = (*mysolution->_SolOld[indVAR[j]])(idof);
          Soli_old[indexVAR[j + dim]][i] = (*mysolution->_SolOld[indVAR[j + dim]])(idof);

          aRhs[indexVAR[j]][i]     = 0.;
          aRhs[indexVAR[j + dim]][i] = 0.;

          //Fixed coordinates (Reference frame)
          vx_hat[j][i] = (*mymsh->_topology->_Sol[j])(idof);
          // displacement dofs
          dofsVAR[j][i] = myLinEqSolver->GetSystemDof(indVAR[j], indexVAR[j], i, iel);
          // velocity dofs
          dofsVAR[j + dim][i] = myLinEqSolver->GetSystemDof(indVAR[j + dim], indexVAR[j + dim], i, iel);
        }
      }

      // pressure dofs
      for (unsigned i = 0; i < nve1; i++) {
        unsigned idof = mymsh->GetSolutionDof(i, iel, SolType[2 * dim]);
        dofsVAR[2 * dim][i] = myLinEqSolver->GetSystemDof(indVAR[2 * dim], indexVAR[2 * dim], i, iel);
        Soli[indexVAR[2 * dim]][i]     = (*mysolution->_Sol[indVAR[2 * dim]])(idof);
        Soli_old[indexVAR[2 * dim]][i] = (*mysolution->_SolOld[indVAR[2 * dim]])(idof);
        aRhs[indexVAR[2 * dim]][i] = 0.;
      }

      // build dof ccomposition
      for (int idim = 0; idim < 2 * dim; idim++) {
        dofsAll.insert(dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end());
      }

      dofsAll.insert(dofsAll.end(), dofsVAR[2 * dim].begin(), dofsVAR[2 * dim].end());

      //if (igrid==gridn || !myel->GetRefinedElementIndex(iel) ) {

      s.new_recording();

      for (unsigned idim = 0; idim < dim; idim++) {
        for (int j = 0; j < nve; j++) {
          vx[idim][j]    = vx_hat[idim][j] + Soli[indexVAR[idim]][j];
          vx_old[idim][j] = vx_hat[idim][j] + Soli_old[indexVAR[idim]][j];
        }
      }

      // Boundary integral
      {
        double tau = 0.;
        vector < adept::adouble> normal(dim, 0);
        vector < double > normal_old(dim, 0);

        // loop on faces
        for (unsigned jface = 0; jface < mymsh->GetElementFaceNumber(iel); jface++) {
          std::vector< double > xx(dim, 0.);

          // look for boundary faces
          if (myel->GetFaceElementIndex(iel, jface) < 0) {

            unsigned int face = -(mymsh->el->GetFaceElementIndex(iel, jface) + 1);

            if (!ml_sol->GetBdcFunction()(xx, "U", tau, face, time) && tau != 0.) {
              double tau_old;
              ml_sol->GetBdcFunction()(xx, "U", tau_old, face, time - dt);
              unsigned nve = mymsh->GetElementFaceDofNumber(iel, jface, SolType2);
              const unsigned felt = mymsh->GetElementFaceType(iel, jface);

              for (unsigned i = 0; i < nve; i++) {
                unsigned int ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);
                unsigned idof = mymsh->GetSolutionDof(ilocal, iel, 2);

                for (unsigned idim = 0; idim < dim; idim++) {
                  vx_face[idim][i]    = (*mymsh->_topology->_Sol[idim])(idof) + Soli[indexVAR[idim]][ilocal];
                  vx_face_old[idim][i] = (*mymsh->_topology->_Sol[idim])(idof) + Soli_old[indexVAR[idim]][ilocal];
                }
              }

              for (unsigned igs = 0; igs < mymsh->_finiteElement[felt][SolType2]->GetGaussPointNumber(); igs++) {
                mymsh->_finiteElement[felt][SolType2]->JacobianSur(vx_face, igs, Weight, phi, gradphi, normal);
                mymsh->_finiteElement[felt][SolType2]->JacobianSur(vx_face_old, igs, Weight_old, phi_old, gradphi_old, normal_old);

                //phi1 =mymsh->_finiteElement[felt][SolType2]->GetPhi(igs);
                // *** phi_i loop ***
                for (unsigned i = 0; i < nve; i++) {
                  adept::adouble value = - phi[i] * tau / rhof * Weight;
                  double value_old = - phi[i] * tau_old / rhof * Weight_old;
                  unsigned int ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);

                  for (unsigned idim = 0; idim < dim; idim++) {
                    if ((!solidmark[ilocal])) {
                      aRhs[indexVAR[dim + idim]][ilocal]   += dt*(theta*value + (1. - theta)*value_old) * normal[idim];
                    }
                    else { //if interface node it goes to solid
                      aRhs[indexVAR[idim]][ilocal]   += dt*(theta*value + (1. - theta)*value_old) * normal[idim];
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

      for (unsigned ig = 0; ig < mymsh->_finiteElement[ielt][SolType2]->GetGaussPointNumber(); ig++) {
        // *** get Jacobian and test function and test function derivatives in the moving frame***
        mymsh->_finiteElement[ielt][SolType2]->Jacobian(vx,  ig, Weight,     phi,     gradphi,     nablaphi);
        mymsh->_finiteElement[ielt][SolType2]->Jacobian(vx_hat, ig, Weight_hat, phi_hat, gradphi_hat, nablaphi_hat);
        mymsh->_finiteElement[ielt][SolType2]->Jacobian(vx_old, ig, Weight_old, phi_old, gradphi_old, nablaphi_old);
        phi1 = mymsh->_finiteElement[ielt][SolType1]->GetPhi(ig);

        if (flag_mat == 2 || iel == mymsh->_elementOffset[iproc]) {
          if (ig == 0) {
            double GaussWeight = mymsh->_finiteElement[ielt][SolType2]->GetGaussWeight(ig);
            area = Weight_hat / GaussWeight;

            if (iel == mymsh->_elementOffset[iproc]) {
              area_elem_first->add(mymsh->processor_id(), area);
              area_elem_first->close();
              rapresentative_area = area_elem_first->l1_norm() / nprocs;
            }
          }

          Weight_nojac = Weight / area * rapresentative_area;
        }

        // ---------------------------------------------------------------------------
        // displacement and velocity
        for (int i = 0; i < 2 * dim; i++) {
          SolVAR[i] = 0.;
          SolVAR_old[i] = 0.;

          for (int j = 0; j < dim; j++) {
            GradSolVAR[i][j] = 0.;
            GradSolVAR_old[i][j] = 0.;

            GradSolhatVAR[i][j] = 0.;
            GradSolhatVAR_old[i][j] = 0.;
          }

          for (unsigned inode = 0; inode < nve; inode++) {
            SolVAR[i] += phi[inode] * Soli[indexVAR[i]][inode];
            SolVAR_old[i] += phi_old[inode] * Soli_old[indexVAR[i]][inode];

            for (int j = 0; j < dim; j++) {
              GradSolVAR[i][j]     += gradphi[inode * dim + j]    * Soli[indexVAR[i]][inode];
              GradSolVAR_old[i][j] += gradphi_old[inode * dim + j] * Soli_old[indexVAR[i]][inode];

              GradSolhatVAR[i][j]     += gradphi_hat[inode * dim + j] * Soli[indexVAR[i]][inode];
              GradSolhatVAR_old[i][j] += gradphi_hat[inode * dim + j] * Soli_old[indexVAR[i]][inode];
            }
          }
        }

        // pressure
        SolVAR[2 * dim] = 0.;
        SolVAR_old[2 * dim] = 0.;

        for (unsigned inode = 0; inode < nve1; inode++) {
          SolVAR[2 * dim]    += phi1[inode] * Soli[indexVAR[2 * dim]][inode];
          SolVAR_old[2 * dim] += phi1[inode] * Soli_old[indexVAR[2 * dim]][inode];
        }

        // Lagrangian mesh velocity at t = time + dt/2
        vector < adept::adouble > meshVel(dim);
        vector < vector < adept::adouble > > GradMeshVel(dim);

        for (unsigned i = 0; i < dim; i++) {
          meshVel[i] = (SolVAR[i] - SolVAR_old[i]) / dt;
          GradMeshVel[i].resize(dim);

          for (unsigned j = 0; j < dim; j++) {
            GradMeshVel[i][j] = (GradSolVAR[i][j] - GradSolVAR_old[i][j]) / dt;
          }
        }

        // ---------------------------------------------------------------------------
        //BEGIN FLUID ASSEMBLY ============
        if (flag_mat == 2) {
          //BEGIN ALE + Momentum (Navier-Stokes)
          {
            for (unsigned i = 0; i < nve; i++) {

              //BEGIN redidual Laplacian ALE map in the reference domain
              adept::adouble LapmapVAR[3] = {0., 0., 0.};

              for (int idim = 0; idim < dim; idim++) {
                for (int jdim = 0; jdim < dim; jdim++) {
                  LapmapVAR[idim] += (GradSolVAR[idim][jdim] + GradSolVAR[jdim][idim]) * gradphi[i * dim + jdim];
                }
              }

              for (int idim = 0; idim < dim; idim++) {
                aRhs[indexVAR[idim]][i] += (!solidmark[i]) * (- LapmapVAR[idim] * Weight_nojac);
              }

              //END redidual Laplacian ALE map in reference domain

              //BEGIN redidual Navier-Stokes in moving domain
              adept::adouble LapvelVAR[3] = {0., 0., 0.};
              adept::adouble LapvelVAR_old[3] = {0., 0., 0.};
              adept::adouble AdvaleVAR[3] = {0., 0., 0.};
              adept::adouble AdvaleVAR_old[3] = {0., 0., 0.};

              for (int idim = 0.; idim < dim; idim++) {
                for (int jdim = 0.; jdim < dim; jdim++) {

                  LapvelVAR[idim]     += (GradSolVAR[dim + idim][jdim] + GradSolVAR[dim + jdim][idim]) * gradphi[i * dim + jdim];
                  LapvelVAR_old[idim] += (GradSolVAR_old[dim + idim][jdim] + GradSolVAR_old[dim + jdim][idim]) * gradphi_old[i * dim + jdim];

                  AdvaleVAR[idim]	+= ((SolVAR[dim + jdim] - meshVel[jdim]) * GradSolVAR[dim + idim][jdim]
                                      + (GradSolVAR[dim + jdim][jdim] - GradMeshVel[jdim][jdim]) * SolVAR[dim + idim]
                                     ) * phi[i];

                  AdvaleVAR_old[idim]	+= ((SolVAR_old[dim + jdim] - meshVel[jdim]) * GradSolVAR_old[dim + idim][jdim]
                                          + (GradSolVAR_old[dim + jdim][jdim] - GradMeshVel[jdim][jdim]) * SolVAR_old[dim + idim]
                                         ) * phi_old[i];
                }
              }

              for (int idim = 0; idim < dim; idim++) {

                adept::adouble timeDerivative = -(SolVAR[dim + idim] * phi[i] * Weight
                                                  - SolVAR_old[dim + idim] * phi_old[i] * Weight_old);

                adept::adouble value =  theta * dt * (
                                          - AdvaleVAR[idim]      	             // advection term
                                          - IRe * LapvelVAR[idim]	             // viscous dissipation
                                          + SolVAR[2 * dim] * gradphi[i * dim + idim] // pressure gradient
                                        ) * Weight;                                // at time t

                adept::adouble value_old = (1. - theta) * dt * (
                                              - AdvaleVAR_old[idim]               	         // advection term
                                              - IRe * LapvelVAR_old[idim]	       	         // viscous dissipation
                                              + SolVAR[2 * dim] * gradphi_old[i * dim + idim]  // pressure gradient
                                            ) * Weight_old;			                 // at time t-dt

                if (!solidmark[i]) {
                  aRhs[indexVAR[dim + idim]][i] += timeDerivative + value + value_old;
                }
                else {
                  aRhs[indexVAR[idim]][i] += timeDerivative + value + value_old;
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
              aRhs[indexVAR[2 * dim]][i] += -(-phi1[i] * div_vel) * Weight;
            }
          }
          //END continuity block ===========================
        }
        //END FLUID ASSEMBLY ============

        //*******************************************************************************************************

        //BEGIN SOLID ASSEMBLY ============
        else {
          //BEGIN build Chauchy Stress in moving domain
          //physical quantity
          adept::adouble J_hat;
          double J_hat_old;
          adept::adouble I_e;
          double I_e_old;
          adept::adouble Cauchy[3][3];
          adept::adouble Cauchy_old[3][3];
          double Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};




          if (solid_model == 0) { // Saint-Venant
            adept::adouble e[3][3];
            double e_old[3][3];

            //computation of the stress tensor
            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                e[i][j]    = 0.5 * (GradSolhatVAR[i][j]    + GradSolhatVAR[j][i]);
                e_old[i][j] = 0.5 * (GradSolhatVAR_old[i][j] + GradSolhatVAR_old[j][i]);
              }
            }

            I_e = 0;
            I_e_old = 0;

            for (int i = 0; i < dim; i++) {
              I_e     += e[i][i];
              I_e_old += e_old[i][i];
            }

            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                //incompressible
                Cauchy[i][j]     = 2 * mus * e[i][j] - 2 * mus * I_e * SolVAR[2 * dim] * Id2th[i][j];
                Cauchy_old[i][j] = 2 * mus * e_old[i][j] - 2 * mus * I_e_old * SolVAR_old[2 * dim] * Id2th[i][j];
                //+(penalty)*lambda*I_e*Id2th[i][j];
              }
            }
          }

          else { // hyperelastic non linear material

            adept::adouble I1_B = 0.;
            adept::adouble I2_B = 0.;

            double I1_B_old = 0.;
            double I2_B_old = 0.;

            adept::adouble F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
            adept::adouble B[3][3];

            double F_old[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
            double B_old[3][3];

            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                F[i][j] += GradSolhatVAR[i][j];
                F_old[i][j] += GradSolhatVAR_old[i][j];
              }
            }

            J_hat =   F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                      - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

            J_hat_old =    F_old[0][0] * F_old[1][1] * F_old[2][2] + F_old[0][1] * F_old[1][2] * F_old[2][0]
                           + F_old[0][2] * F_old[1][0] * F_old[2][1] - F_old[2][0] * F_old[1][1] * F_old[0][2]
                           - F_old[2][1] * F_old[1][2] * F_old[0][0] - F_old[2][2] * F_old[1][0] * F_old[0][1];


            for (int I = 0; I < 3; ++I) {
              for (int J = 0; J < 3; ++J) {
                B[I][J] = 0.;
                B_old[I][J] = 0.;

                for (int K = 0; K < 3; ++K) {
                  //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
                  B[I][J]     += F[I][K] * F[J][K];
                  B_old[I][J] += F_old[I][K] * F_old[J][K];
                }
              }
            }

            if (solid_model <= 4) { // Neo-Hookean
              I1_B = B[0][0] + B[1][1] + B[2][2];
              I1_B_old = B_old[0][0] + B_old[1][1] + B_old[2][2];

              for (int I = 0; I < 3; ++I) {
                for (int J = 0; J < 3; ++J) {
                  if (1 == solid_model) {   //Wood-Bonet J_hat  =1;
                    Cauchy[I][J]     = mus * (B[I][J]     - Id2th[I][J]) - mus / 3.*I1_B     * SolVAR[2 * dim] * Id2th[I][J];
                    Cauchy_old[I][J] = mus * (B_old[I][J] - Id2th[I][J]) - mus / 3.*I1_B_old * SolVAR[2 * dim] * Id2th[I][J];
                  }

// 		    else if ( 2 == solid_model ) Cauchy[I][J] = mus/J_hat*B[I][J]
// 							       -mus/J_hat*SolVAR[2*dim]*Id2th[I][J];    //Wood-Bonet J_hat !=1;
// 		    else if ( 3 == solid_model ) Cauchy[I][J] = mus*(B[I][J] - Id2th[I][J])/J_hat
// 							      + lambda/J_hat*log(J_hat)*Id2th[I][J]; 	//Wood-Bonet penalty
// 		    else if ( 4 == solid_model ) Cauchy[I][J] = mus*(B[I][J] - I1_B*Id2th[I][J]/3.)/pow(J_hat,5./3.)
// 						              + lambda*(J_hat-1.)*Id2th[I][J];  	  //Allan-Bower
//
                }
              }
            }
            else if (5 == solid_model) {  //Mooney-Rivlin
              adept::adouble detB =   B[0][0] * (B[1][1] * B[2][2] - B[2][1] * B[1][2])
                                      - B[0][1] * (B[2][2] * B[1][0] - B[1][2] * B[2][0])
                                      + B[0][2] * (B[1][0] * B[2][1] - B[2][0] * B[1][1]);

              double detB_old =  B_old[0][0] * (B_old[1][1] * B_old[2][2] - B_old[2][1] * B_old[1][2])
                                 - B_old[0][1] * (B_old[2][2] * B_old[1][0] - B_old[1][2] * B_old[2][0])
                                 + B_old[0][2] * (B_old[1][0] * B_old[2][1] - B_old[2][0] * B_old[1][1]);

              adept::adouble invdetB = 1. / detB;
              double invdetB_old = 1. / detB_old;
              adept::adouble invB[3][3];
              double invB_old[3][3];

              invB[0][0] = (B[1][1] * B[2][2] - B[1][2] * B[2][1]) * invdetB;
              invB[1][0] = -(B[0][1] * B[2][2] - B[0][2] * B[2][1]) * invdetB;
              invB[2][0] = (B[0][1] * B[1][2] - B[0][2] * B[1][1]) * invdetB;
              invB[0][1] = -(B[1][0] * B[2][2] - B[1][2] * B[2][0]) * invdetB;
              invB[1][1] = (B[0][0] * B[2][2] - B[0][2] * B[2][0]) * invdetB;
              invB[2][1] = -(B[0][0] * B[1][2] - B[1][0] * B[0][2]) * invdetB;
              invB[0][2] = (B[1][0] * B[2][1] - B[2][0] * B[1][1]) * invdetB;
              invB[1][2] = -(B[0][0] * B[2][1] - B[2][0] * B[0][1]) * invdetB;
              invB[2][2] = (B[0][0] * B[1][1] - B[1][0] * B[0][1]) * invdetB;

              invB_old[0][0] = (B_old[1][1] * B_old[2][2] - B_old[1][2] * B_old[2][1]) * invdetB_old;
              invB_old[1][0] = -(B_old[0][1] * B_old[2][2] - B_old[0][2] * B_old[2][1]) * invdetB_old;
              invB_old[2][0] = (B_old[0][1] * B_old[1][2] - B_old[0][2] * B_old[1][1]) * invdetB_old;
              invB_old[0][1] = -(B_old[1][0] * B_old[2][2] - B_old[1][2] * B_old[2][0]) * invdetB_old;
              invB_old[1][1] = (B_old[0][0] * B_old[2][2] - B_old[0][2] * B_old[2][0]) * invdetB_old;
              invB_old[2][1] = -(B_old[0][0] * B_old[1][2] - B_old[1][0] * B_old[0][2]) * invdetB_old;
              invB_old[0][2] = (B_old[1][0] * B_old[2][1] - B_old[2][0] * B_old[1][1]) * invdetB_old;
              invB_old[1][2] = -(B_old[0][0] * B_old[2][1] - B_old[2][0] * B_old[0][1]) * invdetB_old;
              invB_old[2][2] = (B_old[0][0] * B_old[1][1] - B_old[1][0] * B_old[0][1]) * invdetB_old;


              I1_B = B[0][0] + B[1][1] + B[2][2];
              I2_B = B[0][0] * B[1][1] + B[1][1] * B[2][2] + B[2][2] * B[0][0]
                     - B[0][1] * B[1][0] - B[1][2] * B[2][1] - B[2][0] * B[0][2];

              I1_B_old = B_old[0][0] + B_old[1][1] + B_old[2][2];
              I2_B_old = B_old[0][0] * B_old[1][1] + B_old[1][1] * B_old[2][2] + B_old[2][2] * B_old[0][0]
                         - B_old[0][1] * B_old[1][0] - B_old[1][2] * B_old[2][1] - B_old[2][0] * B_old[0][2];

              double C1 = mus / 3.;
              double C2 = C1 / 2.;

              for (int I = 0; I < 3; ++I) {
                for (int J = 0; J < 3; ++J) {
                  Cauchy[I][J] =  2.*(C1 * B[I][J] - C2 * invB[I][J])
                                  //- (2. / 3.) * (C1 * I1_B - C2 * I2_B) * SolVAR[2 * dim] * Id2th[I][J];
                                  - SolVAR[2 * dim] * Id2th[I][J];

                  Cauchy_old[I][J] =  2.* (C1 * B_old[I][J] - C2 * invB_old[I][J])
                                      //- (2. / 3.) * (C1 * I1_B_old - C2 * I2_B_old) * SolVAR[2 * dim] * Id2th[I][J];
                                      - SolVAR[2 * dim] * Id2th[I][J];
                  
                }
              }

            }
          }

          //END build Chauchy Stress in moving domain

          //BEGIN v=0 + Momentum (Solid)
          {
            for (unsigned i = 0; i < nve; i++) {

              //BEGIN redidual d_t - v = 0 in fixed domain
              for (int idim = 0; idim < dim; idim++) {
                aRhs[indexVAR[dim + idim]][i] +=  - phi[i] * (-SolVAR[idim] + SolVAR_old[idim] +
                                                        dt * (theta*SolVAR[dim + idim] + (1. - theta)*SolVAR_old[dim + idim])
                                                             ) * Weight_hat;
              }

              //END redidual d_t - v = 0 in fixed domain

              //BEGIN redidual Solid Momentum in moving domain
              adept::adouble CauchyDIR[3] = {0., 0., 0.};
              adept::adouble CauchyDIR_old[3] = {0., 0., 0.};

              for (int idim = 0.; idim < dim; idim++) {
                for (int jdim = 0.; jdim < dim; jdim++) {
                  CauchyDIR[idim] += gradphi[i * dim + jdim] * Cauchy[idim][jdim];
                  CauchyDIR_old[idim] += gradphi_old[i * dim + jdim] * Cauchy_old[idim][jdim];
                }
              }

              for (int idim = 0; idim < dim; idim++) {

                adept::adouble timeDerivative = -(rhos * SolVAR[dim + idim] * phi[i] * Weight
                                                  - rhos * SolVAR_old[dim + idim] * phi_old[i] * Weight_old);

                adept::adouble value =  theta * dt * (rhos * phi[i] * _gravity[idim]      // body force
                                                    - CauchyDIR[idim]			  // stress
                                                   ) * Weight;                         // at time t

                adept::adouble value_old =  (1. - theta) * dt * (rhos * phi_old[i] * _gravity[idim]     // body force
                                                        - CauchyDIR_old[idim]			 // stress
                                                       ) * Weight_old;                         // at time t-dt

                aRhs[indexVAR[idim]][i] += timeDerivative + value + value_old;
              }

              //END redidual Solid Momentum in moving domain
            }
          }
          //END v=0 + Momentum (Solid)

          //BEGIN continuity block
          {
            for (unsigned i = 0; i < nve1; i++) {
              if (!penalty) {
                if (0 == solid_model) {
                  aRhs[indexVAR[2 * dim]][i] += -(-phi1[i] * (I_e + (!incompressible) / lambda * SolVAR[2 * dim])) * Weight_hat;
                }
                else if (1 == solid_model || 5 == solid_model) {
                  aRhs[indexVAR[2 * dim]][i] += phi1[i] * (J_hat - 1. + (!incompressible) / lambda * SolVAR[2 * dim]) * Weight_hat;
                }

// 		  else if (2 == solid_model){
// 		    aRhs[indexVAR[2*dim]][i] += -(-phi1[i]*( log(J_hat)/J_hat + (!incompressible)/lambda*SolVAR[2*dim] ) )*Weight_hat;
// 		  }
              }

// 		else if (3 == solid_model || 4 == solid_model){ // pressure = 0 in the solid
// 		  aRhs[indexVAR[2*dim]][i] += -(-phi1[i]*( SolVAR[2*dim] ) )*Weight_hat;
// 		}
            }
          }
          //END continuity block
        }

        //END SOLID ASSEMBLY ============
      }

      //}

      //BEGIN local to global assembly
      //copy adouble aRhs into double Rhs
      for (unsigned i = 0; i < 2 * dim; i++) {
        for (int j = 0; j < nve; j++) {
          Rhs[indexVAR[i]][j] = aRhs[indexVAR[i]][j].value();
        }
      }

      for (unsigned j = 0; j < nve1; j++) {
        Rhs[indexVAR[2 * dim]][j] = aRhs[indexVAR[2 * dim]][j].value();
      }

      for (int i = 0; i < 2 * dim + 1; i++) {
        myRES->add_vector_blocked(Rhs[indexVAR[i]], dofsVAR[i]);
      }

      //Store equations
      for (int i = 0; i < 2 * dim; i++) {
        s.dependent(&aRhs[indexVAR[i]][0], nve);
        s.independent(&Soli[indexVAR[i]][0], nve);
      }

      s.dependent(&aRhs[indexVAR[2 * dim]][0], nve1);
      s.independent(&Soli[indexVAR[2 * dim]][0], nve1);
      s.jacobian(&Jac[0]);
      unsigned nveAll = (2 * dim * nve + nve1);

      for (int inode = 0; inode < nveAll; inode++) {
        for (int jnode = 0; jnode < nveAll; jnode++) {
          KKloc[inode * nveAll + jnode] = -Jac[jnode * nveAll + inode];
        }
      }

      myKK->add_matrix_blocked(KKloc, dofsAll, dofsAll);
      s.clear_independents();
      s.clear_dependents();

      //END local to global assembly

    } //end list of elements loop

    myKK->close();
    myRES->close();

    delete area_elem_first;

    // *************************************
    end_time = clock();
    AssemblyTime += (end_time - start_time);
    // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************
  }


}

#endif
