#ifndef __femus_include_IncompressibleFSIAssemblySupg_hpp__
#define __femus_include_IncompressibleFSIAssemblySupg_hpp__
#endif

#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "MultiLevelSolution.hpp"
#include "adept.h"

#include "OprtrTypeEnum.hpp"

namespace femus
{

  bool meshIsCurrupted = true;
  //double factordxi[3] =  {1.e7, 1.e7, 1.e7};
  double factordxi[3] =  {1.e13, 1.e13, 1.e13};


  void FSIConstrainLeaflet(MultiLevelSolution& mlSol);

  void FSITimeDependentAssemblySupg(MultiLevelProblem & ml_prob)
  {

    clock_t AssemblyTime = 0;
    clock_t start_time, end_time;

    //pointers and references

    //MonolithicFSINonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<MonolithicFSINonLinearImplicitSystem>("Fluid-Structure-Interaction");
    TransientNonlinearImplicitSystem & my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("Fluid-Structure-Interaction");
    const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();

    MultiLevelSolution	* ml_sol		= ml_prob._ml_sol;

    //FSIConstrainLeaflet(*ml_sol);

    Solution	* mysolution		= ml_sol->GetSolutionLevel(level);

    LinearEquationSolver * myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];
    Mesh	*	mymsh		=  ml_prob._ml_msh->GetLevel(level);
    elem	*	myel		=  mymsh->el;
    SparseMatrix	* myKK		=  myLinEqSolver->_KK;
    NumericVector	* myRES		=  myLinEqSolver->_RES;


    bool assembleMatrix = my_nnlin_impl_sys.GetAssembleMatrix();
    // call the adept stack object
    adept::Stack& s = FemusInit::_adeptStack;
    if (assembleMatrix) s.continue_recording();
    else s.pause_recording();

    const unsigned dim = mymsh->GetDimension();
    const unsigned nabla_dim = 3 * (dim - 1);
    const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

    double theta = 0.5;
    //double theta = 1.;

    // local objects
    vector<adept::adouble> SolVAR(2 * dim + 1);
    vector<double> SolVAR_old(2 * dim);

    vector<vector < adept::adouble > > GradSolVAR(2 * dim);
    vector<vector < double > > GradSolVAR_old(2 * dim);

    vector<vector < adept::adouble > > GradSolhatVAR(2 * dim);
    vector<vector < double > > GradSolhatVAR_old(2 * dim);

    vector<vector<adept::adouble> > NablaSolVAR(2 * dim);
    vector<vector < double > > NablaSolVAR_old(2 * dim);

    for (int i = 0; i < 2 * dim; i++) {
      GradSolVAR[i].resize(dim);
      GradSolVAR_old[i].resize(dim);

      GradSolhatVAR[i].resize(dim);
      GradSolhatVAR_old[i].resize(dim);

      NablaSolVAR[i].resize(nabla_dim);
      NablaSolVAR_old[i].resize(nabla_dim);
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

    const double * phi1;

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

    //vector < double > KKloc;
    //KKloc.reserve ( dim * max_size * ( 2 * dim + 1 ) *dim * max_size * ( 2 * dim + 1 ) );

    vector < double > Jac;
    Jac.reserve(dim * max_size * (2 * dim + 1) *dim * max_size * (2 * dim + 1));

    // ------------------------------------------------------------------------
    // Physical parameters
    double rhof	 	= ml_prob.parameters.get<Fluid> ("Fluid").get_density();
    double rhos	 	= ml_prob.parameters.get<Fluid> ("Solid").get_density() / rhof;
    double mu_lame 	= ml_prob.parameters.get<Solid> ("Solid").get_lame_shear_modulus();
    double lambda_lame 	= ml_prob.parameters.get<Solid> ("Solid").get_lame_lambda();
    double mus		= mu_lame / rhof;
    double mu_lame1 	= ml_prob.parameters.get < Solid> ("Solid1").get_lame_shear_modulus();
    double mus1 	= mu_lame1 / rhof;
    double IRe 		= ml_prob.parameters.get<Fluid> ("Fluid").get_IReynolds_number();
    double lambda	= lambda_lame / rhof;
    double betans	= 1.;
    int    solid_model	= ml_prob.parameters.get<Solid> ("Solid").get_physical_model();

    if (solid_model >= 2 && solid_model <= 4) {
      std::cout << "Error! Solid Model " << solid_model << "not implemented\n";
      abort();
    }

    bool incompressible = (0.5 == ml_prob.parameters.get<Solid> ("Solid").get_poisson_coeff()) ? 1 : 0;
    const bool penalty = ml_prob.parameters.get<Solid> ("Solid").get_if_penalty();

    // gravity
    double _gravity[3] = {0., 0., 0.};

    double dt =  my_nnlin_impl_sys.GetIntervalTime();
    double time =  my_nnlin_impl_sys.GetTime();
    // -----------------------------------------------------------------
    // space discretization parameters
    unsigned SolType2 = ml_sol->GetSolutionType(ml_sol->GetIndex("U"));
    unsigned SolType1 = ml_sol->GetSolutionType(ml_sol->GetIndex("PS"));

    // mesh and procs
    unsigned nel    = mymsh->GetNumberOfElements();
    unsigned igrid  = mymsh->GetLevel();
    unsigned iproc  = mymsh->processor_id();

    unsigned indLmbd = ml_sol->GetIndex("lmbd");

    //----------------------------------------------------------------------------------
    //variable-name handling
    const char varname[7][3] = {"DX", "DY", "DZ", "U", "V", "W", "PS"};
    //const char varname1[3][4] = {"DX1", "DY1", "DZ1"};
    vector <unsigned> indexVAR(2 * dim + 1);
    vector <unsigned> indVAR(2 * dim + 1);
    //vector <unsigned> indVAR1(dim);
    vector <unsigned> SolType(2 * dim + 1);

    for (unsigned ivar = 0; ivar < dim; ivar++) {
      indVAR[ivar] = ml_sol->GetIndex(&varname[ivar][0]);
      //indVAR1[ivar] = ml_sol->GetIndex(&varname1[ivar][0]);
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
    NumericVector * area_elem_first;
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

    if (assembleMatrix) myKK->zero();

    // *** element loop ***
    for (int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++) {

      short unsigned ielt = mymsh->GetElementType(iel);
      unsigned nve        = mymsh->GetElementDofNumber(iel, SolType2);
      unsigned nve1       = mymsh->GetElementDofNumber(iel, SolType1);
      int flag_mat        = mymsh->GetElementMaterial(iel);
      unsigned elementGroup = mymsh->GetElementGroup(iel);

      // *******************************************************************************************************

      //initialization of everything is in common fluid and solid

      //Rhs
      for (int i = 0; i < 2 * dim; i++) {
        dofsVAR[i].resize(nve);
        Soli[indexVAR[i]].resize(nve);
        Soli_old[indexVAR[i]].resize(nve);

        aRhs[indexVAR[i]].resize(nve);
        //    Rhs[indexVAR[i]].resize ( nve );
      }

      dofsVAR[2 * dim].resize(nve1);
      Soli[indexVAR[2 * dim]].resize(nve1);
      Soli_old[indexVAR[2 * dim]].resize(nve1);
      aRhs[indexVAR[2 * dim]].resize(nve1);
      //  Rhs[indexVAR[2 * dim]].resize ( nve1 );

      dofsAll.resize(0);

      //  KKloc.resize ( ( 2 * dim * nve + nve1 ) * ( 2 * dim * nve + nve1 ) );
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
        solidmark[i] = mymsh->GetSolidMark(idof);    // to check

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

      if (assembleMatrix) s.new_recording();


      for (int j = 0; j < nve; j++) {
        unsigned idof = mymsh->GetSolutionDof(j, iel, SolType2);
        for (unsigned idim = 0; idim < dim; idim++) {
//           vx[idim][j] = vx_hat[idim][j] + (!meshIsCurrupted) * Soli[indexVAR[idim]][j]
//                         + meshIsCurrupted * (*mysolution->_Sol[indVAR1[idim]])(idof);
	  vx[idim][j] = vx_hat[idim][j] + Soli[indexVAR[idim]][j];
          vx_old[idim][j] = vx_hat[idim][j] + Soli_old[indexVAR[idim]][j];
        }
      }

      //Boundary integral
      {

        vector < adept::adouble> normal(dim, 0);
        vector < double > normal_old(dim, 0);

        // loop on faces
        for (unsigned jface = 0; jface < mymsh->GetElementFaceNumber(iel); jface++) {
          std::vector< double > xx(dim, 0.);

          // look for boundary faces
          if (myel->GetFaceElementIndex(iel, jface) < 0) {

            unsigned int face = - (mymsh->el->GetFaceElementIndex(iel, jface) + 1);
            double tau = 0.;
            double tau_old = 0.;
            if ((!ml_sol->GetBdcFunction()(xx, "PS", tau, face, time) &&
                 !ml_sol->GetBdcFunction()(xx, "PS", tau_old, face, time - dt))
                && (tau != 0. || tau_old != 0.)) {

              unsigned nve = mymsh->GetElementFaceDofNumber(iel, jface, SolType2);
              const unsigned felt = mymsh->GetElementFaceType(iel, jface);

              for (unsigned i = 0; i < nve; i++) {
                unsigned int ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);
                unsigned idof = mymsh->GetSolutionDof(ilocal, iel, 2);

                for (unsigned idim = 0; idim < dim; idim++) {
                  vx_face[idim][i]    = (*mymsh->_topology->_Sol[idim])(idof) + Soli[indexVAR[idim]][ilocal]; //TODO
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
                  double value_old = - phi_old[i] * tau_old / rhof * Weight_old;
                  unsigned int ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);

                  for (unsigned idim = 0; idim < dim; idim++) {
                    if ((!solidmark[ilocal])) {
                      aRhs[indexVAR[dim + idim]][ilocal]   += (theta * value + (1. - theta) * value_old) * normal[idim];
                    }
                    else {   //if interface node it goes to solid
                      aRhs[indexVAR[idim]][ilocal]   += (theta * value + (1. - theta) * value_old) * normal[idim];
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

        if (flag_mat == 2 || flag_mat == 3  || iel == mymsh->_elementOffset[iproc]) {
          if (ig == 0) {
            double GaussWeight = mymsh->_finiteElement[ielt][SolType2]->GetGaussWeight(ig);
            area = Weight_hat / GaussWeight;

            if (iel == mymsh->_elementOffset[iproc]) {
              area_elem_first->set(mymsh->processor_id(), area);
              area_elem_first->close();
              rapresentative_area = area_elem_first->l1_norm() / nprocs;
            }
          }

          Weight_nojac = Weight / area * rapresentative_area;

          //-----------------------------------------------------------------------//
          // for vein_valve mesh
          std::vector<double> xc(dim, 0.);
          //xc[0] = -0.001013; // vein_valve
          //xc[1] = 0.07000;
          xc[0] = -0.000286; // vein_valve_closed
          xc[1] = 0.07000;
          double distance = 0.;
          for (unsigned k = 0; k < dim; k++) {
            distance += (vx_hat[k][ nve - 1] - xc[k]) * (vx_hat[k][nve - 1] - xc[k]);
          }
          distance = sqrt(distance);
          Weight_nojac *= 1. / (1 + 10000 * distance);

          if (elementGroup == 16) Weight_nojac *= 100;
          //-----------------------------------------------------------------------//
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

          for (int j = 0; j < nabla_dim; j++) {
            NablaSolVAR[i][j] = 0.;
            NablaSolVAR_old[i][j] = 0.;
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
            for (int j = 0; j < nabla_dim; j++) {
              NablaSolVAR[i][j]     += nablaphi[inode * nabla_dim + j] * Soli[indexVAR[i]][inode];
              NablaSolVAR_old[i][j] += nablaphi_old[inode * nabla_dim + j] * Soli_old[indexVAR[i]][inode];
            }
          }
        }

        // pressure
        SolVAR[2 * dim] = 0.;

        for (unsigned inode = 0; inode < nve1; inode++) {
          SolVAR[2 * dim]    += phi1[inode] * Soli[indexVAR[2 * dim]][inode];
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
        //BEGIN FLUID and POROUS MEDIA ASSEMBLY ============
        if (flag_mat == 2 || flag_mat == 3) {

          vector < adept::adouble > a(dim);
          vector < adept::adouble > a_old(dim);
          for (int i = 0; i < dim; i++) {
            a[i] = SolVAR[i + dim] - 0.*meshVel[i]; // TODO maybe we subtract meshVel[i] maybe not
            a_old[i] = SolVAR_old[i + dim] - 0.*meshVel[i];
          }

          // speed
          adept::adouble aL2Norm = 0.;
          adept::adouble a_oldL2Norm = 0.;
          for (int i = 0; i < dim; i++) {
            aL2Norm += a[i] * a[i];
            a_oldL2Norm += a_old[i] * a_old[i];
          }
          aL2Norm = sqrt(aL2Norm);
          a_oldL2Norm = sqrt(a_oldL2Norm);

          double sqrtlambdak = (*mysolution->_Sol[indLmbd])(iel);
          adept::adouble tauSupg = 1. / (sqrtlambdak * sqrtlambdak * 4.*IRe);
          adept::adouble tauSupg_old = 1. / (sqrtlambdak * sqrtlambdak * 4.*IRe);
          adept::adouble Rek = aL2Norm / (4.*sqrtlambdak * IRe);
          adept::adouble Rek_old = a_oldL2Norm / (4.*sqrtlambdak * IRe);

          if (Rek > 1.0e-15) {
            adept::adouble xiRek = (Rek >= 1.) ? 1. : Rek;
            tauSupg   = xiRek / (aL2Norm * sqrtlambdak);
          }

          if (Rek_old > 1.0e-15) {
            adept::adouble xiRek_old = (Rek_old >= 1.) ? 1. : Rek_old;
            tauSupg_old   = xiRek_old / (a_oldL2Norm * sqrtlambdak);
          }

          //std::cout << tauSupg.value() <<"  " <<tauSupg_old.value()<<std::endl;


          vector < adept::adouble > phiSupg(nve, 0.);
          vector < adept::adouble > phiSupg_old(nve, 0.);
          //tauSupg=0;
          //tauSupg_old=0;
          for (unsigned i = 0; i < nve; i++) {
            for (unsigned j = 0; j < dim; j++) {
              phiSupg[i] += ((SolVAR[j + dim] - meshVel[j]) * gradphi[i * dim + j]) * tauSupg;    // * (!solidmark[i]) ;
              phiSupg_old[i] += ((SolVAR_old[j + dim] - meshVel[j]) * gradphi_old[i * dim + j]) * tauSupg_old;    // * (!solidmark[i]);
            }
          }


          //BEGIN ALE + Momentum (Navier-Stokes)
          {
            for (unsigned i = 0; i < nve; i++) {

              //BEGIN redidual Laplacian ALE map in the reference domain
              adept::adouble LapmapVAR[3] = {0., 0., 0.};

              for (int idim = 0; idim < dim; idim++) {
                for (int jdim = 0; jdim < dim; jdim++) {
                  //LapmapVAR[idim] += ( GradSolVAR[idim][jdim] ) * gradphi_hat[i * dim + jdim];
                  LapmapVAR[idim] += (GradSolVAR[idim][jdim] + 0.*GradSolVAR[jdim][idim]) * gradphi[i * dim + jdim];
                }
              }

              for (int idim = 0; idim < dim; idim++) {
                aRhs[indexVAR[idim]][i] += (!solidmark[i]) * (- LapmapVAR[idim] * Weight_nojac);
              }

              //END residual Laplacian ALE map in reference domain

              if (flag_mat == 2) {
                //BEGIN residual Navier-Stokes in moving domain
                adept::adouble LapvelVAR[3] = {0., 0., 0.};
                adept::adouble LapvelVAR_old[3] = {0., 0., 0.};
                adept::adouble LapStrong[3] = {0., 0., 0.};
                adept::adouble LapStrong_old[3] = {0., 0., 0.};
                adept::adouble AdvaleVAR[3] = {0., 0., 0.};
                adept::adouble AdvaleVAR_old[3] = {0., 0., 0.};

                for (int idim = 0.; idim < dim; idim++) {
                  for (int jdim = 0.; jdim < dim; jdim++) {
                    unsigned kdim;
                    if (idim == jdim) kdim = jdim;
                    else if (1 == idim + jdim) kdim = dim;    // xy
                    else if (2 == idim + jdim) kdim = dim + 2;   // xz
                    else if (3 == idim + jdim) kdim = dim + 1;   // yz
                    //laplaciano debole
                    LapvelVAR[idim]     += (GradSolVAR[dim + idim][jdim] + GradSolVAR[dim + jdim][idim]) * gradphi[i * dim + jdim];
                    LapvelVAR_old[idim] += (GradSolVAR_old[dim + idim][jdim] + GradSolVAR_old[dim + jdim][idim]) * gradphi_old[i * dim + jdim];
                    //laplaciano strong
                    LapStrong[idim]     += (NablaSolVAR[dim + idim][jdim] + NablaSolVAR[dim + jdim][kdim]) * phiSupg[i];
                    LapStrong_old[idim] += (NablaSolVAR_old[dim + idim][jdim] + NablaSolVAR_old[dim + jdim][kdim]) * phiSupg_old[i];

                    AdvaleVAR[idim]	+= ((SolVAR[dim + jdim] - meshVel[jdim]) * GradSolVAR[dim + idim][jdim]
                                        + (0.*GradSolVAR[dim + jdim][jdim] - GradMeshVel[jdim][jdim]) * SolVAR[dim + idim]
                                        //) * phi[i];
                                       ) * (phi[i] + phiSupg[i]);

                    AdvaleVAR_old[idim]	+= ((SolVAR_old[dim + jdim] - meshVel[jdim]) * GradSolVAR_old[dim + idim][jdim]
                                            + (0.*GradSolVAR_old[dim + jdim][jdim] - GradMeshVel[jdim][jdim]) * SolVAR_old[dim + idim]
                                            //) * phi_old[i];
                                           ) * (phi_old[i] + phiSupg_old[i]);
                  }
                }

                for (int idim = 0; idim < dim; idim++) {

                  adept::adouble timeDerivative = - (SolVAR[dim + idim] * (phi[i] + phiSupg[i]) * Weight
                                                     - SolVAR_old[dim + idim] * (phi_old[i] + phiSupg_old[i]) * Weight_old) / dt;

                  adept::adouble value =  theta * (
                                            - AdvaleVAR[idim]      	             // advection term
                                            - IRe * LapvelVAR[idim]	             // viscous dissipation
                                            + IRe * LapStrong[idim]
                                            + 1. / rhof * SolVAR[2 * dim] * gradphi[i * dim + idim] // pressure gradient
                                          ) * Weight;                                // at time t

                  adept::adouble value_old = (1. - theta) * (
                                               - AdvaleVAR_old[idim]               	         // advection term
                                               - IRe * LapvelVAR_old[idim]	       	         // viscous dissipation
                                               + IRe * LapStrong_old[idim]
                                               + 1. / rhof * SolVAR[2 * dim] * gradphi_old[i * dim + idim]  // pressure gradient
                                             ) * Weight_old;			                 // at time t-dt

                  if (!solidmark[i]) {
                    aRhs[indexVAR[dim + idim]][i] += timeDerivative + value + value_old;
                  }
                  else {
                    aRhs[indexVAR[idim]][i] += timeDerivative + value + value_old;
                  }

                }
                //END redidual Navier-Stokes in moving domain
              }
              else if (flag_mat == 3) {
                //BEGIN redidual Porous Media in moving domain


                adept::adouble speed = 0.;
                adept::adouble speed_old = 0.;
                for (int idim = 0.; idim < dim; idim++) {
                  speed += SolVAR[dim + idim] * SolVAR[dim + idim];
                  speed_old += SolVAR_old[dim + idim] * SolVAR_old[dim + idim];
                }
                double eps = 1.0e-12;
                speed = sqrt(speed + eps);
                speed_old = sqrt(speed_old + eps);
                double DE = 0.;
                if (dim == 2) {
                  DE = 0.0002; // AAA_thrombus_2D
                  //DE = 0.00006; // turek2D
                }
                else if (dim == 3) {
                  DE = 0.000112; // porous3D
                }
                double b = 4188;
                double a = 1452;
                double K = DE * IRe * rhof / b; // alpha = mu/b * De
                double C2 = 2 * a / (rhof * DE);





// 		adept::adouble LapvelVAR[3] = {0., 0., 0.};
// 		adept::adouble LapvelVAR_old[3] = {0., 0., 0.};
// 		adept::adouble AdvaleVAR[3] = {0., 0., 0.};
// 		adept::adouble AdvaleVAR_old[3] = {0., 0., 0.};
//
// 		for(int idim = 0.; idim < dim; idim++) {
// 		  for(int jdim = 0.; jdim < dim; jdim++) {
//
// 		    LapvelVAR[idim]     += (GradSolVAR[dim + idim][jdim] + GradSolVAR[dim + jdim][idim]) * gradphi[i * dim + jdim];
// 		    LapvelVAR_old[idim] += (GradSolVAR_old[dim + idim][jdim] + GradSolVAR_old[dim + jdim][idim]) * gradphi_old[i * dim + jdim];
//
// 		    AdvaleVAR[idim]	+= ((SolVAR[dim + jdim] - meshVel[jdim]) * GradSolVAR[dim + idim][jdim]
// 					+ (GradSolVAR[dim + jdim][jdim] - GradMeshVel[jdim][jdim]) * SolVAR[dim + idim]
// 				      ) * phi[i];
//
// 		    AdvaleVAR_old[idim]	+= ((SolVAR_old[dim + jdim] - meshVel[jdim]) * GradSolVAR_old[dim + idim][jdim]
// 					    + (GradSolVAR_old[dim + jdim][jdim] - GradMeshVel[jdim][jdim]) * SolVAR_old[dim + idim]
// 					  ) * phi_old[i];
// 		  }
// 		}

                for (int idim = 0; idim < dim; idim++) {

// 		  adept::adouble timeDerivative = -(SolVAR[dim + idim] * phi[i] * Weight
//                                                   - SolVAR_old[dim + idim] * phi_old[i] * Weight_old);

                  adept::adouble value =  theta * (
// 					    - AdvaleVAR[idim]      	             // advection term
// 					    - IRe * LapvelVAR[idim]	             // viscous dissipation
                                            //- SolVAR[dim + idim] * ( IRe / K + 0.5 * C2 * speed ) * phi[i]
                                            - (SolVAR[dim + idim] - meshVel[idim]) * (IRe / K + 0.5 * C2 * speed) * (phi[i] + phiSupg[i])
                                            + 1. / rhof * SolVAR[2 * dim] * gradphi[i * dim + idim] // pressure gradient
                                          ) * Weight;                                // at time t

                  adept::adouble value_old = (1. - theta) * (
// 					      - AdvaleVAR_old[idim]               	         // advection term
// 					      - IRe * LapvelVAR_old[idim]	       	         // viscous dissipation
                                               //- SolVAR_old[dim + idim] * ( IRe / K + 0.5 * C2 * speed_old ) * phi_old[i]
                                               - (SolVAR_old[dim + idim] - meshVel[idim]) * (IRe / K + 0.5 * C2 * speed) * (phi_old[i] + phiSupg_old[i])
                                               + 1. / rhof * SolVAR[2 * dim] * gradphi_old[i * dim + idim]  // pressure gradient
                                             ) * Weight_old;			                 // at time t-dt

                  if (!solidmark[i]) {
                    aRhs[indexVAR[dim + idim]][i] += value + value_old;
                  }
                  else {
                    aRhs[indexVAR[idim]][i] += value + value_old;
                  }

                }
                //END redidual Porous Media in moving domain
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
              aRhs[indexVAR[2 * dim]][i] += - (-phi1[i] * div_vel) * Weight;
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




          if (solid_model == 0) {   // Saint-Venant
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
                Cauchy_old[i][j] = 2 * mus * e_old[i][j] - 2 * mus * I_e_old * SolVAR[2 * dim] * Id2th[i][j];
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

            if (solid_model <= 4) {   // Neo-Hookean
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
            else if (5 == solid_model) {   //Mooney-Rivlin
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
              invB[1][0] = - (B[0][1] * B[2][2] - B[0][2] * B[2][1]) * invdetB;
              invB[2][0] = (B[0][1] * B[1][2] - B[0][2] * B[1][1]) * invdetB;
              invB[0][1] = - (B[1][0] * B[2][2] - B[1][2] * B[2][0]) * invdetB;
              invB[1][1] = (B[0][0] * B[2][2] - B[0][2] * B[2][0]) * invdetB;
              invB[2][1] = - (B[0][0] * B[1][2] - B[1][0] * B[0][2]) * invdetB;
              invB[0][2] = (B[1][0] * B[2][1] - B[2][0] * B[1][1]) * invdetB;
              invB[1][2] = - (B[0][0] * B[2][1] - B[2][0] * B[0][1]) * invdetB;
              invB[2][2] = (B[0][0] * B[1][1] - B[1][0] * B[0][1]) * invdetB;

              invB_old[0][0] = (B_old[1][1] * B_old[2][2] - B_old[1][2] * B_old[2][1]) * invdetB_old;
              invB_old[1][0] = - (B_old[0][1] * B_old[2][2] - B_old[0][2] * B_old[2][1]) * invdetB_old;
              invB_old[2][0] = (B_old[0][1] * B_old[1][2] - B_old[0][2] * B_old[1][1]) * invdetB_old;
              invB_old[0][1] = - (B_old[1][0] * B_old[2][2] - B_old[1][2] * B_old[2][0]) * invdetB_old;
              invB_old[1][1] = (B_old[0][0] * B_old[2][2] - B_old[0][2] * B_old[2][0]) * invdetB_old;
              invB_old[2][1] = - (B_old[0][0] * B_old[1][2] - B_old[1][0] * B_old[0][2]) * invdetB_old;
              invB_old[0][2] = (B_old[1][0] * B_old[2][1] - B_old[2][0] * B_old[1][1]) * invdetB_old;
              invB_old[1][2] = - (B_old[0][0] * B_old[2][1] - B_old[2][0] * B_old[0][1]) * invdetB_old;
              invB_old[2][2] = (B_old[0][0] * B_old[1][1] - B_old[1][0] * B_old[0][1]) * invdetB_old;


              I1_B = B[0][0] + B[1][1] + B[2][2];
              I2_B = B[0][0] * B[1][1] + B[1][1] * B[2][2] + B[2][2] * B[0][0]
                     - B[0][1] * B[1][0] - B[1][2] * B[2][1] - B[2][0] * B[0][2];

              I1_B_old = B_old[0][0] + B_old[1][1] + B_old[2][2];
              I2_B_old = B_old[0][0] * B_old[1][1] + B_old[1][1] * B_old[2][2] + B_old[2][2] * B_old[0][0]
                         - B_old[0][1] * B_old[1][0] - B_old[1][2] * B_old[2][1] - B_old[2][0] * B_old[0][2];

              //double C1 = mus / 3.;
              double C1 = (elementGroup == 15) ? mus1 / 3. : mus / 3.;
              double C2 = C1 / 2.;

              for (int I = 0; I < 3; ++I) {
                for (int J = 0; J < 3; ++J) {
                  Cauchy[I][J] =  2.* (C1 * B[I][J] - C2 * invB[I][J])
                                  //- (2. / 3.) * (C1 * I1_B - C2 * I2_B) * SolVAR[2 * dim] * Id2th[I][J];
                                  - 1. / rhof * SolVAR[2 * dim] * Id2th[I][J];

                  Cauchy_old[I][J] =  2.* (C1 * B_old[I][J] - C2 * invB_old[I][J])
                                      //- (2. / 3.) * (C1 * I1_B_old - C2 * I2_B_old) * SolVAR[2 * dim] * Id2th[I][J];
                                      - 1. / rhof * SolVAR[2 * dim] * Id2th[I][J];

                }
              }

            }
          }

          //END build Cauchy Stress in moving domain

          //BEGIN v=0 + Momentum (Solid)
          {
            for (unsigned i = 0; i < nve; i++) {

              //BEGIN redidual d_t - v = 0 in fixed domain
              for (int idim = 0; idim < dim; idim++) {
                aRhs[indexVAR[dim + idim]][i] +=  - phi[i] * ((-SolVAR[idim] + SolVAR_old[idim]) / dt +
                                                  (theta * SolVAR[dim + idim] + (1. - theta) * SolVAR_old[dim + idim])
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

                adept::adouble timeDerivative = - (rhos * SolVAR[dim + idim] * phi[i] * Weight
                                                   - rhos * SolVAR_old[dim + idim] * phi_old[i] * Weight_old) / dt;

                adept::adouble value =  theta * (rhos * phi[i] * _gravity[idim]      // body force
                                                 - CauchyDIR[idim]			  // stress
                                                ) * Weight;                         // at time t

                adept::adouble value_old = (1. - theta) * (rhos * phi_old[i] * _gravity[idim]      // body force
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
                  aRhs[indexVAR[2 * dim]][i] += - (-phi1[i] * (I_e + (!incompressible) / lambda * SolVAR[2 * dim])) * Weight_hat;
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

      if (assembleMatrix) {
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




//       //BEGIN local to global assembly
//       //copy adouble aRhs into double Rhs
//       for ( unsigned i = 0; i < 2 * dim; i++ ) {
//         for ( int j = 0; j < nve; j++ ) {
//           Rhs[indexVAR[i]][j] = aRhs[indexVAR[i]][j].value();
//         }
//       }
//
//       for ( unsigned j = 0; j < nve1; j++ ) {
//         Rhs[indexVAR[2 * dim]][j] = aRhs[indexVAR[2 * dim]][j].value();
//       }
//
//       for ( int i = 0; i < 2 * dim + 1; i++ ) {
//         myRES->add_vector_blocked ( Rhs[indexVAR[i]], dofsVAR[i] );
//       }
//
//       if ( assembleMatrix ) {
//         //Store equations
//         for ( int i = 0; i < 2 * dim; i++ ) {
//           s.dependent ( &aRhs[indexVAR[i]][0], nve );
//           s.independent ( &Soli[indexVAR[i]][0], nve );
//         }
//
//         s.dependent ( &aRhs[indexVAR[2 * dim]][0], nve1 );
//         s.independent ( &Soli[indexVAR[2 * dim]][0], nve1 );
//         s.jacobian ( &Jac[0] );
//         unsigned nveAll = ( 2 * dim * nve + nve1 );
//
//         for ( int inode = 0; inode < nveAll; inode++ ) {
//           for ( int jnode = 0; jnode < nveAll; jnode++ ) {
//             KKloc[inode * nveAll + jnode] = -Jac[jnode * nveAll + inode];
//           }
//         }
//
//
//         myKK->add_matrix_blocked ( KKloc, dofsAll, dofsAll );
//         s.clear_independents();
//         s.clear_dependents();
//
//
//       }
//       //END local to global assembly
//    } //end list of elements loop

    if (assembleMatrix) myKK->close();

    myRES->close();

    delete area_elem_first;

    // *************************************
    end_time = clock();
    AssemblyTime += (end_time - start_time);
    // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************
  }






  void FSITimeDependentAssemblySupgNew(MultiLevelProblem & ml_prob)
  {

    clock_t AssemblyTime = 0;
    clock_t start_time, end_time;

    bool auxDisp = true;
    unsigned nBlocks = (auxDisp) ? 3 : 2;

    //meshIsCurrupted = true;
    //pointers and references

    //MonolithicFSINonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<MonolithicFSINonLinearImplicitSystem>("Fluid-Structure-Interaction");
    TransientNonlinearImplicitSystem & my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("Fluid-Structure-Interaction");
    const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();

    MultiLevelSolution	* ml_sol		= ml_prob._ml_sol;

    //FSIConstrainLeaflet(*ml_sol);

    Solution	* mysolution		= ml_sol->GetSolutionLevel(level);

    LinearEquationSolver * myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];
    Mesh	*	mymsh		=  ml_prob._ml_msh->GetLevel(level);
    elem	*	myel		=  mymsh->el;
    SparseMatrix	* myKK		=  myLinEqSolver->_KK;
    NumericVector	* myRES		=  myLinEqSolver->_RES;


    bool assembleMatrix = my_nnlin_impl_sys.GetAssembleMatrix();
    // call the adept stack object
    adept::Stack& s = FemusInit::_adeptStack;
    if (assembleMatrix) s.continue_recording();
    else s.pause_recording();

    const unsigned dim = mymsh->GetDimension();
    const unsigned nabla_dim = 3 * (dim - 1);
    const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

    double theta = 0.5;
    //double theta = 1.;

    // local objects
    vector<adept::adouble> SolVAR(nBlocks * dim + 1);
    vector<double> SolVAR_old(nBlocks * dim);

    vector<vector < adept::adouble > > GradSolVAR(nBlocks * dim);
    vector<vector < double > > GradSolVAR_old(nBlocks * dim);

    vector<vector < adept::adouble > > GradSolhatVAR(nBlocks * dim);
    vector<vector < double > > GradSolhatVAR_old(nBlocks * dim);

    vector<vector<adept::adouble> > NablaSolVAR(nBlocks * dim);
    vector<vector < double > > NablaSolVAR_old(nBlocks * dim);

    for (int i = 0; i < nBlocks * dim; i++) {
      GradSolVAR[i].resize(dim);
      GradSolVAR_old[i].resize(dim);

      GradSolhatVAR[i].resize(dim);
      GradSolhatVAR_old[i].resize(dim);

      NablaSolVAR[i].resize(nabla_dim);
      NablaSolVAR_old[i].resize(nabla_dim);
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

    const double * phi1;

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

    vector< vector< adept::adouble > > Soli(nBlocks * dim + 1);
    vector< vector< double > > Soli_old(nBlocks * dim + 1);
    vector< vector< int > > dofsVAR(nBlocks * dim + 1);

    for (int i = 0; i < nBlocks * dim + 1; i++) {
      Soli[i].reserve(max_size);
      Soli_old[i].reserve(max_size);
      dofsVAR[i].reserve(max_size);
    }

    vector< vector< double > > Rhs(nBlocks * dim + 1);
    vector< vector< adept::adouble > > aRhs(nBlocks * dim + 1);

    for (int i = 0; i < nBlocks * dim + 1; i++) {
      aRhs[i].reserve(max_size);
      Rhs[i].reserve(max_size);
    }


    vector < int > dofsAll;
    dofsAll.reserve(max_size * (nBlocks * dim + 1));

    vector < double > Jac;
    Jac.reserve(dim * max_size * (nBlocks * dim + 1) *dim * max_size * (nBlocks * dim + 1));

    // ------------------------------------------------------------------------
    // Physical parameters
    double rhof	 	= ml_prob.parameters.get<Fluid> ("Fluid").get_density();
    double rhos	 	= ml_prob.parameters.get<Fluid> ("Solid").get_density() / rhof;
    double mu_lame 	= ml_prob.parameters.get<Solid> ("Solid").get_lame_shear_modulus();
    double lambda_lame 	= ml_prob.parameters.get<Solid> ("Solid").get_lame_lambda();
    double mus		= mu_lame / rhof;
    double mu_lame1 	= ml_prob.parameters.get < Solid> ("Solid1").get_lame_shear_modulus();
    double mus1 	= mu_lame1 / rhof;
    double IRe 		= ml_prob.parameters.get<Fluid> ("Fluid").get_IReynolds_number();
    double lambda	= lambda_lame / rhof;
    double betans	= 1.;
    int    solid_model	= ml_prob.parameters.get<Solid> ("Solid").get_physical_model();

    if (solid_model >= 2 && solid_model <= 4) {
      std::cout << "Error! Solid Model " << solid_model << "not implemented\n";
      abort();
    }

    bool incompressible = (0.5 == ml_prob.parameters.get<Solid> ("Solid").get_poisson_coeff()) ? 1 : 0;
    const bool penalty = ml_prob.parameters.get<Solid> ("Solid").get_if_penalty();

    // gravity
    double _gravity[3] = {0., 0., 0.};

    double dt =  my_nnlin_impl_sys.GetIntervalTime();
    double time =  my_nnlin_impl_sys.GetTime();
    // -----------------------------------------------------------------
    // space discretization parameters
    unsigned SolType2 = ml_sol->GetSolutionType(ml_sol->GetIndex("U"));
    unsigned SolType1 = ml_sol->GetSolutionType(ml_sol->GetIndex("PS"));

    // mesh and procs
    unsigned nel    = mymsh->GetNumberOfElements();
    unsigned igrid  = mymsh->GetLevel();
    unsigned iproc  = mymsh->processor_id();

    unsigned indLmbd = ml_sol->GetIndex("lmbd");

    //----------------------------------------------------------------------------------
    //variable-name handling
    const char varname[10][4] = {"DX", "DY", "DZ", "U", "V", "W", "DX1", "DY1", "DZ1", "PS"};

    vector <unsigned> indexVAR(nBlocks * dim + 1);
    vector <unsigned> indVAR(nBlocks * dim + 1);
    vector <unsigned> SolType(nBlocks * dim + 1);

    for (unsigned ivar = 0; ivar < dim; ivar++) {
      for (unsigned k = 0; k < nBlocks; k++) {
        indVAR[ivar + k * dim] = ml_sol->GetIndex(&varname[ivar + k * 3][0]);
        SolType[ivar + k * dim] = ml_sol->GetSolutionType(&varname[ivar + k * 3][0]);
        indexVAR[ivar + k * dim] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar + k * 3][0]);
      }
    }

    indexVAR[nBlocks * dim] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[9][0]);
    indVAR[nBlocks * dim] = ml_sol->GetIndex(&varname[9][0]);
    SolType[nBlocks * dim] = ml_sol->GetSolutionType(&varname[9][0]);
    //----------------------------------------------------------------------------------

    int nprocs = mymsh->n_processors();
    NumericVector * area_elem_first;
    area_elem_first = NumericVector::build().release();

    if (nprocs == 1) {
      area_elem_first->init(nprocs, 1, false, SERIAL);
    }
    else {
      area_elem_first->init(nprocs, 1, false, PARALLEL);
    }


    double rapresentative_area = 1.;

    start_time = clock();

    if (assembleMatrix) myKK->zero();

    NumericVector* setIfCorrupted;
    setIfCorrupted = NumericVector::build().release();
    setIfCorrupted->init(mymsh->n_processors(), 1 , false, AUTOMATIC);


  begin:

    area_elem_first->zero();
    setIfCorrupted->zero();

    // *** element loop ***
    for (int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++) {

      short unsigned ielt = mymsh->GetElementType(iel);
      unsigned nve        = mymsh->GetElementDofNumber(iel, SolType2);
      unsigned nve1       = mymsh->GetElementDofNumber(iel, SolType1);
      int flag_mat        = mymsh->GetElementMaterial(iel);
      unsigned elementGroup = mymsh->GetElementGroup(iel);

      // *******************************************************************************************************

      //initialization of everything is in common fluid and solid

      //Rhs
      for (int i = 0; i < nBlocks * dim; i++) {
        dofsVAR[i].resize(nve);
        Soli[indexVAR[i]].resize(nve);
        Soli_old[indexVAR[i]].resize(nve);
        aRhs[indexVAR[i]].resize(nve);
      }

      dofsVAR[nBlocks * dim].resize(nve1);
      Soli[indexVAR[nBlocks * dim]].resize(nve1);
      Soli_old[indexVAR[nBlocks * dim]].resize(nve1);
      aRhs[indexVAR[nBlocks * dim]].resize(nve1);

      dofsAll.resize(0);

      Jac.resize((nBlocks * dim * nve + nve1) * (nBlocks * dim * nve + nve1));

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
        solidmark[i] = mymsh->GetSolidMark(idof);    // to check

        for (int j = 0; j < dim; j++) {
          for (unsigned k = 0; k < nBlocks; k++) {
            Soli[indexVAR[j + k * dim]][i] = (*mysolution->_Sol[indVAR[j + k * dim]])(idof);
            Soli_old[indexVAR[j + k * dim]][i] = (*mysolution->_SolOld[indVAR[j + k * dim]])(idof);
            aRhs[indexVAR[j + k * dim]][i] = 0.;
            dofsVAR[j + k * dim][i] = myLinEqSolver->GetSystemDof(indVAR[j + k * dim], indexVAR[j + k * dim], i, iel);
          }
          vx_hat[j][i] = (*mymsh->_topology->_Sol[j])(idof);
        }
      }

      // pressure dofs
      for (unsigned i = 0; i < nve1; i++) {
        unsigned idof = mymsh->GetSolutionDof(i, iel, SolType[nBlocks * dim]);
        dofsVAR[nBlocks * dim][i] = myLinEqSolver->GetSystemDof(indVAR[nBlocks * dim], indexVAR[nBlocks * dim], i, iel);
        Soli[indexVAR[nBlocks * dim]][i]     = (*mysolution->_Sol[indVAR[nBlocks * dim]])(idof);
        Soli_old[indexVAR[nBlocks * dim]][i] = (*mysolution->_SolOld[indVAR[nBlocks * dim]])(idof);
        aRhs[indexVAR[nBlocks * dim]][i] = 0.;
      }

      // build dof ccomposition
      for (int idim = 0; idim < nBlocks * dim; idim++) {
        dofsAll.insert(dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end());
      }
      dofsAll.insert(dofsAll.end(), dofsVAR[nBlocks * dim].begin(), dofsVAR[nBlocks * dim].end());

      if (assembleMatrix) s.new_recording();


      for (int j = 0; j < nve; j++) {
        unsigned idof = mymsh->GetSolutionDof(j, iel, SolType2);
        for (unsigned idim = 0; idim < dim; idim++) {
          vx[idim][j] = vx_hat[idim][j] +  Soli[indexVAR[idim + meshIsCurrupted * 2 * dim]][j];
          vx_old[idim][j] = vx_hat[idim][j] + Soli_old[indexVAR[idim + meshIsCurrupted * 2 * dim]][j];
        }
      }

      //Boundary integral
      {

        vector < adept::adouble> normal(dim, 0);
        vector < double > normal_old(dim, 0);

        // loop on faces
        for (unsigned jface = 0; jface < mymsh->GetElementFaceNumber(iel); jface++) {
          std::vector< double > xx(dim, 0.);

          // look for boundary faces
          if (myel->GetFaceElementIndex(iel, jface) < 0) {

            unsigned int face = - (mymsh->el->GetFaceElementIndex(iel, jface) + 1);
            double tau = 0.;
            double tau_old = 0.;
            if ((!ml_sol->GetBdcFunction()(xx, "PS", tau, face, time) &&
                 !ml_sol->GetBdcFunction()(xx, "PS", tau_old, face, time - dt))
                && (tau != 0. || tau_old != 0.)) {

              unsigned nve = mymsh->GetElementFaceDofNumber(iel, jface, SolType2);
              const unsigned felt = mymsh->GetElementFaceType(iel, jface);

              for (unsigned i = 0; i < nve; i++) {
                unsigned int ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);
                unsigned idof = mymsh->GetSolutionDof(ilocal, iel, 2);

                for (unsigned idim = 0; idim < dim; idim++) {
                  vx_face[idim][i]    = (*mymsh->_topology->_Sol[idim])(idof) + Soli[indexVAR[idim + meshIsCurrupted * 2 * dim]][ilocal]; //TODO
                  vx_face_old[idim][i] = (*mymsh->_topology->_Sol[idim])(idof) + Soli_old[indexVAR[idim + meshIsCurrupted * 2 * dim]][ilocal];
                }
              }

              for (unsigned igs = 0; igs < mymsh->_finiteElement[felt][SolType2]->GetGaussPointNumber(); igs++) {
                mymsh->_finiteElement[felt][SolType2]->JacobianSur(vx_face, igs, Weight, phi, gradphi, normal);
                mymsh->_finiteElement[felt][SolType2]->JacobianSur(vx_face_old, igs, Weight_old, phi_old, gradphi_old, normal_old);

                //phi1 =mymsh->_finiteElement[felt][SolType2]->GetPhi(igs);
                // *** phi_i loop ***
                for (unsigned i = 0; i < nve; i++) {
                  adept::adouble value = - phi[i] * tau / rhof * Weight;
                  double value_old = - phi_old[i] * tau_old / rhof * Weight_old;
                  unsigned int ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);

                  for (unsigned idim = 0; idim < dim; idim++) {
                    if ((!solidmark[ilocal])) {
                      aRhs[indexVAR[dim + idim]][ilocal]   += (theta * value + (1. - theta) * value_old) * normal[idim];
                    }
                    else {   //if interface node it goes to solid
                      aRhs[indexVAR[idim]][ilocal]   += (theta * value + (1. - theta) * value_old) * normal[idim];
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

        if ( Weight.value() < 0. ) setIfCorrupted->set(iproc, 1.);

        if (flag_mat == 2 || flag_mat == 3  || iel == mymsh->_elementOffset[iproc]) {
          if (ig == 0) {
            double GaussWeight = mymsh->_finiteElement[ielt][SolType2]->GetGaussWeight(ig);
            area = Weight_hat / GaussWeight;

            if (iel == mymsh->_elementOffset[iproc]) {
              area_elem_first->add(mymsh->processor_id(), area);
              area_elem_first->close();
              rapresentative_area = area_elem_first->l1_norm() / nprocs;
            }
          }

          Weight_nojac = Weight;// / area * 1.e-10;// * rapresentative_area;

          //-----------------------------------------------------------------------//
          // for vein_valve mesh
          std::vector<double> xc(dim, 0.);
          //xc[0] = -0.001013; // vein_valve
          //xc[1] = 0.07000;
//           xc[0] = -0.000286; // vein_valve_closed
//           xc[1] = 0.07000;
//
          if (dim == 2) {
            xc[0] = -0.00025; // vein_valve_closed
            xc[1] = 0.07000;
          }
// 	  else if (dim == 3){
// 	    xc[0] = -0.008; // vein_valve_closed
// 	    xc[1] = 0.0;
// 	  }
          else if (dim == 3) {
            xc[0] = 0.0015; // vein_valve_closed
            xc[1] = 0.0;
            xc[2] = 0.01;
          }

          double distance = 0.;
          for (unsigned k = 0; k < dim; k++) {
            distance += (vx_hat[k][ nve - 1] - xc[k]) * (vx_hat[k][nve - 1] - xc[k]);
          }
          distance = sqrt(distance);
          Weight_nojac *= 1. / (1 + 10000 * distance);

          if (elementGroup == 16) Weight_nojac *= 100;
          //-----------------------------------------------------------------------//
        }

        // ---------------------------------------------------------------------------
        // displacement and velocity
        for (int i = 0; i < nBlocks * dim; i++) {
          SolVAR[i] = 0.;
          SolVAR_old[i] = 0.;

          for (int j = 0; j < dim; j++) {
            GradSolVAR[i][j] = 0.;
            GradSolVAR_old[i][j] = 0.;

            GradSolhatVAR[i][j] = 0.;
            GradSolhatVAR_old[i][j] = 0.;
          }

          for (int j = 0; j < nabla_dim; j++) {
            NablaSolVAR[i][j] = 0.;
            NablaSolVAR_old[i][j] = 0.;
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
            for (int j = 0; j < nabla_dim; j++) {
              NablaSolVAR[i][j]     += nablaphi[inode * nabla_dim + j] * Soli[indexVAR[i]][inode];
              NablaSolVAR_old[i][j] += nablaphi_old[inode * nabla_dim + j] * Soli_old[indexVAR[i]][inode];
            }
          }
        }

        // pressure
        SolVAR[nBlocks * dim] = 0.;

        for (unsigned inode = 0; inode < nve1; inode++) {
          SolVAR[nBlocks * dim]    += phi1[inode] * Soli[indexVAR[nBlocks * dim]][inode];
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
        //BEGIN FLUID and POROUS MEDIA ASSEMBLY ============
        if (flag_mat == 2 || flag_mat == 3) {

          vector < adept::adouble > a(dim);
          vector < adept::adouble > a_old(dim);
          for (int i = 0; i < dim; i++) {
            a[i] = SolVAR[i + dim] - 0.*meshVel[i]; // TODO maybe we subtract meshVel[i] maybe not
            a_old[i] = SolVAR_old[i + dim] - 0.*meshVel[i];
          }

          // speed
          adept::adouble aL2Norm = 0.;
          adept::adouble a_oldL2Norm = 0.;
          for (int i = 0; i < dim; i++) {
            aL2Norm += a[i] * a[i];
            a_oldL2Norm += a_old[i] * a_old[i];
          }
          aL2Norm = sqrt(aL2Norm);
          a_oldL2Norm = sqrt(a_oldL2Norm);

          double sqrtlambdak = (*mysolution->_Sol[indLmbd])(iel);
          adept::adouble tauSupg = 1. / (sqrtlambdak * sqrtlambdak * 4.*IRe);
          adept::adouble tauSupg_old = 1. / (sqrtlambdak * sqrtlambdak * 4.*IRe);
          adept::adouble Rek = aL2Norm / (4.*sqrtlambdak * IRe);
          adept::adouble Rek_old = a_oldL2Norm / (4.*sqrtlambdak * IRe);

          if (Rek > 1.0e-15) {
            adept::adouble xiRek = (Rek >= 1.) ? 1. : Rek;
            tauSupg   = xiRek / (aL2Norm * sqrtlambdak);
          }

          if (Rek_old > 1.0e-15) {
            adept::adouble xiRek_old = (Rek_old >= 1.) ? 1. : Rek_old;
            tauSupg_old   = xiRek_old / (a_oldL2Norm * sqrtlambdak);
          }

          //std::cout << tauSupg.value() <<"  " <<tauSupg_old.value()<<std::endl;


          vector < adept::adouble > phiSupg(nve, 0.);
          vector < adept::adouble > phiSupg_old(nve, 0.);
          //tauSupg=0;
          //tauSupg_old=0;
          for (unsigned i = 0; i < nve; i++) {
            for (unsigned j = 0; j < dim; j++) {
              phiSupg[i] += ((SolVAR[j + dim] - meshVel[j]) * gradphi[i * dim + j]) * tauSupg;    // * (!solidmark[i]) ;
              phiSupg_old[i] += ((SolVAR_old[j + dim] - meshVel[j]) * gradphi_old[i * dim + j]) * tauSupg_old;    // * (!solidmark[i]);
            }
          }


          //BEGIN ALE + Momentum (Navier-Stokes)
          {
            for (unsigned i = 0; i < nve; i++) {

              //BEGIN redidual Laplacian ALE map in the reference domain


              if (!solidmark[i]) {
                adept::adouble LapmapVAR[3] = {0., 0., 0.};
                for (int idim = 0; idim < dim; idim++) {
                  for (int jdim = 0; jdim < dim; jdim++) {
                    //LapmapVAR[idim] += ( GradSolVAR[idim][jdim] ) * gradphi_hat[i * dim + jdim];
                    LapmapVAR[idim] += (GradSolVAR[idim][jdim] + 0.*GradSolVAR[jdim][idim]) * gradphi[i * dim + jdim];
                  }
                }

                for (int idim = 0; idim < dim; idim++) {
                  aRhs[indexVAR[idim]][i] +=  (- LapmapVAR[idim] * Weight_nojac);
                }
              }

              //END residual Laplacian ALE map in reference domain

              if (flag_mat == 2) {
                //BEGIN residual Navier-Stokes in moving domain
                adept::adouble LapvelVAR[3] = {0., 0., 0.};
                adept::adouble LapvelVAR_old[3] = {0., 0., 0.};
                adept::adouble LapStrong[3] = {0., 0., 0.};
                adept::adouble LapStrong_old[3] = {0., 0., 0.};
                adept::adouble AdvaleVAR[3] = {0., 0., 0.};
                adept::adouble AdvaleVAR_old[3] = {0., 0., 0.};

                for (int idim = 0.; idim < dim; idim++) {
                  for (int jdim = 0.; jdim < dim; jdim++) {
                    unsigned kdim;
                    if (idim == jdim) kdim = jdim;
                    else if (1 == idim + jdim) kdim = dim;    // xy
                    else if (2 == idim + jdim) kdim = dim + 2;   // xz
                    else if (3 == idim + jdim) kdim = dim + 1;   // yz
                    //laplaciano debole
                    LapvelVAR[idim]     += (GradSolVAR[dim + idim][jdim] + GradSolVAR[dim + jdim][idim]) * gradphi[i * dim + jdim];
                    LapvelVAR_old[idim] += (GradSolVAR_old[dim + idim][jdim] + GradSolVAR_old[dim + jdim][idim]) * gradphi_old[i * dim + jdim];
                    //laplaciano strong
                    LapStrong[idim]     += (NablaSolVAR[dim + idim][jdim] + NablaSolVAR[dim + jdim][kdim]) * phiSupg[i];
                    LapStrong_old[idim] += (NablaSolVAR_old[dim + idim][jdim] + NablaSolVAR_old[dim + jdim][kdim]) * phiSupg_old[i];

                    AdvaleVAR[idim]	+= ((SolVAR[dim + jdim] - meshVel[jdim]) * GradSolVAR[dim + idim][jdim]
                                        + (0.*GradSolVAR[dim + jdim][jdim] - GradMeshVel[jdim][jdim]) * SolVAR[dim + idim]
                                        //) * phi[i];
                                       ) * (phi[i] + phiSupg[i]);

                    AdvaleVAR_old[idim]	+= ((SolVAR_old[dim + jdim] - meshVel[jdim]) * GradSolVAR_old[dim + idim][jdim]
                                            + (0.*GradSolVAR_old[dim + jdim][jdim] - GradMeshVel[jdim][jdim]) * SolVAR_old[dim + idim]
                                            //) * phi_old[i];
                                           ) * (phi_old[i] + phiSupg_old[i]);
                  }
                }

                for (int idim = 0; idim < dim; idim++) {

                  adept::adouble timeDerivative = - (SolVAR[dim + idim] * (phi[i] + phiSupg[i]) * Weight
                                                     - SolVAR_old[dim + idim] * (phi_old[i] + phiSupg_old[i]) * Weight_old) / dt;

                  adept::adouble value =  theta * (
                                            - AdvaleVAR[idim]      	             // advection term
                                            - IRe * LapvelVAR[idim]	             // viscous dissipation
                                            + IRe * LapStrong[idim]
                                            + 1. / rhof * SolVAR[nBlocks * dim] * gradphi[i * dim + idim] // pressure gradient
                                          ) * Weight;                                // at time t

                  adept::adouble value_old = (1. - theta) * (
                                               - AdvaleVAR_old[idim]               	         // advection term
                                               - IRe * LapvelVAR_old[idim]	       	         // viscous dissipation
                                               + IRe * LapStrong_old[idim]
                                               + 1. / rhof * SolVAR[nBlocks * dim] * gradphi_old[i * dim + idim]  // pressure gradient
                                             ) * Weight_old;			                 // at time t-dt

                  if (!solidmark[i]) {
                    aRhs[indexVAR[dim + idim]][i] += timeDerivative + value + value_old;
                  }
                  else {
                    aRhs[indexVAR[idim]][i] += timeDerivative + value + value_old;
                  }

                }
                //END redidual Navier-Stokes in moving domain
              }
              else if (flag_mat == 3) {
                //BEGIN redidual Porous Media in moving domain


                adept::adouble speed = 0.;
                adept::adouble speed_old = 0.;
                for (int idim = 0.; idim < dim; idim++) {
                  speed += SolVAR[dim + idim] * SolVAR[dim + idim];
                  speed_old += SolVAR_old[dim + idim] * SolVAR_old[dim + idim];
                }
                double eps = 1.0e-12;
                speed = sqrt(speed + eps);
                speed_old = sqrt(speed_old + eps);
                double DE = 0.;
                if (dim == 2) {
                  DE = 0.0002; // AAA_thrombus_2D
                  //DE = 0.00006; // turek2D
                }
                else if (dim == 3) {
                  DE = 0.000112; // porous3D
                }
                double b = 4188;
                double a = 1452;
                double K = DE * IRe * rhof / b; // alpha = mu/b * De
                double C2 = 2 * a / (rhof * DE);





// 		adept::adouble LapvelVAR[3] = {0., 0., 0.};
// 		adept::adouble LapvelVAR_old[3] = {0., 0., 0.};
// 		adept::adouble AdvaleVAR[3] = {0., 0., 0.};
// 		adept::adouble AdvaleVAR_old[3] = {0., 0., 0.};
//
// 		for(int idim = 0.; idim < dim; idim++) {
// 		  for(int jdim = 0.; jdim < dim; jdim++) {
//
// 		    LapvelVAR[idim]     += (GradSolVAR[dim + idim][jdim] + GradSolVAR[dim + jdim][idim]) * gradphi[i * dim + jdim];
// 		    LapvelVAR_old[idim] += (GradSolVAR_old[dim + idim][jdim] + GradSolVAR_old[dim + jdim][idim]) * gradphi_old[i * dim + jdim];
//
// 		    AdvaleVAR[idim]	+= ((SolVAR[dim + jdim] - meshVel[jdim]) * GradSolVAR[dim + idim][jdim]
// 					+ (GradSolVAR[dim + jdim][jdim] - GradMeshVel[jdim][jdim]) * SolVAR[dim + idim]
// 				      ) * phi[i];
//
// 		    AdvaleVAR_old[idim]	+= ((SolVAR_old[dim + jdim] - meshVel[jdim]) * GradSolVAR_old[dim + idim][jdim]
// 					    + (GradSolVAR_old[dim + jdim][jdim] - GradMeshVel[jdim][jdim]) * SolVAR_old[dim + idim]
// 					  ) * phi_old[i];
// 		  }
// 		}

                for (int idim = 0; idim < dim; idim++) {

// 		  adept::adouble timeDerivative = -(SolVAR[dim + idim] * phi[i] * Weight
//                                                   - SolVAR_old[dim + idim] * phi_old[i] * Weight_old);

                  adept::adouble value =  theta * (
// 					    - AdvaleVAR[idim]      	             // advection term
// 					    - IRe * LapvelVAR[idim]	             // viscous dissipation
                                            //- SolVAR[dim + idim] * ( IRe / K + 0.5 * C2 * speed ) * phi[i]
                                            - (SolVAR[dim + idim] - meshVel[idim]) * (IRe / K + 0.5 * C2 * speed) * (phi[i] + phiSupg[i])
                                            + 1. / rhof * SolVAR[nBlocks * dim] * gradphi[i * dim + idim] // pressure gradient
                                          ) * Weight;                                // at time t

                  adept::adouble value_old = (1. - theta) * (
// 					      - AdvaleVAR_old[idim]               	         // advection term
// 					      - IRe * LapvelVAR_old[idim]	       	         // viscous dissipation
                                               //- SolVAR_old[dim + idim] * ( IRe / K + 0.5 * C2 * speed_old ) * phi_old[i]
                                               - (SolVAR_old[dim + idim] - meshVel[idim]) * (IRe / K + 0.5 * C2 * speed) * (phi_old[i] + phiSupg_old[i])
                                               + 1. / rhof * SolVAR[nBlocks * dim] * gradphi_old[i * dim + idim]  // pressure gradient
                                             ) * Weight_old;			                 // at time t-dt

                  if (!solidmark[i]) {
                    aRhs[indexVAR[dim + idim]][i] += value + value_old;
                  }
                  else {
                    aRhs[indexVAR[idim]][i] += value + value_old;
                  }

                }
                //END redidual Porous Media in moving domain
              }

              if (auxDisp) {
                //BEGIN auxiliary displacement equation
                adept::adouble LapAuxVAR[3] = {0., 0., 0.};

                for (int idim = 0; idim < dim; idim++) {
                  for (int jdim = 0; jdim < dim; jdim++) {
                    LapAuxVAR[idim] += ( GradSolVAR[idim + 2 * dim][jdim] +  0. * GradSolVAR[jdim + 2 * dim][idim] ) * gradphi[i * dim + jdim];
                  }
                }

                for (int idim = 0; idim < dim; idim++) {
                  aRhs[indexVAR[idim + 2 * dim]][i] += (  meshIsCurrupted * 1.0e-1 * LapAuxVAR[idim] +  0 * ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi[i] ) * Weight_nojac;
                }
//                 for (int idim = 1; idim < dim; idim++) {
//                   aRhs[indexVAR[idim + 2 * dim]][i] += (  meshIsCurrupted * 1.0e-1 * LapAuxVAR[idim] +  0 * ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi[i] ) * Weight_nojac;
//                 }
//                 for (int idim = 1; idim < dim; idim++) {
//                   aRhs[indexVAR[idim + 2 * dim]][i] += (  meshIsCurrupted * 1.0e-1 * LapAuxVAR[idim] +  0 * ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi[i] ) * Weight_nojac;
//                 }

//                 adept::adouble LapAuxVAR[3] = {0., 0., 0.};
//
//                 for (int idim = 0; idim < dim; idim++) {
//                   for (int jdim = 0; jdim < dim; jdim++) {
//                      LapAuxVAR[idim] += GradSolhatVAR[idim + 2 * dim][jdim] * gradphi_hat[i * dim + jdim];
//                   }
//                 }
//
//                 for (int idim = 0; idim < dim; idim++) {
//                   aRhs[indexVAR[idim + 2 * dim]][i] += (  1.0e-10 * LapAuxVAR[idim] +  ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi_hat[i] ) * Weight_hat;
//                 }


                //END auxiliary displacement equation
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
              aRhs[indexVAR[nBlocks * dim]][i] += - (-phi1[i] * div_vel) * Weight;
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
          adept::adouble J_hat1;
          double J_hat_old;
          adept::adouble I_e;
          double I_e_old;
          adept::adouble Cauchy[3][3];
          adept::adouble Cauchy_old[3][3];
          double Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};




          if (solid_model == 0) {   // Saint-Venant
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
                Cauchy[i][j]     = 2 * mus * e[i][j] - 2 * mus * I_e * SolVAR[nBlocks * dim] * Id2th[i][j];
                Cauchy_old[i][j] = 2 * mus * e_old[i][j] - 2 * mus * I_e_old * SolVAR[nBlocks * dim] * Id2th[i][j];
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
            adept::adouble F1[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
            adept::adouble B[3][3];

            double F_old[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
            double B_old[3][3];

            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                F[i][j] += GradSolhatVAR[i][j];
                F1[i][j] += GradSolhatVAR[i + 2 * dim][j];
                F_old[i][j] += GradSolhatVAR_old[i][j];
              }
            }

            J_hat =   F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                      - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

            J_hat1 =   F1[0][0] * F1[1][1] * F1[2][2] + F1[0][1] * F1[1][2] * F1[2][0] + F1[0][2] * F1[1][0] * F1[2][1]
                       - F1[2][0] * F1[1][1] * F1[0][2] - F1[2][1] * F1[1][2] * F1[0][0] - F1[2][2] * F1[1][0] * F1[0][1];

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

            if (solid_model <= 4) {   // Neo-Hookean
              I1_B = B[0][0] + B[1][1] + B[2][2];
              I1_B_old = B_old[0][0] + B_old[1][1] + B_old[2][2];

              for (int I = 0; I < 3; ++I) {
                for (int J = 0; J < 3; ++J) {
                  if (1 == solid_model) {   //Wood-Bonet J_hat  =1;
                    Cauchy[I][J]     = mus * (B[I][J]     - Id2th[I][J]) - mus / 3.*I1_B     * SolVAR[nBlocks * dim] * Id2th[I][J];
                    Cauchy_old[I][J] = mus * (B_old[I][J] - Id2th[I][J]) - mus / 3.*I1_B_old * SolVAR[nBlocks * dim] * Id2th[I][J];
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
            else if (5 == solid_model) {   //Mooney-Rivlin
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
              invB[1][0] = - (B[0][1] * B[2][2] - B[0][2] * B[2][1]) * invdetB;
              invB[2][0] = (B[0][1] * B[1][2] - B[0][2] * B[1][1]) * invdetB;
              invB[0][1] = - (B[1][0] * B[2][2] - B[1][2] * B[2][0]) * invdetB;
              invB[1][1] = (B[0][0] * B[2][2] - B[0][2] * B[2][0]) * invdetB;
              invB[2][1] = - (B[0][0] * B[1][2] - B[1][0] * B[0][2]) * invdetB;
              invB[0][2] = (B[1][0] * B[2][1] - B[2][0] * B[1][1]) * invdetB;
              invB[1][2] = - (B[0][0] * B[2][1] - B[2][0] * B[0][1]) * invdetB;
              invB[2][2] = (B[0][0] * B[1][1] - B[1][0] * B[0][1]) * invdetB;

              invB_old[0][0] = (B_old[1][1] * B_old[2][2] - B_old[1][2] * B_old[2][1]) * invdetB_old;
              invB_old[1][0] = - (B_old[0][1] * B_old[2][2] - B_old[0][2] * B_old[2][1]) * invdetB_old;
              invB_old[2][0] = (B_old[0][1] * B_old[1][2] - B_old[0][2] * B_old[1][1]) * invdetB_old;
              invB_old[0][1] = - (B_old[1][0] * B_old[2][2] - B_old[1][2] * B_old[2][0]) * invdetB_old;
              invB_old[1][1] = (B_old[0][0] * B_old[2][2] - B_old[0][2] * B_old[2][0]) * invdetB_old;
              invB_old[2][1] = - (B_old[0][0] * B_old[1][2] - B_old[1][0] * B_old[0][2]) * invdetB_old;
              invB_old[0][2] = (B_old[1][0] * B_old[2][1] - B_old[2][0] * B_old[1][1]) * invdetB_old;
              invB_old[1][2] = - (B_old[0][0] * B_old[2][1] - B_old[2][0] * B_old[0][1]) * invdetB_old;
              invB_old[2][2] = (B_old[0][0] * B_old[1][1] - B_old[1][0] * B_old[0][1]) * invdetB_old;


              I1_B = B[0][0] + B[1][1] + B[2][2];
              I2_B = B[0][0] * B[1][1] + B[1][1] * B[2][2] + B[2][2] * B[0][0]
                     - B[0][1] * B[1][0] - B[1][2] * B[2][1] - B[2][0] * B[0][2];

              I1_B_old = B_old[0][0] + B_old[1][1] + B_old[2][2];
              I2_B_old = B_old[0][0] * B_old[1][1] + B_old[1][1] * B_old[2][2] + B_old[2][2] * B_old[0][0]
                         - B_old[0][1] * B_old[1][0] - B_old[1][2] * B_old[2][1] - B_old[2][0] * B_old[0][2];

              //double C1 = mus / 3.;
              double C1 = (elementGroup == 15) ? mus1 / 3. : mus / 3.;
              double C2 = C1 / 2.;

              for (int I = 0; I < 3; ++I) {
                for (int J = 0; J < 3; ++J) {
                  Cauchy[I][J] =  2.* (C1 * B[I][J] - C2 * invB[I][J])
                                  //- (2. / 3.) * (C1 * I1_B - C2 * I2_B) * SolVAR[nBlocks * dim] * Id2th[I][J];
                                  - 1. / rhof * SolVAR[nBlocks * dim] * Id2th[I][J];

                  Cauchy_old[I][J] =  2.* (C1 * B_old[I][J] - C2 * invB_old[I][J])
                                      //- (2. / 3.) * (C1 * I1_B_old - C2 * I2_B_old) * SolVAR[nBlocks * dim] * Id2th[I][J];
                                      - 1. / rhof * SolVAR[nBlocks * dim] * Id2th[I][J];

                }
              }

            }
          }

          //END build Cauchy Stress in moving domain

          //BEGIN v=0 + Momentum (Solid)
          {
            for (unsigned i = 0; i < nve; i++) {

              //BEGIN redidual d_t - v = 0 in fixed domain
              for (int idim = 0; idim < dim; idim++) {
                aRhs[indexVAR[dim + idim]][i] +=  - phi[i] * ((-SolVAR[idim] + SolVAR_old[idim]) / dt +
                                                  (theta * SolVAR[dim + idim] + (1. - theta) * SolVAR_old[dim + idim])
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

                adept::adouble timeDerivative = - (rhos * SolVAR[dim + idim] * phi[i] * Weight
                                                   - rhos * SolVAR_old[dim + idim] * phi_old[i] * Weight_old) / dt;

                adept::adouble value =  theta * (rhos * phi[i] * _gravity[idim]      // body force
                                                 - CauchyDIR[idim]			  // stress
                                                ) * Weight;                         // at time t

                adept::adouble value_old = (1. - theta) * (rhos * phi_old[i] * _gravity[idim]      // body force
                                           - CauchyDIR_old[idim]			 // stress
                                                          ) * Weight_old;                         // at time t-dt

                aRhs[indexVAR[idim]][i] += timeDerivative + value + value_old;
              }

              //END redidual Solid Momentum in moving domain


              if (auxDisp) {
                //BEGIN auxiliary displacement equation
                adept::adouble LapAuxVAR[3] = {0., 0., 0.};

                for (int idim = 0; idim < dim; idim++) {
                  for (int jdim = 0; jdim < dim; jdim++) {
                    LapAuxVAR[idim] += ( GradSolVAR[idim + 2 * dim][jdim] + 0 * GradSolVAR[jdim + 2 * dim][idim] ) * gradphi[i * dim + jdim];
                  }
                }

//                 for (int idim = 0; idim < dim; idim++) {
//                   aRhs[indexVAR[idim + 2 * dim]][i] += (  1.0e-8 * LapAuxVAR[idim] +  ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi_hat[i] ) * Weight;
//                 }

                for (int idim = 0; idim < dim; idim++) {
                  aRhs[indexVAR[idim + 2 * dim]][i] += (  meshIsCurrupted * LapAuxVAR[idim] +  factordxi[idim] * ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi[i] ) * Weight;
                }
//                 for (int idim = 1; idim < dim; idim++) {
//                   aRhs[indexVAR[idim + 2 * dim]][i] += (  meshIsCurrupted * LapAuxVAR[idim] +  factordy * ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi[i] ) * Weight;
//                 }




//                 adept::adouble LapAuxVAR[3] = {0., 0., 0.};
//
//                 for (int idim = 0; idim < dim; idim++) {
//                   for (int jdim = 0; jdim < dim; jdim++) {
//                      LapAuxVAR[idim] += GradSolhatVAR[idim + 2 * dim][jdim] * gradphi_hat[i * dim + jdim];
//                   }
//                 }
//
//                 for (int idim = 0; idim < dim; idim++) {
//                   aRhs[indexVAR[idim + 2 * dim]][i] += (  1.0e-8 * LapAuxVAR[idim] +  ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi_hat[i] ) * Weight_hat;
//                 }
                //END auxiliary displacement equation
              }

            }
          }
          //END v=0 + Momentum (Solid)

          //BEGIN continuity block
          {
            for (unsigned i = 0; i < nve1; i++) {
              if (!penalty) {
                if (0 == solid_model) {
                  aRhs[indexVAR[nBlocks * dim]][i] += - (-phi1[i] * (I_e + (!incompressible) / lambda * SolVAR[nBlocks * dim])) * Weight_hat;
                }
                else if (1 == solid_model || 5 == solid_model) {
                  aRhs[indexVAR[nBlocks * dim]][i] += phi1[i] * ( J_hat  -  1. +
                                                      (!incompressible) / lambda * SolVAR[nBlocks * dim]) * Weight_hat;
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
      for (unsigned i = 0; i < nBlocks * dim; i++) {
        Rhs[indexVAR[i]].resize(nve);

        for (int j = 0; j < nve; j++) {
          Rhs[indexVAR[i]][j] = -aRhs[indexVAR[i]][j].value();
        }
      }

      Rhs[indexVAR[nBlocks * dim]].resize(nve1);

      for (unsigned j = 0; j < nve1; j++) {
        Rhs[indexVAR[nBlocks * dim]][j] = -aRhs[indexVAR[nBlocks * dim]][j].value();
      }

      for (int i = 0; i < nBlocks * dim + 1; i++) {
        myRES->add_vector_blocked(Rhs[indexVAR[i]], dofsVAR[i]);
      }

      if (assembleMatrix) {
        //Store equations
        for (int i = 0; i < nBlocks * dim; i++) {
          s.dependent(&aRhs[indexVAR[i]][0], nve);
          s.independent(&Soli[indexVAR[i]][0], nve);
        }

        s.dependent(&aRhs[indexVAR[nBlocks * dim]][0], nve1);
        s.independent(&Soli[indexVAR[nBlocks * dim]][0], nve1);

        Jac.resize((nBlocks * dim * nve + nve1) * (nBlocks * dim * nve + nve1));

        s.jacobian(&Jac[0], true);

        myKK->add_matrix_blocked(Jac, dofsAll, dofsAll);
        s.clear_independents();
        s.clear_dependents();

        //END local to global assembly
      }
    } //end list of elements loop

    if (assembleMatrix) myKK->close();

    myRES->close();

    setIfCorrupted->close();
    double setIfCorruptedNorm = setIfCorrupted->l1_norm();

//     std::cout << "I am in Assembly and I belived the mesh is ";
//
//     if (!meshIsCurrupted) {
//       std::cout << " not corrupted";
//       if (setIfCorruptedNorm > 0) {
//         meshIsCurrupted = true;
//         std::cout << ", but I am wrong!!!!!!!!!!!!!!!!!!!!!!!!";
// 	if (assembleMatrix) myKK->zero();
// 	myRES->zero();
// 	mysolution->ResetSolutionToOldSolution();
//         goto begin;
//       }
//       std::cout << std::endl;
//     }
//     else {
//       std::cout << " corrupted \n";
//       if (setIfCorruptedNorm > 0) {
// 	std::cout << ", but I am wrong!!!!!!!!!!!!!!!!!!!!!!!!";
//         //play with the factor TODO
//       }
//     }

    meshIsCurrupted = true;

    delete setIfCorrupted;
    delete area_elem_first;

    // *************************************
    end_time = clock();
    AssemblyTime += (end_time - start_time);
    // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************



  }

//****************************************************************************************

  void SetLambda(MultiLevelSolution &mlSol, const unsigned &level, const  FEOrder &order, Operator operatorType)
  {

    //FSIConstrainLeaflet(mlSol);

    unsigned SolType;
    if (order < FIRST || order > SECOND) {
      std::cout << "Wong Solution Order" << std::endl;
      exit(0);
    }
    else if (order == FIRST) SolType = 0;
    else if (order == SERENDIPITY) SolType = 1;
    else if (order == SECOND) SolType = 2;



    clock_t GetLambdaTime = 0;
    clock_t start_time, end_time;
    start_time = clock();

    adept::Stack & adeptStack = FemusInit::_adeptStack;

    Solution *mysolution = mlSol.GetSolutionLevel(level);
    Mesh *mymsh	=  mlSol._mlMesh->GetLevel(level);
    elem *myel	=  mymsh->el;


    unsigned indLmbd = mlSol.GetIndex("lmbd");

    const unsigned geoDim = mymsh->GetDimension();
    const unsigned nablaGoeDim = (3 * (geoDim - 1) + !(geoDim - 1));
    const unsigned max_size = static_cast< unsigned >(ceil(pow(3, geoDim)));


    const char varname[3][3] = {"DX", "DY", "DZ"};
    vector <unsigned> indVAR(geoDim);

    //const char varname1[3][4] = {"DX1", "DY1", "DZ1"};
    //vector <unsigned> indVAR1(geoDim);

    for (unsigned ivar = 0; ivar < geoDim; ivar++) {
      indVAR[ivar] = mlSol.GetIndex(&varname[ivar][0]);
      //indVAR1[ivar] = mlSol.GetIndex(&varname1[ivar][0]);
    }


    bool diffusion, elasticity;
    if (operatorType == DIFFUSION) {
      diffusion  = true;
      elasticity = false;
    }
    if (operatorType == ELASTICITY) {
      diffusion  = false;
      elasticity = true;
    }
    else {
      std::cout << "wrong operator name in SetLambda\n"
                << "valid options are diffusion or elasicity\n";
      abort();
    }

    unsigned varDim = geoDim * elasticity + diffusion;

    // local objects
    vector<vector<adept::adouble> > GradSolVAR(varDim);
    vector<vector<adept::adouble> > NablaSolVAR(varDim);

    for (int ivar = 0; ivar < varDim; ivar++) {
      GradSolVAR[ivar].resize(geoDim);
      NablaSolVAR[ivar].resize(nablaGoeDim);
    }

    vector <double > phi;
    vector <adept::adouble> gradphi;
    vector <adept::adouble> nablaphi;
    adept::adouble Weight;

    phi.reserve(max_size);
    gradphi.reserve(max_size * geoDim);
    nablaphi.reserve(max_size * nablaGoeDim);

    vector <vector < adept::adouble> > vx(geoDim);
    for (int ivar = 0; ivar < geoDim; ivar++) {
      vx[ivar].reserve(max_size);
    }
    unsigned SolTypeVx = 2.;

    vector< vector< adept::adouble > > Soli(varDim);
    vector< vector< adept::adouble > > aRhs(varDim);
    vector< vector< adept::adouble > > aLhs(varDim);
    for (int ivar = 0; ivar < varDim; ivar++) {
      Soli[ivar].reserve(max_size);
      aRhs[ivar].reserve(max_size);
      aLhs[ivar].reserve(max_size);
    }
    vector < double > K;
    K.reserve((max_size * varDim) * (max_size * varDim));
    vector < double > M;
    M.reserve((max_size * varDim) * (max_size * varDim));

    // mesh and procs
    unsigned nel    = mymsh->GetNumberOfElements();
    unsigned iproc  = mymsh->processor_id();

    // *** element loop ***
    for (int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++) {

      unsigned kel        = iel;
      short unsigned kelt = mymsh->GetElementType(kel);
      unsigned nve        = mymsh->GetElementDofNumber(kel, SolType) - 1;
      unsigned nveVx      = mymsh->GetElementDofNumber(kel, SolTypeVx);

      // -------------- resize --------------
      for (int ivar = 0; ivar < varDim; ivar++) {
        Soli[ivar].resize(nve);
        aRhs[ivar].resize(nve);
        aLhs[ivar].resize(nve);
      }

      M.resize((varDim * nve) * (varDim * nve));
      K.resize((varDim * nve) * (varDim * nve));
      // ------------------------------------

      // ------------ get coordinates -------
      for (int i = 0; i < geoDim; i++) {
        vx[i].resize(nveVx);
      }
      for (unsigned i = 0; i < nveVx; i++) {
        unsigned inodeVx_Metis = mymsh->GetSolutionDof(i, iel, SolTypeVx);
        for (int j = 0; j < geoDim; j++) {
          //coordinates
//           vx[j][i] = (*mymsh->_topology->_Sol[j])(inodeVx_Metis) +
//                      (!meshIsCurrupted) * (*mysolution->_Sol[indVAR[j]])(inodeVx_Metis) +
//                      meshIsCurrupted * (*mysolution->_Sol[indVAR1[j]])(inodeVx_Metis);
          vx[j][i] = (*mymsh->_topology->_Sol[j])(inodeVx_Metis) + (*mysolution->_Sol[indVAR[j]])(inodeVx_Metis);
        }
      }
      // ------------------------------------

      // ------------ init ------------------
      for (unsigned i = 0; i < nve; i++) {
        for (int ivar = 0; ivar < varDim; ivar++) {
          Soli[ivar][i] = 1.;
          aRhs[ivar][i] = 0.;
          aLhs[ivar][i] = 0.;
        }
      }
      // ------------------------------------

      adeptStack.new_recording();
      double hk = 1.;
      for (unsigned ig = 0; ig < mymsh->_finiteElement[kelt][SolType]->GetGaussPointNumber(); ig++) {
        // *** get Jacobian and test function and test function derivatives in the moving frame***
        mymsh->_finiteElement[kelt][SolType]->Jacobian(vx, ig, Weight, phi, gradphi, nablaphi);
        if (ig == 0) {
          double referenceElementScale[6] = {8., 1. / 6., 1., 4., 1., 2.};
          double GaussWeight = mymsh->_finiteElement[kelt][SolType]->GetGaussWeight(ig);
          double area = referenceElementScale[kelt] * Weight.value() / GaussWeight;
          hk = pow(area, 1. / geoDim);
          //cout<<hk<<endl;
          if (0 == SolType) break;
        }

        for (int ivar = 0; ivar < varDim; ivar++) {
          for (int jvar = 0; jvar < geoDim; jvar++) {
            GradSolVAR[ivar][jvar] = 0.;
          }
          for (int jvar = 0; jvar < nablaGoeDim; jvar++) {
            NablaSolVAR[ivar][jvar] = 0.;
          }
          for (unsigned inode = 0; inode < nve; inode++) {
            adept::adouble soli = Soli[ivar][inode];
            for (int jvar = 0; jvar < geoDim; jvar++) {
              GradSolVAR[ivar][jvar] += gradphi[inode * geoDim + jvar] * soli;
            }
            for (int jvar = 0; jvar < nablaGoeDim; jvar++) {
              NablaSolVAR[ivar][jvar] += nablaphi[inode * nablaGoeDim + jvar] * soli;
            }
          }
        }


        vector < adept::adouble > divGradSol(varDim, 0.);
        for (unsigned ivar = 0; ivar < varDim; ivar++) {
          for (unsigned jvar = 0; jvar < geoDim; jvar++) {
            if (diffusion) {
              divGradSol[ivar] += NablaSolVAR[ivar][jvar];
            }
            else if (elasticity) {
              unsigned kvar;
              if (ivar == jvar) kvar = jvar;
              else if (1 == ivar + jvar) kvar = geoDim;    // xy
              else if (2 == ivar + jvar) kvar = geoDim + 2;   // xz
              else if (3 == ivar + jvar) kvar = geoDim + 1;   // yz
              divGradSol[ivar]   += 0.5 * (NablaSolVAR[ivar][jvar] + NablaSolVAR[jvar][kvar]);
            }
          }
        }

        //BEGIN local assembly
        for (unsigned i = 0; i < nve; i++) {
          for (unsigned ivar = 0; ivar < varDim; ivar++) {
            for (unsigned jvar = 0; jvar < geoDim; jvar++) {
              aRhs[ivar][i] += gradphi[i * geoDim + jvar] * (GradSolVAR[ivar][jvar]) * Weight;
              if (diffusion) {
                aLhs[ivar][i] +=  divGradSol[ivar] * nablaphi[i * nablaGoeDim + jvar] * Weight;
                //aRhs[ivar][i] += gradphi[i*geoDim+jvar]*(GradSolVAR[ivar][jvar]) * Weight;
              }
              else if (elasticity) {
                unsigned kvar;
                if (ivar == jvar) kvar = jvar;
                else if (1 == ivar + jvar) kvar = geoDim;    // xy
                else if (2 == ivar + jvar) kvar = geoDim + 2;   // xz
                else if (3 == ivar + jvar) kvar = geoDim + 1;   // yz
                aLhs[ivar][i] +=  divGradSol[ivar] * 0.5 * nablaphi[i * nablaGoeDim + jvar] * Weight;
                aLhs[jvar][i] +=  divGradSol[ivar] * 0.5 * nablaphi[i * nablaGoeDim + kvar] * Weight;
                //aRhs[ivar][i] += 0.5*gradphi[i*geoDim+jvar]*0.5*(GradSolVAR[ivar][jvar]+GradSolVAR[jvar][ivar]) * Weight;
                //aRhs[jvar][i] += 0.5*gradphi[i*geoDim+ivar]*0.5*(GradSolVAR[ivar][jvar]+GradSolVAR[jvar][ivar]) * Weight;
              }
            }
          }
        }
        //END local assembly
      }
      //std::cout<<hk<<std::endl;
      double lambdak = 6. / (hk * hk);   //if SolType is linear

      if (SolType == 1 || SolType == 2) {   // only if solType is quadratic or biquadratic
        for (int ivar = 0; ivar < varDim; ivar++) {
          adeptStack.independent(&Soli[ivar][0], nve);
        }

        //Store RHS in M
        for (int ivar = 0; ivar < varDim; ivar++) {
          adeptStack.dependent(&aRhs[ivar][0], nve);
        }
        adeptStack.jacobian(&M[0]);
        adeptStack.clear_dependents();

        //Store LHS in K
        for (int ivar = 0; ivar < varDim; ivar++) {
          adeptStack.dependent(&aLhs[ivar][0], nve);
        }
        adeptStack.jacobian(&K[0]);
        adeptStack.clear_dependents();

        adeptStack.clear_independents();

        unsigned matSize = nve * varDim;

        int remove[6][3] = {{}, {}, {}, {0, 1, 2}, {0, 2, 2}, {}};

        unsigned indSize = matSize;// - remove[kelt][SolType]*elasticity;

//       for(int i=0;i<indSize;i++){
// 	for(int j=0;j<indSize;j++){
// 	  cout<<K[matSize*i+j]<<" ";
// 	}
// 	cout<<endl;
//       }
//       cout<<endl;
//       for(int i=0;i<indSize;i++){
// 	for(int j=0;j<indSize;j++){
// 	  cout<<M[matSize*i+j]<<" ";
// 	}
// 	cout<<endl;
//      }

        // LU = M factorization
        for (int k = 0 ; k < indSize - 1 ; k++) {
          for (int i = k + 1 ; i < indSize ; i++) {
            M[i * matSize + k] /= M[k * matSize + k];
            for (int j = k + 1 ; j < indSize ; j++) {
              M[i * matSize + j] -=  M[i * matSize + k] * M[k * matSize + j] ;
            }
          }
        }

        // Power Method for the largest eigenvalue of K x = lambda LU x :
        // iteration step:
        // y = U^(-1) L^(-1) K x
        // lambda= phi(y)/phi(x)
        // y = y / l2norm(y)
        // x = y


        vector < double > x(matSize, 1.);
        vector < double > y(matSize);

        double phik = x[0] + x[1];
        lambdak = 1.;
        double error = 1.;
        while (error > 1.0e-10) {
          double phikm1 = phik;
          double lambdakm1 = lambdak;

          // y = K x
          for (int i = 0; i < indSize; i++) {
            y[i] = 0.;
            for (int j = 0; j < indSize; j++) {
              y[i] += K[ i * matSize + j ] * x[j];
            }
          }

          // y = L^(-1) y
          for (int i = 0; i < indSize; i++) {
            for (int j = 0; j < i; j++) {
              y[i] -= M[i * matSize + j] * y[j];
            }
          }

          // x <--  y = U^(-1) y
          double l2norm = 0.;
          for (int i = indSize - 1; i >= 0; i--) {
            x[i] = y[i];
            for (int j = i + 1; j < indSize; j++) {
              x[i] -= M[ i * matSize + j] * x[j];
            }
            x[i] /= M[i * matSize + i];
            l2norm += x[i] * x[i];
          }
          l2norm = sqrt(l2norm);

          phik = (x[0] + x[1]);
          lambdak =  phik / phikm1;

          for (int i = 0; i < indSize; i++) {
            x[i] /= l2norm;
          }
          phik /= l2norm;
          error = fabs((lambdak - lambdakm1) / lambdak);
        }
      }

      //std::cout << lambdak*hk*hk << std::endl;
      mysolution->_Sol[indLmbd]->set(iel, sqrt(lambdak));
      //abort();
    } //end list of elements loop


    mysolution->_Sol[indLmbd]->close();
    // *************************************
    end_time = clock();
    GetLambdaTime += (end_time - start_time);

    std::cout << "GetLambda Time = " << GetLambdaTime / CLOCKS_PER_SEC << std::endl;
    //abort();
  }



  void SetLambdaNew(MultiLevelSolution &mlSol, const unsigned &level, const  FEOrder &order, Operator operatorType)
  {

    //FSIConstrainLeaflet(mlSol);

    unsigned SolType;
    if (order < FIRST || order > SECOND) {
      std::cout << "Wong Solution Order" << std::endl;
      exit(0);
    }
    else if (order == FIRST) SolType = 0;
    else if (order == SERENDIPITY) SolType = 1;
    else if (order == SECOND) SolType = 2;



    clock_t GetLambdaTime = 0;
    clock_t start_time, end_time;
    start_time = clock();

    adept::Stack & adeptStack = FemusInit::_adeptStack;

    Solution *mysolution = mlSol.GetSolutionLevel(level);
    Mesh *mymsh	=  mlSol._mlMesh->GetLevel(level);
    elem *myel	=  mymsh->el;

    unsigned indLmbd = mlSol.GetIndex("lmbd");

    const unsigned geoDim = mymsh->GetDimension();
    const unsigned nablaGoeDim = (3 * (geoDim - 1) + !(geoDim - 1));
    const unsigned max_size = static_cast< unsigned >(ceil(pow(3, geoDim)));


    const char varname[3][4] = {"DX", "DY", "DZ"};
    vector <unsigned> indVAR(geoDim);
    const char varname1[3][4] = {"DX1", "DY1", "DZ1"};
    vector <unsigned> indVAR1(geoDim);


    for (unsigned ivar = 0; ivar < geoDim; ivar++) {
      indVAR[ivar] = mlSol.GetIndex(&varname[ivar][0]);
      indVAR1[ivar] = mlSol.GetIndex(&varname1[ivar][0]);
    }


    bool diffusion, elasticity;
    if (operatorType == DIFFUSION) {
      diffusion  = true;
      elasticity = false;
    }
    if (operatorType == ELASTICITY) {
      diffusion  = false;
      elasticity = true;
    }
    else {
      std::cout << "wrong operator name in SetLambda\n"
                << "valid options are diffusion or elasicity\n";
      abort();
    }

    unsigned varDim = geoDim * elasticity + diffusion;

    // local objects
    vector<vector<adept::adouble> > GradSolVAR(varDim);
    vector<vector<adept::adouble> > NablaSolVAR(varDim);

    for (int ivar = 0; ivar < varDim; ivar++) {
      GradSolVAR[ivar].resize(geoDim);
      NablaSolVAR[ivar].resize(nablaGoeDim);
    }

    vector <double > phi;
    vector <adept::adouble> gradphi;
    vector <adept::adouble> nablaphi;
    adept::adouble Weight;

    phi.reserve(max_size);
    gradphi.reserve(max_size * geoDim);
    nablaphi.reserve(max_size * nablaGoeDim);

    vector <vector < adept::adouble> > vx(geoDim);
    vector <vector < adept::adouble> > vx1(geoDim);
    for (int ivar = 0; ivar < geoDim; ivar++) {
      vx[ivar].reserve(max_size);
      vx1[ivar].reserve(max_size);
    }
    unsigned SolTypeVx = 2.;

    vector< vector< adept::adouble > > Soli(varDim);
    vector< vector< adept::adouble > > aRhs(varDim);
    vector< vector< adept::adouble > > aLhs(varDim);
    for (int ivar = 0; ivar < varDim; ivar++) {
      Soli[ivar].reserve(max_size);
      aRhs[ivar].reserve(max_size);
      aLhs[ivar].reserve(max_size);
    }
    vector < double > K;
    K.reserve((max_size * varDim) * (max_size * varDim));
    vector < double > M;
    M.reserve((max_size * varDim) * (max_size * varDim));

    // mesh and procs
    unsigned nel    = mymsh->GetNumberOfElements();
    unsigned iproc  = mymsh->processor_id();

    NumericVector* setIfCorrupted;
    setIfCorrupted = NumericVector::build().release();
    setIfCorrupted->init(mymsh->n_processors(), 1 , false, AUTOMATIC);
    setIfCorrupted->set(iproc, 0.);

    // *** element loop ***
    for (int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++) {

      unsigned kel        = iel;
      short unsigned kelt = mymsh->GetElementType(kel);
      unsigned nve        = mymsh->GetElementDofNumber(kel, SolType) - 1;
      unsigned nveVx      = mymsh->GetElementDofNumber(kel, SolTypeVx);

      // -------------- resize --------------
      for (int ivar = 0; ivar < varDim; ivar++) {
        Soli[ivar].resize(nve);
        aRhs[ivar].resize(nve);
        aLhs[ivar].resize(nve);
      }

      M.resize((varDim * nve) * (varDim * nve));
      K.resize((varDim * nve) * (varDim * nve));
      // ------------------------------------

      // ------------ get coordinates -------
      for (int i = 0; i < geoDim; i++) {
        vx[i].resize(nveVx);
        vx1[i].resize(nveVx);
      }
      for (unsigned i = 0; i < nveVx; i++) {
        unsigned inodeVx_Metis = mymsh->GetSolutionDof(i, iel, SolTypeVx);
        for (int j = 0; j < geoDim; j++) {
          //coordinates

          vx[j][i] = (*mymsh->_topology->_Sol[j])(inodeVx_Metis) + (*mysolution->_Sol[indVAR[j]])(inodeVx_Metis);

          if (meshIsCurrupted) {
            vx1[j][i] = (*mymsh->_topology->_Sol[j])(inodeVx_Metis) + (*mysolution->_Sol[indVAR1[j]])(inodeVx_Metis);
          }
        }
      }
      // ------------------------------------

      // ------------ init ------------------
      for (unsigned i = 0; i < nve; i++) {
        for (int ivar = 0; ivar < varDim; ivar++) {
          Soli[ivar][i] = 1.;
          aRhs[ivar][i] = 0.;
          aLhs[ivar][i] = 0.;
        }
      }
      // ------------------------------------

      adeptStack.new_recording();
      double hk = 1.;
      for (unsigned ig = 0; ig < mymsh->_finiteElement[kelt][SolType]->GetGaussPointNumber(); ig++) {
        // *** get Jacobian and test function and test function derivatives in the moving frame***
        mymsh->_finiteElement[kelt][SolType]->Jacobian(vx, ig, Weight, phi, gradphi, nablaphi);
        if (!meshIsCurrupted) {
          //check if current mesh is really not corrupted
          if (Weight < 0.) setIfCorrupted->set(iproc, 1.);
        }
        else {
          // check if the current mesh is really corrupted
          if ( Weight < 0.) setIfCorrupted->set(iproc, 1.);
          mymsh->_finiteElement[kelt][SolType]->Jacobian(vx1, ig, Weight, phi, gradphi, nablaphi);
        }

        if (ig == 0) {
          double referenceElementScale[6] = {8., 1. / 6., 1., 4., 1., 2.};
          double GaussWeight = mymsh->_finiteElement[kelt][SolType]->GetGaussWeight(ig);
          double area = referenceElementScale[kelt] * Weight.value() / GaussWeight;
          hk = pow(area, 1. / geoDim);
          //cout<<hk<<endl;
          if (0 == SolType) break;
        }

        for (int ivar = 0; ivar < varDim; ivar++) {
          for (int jvar = 0; jvar < geoDim; jvar++) {
            GradSolVAR[ivar][jvar] = 0.;
          }
          for (int jvar = 0; jvar < nablaGoeDim; jvar++) {
            NablaSolVAR[ivar][jvar] = 0.;
          }
          for (unsigned inode = 0; inode < nve; inode++) {
            adept::adouble soli = Soli[ivar][inode];
            for (int jvar = 0; jvar < geoDim; jvar++) {
              GradSolVAR[ivar][jvar] += gradphi[inode * geoDim + jvar] * soli;
            }
            for (int jvar = 0; jvar < nablaGoeDim; jvar++) {
              NablaSolVAR[ivar][jvar] += nablaphi[inode * nablaGoeDim + jvar] * soli;
            }
          }
        }


        vector < adept::adouble > divGradSol(varDim, 0.);
        for (unsigned ivar = 0; ivar < varDim; ivar++) {
          for (unsigned jvar = 0; jvar < geoDim; jvar++) {
            if (diffusion) {
              divGradSol[ivar] += NablaSolVAR[ivar][jvar];
            }
            else if (elasticity) {
              unsigned kvar;
              if (ivar == jvar) kvar = jvar;
              else if (1 == ivar + jvar) kvar = geoDim;    // xy
              else if (2 == ivar + jvar) kvar = geoDim + 2;   // xz
              else if (3 == ivar + jvar) kvar = geoDim + 1;   // yz
              divGradSol[ivar]   += 0.5 * (NablaSolVAR[ivar][jvar] + NablaSolVAR[jvar][kvar]);
            }
          }
        }

        //BEGIN local assembly
        for (unsigned i = 0; i < nve; i++) {
          for (unsigned ivar = 0; ivar < varDim; ivar++) {
            for (unsigned jvar = 0; jvar < geoDim; jvar++) {
              aRhs[ivar][i] += gradphi[i * geoDim + jvar] * (GradSolVAR[ivar][jvar]) * Weight;
              if (diffusion) {
                aLhs[ivar][i] +=  divGradSol[ivar] * nablaphi[i * nablaGoeDim + jvar] * Weight;
                //aRhs[ivar][i] += gradphi[i*geoDim+jvar]*(GradSolVAR[ivar][jvar]) * Weight;
              }
              else if (elasticity) {
                unsigned kvar;
                if (ivar == jvar) kvar = jvar;
                else if (1 == ivar + jvar) kvar = geoDim;    // xy
                else if (2 == ivar + jvar) kvar = geoDim + 2;   // xz
                else if (3 == ivar + jvar) kvar = geoDim + 1;   // yz
                aLhs[ivar][i] +=  divGradSol[ivar] * 0.5 * nablaphi[i * nablaGoeDim + jvar] * Weight;
                aLhs[jvar][i] +=  divGradSol[ivar] * 0.5 * nablaphi[i * nablaGoeDim + kvar] * Weight;
                //aRhs[ivar][i] += 0.5*gradphi[i*geoDim+jvar]*0.5*(GradSolVAR[ivar][jvar]+GradSolVAR[jvar][ivar]) * Weight;
                //aRhs[jvar][i] += 0.5*gradphi[i*geoDim+ivar]*0.5*(GradSolVAR[ivar][jvar]+GradSolVAR[jvar][ivar]) * Weight;
              }
            }
          }
        }
        //END local assembly
      }
      //std::cout<<hk<<std::endl;
      double lambdak = 6. / (hk * hk);   //if SolType is linear

      if (SolType == 1 || SolType == 2) {   // only if solType is quadratic or biquadratic
        for (int ivar = 0; ivar < varDim; ivar++) {
          adeptStack.independent(&Soli[ivar][0], nve);
        }

        //Store RHS in M
        for (int ivar = 0; ivar < varDim; ivar++) {
          adeptStack.dependent(&aRhs[ivar][0], nve);
        }
        adeptStack.jacobian(&M[0]);
        adeptStack.clear_dependents();

        //Store LHS in K
        for (int ivar = 0; ivar < varDim; ivar++) {
          adeptStack.dependent(&aLhs[ivar][0], nve);
        }
        adeptStack.jacobian(&K[0]);
        adeptStack.clear_dependents();

        adeptStack.clear_independents();

        unsigned matSize = nve * varDim;

        int remove[6][3] = {{}, {}, {}, {0, 1, 2}, {0, 2, 2}, {}};

        unsigned indSize = matSize;// - remove[kelt][SolType]*elasticity;

//       for(int i=0;i<indSize;i++){
// 	for(int j=0;j<indSize;j++){
// 	  cout<<K[matSize*i+j]<<" ";
// 	}
// 	cout<<endl;
//       }
//       cout<<endl;
//       for(int i=0;i<indSize;i++){
// 	for(int j=0;j<indSize;j++){
// 	  cout<<M[matSize*i+j]<<" ";
// 	}
// 	cout<<endl;
//      }

        // LU = M factorization
        for (int k = 0 ; k < indSize - 1 ; k++) {
          for (int i = k + 1 ; i < indSize ; i++) {
            M[i * matSize + k] /= M[k * matSize + k];
            for (int j = k + 1 ; j < indSize ; j++) {
              M[i * matSize + j] -=  M[i * matSize + k] * M[k * matSize + j] ;
            }
          }
        }

        // Power Method for the largest eigenvalue of K x = lambda LU x :
        // iteration step:
        // y = U^(-1) L^(-1) K x
        // lambda= phi(y)/phi(x)
        // y = y / l2norm(y)
        // x = y


        vector < double > x(matSize, 1.);
        vector < double > y(matSize);

        double phik = x[0] + x[1];
        lambdak = 1.;
        double error = 1.;
        while (error > 1.0e-10) {
          double phikm1 = phik;
          double lambdakm1 = lambdak;

          // y = K x
          for (int i = 0; i < indSize; i++) {
            y[i] = 0.;
            for (int j = 0; j < indSize; j++) {
              y[i] += K[ i * matSize + j ] * x[j];
            }
          }

          // y = L^(-1) y
          for (int i = 0; i < indSize; i++) {
            for (int j = 0; j < i; j++) {
              y[i] -= M[i * matSize + j] * y[j];
            }
          }

          // x <--  y = U^(-1) y
          double l2norm = 0.;
          for (int i = indSize - 1; i >= 0; i--) {
            x[i] = y[i];
            for (int j = i + 1; j < indSize; j++) {
              x[i] -= M[ i * matSize + j] * x[j];
            }
            x[i] /= M[i * matSize + i];
            l2norm += x[i] * x[i];
          }
          l2norm = sqrt(l2norm);

          phik = (x[0] + x[1]);
          lambdak =  phik / phikm1;

          for (int i = 0; i < indSize; i++) {
            x[i] /= l2norm;
          }
          phik /= l2norm;
          error = fabs((lambdak - lambdakm1) / lambdak);
        }
      }

      //std::cout << lambdak*hk*hk << std::endl;
      mysolution->_Sol[indLmbd]->set(iel, sqrt(lambdak));
      //abort();
    } //end list of elements loop


    mysolution->_Sol[indLmbd]->close();
    // *************************************
    end_time = clock();
    GetLambdaTime += (end_time - start_time);

    std::cout << "GetLambda Time = " << GetLambdaTime / CLOCKS_PER_SEC << std::endl;


    setIfCorrupted->close();
    double setIfCorruptedNorm = setIfCorrupted->l1_norm();


//     std::cout << "I am in Set Lambda and I belived the mesh is ";
//     if (!meshIsCurrupted) {
//       std::cout << " not corrupted";
//       if (setIfCorruptedNorm > 0) {
//         meshIsCurrupted = true;
//         std::cout << ", but I am wrong!!!!!!!!!!!!!!!!!!!!!!!!";
//       }
//       std::cout << std::endl;
//     }
//     else {
//       std::cout << " corrupted";
//       if (setIfCorruptedNorm < 1.0e-10) {
//         meshIsCurrupted = false;
//         std::cout << ", but I am wrong!!!!!!!!!!!!!!!!!!!!!!!!";
//       }
//       std::cout << std::endl;
//     }

    meshIsCurrupted = true;

    delete setIfCorrupted;


    //abort();
  }



  const double leaflet[129][2] = {
    {   -0.00402782    ,  0.0601179     },  {   -0.00391292    ,	0.0601527     },  {   -0.00379967    ,	0.0601912     },  {   -0.00368805    ,	0.0602334     },
    {   -0.00357807    ,	0.0602794     },  {   -0.00346973    ,	0.0603291     },  {   -0.00336303    ,	0.0603826     },  {   -0.00325798    ,	0.0604397     },
    {   -0.00315456    ,	0.0605006     },  {   -0.00305329    ,	0.0605649     },  {   -0.00295415    ,	0.0606319     },  {   -0.00285713    ,	0.0607016     },
    {   -0.00276223    ,	0.060774      },  {   -0.00266945    ,	0.060849      },  {   -0.00257879    ,	0.0609268     },  {   -0.00249026    ,	0.0610072     },
    {   -0.00240385    ,	0.0610904     },  {   -0.00231958    ,	0.0611756     },  {   -0.00223735    ,	0.0612626     },  {   -0.00215715    ,	0.0613513     },
    {   -0.00207899    ,	0.0614417     },  {   -0.00200286    ,	0.0615339     },  {   -0.00192876    ,	0.0616277     },  {   -0.0018567     ,  0.0617232     },
    {   -0.00178666    ,	0.0618205     },  {   -0.00171853    ,	0.061919      },  {   -0.00165215    ,	0.0620186     },  {   -0.00158754    ,	0.0621193     },
    {   -0.00152469    ,	0.062221      },  {   -0.00146361    ,	0.0623238     },  {   -0.00140429    ,	0.0624276     },  {   -0.00134673    ,	0.0625325     },
    {   -0.00129093    ,	0.0626385     },  {   -0.00123672    ,	0.0627453     },  {   -0.00118402    ,	0.0628528     },  {   -0.00113285    ,	0.062961      },
    {   -0.00108319    ,	0.0630698     },  {   -0.00103505    ,  0.0631794     },  {   -0.000988434   ,	0.0632896     },  {   -0.000943336   ,	0.0634005     },
    {   -0.000899757   ,	0.063512      },  {   -0.000857621   ,	0.0636239     },  {   -0.000816835   ,	0.0637363     },  {   -0.000777398   ,	0.0638492     },
    {   -0.000739311   ,	0.0639624     },  {   -0.000702573   ,	0.0640761     },  {   -0.000667185   ,	0.0641903     },  {   -0.000633146   ,	0.0643049     },
    {   -0.000600456   ,	0.0644199     },  {   -0.000568998   ,	0.0645354     },  {   -0.000538812   ,	0.0646513     },  {   -0.000509899   ,	0.0647674     },
    {   -0.000482258   ,	0.0648837     },  {   -0.00045589    ,  0.0650004     },  {   -0.000430794   ,	0.0651174     },  {   -0.000406971   ,	0.0652347     },
    {   -0.00038442    ,  0.0653523     },  {   -0.000363115   ,	0.06547       },  {   -0.000343096   ,	0.0655879     },  {   -0.000324363   ,	0.0657059     },
    {   -0.000306917   ,	0.0658242     },  {   -0.000290758   ,	0.0659426     },  {   -0.000275886   ,	0.0660612     },  {   -0.000262299   ,	0.06618       },
    {   -0.00025	     ,  0.066299      },  {   -0.000239534   ,	0.0664144     },  {   -0.000230166   ,	0.0665298     },  {   -0.000221896   ,	0.0666454     },
    {   -0.000214725   ,	0.0667609     },  {   -0.000208652   ,	0.0668766     },  {   -0.000203678   ,	0.0669923     },  {   -0.000199802   ,	0.067108      },
    {   -0.000197024   ,	0.0672238     },  {   -0.000195205   ,	0.0673395     },  {   -0.000194055   ,	0.0674552     },  {   -0.000193574   ,	0.0675708     },
    {   -0.000193763   ,	0.0676864     },  {   -0.000194621   ,	0.0678021     },  {   -0.000196148   ,	0.0679177     },  {   -0.000198345   ,	0.0680334     },
    {   -0.000201211   ,	0.068149      },  {   -0.000204526   ,	0.0682648     },  {   -0.000207968   ,	0.0683806     },  {   -0.000211535   ,	0.0684963     },
    {   -0.000215229   ,	0.0686121     },  {   -0.000219049   ,	0.0687279     },  {   -0.000222994   ,	0.0688437     },  {   -0.000227066   ,	0.0689594     },
    {   -0.000231264   ,	0.0690752     },  {   -0.000235294   ,	0.0691908     },  {   -0.000238841   ,	0.0693063     },  {   -0.000241906   ,	0.0694219     },
    {   -0.000244489   ,	0.0695375     },  {   -0.00024659    ,  0.0696531     },  {   -0.000248209   ,	0.0697687     },  {   -0.000249345   ,	0.0698844     },
    {   -0.00025       ,  0.07          },  {   -0.000255613   ,	0.0701159     },  {   -0.000261276   ,	0.0702317     },  {   -0.00026699    ,  0.0703476     },
    {   -0.000272754   ,	0.0704635     },  {   -0.000278569   ,	0.0705793     },  {   -0.000284435   ,	0.0706952     },  {   -0.00029035    ,  0.070811      },
    {   -0.000296317   ,	0.0709269     },  {   -0.000302377   ,	0.0710426     },  {   -0.000308591   ,	0.0711583     },  {   -0.000314961   ,	0.0712741     },
    {   -0.000321485   ,	0.0713898     },  {   -0.000328164   ,	0.0715055     },  {   -0.000334998   ,	0.0716212     },  {   -0.000341987   ,	0.0717368     },
    {   -0.000349131   ,	0.0718525     },  {   -0.000356484   ,	0.0719682     },  {   -0.000364109   ,	0.072084      },  {   -0.000372005   ,	0.0721997     },
    {   -0.000380173   ,	0.0723154     },  {   -0.000388611   ,	0.072431      },  {   -0.000397321   ,	0.0725467     },  {   -0.000406302   ,	0.0726623     },
    {   -0.000415554   ,	0.0727779     },  {   -0.000425128   ,	0.0728934     },  {   -0.00043511    ,  0.0730089     },  {   -0.0004455     ,  0.0731244     },
    {   -0.000456299   ,	0.0732398     },  {   -0.000467505   ,	0.0733551     },  {   -0.000479121   ,	0.0734705     },  {   -0.000491144   ,	0.0735857     },
    {   -0.000503576   ,	0.073701      }
  };

  bool MeshIsCurrupted(const std::vector < double >& x, const std::vector < double >& dX)
  {
    bool movedNode = false;
    double epsilon = .0; // maximum leaflet colsing
    double  delta = 7.1e-05; //leaflet thicknes;
    unsigned N = 128;
    if (x[1] >= leaflet[0][1] && x[1] <= leaflet[N][1]) {
      unsigned i0 = 0;
      unsigned i1 = N;
      while (i1 - i0 > 1) {
        unsigned i2 = i0 + (i1 - i0) / 2;
        if (leaflet[i2][1] > x[1]) {
          i1 = i2;
        }
        else {
          i0 = i2;
        }
      }

      double s = (x[1] - leaflet[i0][1]) / (leaflet[i1][1] - leaflet[i0][1]);
      double xl = leaflet[i0][0] * (1. - s) + leaflet[i1][0] * s;

      if (x[0] - xl >= 0) { //on the right of the lealflet
        if (x[0] + dX[0] > -epsilon * x[0] / xl) {
          //dX[0] = -epsilon * x[0] / xl - x[0] ;
          movedNode = true;
        }
      }
//       else if(x[0] - xl >= -delta) {  //inside the lealflet
//         if(x[0] + dX[0] > -epsilon + (x[0] - xl)) {
//           dX[0] = -epsilon - (xl - x[0]) - x[0];
//           movedNode = true;
//         }
//       }
//       else { //on the left of the lealflet
//         if(x[0] + dX[0] > -epsilon - delta - 3.* epsilon * (xl - delta - x[0]) / 0.001) {
//           dX[0] = -epsilon - delta - 3.* epsilon * (xl - delta - x[0]) / 0.001 - x[0];
//         }
//       }
    }
    return movedNode;
  }


  unsigned counter = 0;

  void FSIConstrainLeaflet(MultiLevelSolution& mlSol)
  {

    unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

    Solution* sol  = mlSol.GetSolutionLevel(level);
    Mesh* msh = mlSol._mlMesh->GetLevel(level);

    const unsigned dim = msh->GetDimension();

    //----------------------------------------------------------------------------------
    //variable-name handling
    const char varname[3][3] = {"DX", "DY", "DZ"};
    const char varname1[3][4] = {"DX1", "DY1", "DZ1"};
    vector <unsigned> indVAR(dim);
    vector <unsigned> indVAR1(dim);
    unsigned SolType;

    for (unsigned k = 0; k < dim; k++) {
      indVAR[k] = mlSol.GetIndex(&varname[k][0]);
      indVAR1[k] = mlSol.GetIndex(&varname1[k][0]);
      SolType = mlSol.GetSolutionType(&varname[k][0]);
    }

    std::vector < double > x(dim);
    std::vector < double > dx(dim);
    std::vector < double > dxOld(dim);

    unsigned iproc  = msh->processor_id();

    int nprocs = msh->n_processors();
    NumericVector* meshMoved;
    NumericVector* meshOldMoved;
    meshMoved = NumericVector::build().release();
    meshOldMoved = NumericVector::build().release();

    if (nprocs == 1) {
      meshMoved->init(nprocs, 1, false, SERIAL);
      meshOldMoved->init(nprocs, 1, false, SERIAL);
    }
    else {
      meshMoved->init(nprocs, 1, false, PARALLEL);
      meshOldMoved->init(nprocs, 1, false, PARALLEL);
    }
    meshMoved->zero();
    meshOldMoved->zero();

    for (unsigned idof = msh->_dofOffset[SolType][iproc]; idof < msh->_dofOffset[SolType][iproc + 1]; idof++) {
      for (unsigned k = 0; k < dim; k++) {
        x[k] = (*msh->_topology->_Sol[k])(idof);
        dx[k] = (*sol->_Sol[indVAR[k]])(idof);
        dxOld[k] = (*sol->_SolOld[indVAR[k]])(idof);
      }
      bool movedNode;
      movedNode = MeshIsCurrupted(x, dx);
      if (movedNode) {
        //for(unsigned k = 0; k < dim; k++) {
        //sol->_Sol[indVAR[k]]->set(idof, dx[k]);
        //}
        meshMoved->set(iproc, 1.);
      }
      movedNode = MeshIsCurrupted(x, dxOld);
      if (movedNode) {
//         for(unsigned k = 0; k < dim; k++) {
//           sol->_SolOld[indVAR[k]]->set(idof, dxOld[k]);
//         }
        meshOldMoved->set(iproc, 1.);
      }
    }

    meshMoved->close();
    meshOldMoved->close();

    if (meshMoved->l1_norm() < 0.5) {
      meshIsCurrupted = false;
      for (unsigned k = 0; k < dim; k++) {
        *(sol->_Sol[indVAR1[k]]) = *(sol->_Sol[indVAR[k]]);
      }
    }
    else {
      std::cout << "Warning the Mesh Constrained Function has been called for DX" << std::endl;
      for (unsigned k = 0; k < dim; k++) {
        meshIsCurrupted = true;
        //*(sol->_Sol[indVAR[k]]) = *(sol->_Sol[indVAR1[k]]);
      }
    }


    if (meshOldMoved->l1_norm() > 0.5) {
      std::cout << "Warning the Mesh Constrained Function has been called for DXOld" << std::endl;
      for (unsigned k = 0; k < dim; k++) {
        *(sol->_SolOld[indVAR[k]]) = *(sol->_Sol[indVAR1[k]]);
      }
    }


    //if( meshMoved->l1_norm() ){
    std::vector<std::string> print_vars;
    print_vars.push_back("All");

    mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 1000 + counter);

    counter++;
    //}

  }





  void FSITimeDependentAssemblySupgNew2(MultiLevelProblem & ml_prob)
  {

    clock_t AssemblyTime = 0;
    clock_t start_time, end_time;

    bool auxDisp = true;
    unsigned nBlocks = (auxDisp) ? 3 : 2;

    //meshIsCurrupted = true;
    //pointers and references

    //MonolithicFSINonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<MonolithicFSINonLinearImplicitSystem>("Fluid-Structure-Interaction");
    TransientNonlinearImplicitSystem & my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("Fluid-Structure-Interaction");
    const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();

    MultiLevelSolution	* ml_sol		= ml_prob._ml_sol;

    //FSIConstrainLeaflet(*ml_sol);

    Solution	* mysolution		= ml_sol->GetSolutionLevel(level);

    LinearEquationSolver * myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];
    Mesh	*	mymsh		=  ml_prob._ml_msh->GetLevel(level);
    elem	*	myel		=  mymsh->el;
    SparseMatrix	* myKK		=  myLinEqSolver->_KK;
    NumericVector	* myRES		=  myLinEqSolver->_RES;


    bool assembleMatrix = my_nnlin_impl_sys.GetAssembleMatrix();
    // call the adept stack object
    adept::Stack& s = FemusInit::_adeptStack;
    if (assembleMatrix) s.continue_recording();
    else s.pause_recording();

    const unsigned dim = mymsh->GetDimension();
    const unsigned nabla_dim = 3 * (dim - 1);
    const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

    double theta = 0.5;
    //double theta = 1.;

    // local objects
    vector<adept::adouble> SolVAR(nBlocks * dim + 1);
    vector<double> SolVAR_old(nBlocks * dim);
    vector<double> meshVel_old(dim);

    vector<vector < adept::adouble > > GradSolVAR(nBlocks * dim);
    vector<vector < double > > GradSolVAR_old(nBlocks * dim);
    vector<vector < double > > GradMeshVel_old(dim);

    vector<vector < adept::adouble > > GradSolhatVAR(nBlocks * dim);
    vector<vector < double > > GradSolhatVAR_old(nBlocks * dim);

    vector<vector<adept::adouble> > NablaSolVAR(nBlocks * dim);
    vector<vector < double > > NablaSolVAR_old(nBlocks * dim);

    for (int i = 0; i < nBlocks * dim; i++) {
      GradSolVAR[i].resize(dim);
      GradSolVAR_old[i].resize(dim);

      GradSolhatVAR[i].resize(dim);
      GradSolhatVAR_old[i].resize(dim);

      NablaSolVAR[i].resize(nabla_dim);
      NablaSolVAR_old[i].resize(nabla_dim);
    }
    for (int i = 0; i < dim; i++) {
      GradMeshVel_old[i].resize(dim);
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

    const double * phi1;

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

    vector< vector< adept::adouble > > Soli(nBlocks * dim + 1);
    vector< vector< double > > Soli_old(nBlocks * dim + 1);
    vector< vector< int > > dofsVAR(nBlocks * dim + 1);
    vector< vector< double > > meshVelOldNode(dim);



    for (int i = 0; i < nBlocks * dim + 1; i++) {
      Soli[i].reserve(max_size);
      Soli_old[i].reserve(max_size);
      dofsVAR[i].reserve(max_size);
    }
    for (int i = 0; i < dim; i++) {
      meshVelOldNode[i].reserve(max_size);
    }

    vector< vector< double > > Rhs(nBlocks * dim + 1);
    vector< vector< adept::adouble > > aRhs(nBlocks * dim + 1);

    for (int i = 0; i < nBlocks * dim + 1; i++) {
      aRhs[i].reserve(max_size);
      Rhs[i].reserve(max_size);
    }


    vector < int > dofsAll;
    dofsAll.reserve(max_size * (nBlocks * dim + 1));

    vector < double > Jac;
    Jac.reserve(dim * max_size * (nBlocks * dim + 1) *dim * max_size * (nBlocks * dim + 1));

    // ------------------------------------------------------------------------
    // Physical parameters
    double rhof	 	= ml_prob.parameters.get<Fluid> ("Fluid").get_density();
    double rhos	 	= ml_prob.parameters.get<Fluid> ("Solid").get_density() / rhof;
    double mu_lame 	= ml_prob.parameters.get<Solid> ("Solid").get_lame_shear_modulus();
    double lambda_lame 	= ml_prob.parameters.get<Solid> ("Solid").get_lame_lambda();
    double mus		= mu_lame / rhof;
    double mu_lame1 	= ml_prob.parameters.get < Solid> ("Solid1").get_lame_shear_modulus();
    double mus1 	= mu_lame1 / rhof;
    double IRe 		= ml_prob.parameters.get<Fluid> ("Fluid").get_IReynolds_number();
    double lambda	= lambda_lame / rhof;
    double betans	= 1.;
    int    solid_model	= ml_prob.parameters.get<Solid> ("Solid").get_physical_model();

    if (solid_model >= 2 && solid_model <= 4) {
      std::cout << "Error! Solid Model " << solid_model << "not implemented\n";
      abort();
    }

    bool incompressible = (0.5 == ml_prob.parameters.get<Solid> ("Solid").get_poisson_coeff()) ? 1 : 0;
    const bool penalty = ml_prob.parameters.get<Solid> ("Solid").get_if_penalty();

    // gravity
    double _gravity[3] = {0., 0., 0.};

    double dt =  my_nnlin_impl_sys.GetIntervalTime();
    double time =  my_nnlin_impl_sys.GetTime();
    // -----------------------------------------------------------------
    // space discretization parameters
    unsigned SolType2 = ml_sol->GetSolutionType(ml_sol->GetIndex("U"));
    unsigned SolType1 = ml_sol->GetSolutionType(ml_sol->GetIndex("PS"));

    // mesh and procs
    unsigned nel    = mymsh->GetNumberOfElements();
    unsigned igrid  = mymsh->GetLevel();
    unsigned iproc  = mymsh->processor_id();

    unsigned indLmbd = ml_sol->GetIndex("lmbd");

    //----------------------------------------------------------------------------------
    //variable-name handling
    const char varname[10][4] = {"DX", "DY", "DZ", "U", "V", "W", "DX1", "DY1", "DZ1", "PS"};

    const char varname2[10][4] = {"Um", "Vm", "Wm"};

    vector <unsigned> indexVAR(nBlocks * dim + 1);
    vector <unsigned> indVAR(nBlocks * dim + 1);
    vector <unsigned> SolType(nBlocks * dim + 1);


    vector <unsigned> indVAR2(dim);

    for (unsigned ivar = 0; ivar < dim; ivar++) {
      for (unsigned k = 0; k < nBlocks; k++) {
        indVAR[ivar + k * dim] = ml_sol->GetIndex(&varname[ivar + k * 3][0]);
        SolType[ivar + k * dim] = ml_sol->GetSolutionType(&varname[ivar + k * 3][0]);
        indexVAR[ivar + k * dim] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar + k * 3][0]);
      }
      indVAR2[ivar] = ml_sol->GetIndex(&varname2[ivar][0]);
    }

    indexVAR[nBlocks * dim] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[9][0]);
    indVAR[nBlocks * dim] = ml_sol->GetIndex(&varname[9][0]);
    SolType[nBlocks * dim] = ml_sol->GetSolutionType(&varname[9][0]);
    //----------------------------------------------------------------------------------

    int nprocs = mymsh->n_processors();
    NumericVector * area_elem_first;
    area_elem_first = NumericVector::build().release();

    if (nprocs == 1) {
      area_elem_first->init(nprocs, 1, false, SERIAL);
    }
    else {
      area_elem_first->init(nprocs, 1, false, PARALLEL);
    }


    double rapresentative_area = 1.;

    start_time = clock();

    if (assembleMatrix) myKK->zero();

    NumericVector* setIfCorrupted;
    setIfCorrupted = NumericVector::build().release();
    setIfCorrupted->init(mymsh->n_processors(), 1 , false, AUTOMATIC);


  begin:

    area_elem_first->zero();
    setIfCorrupted->zero();

    // *** element loop ***
    for (int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++) {

      short unsigned ielt = mymsh->GetElementType(iel);
      unsigned nve        = mymsh->GetElementDofNumber(iel, SolType2);
      unsigned nve1       = mymsh->GetElementDofNumber(iel, SolType1);
      int flag_mat        = mymsh->GetElementMaterial(iel);
      unsigned elementGroup = mymsh->GetElementGroup(iel);

      // *******************************************************************************************************

      //initialization of everything is in common fluid and solid

      //Rhs
      for (int i = 0; i < nBlocks * dim; i++) {
        dofsVAR[i].resize(nve);
        Soli[indexVAR[i]].resize(nve);
        Soli_old[indexVAR[i]].resize(nve);
        aRhs[indexVAR[i]].resize(nve);
      }
      for (int i = 0; i < dim; i++) {
        meshVelOldNode[i].resize(nve);
      }

      dofsVAR[nBlocks * dim].resize(nve1);
      Soli[indexVAR[nBlocks * dim]].resize(nve1);
      Soli_old[indexVAR[nBlocks * dim]].resize(nve1);
      aRhs[indexVAR[nBlocks * dim]].resize(nve1);

      dofsAll.resize(0);

      Jac.resize((nBlocks * dim * nve + nve1) * (nBlocks * dim * nve + nve1));

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
        solidmark[i] = mymsh->GetSolidMark(idof);    // to check

        for (int j = 0; j < dim; j++) {
          for (unsigned k = 0; k < nBlocks; k++) {
            Soli[indexVAR[j + k * dim]][i] = (*mysolution->_Sol[indVAR[j + k * dim]])(idof);
            Soli_old[indexVAR[j + k * dim]][i] = (*mysolution->_SolOld[indVAR[j + k * dim]])(idof);
            aRhs[indexVAR[j + k * dim]][i] = 0.;
            dofsVAR[j + k * dim][i] = myLinEqSolver->GetSystemDof(indVAR[j + k * dim], indexVAR[j + k * dim], i, iel);
          }
          meshVelOldNode[j][i] = (*mysolution->_SolOld[indVAR2[j]])(idof);
          vx_hat[j][i] = (*mymsh->_topology->_Sol[j])(idof);
        }
      }

      // pressure dofs
      for (unsigned i = 0; i < nve1; i++) {
        unsigned idof = mymsh->GetSolutionDof(i, iel, SolType[nBlocks * dim]);
        dofsVAR[nBlocks * dim][i] = myLinEqSolver->GetSystemDof(indVAR[nBlocks * dim], indexVAR[nBlocks * dim], i, iel);
        Soli[indexVAR[nBlocks * dim]][i]     = (*mysolution->_Sol[indVAR[nBlocks * dim]])(idof);
        Soli_old[indexVAR[nBlocks * dim]][i] = (*mysolution->_SolOld[indVAR[nBlocks * dim]])(idof);
        aRhs[indexVAR[nBlocks * dim]][i] = 0.;
      }

      // build dof ccomposition
      for (int idim = 0; idim < nBlocks * dim; idim++) {
        dofsAll.insert(dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end());
      }
      dofsAll.insert(dofsAll.end(), dofsVAR[nBlocks * dim].begin(), dofsVAR[nBlocks * dim].end());

      if (assembleMatrix) s.new_recording();


      for (int j = 0; j < nve; j++) {
        unsigned idof = mymsh->GetSolutionDof(j, iel, SolType2);
        for (unsigned idim = 0; idim < dim; idim++) {
          vx[idim][j] = vx_hat[idim][j] +  Soli[indexVAR[idim + meshIsCurrupted * 2 * dim]][j];
          vx_old[idim][j] = vx_hat[idim][j] + Soli_old[indexVAR[idim + meshIsCurrupted * 2 * dim]][j];
        }
      }

      //Boundary integral
      {

        vector < adept::adouble> normal(dim, 0);
        vector < double > normal_old(dim, 0);

        // loop on faces
        for (unsigned jface = 0; jface < mymsh->GetElementFaceNumber(iel); jface++) {
          std::vector< double > xx(dim, 0.);

          // look for boundary faces
          if (myel->GetFaceElementIndex(iel, jface) < 0) {

            unsigned int face = - (mymsh->el->GetFaceElementIndex(iel, jface) + 1);
            double tau = 0.;
            double tau_old = 0.;
            if ((!ml_sol->GetBdcFunction()(xx, "PS", tau, face, time) &&
                 !ml_sol->GetBdcFunction()(xx, "PS", tau_old, face, time - dt))
                && (tau != 0. || tau_old != 0.)) {

              unsigned nve = mymsh->GetElementFaceDofNumber(iel, jface, SolType2);
              const unsigned felt = mymsh->GetElementFaceType(iel, jface);

              for (unsigned i = 0; i < nve; i++) {
                unsigned int ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);
                unsigned idof = mymsh->GetSolutionDof(ilocal, iel, 2);

                for (unsigned idim = 0; idim < dim; idim++) {
                  vx_face[idim][i]    = (*mymsh->_topology->_Sol[idim])(idof) + Soli[indexVAR[idim + meshIsCurrupted * 2 * dim]][ilocal]; //TODO
                  vx_face_old[idim][i] = (*mymsh->_topology->_Sol[idim])(idof) + Soli_old[indexVAR[idim + meshIsCurrupted * 2 * dim]][ilocal];
                }
              }

              for (unsigned igs = 0; igs < mymsh->_finiteElement[felt][SolType2]->GetGaussPointNumber(); igs++) {
                mymsh->_finiteElement[felt][SolType2]->JacobianSur(vx_face, igs, Weight, phi, gradphi, normal);
                mymsh->_finiteElement[felt][SolType2]->JacobianSur(vx_face_old, igs, Weight_old, phi_old, gradphi_old, normal_old);

                //phi1 =mymsh->_finiteElement[felt][SolType2]->GetPhi(igs);
                // *** phi_i loop ***
                for (unsigned i = 0; i < nve; i++) {
                  adept::adouble value = - phi[i] * tau / rhof * Weight;
                  double value_old = - phi_old[i] * tau_old / rhof * Weight_old;
                  unsigned int ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);

                  for (unsigned idim = 0; idim < dim; idim++) {
                    if ((!solidmark[ilocal])) {
                      aRhs[indexVAR[dim + idim]][ilocal]   += (theta * value + (1. - theta) * value_old) * normal[idim];
                    }
                    else {   //if interface node it goes to solid
                      aRhs[indexVAR[idim]][ilocal]   += (theta * value + (1. - theta) * value_old) * normal[idim];
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

        if ( Weight.value() < 0. ) setIfCorrupted->set(iproc, 1.);

        if (flag_mat == 2 || flag_mat == 3  || iel == mymsh->_elementOffset[iproc]) {
          if (ig == 0) {
            double GaussWeight = mymsh->_finiteElement[ielt][SolType2]->GetGaussWeight(ig);
            area = Weight_hat / GaussWeight;

            if (iel == mymsh->_elementOffset[iproc]) {
              area_elem_first->add(mymsh->processor_id(), area);
              area_elem_first->close();
              rapresentative_area = area_elem_first->l1_norm() / nprocs;
            }
          }

          Weight_nojac = Weight;// / area * 1.e-10;// * rapresentative_area;

          //-----------------------------------------------------------------------//
          // for vein_valve mesh
          std::vector<double> xc(dim, 0.);
          //xc[0] = -0.001013; // vein_valve
          //xc[1] = 0.07000;
//           xc[0] = -0.000286; // vein_valve_closed
//           xc[1] = 0.07000;
//
//           if (dim == 2) {
//             xc[0] = -0.00025; // vein_valve_closed in valve2.neu
//             xc[1] = 0.07000;
//           }
//           if (dim == 2) {
//             xc[0] = -0.000220; // vein_valve_closed valve2_corta.neu
//             xc[1] = 0.068718;
//           }
           if (dim == 2) {
            xc[0] = -0.000193; // vein_valve_closed valve2_corta2.neu
            xc[1] = 0.067541;
          }
// 	  else if (dim == 3){
// 	    xc[0] = -0.008; // vein_valve_closed
// 	    xc[1] = 0.0;
// 	  }
          else if (dim == 3) {
            xc[0] = 0.0015; // vein_valve_closed
            xc[1] = 0.0;
            xc[2] = 0.01;
          }

          double distance = 0.;
          for (unsigned k = 0; k < dim; k++) {
            distance += (vx_hat[k][ nve - 1] - xc[k]) * (vx_hat[k][nve - 1] - xc[k]);
          }
          distance = sqrt(distance);
          Weight_nojac *= 1. / (1 + 10000 * distance);

          if (elementGroup == 16) Weight_nojac *= 100;
          //-----------------------------------------------------------------------//
        }

        // ---------------------------------------------------------------------------
        // displacement and velocity
        
        for (int i = 0; i < nBlocks * dim; i++) {
          SolVAR[i] = 0.;
          SolVAR_old[i] = 0.;

          for (int j = 0; j < dim; j++) {
            GradSolVAR[i][j] = 0.;
            GradSolVAR_old[i][j] = 0.;

            GradSolhatVAR[i][j] = 0.;
            GradSolhatVAR_old[i][j] = 0.;
          }

          for (int j = 0; j < nabla_dim; j++) {
            NablaSolVAR[i][j] = 0.;
            NablaSolVAR_old[i][j] = 0.;
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
            for (int j = 0; j < nabla_dim; j++) {
              NablaSolVAR[i][j]     += nablaphi[inode * nabla_dim + j] * Soli[indexVAR[i]][inode];
              NablaSolVAR_old[i][j] += nablaphi_old[inode * nabla_dim + j] * Soli_old[indexVAR[i]][inode];
            }
          }
        }

        // pressure
        SolVAR[nBlocks * dim] = 0.;
        for (unsigned inode = 0; inode < nve1; inode++) {
          SolVAR[nBlocks * dim]    += phi1[inode] * Soli[indexVAR[nBlocks * dim]][inode];
        }

        vector < adept::adouble > vx_ig(dim);
	vector < double > vxOld_ig(dim);
        
        // mesh velocity OLD
        for (int i = 0; i < dim; i++) {
	  vx_ig[i]=0.;
	  vxOld_ig[i]=0.;
          meshVel_old[i] = 0.;
          for (int j = 0; j < dim; j++) {
            GradMeshVel_old[i][j] = 0.;
          }
          for (unsigned inode = 0; inode < nve; inode++) {
	    vx_ig[i]+=phi[inode] * vx[i][inode];
	    vxOld_ig[i]+=phi[inode] * vx_old[i][inode];
            meshVel_old[i] += phi[inode] * meshVelOldNode[i][inode];
            for (int j = 0; j < dim; j++) {
              GradMeshVel_old[i][j] += gradphi_old[inode * dim + j] * meshVelOldNode[i][inode];//TODO
            }
          }
        }
        
        // mesh velocity NEW
        vector < adept::adouble > meshVel(dim);
        vector < vector < adept::adouble > > GradMeshVel(dim);
	unsigned dBlock = ( meshIsCurrupted ) ? (nBlocks-1)*dim : 0;
        for (unsigned i = 0; i < dim; i++) {
          meshVel[i] = 2./dt * (SolVAR[dBlock + i] - SolVAR_old[dBlock + i]) - meshVel_old[i];
          GradMeshVel[i].resize(dim);
          for (unsigned j = 0; j < dim; j++) {	    
	    GradMeshVel[i][j] = 2./dt * (GradSolVAR[dBlock + i][j] - GradSolVAR_old[dBlock + i][j]) - GradMeshVel_old[i][j];
          }
        }

        // ---------------------------------------------------------------------------
        //BEGIN FLUID and POROUS MEDIA ASSEMBLY ============
        if (flag_mat == 2 || flag_mat == 3) {

          vector < adept::adouble > a(dim);
          vector < double > a_old(dim);
          for (int i = 0; i < dim; i++) {
            a[i] = SolVAR[i + dim];// - 0.*meshVel[i]; // TODO maybe we subtract meshVel[i] maybe not
            a_old[i] = SolVAR_old[i + dim];// - 0.*meshVel[i];
          }

          // speed
          adept::adouble aL2Norm = 0.;
          double a_oldL2Norm = 0.;
          for (int i = 0; i < dim; i++) {
            aL2Norm += a[i] * a[i];
            a_oldL2Norm += a_old[i] * a_old[i];
          }
          aL2Norm = sqrt(aL2Norm);
          a_oldL2Norm = sqrt(a_oldL2Norm);

          double sqrtlambdak = (*mysolution->_Sol[indLmbd])(iel);
          adept::adouble tauSupg = 1. / (sqrtlambdak * sqrtlambdak * 4.*IRe);
          double tauSupg_old = 1. / (sqrtlambdak * sqrtlambdak * 4.*IRe);
          adept::adouble Rek = aL2Norm / (4.*sqrtlambdak * IRe);
          double Rek_old = a_oldL2Norm / (4.*sqrtlambdak * IRe);

          if (Rek > 1.0e-15) {
            adept::adouble xiRek = (Rek >= 1.) ? 1. : Rek;
            tauSupg   = xiRek / (aL2Norm * sqrtlambdak);
          }

          if (Rek_old > 1.0e-15) {
            double xiRek_old = (Rek_old >= 1.) ? 1. : Rek_old;
            tauSupg_old   = xiRek_old / (a_oldL2Norm * sqrtlambdak);
          }

          //std::cout << tauSupg.value() <<"  " <<tauSupg_old.value()<<std::endl;


          vector < adept::adouble > phiSupg(nve, 0.);
          vector < double > phiSupg_old(nve, 0.);
          //tauSupg=0;
          //tauSupg_old=0;
          for (unsigned i = 0; i < nve; i++) {
            for (unsigned j = 0; j < dim; j++) {
              phiSupg[i] += ((SolVAR[j + dim] - meshVel[j]) * gradphi[i * dim + j]) * tauSupg;    // * (!solidmark[i]) ;
              phiSupg_old[i] += ((SolVAR_old[j + dim] - meshVel_old[j]) * gradphi_old[i * dim + j]) * tauSupg_old;    // * (!solidmark[i]);
            }
          }


          //BEGIN ALE + Momentum (Navier-Stokes)
          {
            for (unsigned i = 0; i < nve; i++) {

              //BEGIN redidual Laplacian ALE map in the reference domain


              if (!solidmark[i]) {
                adept::adouble LapmapVAR[3] = {0., 0., 0.};
                for (int idim = 0; idim < dim; idim++) {
                  for (int jdim = 0; jdim < dim; jdim++) {
                    //LapmapVAR[idim] += ( GradSolVAR[idim][jdim] ) * gradphi_hat[i * dim + jdim];
                    LapmapVAR[idim] += (GradSolVAR[idim][jdim] + 0.*GradSolVAR[jdim][idim]) * gradphi[i * dim + jdim];
                  }
                }

                for (int idim = 0; idim < dim; idim++) {
                  aRhs[indexVAR[idim]][i] +=  (- LapmapVAR[idim] * Weight_nojac);
                }
              }

              //END residual Laplacian ALE map in reference domain

              if (flag_mat == 2) {
                //BEGIN residual Navier-Stokes in moving domain
                adept::adouble Lapdisp[3] = {0., 0., 0.};
		adept::adouble Lapdisp_old[3] = {0., 0., 0.};
		adept::adouble LapvelVAR[3] = {0., 0., 0.};
	        double LapvelVAR_old[3] = {0., 0., 0.};
                adept::adouble LapStrong[3] = {0., 0., 0.};
                double LapStrong_old[3] = {0., 0., 0.};
                adept::adouble AdvaleVAR[3] = {0., 0., 0.};
                double AdvaleVAR_old[3] = {0., 0., 0.};

                for (int idim = 0.; idim < dim; idim++) {
                  for (int jdim = 0.; jdim < dim; jdim++) {
                    unsigned kdim;
                    if (idim == jdim) kdim = jdim;
                    else if (1 == idim + jdim) kdim = dim;    // xy
                    else if (2 == idim + jdim) kdim = dim + 2;   // xz
                    else if (3 == idim + jdim) kdim = dim + 1;   // yz
                    //laplaciano debole
                    LapvelVAR[idim]     += (GradSolVAR[dim + idim][jdim] + GradSolVAR[dim + jdim][idim]) * gradphi[i * dim + jdim];
                    LapvelVAR_old[idim] += (GradSolVAR_old[dim + idim][jdim] + GradSolVAR_old[dim + jdim][idim]) * gradphi_old[i * dim + jdim];
		    
		    Lapdisp[idim]     += (GradSolVAR[idim][jdim] + GradSolVAR[jdim][idim]) * gradphi[i * dim + jdim];
                    Lapdisp_old[idim] += (GradSolVAR_old[idim][jdim] + GradSolVAR_old[jdim][idim]) * gradphi_old[i * dim + jdim];
		    
                    //laplaciano strong
                    LapStrong[idim]     += (NablaSolVAR[dim + idim][jdim] + NablaSolVAR[dim + jdim][kdim]) * phiSupg[i];
                    LapStrong_old[idim] += (NablaSolVAR_old[dim + idim][jdim] + NablaSolVAR_old[dim + jdim][kdim]) * phiSupg_old[i];

                    AdvaleVAR[idim]	+= ((SolVAR[dim + jdim] - meshVel[jdim]) * GradSolVAR[dim + idim][jdim]
					    + (GradSolVAR[dim + jdim][jdim] - GradMeshVel[jdim][jdim]) * SolVAR[dim + idim]
					   ) * (phi[i] + phiSupg[i]);

                    AdvaleVAR_old[idim]	+= ((SolVAR_old[dim + jdim] - meshVel_old[jdim]) * GradSolVAR_old[dim + idim][jdim]
                                            + (GradSolVAR_old[dim + jdim][jdim] - GradMeshVel_old[jdim][jdim]) * SolVAR_old[dim + idim]
                                           ) * (phi_old[i] + phiSupg_old[i]);
                  }
                }

                for (int idim = 0; idim < dim; idim++) {

                  adept::adouble timeDerivative = - (SolVAR[dim + idim] * (phi[i] + phiSupg[i]) * Weight
                                                     - SolVAR_old[dim + idim] * (phi_old[i] + phiSupg_old[i]) * Weight_old) / dt;

                  adept::adouble value =  theta * (
                                            - AdvaleVAR[idim]      	             // advection term
                                            - IRe * LapvelVAR[idim]	             // viscous dissipation
                                            - (idim == 0) * Lapdisp[idim] * 1.e-3 * exp( (vx_ig[0] + 2.5e-5)/0.0001)*(dim==2)
					    - (idim == 0) * Lapdisp[idim] * 1.e-3 * exp(-(vx_ig[0] - 1e-4)/0.00004)*(dim==3)
                                            + IRe * LapStrong[idim]
                                            + 1. / rhof * SolVAR[nBlocks * dim] * gradphi[i * dim + idim] // pressure gradient
                                          ) * Weight;                                // at time t

                  adept::adouble value_old = (1. - theta) * (
                                               - AdvaleVAR_old[idim]               	         // advection term
                                               - IRe * LapvelVAR_old[idim]	       	         // viscous dissipation
                                               - (idim == 0) * Lapdisp_old[idim] * 1.e-3 * exp( (vxOld_ig[0] + 2.5e-5)/0.0001)*(dim==2)
					       - (idim == 0) * Lapdisp_old[idim] * 1.e-3 * exp(-(vxOld_ig[0] - 1.e-4)/0.0004)*(dim==3)
                                               + IRe * LapStrong_old[idim]
                                               + 1. / rhof * SolVAR[nBlocks * dim] * gradphi_old[i * dim + idim]  // pressure gradient
                                             ) * Weight_old;			                 // at time t-dt

                  if (!solidmark[i]) {
                    aRhs[indexVAR[dim + idim]][i] += timeDerivative + value + value_old;
                  }
                  else {
                    aRhs[indexVAR[idim]][i] += timeDerivative + value + value_old;
                  }

                }
                //END redidual Navier-Stokes in moving domain
              }
              else if (flag_mat == 3) {
                //BEGIN redidual Porous Media in moving domain


                adept::adouble speed = 0.;
                adept::adouble speed_old = 0.;
                for (int idim = 0.; idim < dim; idim++) {
                  speed += SolVAR[dim + idim] * SolVAR[dim + idim];
                  speed_old += SolVAR_old[dim + idim] * SolVAR_old[dim + idim];
                }
                double eps = 1.0e-12;
                speed = sqrt(speed + eps);
                speed_old = sqrt(speed_old + eps);
                double DE = 0.;
                if (dim == 2) {
                  DE = 0.0002; // AAA_thrombus_2D
                  //DE = 0.00006; // turek2D
                }
                else if (dim == 3) {
                  DE = 0.000112; // porous3D
                }
                double b = 4188;
                double a = 1452;
                double K = DE * IRe * rhof / b; // alpha = mu/b * De
                double C2 = 2 * a / (rhof * DE);





// 		adept::adouble LapvelVAR[3] = {0., 0., 0.};
// 		adept::adouble LapvelVAR_old[3] = {0., 0., 0.};
// 		adept::adouble AdvaleVAR[3] = {0., 0., 0.};
// 		adept::adouble AdvaleVAR_old[3] = {0., 0., 0.};
//
// 		for(int idim = 0.; idim < dim; idim++) {
// 		  for(int jdim = 0.; jdim < dim; jdim++) {
//
// 		    LapvelVAR[idim]     += (GradSolVAR[dim + idim][jdim] + GradSolVAR[dim + jdim][idim]) * gradphi[i * dim + jdim];
// 		    LapvelVAR_old[idim] += (GradSolVAR_old[dim + idim][jdim] + GradSolVAR_old[dim + jdim][idim]) * gradphi_old[i * dim + jdim];
//
// 		    AdvaleVAR[idim]	+= ((SolVAR[dim + jdim] - meshVel[jdim]) * GradSolVAR[dim + idim][jdim]
// 					+ (GradSolVAR[dim + jdim][jdim] - GradMeshVel[jdim][jdim]) * SolVAR[dim + idim]
// 				      ) * phi[i];
//
// 		    AdvaleVAR_old[idim]	+= ((SolVAR_old[dim + jdim] - meshVel[jdim]) * GradSolVAR_old[dim + idim][jdim]
// 					    + (GradSolVAR_old[dim + jdim][jdim] - GradMeshVel[jdim][jdim]) * SolVAR_old[dim + idim]
// 					  ) * phi_old[i];
// 		  }
// 		}

                for (int idim = 0; idim < dim; idim++) {

// 		  adept::adouble timeDerivative = -(SolVAR[dim + idim] * phi[i] * Weight
//                                                   - SolVAR_old[dim + idim] * phi_old[i] * Weight_old);

                  adept::adouble value =  theta * (
// 					    - AdvaleVAR[idim]      	             // advection term
// 					    - IRe * LapvelVAR[idim]	             // viscous dissipation
                                            //- SolVAR[dim + idim] * ( IRe / K + 0.5 * C2 * speed ) * phi[i]
                                            - (SolVAR[dim + idim] - meshVel[idim]) * (IRe / K + 0.5 * C2 * speed) * (phi[i] + phiSupg[i])
                                            + 1. / rhof * SolVAR[nBlocks * dim] * gradphi[i * dim + idim] // pressure gradient
                                          ) * Weight;                                // at time t

                  adept::adouble value_old = (1. - theta) * (
// 					      - AdvaleVAR_old[idim]               	         // advection term
// 					      - IRe * LapvelVAR_old[idim]	       	         // viscous dissipation
                                               //- SolVAR_old[dim + idim] * ( IRe / K + 0.5 * C2 * speed_old ) * phi_old[i]
                                               - (SolVAR_old[dim + idim] - meshVel_old[idim]) * (IRe / K + 0.5 * C2 * speed) * (phi_old[i] + phiSupg_old[i])
                                               + 1. / rhof * SolVAR[nBlocks * dim] * gradphi_old[i * dim + idim]  // pressure gradient
                                             ) * Weight_old;			                 // at time t-dt

                  if (!solidmark[i]) {
                    aRhs[indexVAR[dim + idim]][i] += value + value_old;
                  }
                  else {
                    aRhs[indexVAR[idim]][i] += value + value_old;
                  }

                }
                //END redidual Porous Media in moving domain
              }

              if (auxDisp) {
                //BEGIN auxiliary displacement equation
                adept::adouble LapAuxVAR[3] = {0., 0., 0.};

                for (int idim = 0; idim < dim; idim++) {
                  for (int jdim = 0; jdim < dim; jdim++) {
                    LapAuxVAR[idim] += ( GradSolVAR[idim + 2 * dim][jdim] +  0. * GradSolVAR[jdim + 2 * dim][idim] ) * gradphi[i * dim + jdim];
                  }
                }

                for (int idim = 0; idim < dim; idim++) {
                  aRhs[indexVAR[idim + 2 * dim]][i] += (  meshIsCurrupted * 1.0e-1 * LapAuxVAR[idim] +  0 * ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi[i] ) * Weight_nojac;
                }
//                 for (int idim = 1; idim < dim; idim++) {
//                   aRhs[indexVAR[idim + 2 * dim]][i] += (  meshIsCurrupted * 1.0e-1 * LapAuxVAR[idim] +  0 * ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi[i] ) * Weight_nojac;
//                 }
//                 for (int idim = 1; idim < dim; idim++) {
//                   aRhs[indexVAR[idim + 2 * dim]][i] += (  meshIsCurrupted * 1.0e-1 * LapAuxVAR[idim] +  0 * ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi[i] ) * Weight_nojac;
//                 }

//                 adept::adouble LapAuxVAR[3] = {0., 0., 0.};
//
//                 for (int idim = 0; idim < dim; idim++) {
//                   for (int jdim = 0; jdim < dim; jdim++) {
//                      LapAuxVAR[idim] += GradSolhatVAR[idim + 2 * dim][jdim] * gradphi_hat[i * dim + jdim];
//                   }
//                 }
//
//                 for (int idim = 0; idim < dim; idim++) {
//                   aRhs[indexVAR[idim + 2 * dim]][i] += (  1.0e-10 * LapAuxVAR[idim] +  ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi_hat[i] ) * Weight_hat;
//                 }


                //END auxiliary displacement equation
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
              aRhs[indexVAR[nBlocks * dim]][i] += - (-phi1[i] * div_vel) * Weight;
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
          adept::adouble J_hat1;
          double J_hat_old;
          adept::adouble I_e;
          double I_e_old;
          adept::adouble Cauchy[3][3];
          adept::adouble Cauchy_old[3][3];
          double Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};




          if (solid_model == 0) {   // Saint-Venant
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
                Cauchy[i][j]     = 2 * mus * e[i][j] - 2 * mus * I_e * SolVAR[nBlocks * dim] * Id2th[i][j];
                Cauchy_old[i][j] = 2 * mus * e_old[i][j] - 2 * mus * I_e_old * SolVAR[nBlocks * dim] * Id2th[i][j];
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
            adept::adouble F1[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
            adept::adouble B[3][3];

            double F_old[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
            double B_old[3][3];

            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                F[i][j] += GradSolhatVAR[i][j];
                F1[i][j] += GradSolhatVAR[i + 2 * dim][j];
                F_old[i][j] += GradSolhatVAR_old[i][j];
              }
            }

            J_hat =   F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                      - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

            J_hat1 =   F1[0][0] * F1[1][1] * F1[2][2] + F1[0][1] * F1[1][2] * F1[2][0] + F1[0][2] * F1[1][0] * F1[2][1]
                       - F1[2][0] * F1[1][1] * F1[0][2] - F1[2][1] * F1[1][2] * F1[0][0] - F1[2][2] * F1[1][0] * F1[0][1];

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

            if (solid_model <= 4) {   // Neo-Hookean
              I1_B = B[0][0] + B[1][1] + B[2][2];
              I1_B_old = B_old[0][0] + B_old[1][1] + B_old[2][2];

              for (int I = 0; I < 3; ++I) {
                for (int J = 0; J < 3; ++J) {
                  if (1 == solid_model) {   //Wood-Bonet J_hat  =1;
                    Cauchy[I][J]     = mus * (B[I][J]     - Id2th[I][J]) - mus / 3.*I1_B     * SolVAR[nBlocks * dim] * Id2th[I][J];
                    Cauchy_old[I][J] = mus * (B_old[I][J] - Id2th[I][J]) - mus / 3.*I1_B_old * SolVAR[nBlocks * dim] * Id2th[I][J];
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
            else if (5 == solid_model) {   //Mooney-Rivlin
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
              invB[1][0] = - (B[0][1] * B[2][2] - B[0][2] * B[2][1]) * invdetB;
              invB[2][0] = (B[0][1] * B[1][2] - B[0][2] * B[1][1]) * invdetB;
              invB[0][1] = - (B[1][0] * B[2][2] - B[1][2] * B[2][0]) * invdetB;
              invB[1][1] = (B[0][0] * B[2][2] - B[0][2] * B[2][0]) * invdetB;
              invB[2][1] = - (B[0][0] * B[1][2] - B[1][0] * B[0][2]) * invdetB;
              invB[0][2] = (B[1][0] * B[2][1] - B[2][0] * B[1][1]) * invdetB;
              invB[1][2] = - (B[0][0] * B[2][1] - B[2][0] * B[0][1]) * invdetB;
              invB[2][2] = (B[0][0] * B[1][1] - B[1][0] * B[0][1]) * invdetB;

              invB_old[0][0] = (B_old[1][1] * B_old[2][2] - B_old[1][2] * B_old[2][1]) * invdetB_old;
              invB_old[1][0] = - (B_old[0][1] * B_old[2][2] - B_old[0][2] * B_old[2][1]) * invdetB_old;
              invB_old[2][0] = (B_old[0][1] * B_old[1][2] - B_old[0][2] * B_old[1][1]) * invdetB_old;
              invB_old[0][1] = - (B_old[1][0] * B_old[2][2] - B_old[1][2] * B_old[2][0]) * invdetB_old;
              invB_old[1][1] = (B_old[0][0] * B_old[2][2] - B_old[0][2] * B_old[2][0]) * invdetB_old;
              invB_old[2][1] = - (B_old[0][0] * B_old[1][2] - B_old[1][0] * B_old[0][2]) * invdetB_old;
              invB_old[0][2] = (B_old[1][0] * B_old[2][1] - B_old[2][0] * B_old[1][1]) * invdetB_old;
              invB_old[1][2] = - (B_old[0][0] * B_old[2][1] - B_old[2][0] * B_old[0][1]) * invdetB_old;
              invB_old[2][2] = (B_old[0][0] * B_old[1][1] - B_old[1][0] * B_old[0][1]) * invdetB_old;


              I1_B = B[0][0] + B[1][1] + B[2][2];
              I2_B = B[0][0] * B[1][1] + B[1][1] * B[2][2] + B[2][2] * B[0][0]
                     - B[0][1] * B[1][0] - B[1][2] * B[2][1] - B[2][0] * B[0][2];

              I1_B_old = B_old[0][0] + B_old[1][1] + B_old[2][2];
              I2_B_old = B_old[0][0] * B_old[1][1] + B_old[1][1] * B_old[2][2] + B_old[2][2] * B_old[0][0]
                         - B_old[0][1] * B_old[1][0] - B_old[1][2] * B_old[2][1] - B_old[2][0] * B_old[0][2];

              //double C1 = mus / 3.;
              double C1 = (elementGroup == 15) ? mus1 / 3. : mus / 3.;
              double C2 = C1 / 2.;

              for (int I = 0; I < 3; ++I) {
                for (int J = 0; J < 3; ++J) {
                  Cauchy[I][J] =  2.* (C1 * B[I][J] - C2 * invB[I][J])
                                  //- (2. / 3.) * (C1 * I1_B - C2 * I2_B) * SolVAR[nBlocks * dim] * Id2th[I][J];
                                  - 1. / rhof * SolVAR[nBlocks * dim] * Id2th[I][J];

                  Cauchy_old[I][J] =  2.* (C1 * B_old[I][J] - C2 * invB_old[I][J])
                                      //- (2. / 3.) * (C1 * I1_B_old - C2 * I2_B_old) * SolVAR[nBlocks * dim] * Id2th[I][J];
                                      - 1. / rhof * SolVAR[nBlocks * dim] * Id2th[I][J];

                }
              }

            }
          }

          //END build Cauchy Stress in moving domain

          //BEGIN v=0 + Momentum (Solid)
          {
            for (unsigned i = 0; i < nve; i++) {

              //BEGIN redidual d_t - v = 0 in fixed domain
              for (int idim = 0; idim < dim; idim++) {
                aRhs[indexVAR[dim + idim]][i] +=  - phi[i] * ((-SolVAR[idim] + SolVAR_old[idim]) / dt +
                                                  (theta * SolVAR[dim + idim] + (1. - theta) * SolVAR_old[dim + idim])
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

                adept::adouble timeDerivative = - (rhos * SolVAR[dim + idim] * phi[i] * Weight
                                                   - rhos * SolVAR_old[dim + idim] * phi_old[i] * Weight_old) / dt;

                adept::adouble value =  theta * (rhos * phi[i] * _gravity[idim]      // body force
                                                 - CauchyDIR[idim]			  // stress
                                                ) * Weight;                         // at time t

                adept::adouble value_old = (1. - theta) * (rhos * phi_old[i] * _gravity[idim]      // body force
                                           - CauchyDIR_old[idim]			 // stress
                                                          ) * Weight_old;                         // at time t-dt

                aRhs[indexVAR[idim]][i] += timeDerivative + value + value_old;
              }

              //END redidual Solid Momentum in moving domain


              if (auxDisp) {
                //BEGIN auxiliary displacement equation
                adept::adouble LapAuxVAR[3] = {0., 0., 0.};

                for (int idim = 0; idim < dim; idim++) {
                  for (int jdim = 0; jdim < dim; jdim++) {
                    LapAuxVAR[idim] += ( GradSolVAR[idim + 2 * dim][jdim] + 0 * GradSolVAR[jdim + 2 * dim][idim] ) * gradphi[i * dim + jdim];
                  }
                }

//                 for (int idim = 0; idim < dim; idim++) {
//                   aRhs[indexVAR[idim + 2 * dim]][i] += (  1.0e-8 * LapAuxVAR[idim] +  ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi_hat[i] ) * Weight;
//                 }

                for (int idim = 0; idim < dim; idim++) {
                  aRhs[indexVAR[idim + 2 * dim]][i] += (  meshIsCurrupted * LapAuxVAR[idim] +  factordxi[idim] * ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi[i] ) * Weight;
                }
//                 for (int idim = 1; idim < dim; idim++) {
//                   aRhs[indexVAR[idim + 2 * dim]][i] += (  meshIsCurrupted * LapAuxVAR[idim] +  factordy * ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi[i] ) * Weight;
//                 }




//                 adept::adouble LapAuxVAR[3] = {0., 0., 0.};
//
//                 for (int idim = 0; idim < dim; idim++) {
//                   for (int jdim = 0; jdim < dim; jdim++) {
//                      LapAuxVAR[idim] += GradSolhatVAR[idim + 2 * dim][jdim] * gradphi_hat[i * dim + jdim];
//                   }
//                 }
//
//                 for (int idim = 0; idim < dim; idim++) {
//                   aRhs[indexVAR[idim + 2 * dim]][i] += (  1.0e-8 * LapAuxVAR[idim] +  ( SolVAR[idim + 2 * dim] - SolVAR[idim] ) * phi_hat[i] ) * Weight_hat;
//                 }
                //END auxiliary displacement equation
              }

            }
          }
          //END v=0 + Momentum (Solid)

          //BEGIN continuity block
          {
            for (unsigned i = 0; i < nve1; i++) {
              if (!penalty) {
                if (0 == solid_model) {
                  aRhs[indexVAR[nBlocks * dim]][i] += - (-phi1[i] * (I_e + (!incompressible) / lambda * SolVAR[nBlocks * dim])) * Weight_hat;
                }
                else if (1 == solid_model || 5 == solid_model) {
                  aRhs[indexVAR[nBlocks * dim]][i] += phi1[i] * ( J_hat  -  1. +
                                                      (!incompressible) / lambda * SolVAR[nBlocks * dim]) * Weight_hat;
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
      for (unsigned i = 0; i < nBlocks * dim; i++) {
        Rhs[indexVAR[i]].resize(nve);

        for (int j = 0; j < nve; j++) {
          Rhs[indexVAR[i]][j] = -aRhs[indexVAR[i]][j].value();
        }
      }

      Rhs[indexVAR[nBlocks * dim]].resize(nve1);

      for (unsigned j = 0; j < nve1; j++) {
        Rhs[indexVAR[nBlocks * dim]][j] = -aRhs[indexVAR[nBlocks * dim]][j].value();
      }

      for (int i = 0; i < nBlocks * dim + 1; i++) {
        myRES->add_vector_blocked(Rhs[indexVAR[i]], dofsVAR[i]);
      }

      if (assembleMatrix) {
        //Store equations
        for (int i = 0; i < nBlocks * dim; i++) {
          s.dependent(&aRhs[indexVAR[i]][0], nve);
          s.independent(&Soli[indexVAR[i]][0], nve);
        }

        s.dependent(&aRhs[indexVAR[nBlocks * dim]][0], nve1);
        s.independent(&Soli[indexVAR[nBlocks * dim]][0], nve1);

        Jac.resize((nBlocks * dim * nve + nve1) * (nBlocks * dim * nve + nve1));

        s.jacobian(&Jac[0], true);

        myKK->add_matrix_blocked(Jac, dofsAll, dofsAll);
        s.clear_independents();
        s.clear_dependents();

        //END local to global assembly
      }
    } //end list of elements loop

    if (assembleMatrix) myKK->close();

    myRES->close();

    setIfCorrupted->close();
    double setIfCorruptedNorm = setIfCorrupted->l1_norm();

//     std::cout << "I am in Assembly and I belived the mesh is ";
//
//     if (!meshIsCurrupted) {
//       std::cout << " not corrupted";
//       if (setIfCorruptedNorm > 0) {
//         meshIsCurrupted = true;
//         std::cout << ", but I am wrong!!!!!!!!!!!!!!!!!!!!!!!!";
// 	if (assembleMatrix) myKK->zero();
// 	myRES->zero();
// 	mysolution->ResetSolutionToOldSolution();
//         goto begin;
//       }
//       std::cout << std::endl;
//     }
//     else {
//       std::cout << " corrupted \n";
//       if (setIfCorruptedNorm > 0) {
// 	std::cout << ", but I am wrong!!!!!!!!!!!!!!!!!!!!!!!!";
//         //play with the factor TODO
//       }
//     }

    meshIsCurrupted = true;

    delete setIfCorrupted;
    delete area_elem_first;

    // *************************************
    end_time = clock();
    AssemblyTime += (end_time - start_time);
    // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************



  }
  
  void StoreMeshVelocity(MultiLevelProblem & ml_prob){
  
  TransientNonlinearImplicitSystem & my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("Fluid-Structure-Interaction");
  const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();
  
  MultiLevelSolution *ml_sol = ml_prob._ml_sol;
  
  Solution *solution = ml_sol->GetSolutionLevel(level);
  Mesh *msh = ml_prob._ml_msh->GetLevel(level);
  
  //const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();
  const unsigned dim = msh->GetDimension();
  const char varname[9][4] = {"DX", "U", "Um", "DY", "V", "Vm", "DZ", "W", "Wm"};
  const char varname1[9][4] = {"DX1", "U", "Um", "DY1", "V", "Vm", "DZ1", "W", "Wm"};

  vector <unsigned> indVAR(3 * dim);
  
  for (unsigned ivar = 0; ivar < 3 * dim; ivar++) {
    if(!meshIsCurrupted){
      indVAR[ivar] = ml_sol->GetIndex(&varname[ivar][0]);
    }
    else{
      indVAR[ivar] = ml_sol->GetIndex(&varname1[ivar][0]);
    }
  }
  
  double dt =  my_nnlin_impl_sys.GetIntervalTime();
  
  
  int  iproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    
  for(unsigned idim = 0; idim < dim; idim++) {
    unsigned ivar = idim * 3; 
    for(unsigned jdof = msh->_dofOffset[2][iproc]; jdof < msh->_dofOffset[2][iproc + 1]; jdof++) {
      bool solidmark = msh->GetSolidMark(jdof);
      double vnew;
      if(solidmark != solidmark){
	vnew = (*solution->_Sol[indVAR[ivar + 1]])(jdof); // solid: mesh velocity equals solid velocity
      }
      else{
	double unew = (*solution->_Sol[indVAR[ivar]])(jdof);
	double uold = (*solution->_SolOld[indVAR[ivar]])(jdof);
	double vold = (*solution->_SolOld[indVAR[ivar + 2]])(jdof);
	vnew = 2./dt*(unew-uold) - vold;
	
      }
      solution->_Sol[indVAR[ivar + 2]]->set(jdof,vnew);
    }
    solution->_Sol[indVAR[ivar + 2]]->close();
  }
  
  }
  
}




