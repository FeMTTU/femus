#ifndef __femus_include_IncompressibleFSIAssemblySupg_hpp__
#define __femus_include_IncompressibleFSIAssemblySupg_hpp__
#endif

#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "adept.h"

#include "OprtrTypeEnum.hpp"

double Gamma = 0.5;
double beta = 0.25;

namespace femus
{

  double sm[3][4] = {
    {0.5, 1.},
    {0.211324865405, 0.711324865405, 1.},
    {0.112701665379, 0.5, 0.887298334621, 1.}
  };

  double thetam[3][4] = {
    { 1. , 0.},
    { 0.5, 0.5, 0.},
    { 0.277777777778, 0.444444444444, 0.277777777778, 0.}
  };

  unsigned TimeIntPoints = 3;

  double* s = sm[TimeIntPoints - 1];
  double* theta = thetam[TimeIntPoints - 1];

  bool twoPressure = true;

  bool valve = true;

  double Kslip = 0.;

  void FSITimeDependentAssemblySupgNew2 ( MultiLevelProblem& ml_prob )
  {

    clock_t AdeptTime = 0;
    clock_t start_time;

    unsigned nBlocks = 2;
    unsigned nP = ( twoPressure ) ? 2 : 1;
    unsigned jP = ( twoPressure ) ? 1 : 0;

    //pointers and references
    TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ( "Fluid-Structure-Interaction" );
    const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();

    MultiLevelSolution* ml_sol = ml_prob._ml_sol;

    Solution* mysolution = ml_sol->GetSolutionLevel ( level );

    LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];
    Mesh* mymsh =  ml_prob._ml_msh->GetLevel ( level );
    elem* myel =  mymsh->el;
    SparseMatrix* myKK = myLinEqSolver->_KK;
    NumericVector* myRES = myLinEqSolver->_RES;


    bool assembleMatrix = my_nnlin_impl_sys.GetAssembleMatrix();
    // call the adept stack object
    adept::Stack& stack = FemusInit::_adeptStack;

    if ( assembleMatrix ) stack.continue_recording();
    else stack.pause_recording();

    const unsigned dim = mymsh->GetDimension();
    const unsigned nabla_dim = 3 * ( dim - 1 );
    const unsigned max_size = static_cast< unsigned > ( ceil ( pow ( 3, dim ) ) );

    // local objects
    vector<adept::adouble> SolVAR ( nBlocks * dim + nP );
    vector<double> SolVAROld ( nBlocks * dim );
    vector<adept::adouble> SolVARNew ( nBlocks * dim );

    vector < double > vxOld_ig ( dim );
    vector < adept::adouble > vxNew_ig ( dim );
    vector < adept::adouble > vx_ig ( dim );

    vector<vector < adept::adouble > > GradSolVAR ( nBlocks * dim );
    vector<vector < adept::adouble > > GradSolhatVAR ( dim );
    vector<vector<adept::adouble> > NablaSolVAR ( nBlocks * dim );

    vector<double> meshVelOld ( dim );
    vector<double> VelOldest ( dim );
    vector < adept::adouble > meshVel ( dim );

    for ( int i = 0; i < nBlocks * dim; i++ ) {
      GradSolVAR[i].resize ( dim );
      if ( i < dim ) GradSolhatVAR[i].resize ( dim );
      NablaSolVAR[i].resize ( nabla_dim );
    }

    vector < bool> solidmark;
    vector < double > phi;
    vector < double > phi_hat;
    vector < adept::adouble> gradphi;
    vector < double > gradphi_hat;
    vector < adept::adouble> nablaphi;
    vector < double > nablaphi_hat;

    //vector < double > *nullDoublePointer = NULL;
    //vector < adept::adouble > *nullAdoublePointer = NULL;

    const double* phi1;
    adept::adouble Weight = 0.;
    double Weight_hat = 0.;

    vector <vector < double> > vx_hat ( dim );
    vector <vector < double> > vxOld ( dim );
    vector <vector < adept::adouble> > vxNew ( dim );
    vector <vector < adept::adouble> > vx ( dim );
    vector <vector < adept::adouble > > vx_face ( dim );

    vector< vector< adept::adouble > > Soli ( nBlocks * dim + nP );
    vector< vector< double > > Soli_old ( nBlocks * dim );
    vector< vector< double > > aold ( dim );
    vector< vector< int > > dofsVAR ( nBlocks * dim + nP );
    vector< vector< double > > meshVelOldNode ( dim );

    vector< vector< double > > Rhs ( nBlocks * dim + nP );
    vector< vector< adept::adouble > > aRhs ( nBlocks * dim + nP );

    vector < int > dofsAll;
    vector < double > Jac;

    // ------------------------------------------------------------------------
    // Physical parameters
    double rhof	 	= ml_prob.parameters.get<Fluid> ( "Fluid" ).get_density();
    double rhos	 	= ml_prob.parameters.get<Fluid> ( "Solid" ).get_density() / rhof;
    double mu_lame 	= ml_prob.parameters.get<Solid> ( "Solid" ).get_lame_shear_modulus();
    double lambda_lame 	= ml_prob.parameters.get<Solid> ( "Solid" ).get_lame_lambda();
    double lambda_lame1	= ml_prob.parameters.get<Solid> ( "Solid1" ).get_lame_lambda();
    double mus		= mu_lame / rhof;
    double mu_lame1 	= ml_prob.parameters.get < Solid> ( "Solid1" ).get_lame_shear_modulus();
    double mus1 	= mu_lame1 / rhof;
    double IRe 		= ml_prob.parameters.get<Fluid> ( "Fluid" ).get_IReynolds_number();
    double lambda	= lambda_lame / rhof;
    double lambda1	= lambda_lame1 / rhof;
    double betans	= 1.;
    int    solid_model	= ml_prob.parameters.get<Solid> ( "Solid" ).get_physical_model();

    double bulkModule = lambda + 2./3. * mus;
    double bulkModule1 = lambda1 + 2./3. * mus1;
    
    bool incompressible = (0.5  ==  ml_prob.parameters.get < Solid>("Solid").get_poisson_coeff()) ? 1 : 0;
    bool incompressible1 = (0.5  ==  ml_prob.parameters.get < Solid>("Solid1").get_poisson_coeff()) ? 1 : 0;
    
    if ( solid_model >= 2 && solid_model <= 4 ) {
      std::cout << "Error! Solid Model " << solid_model << "not implemented\n";
      abort();
    }

//     bool incompressible = ( 0.5 == ml_prob.parameters.get<Solid> ( "Solid" ).get_poisson_coeff() ) ? 1 : 0;
//     bool incompressible1 = (0.5  ==  ml_prob.parameters.get < Solid>("Solid1").get_poisson_coeff()) ? 1 : 0;
    //const bool penalty = ml_prob.parameters.get<Solid> ( "Solid" ).get_if_penalty();

    // gravity
    double _gravity[3] = {0., 0., 0.};

    double dt =  my_nnlin_impl_sys.GetIntervalTime();
    double time =  my_nnlin_impl_sys.GetTime();
    // -----------------------------------------------------------------
    // space discretization parameters
    unsigned SolType2 = ml_sol->GetSolutionType ( ml_sol->GetIndex ( "U" ) );
    unsigned SolType1 = ml_sol->GetSolutionType ( ml_sol->GetIndex ( "PS" ) );

    // mesh and procs
    unsigned nel    = mymsh->GetNumberOfElements();
    unsigned igrid  = mymsh->GetLevel();
    unsigned iproc  = mymsh->processor_id();

    unsigned indLmbd = ml_sol->GetIndex ( "lmbd" );

    //----------------------------------------------------------------------------------
    //variable-name handling
    const char varname[8][3] = {"DX", "DY", "DZ", "U", "V", "W", "PS", "PF"};
    const char varname2[10][4] = {"Um", "Vm", "Wm"};
    const char varname3[3][4] = {"AX", "AY", "AZ"};

    vector <unsigned> indexVAR ( nBlocks * dim + nP );
    vector <unsigned> indVAR ( nBlocks * dim + nP );
    vector <unsigned> SolType ( nBlocks * dim + nP );

    vector <unsigned> indVAR2 ( dim );
    vector <unsigned> indVAR3 ( dim );

    for ( unsigned ivar = 0; ivar < dim; ivar++ ) {
      for ( unsigned k = 0; k < nBlocks; k++ ) {
        indVAR[ivar + k * dim] = ml_sol->GetIndex ( &varname[ivar + k * 3][0] );
        SolType[ivar + k * dim] = ml_sol->GetSolutionType ( &varname[ivar + k * 3][0] );
        indexVAR[ivar + k * dim] = my_nnlin_impl_sys.GetSolPdeIndex ( &varname[ivar + k * 3][0] );
      }

      indVAR2[ivar] = ml_sol->GetIndex ( &varname2[ivar][0] );
      indVAR3[ivar] = ml_sol->GetIndex ( &varname3[ivar][0] );
    }

    for ( unsigned ivar = 0; ivar < nP; ivar++ ) {
      indexVAR[nBlocks * dim + ivar] = my_nnlin_impl_sys.GetSolPdeIndex ( &varname[6 + ivar][0] );
      indVAR[nBlocks * dim + ivar] = ml_sol->GetIndex ( &varname[6 + ivar][0] );
      SolType[nBlocks * dim + ivar] = ml_sol->GetSolutionType ( &varname[6 + ivar][0] );
    }

    //----------------------------------------------------------------------------------

    if ( assembleMatrix ) myKK->zero();


    int nprocs = mymsh->n_processors();
    NumericVector * area_elem_first;
    area_elem_first = NumericVector::build().release();
    if ( nprocs == 1 ) {
      area_elem_first->init ( nprocs, 1, false, SERIAL );
    }
    else {
      area_elem_first->init ( nprocs, 1, false, PARALLEL );
    }
    area_elem_first->zero();
    double rapresentative_area = 1.;

    // *** element loop ***
    for ( int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++ ) {

      short unsigned ielt = mymsh->GetElementType ( iel );
      unsigned nve        = mymsh->GetElementDofNumber ( iel, SolType2 );
      unsigned nve1       = mymsh->GetElementDofNumber ( iel, SolType1 );
      int flag_mat        = mymsh->GetElementMaterial ( iel );
      unsigned elementGroup = mymsh->GetElementGroup ( iel );

      // *******************************************************************************************************

      // displacement, velocity, mesh dof values

      solidmark.resize ( nve );
      for ( int i = 0; i < nBlocks * dim; i++ ) {
        dofsVAR[i].resize ( nve );
        Soli[indexVAR[i]].resize ( nve );
        Soli_old[indexVAR[i]].resize ( nve );
        aRhs[indexVAR[i]].resize ( nve );
      }

      for ( int i = 0; i < dim; i++ ) {
        meshVelOldNode[i].resize ( nve );
        vx_hat[i].resize ( nve );
        vx[i].resize ( nve );
        vxOld[i].resize ( nve );
        vxNew[i].resize ( nve );
	aold[i].resize ( nve );
      }

      for ( unsigned i = 0; i < nve; i++ ) {
        unsigned idof = mymsh->GetSolutionDof ( i, iel, SolType2 );
        solidmark[i] = mymsh->GetSolidMark ( idof ); // flag to know if the node "idof" lays on the fluid-solid interface
        for ( int j = 0; j < dim; j++ ) {
          for ( unsigned k = 0; k < nBlocks; k++ ) {
            dofsVAR[j + k * dim][i] = myLinEqSolver->GetSystemDof ( indVAR[j + k * dim], indexVAR[j + k * dim], i, iel );
            Soli[indexVAR[j + k * dim]][i] = ( *mysolution->_Sol[indVAR[j + k * dim]] ) ( idof );
            Soli_old[indexVAR[j + k * dim]][i] = ( *mysolution->_SolOld[indVAR[j + k * dim]] ) ( idof );
            aRhs[indexVAR[j + k * dim]][i] = 0.;
          }
          meshVelOldNode[j][i] = ( *mysolution->_Sol[indVAR2[j]] ) ( idof );
          vx_hat[j][i] = ( *mymsh->_topology->_Sol[j] ) ( idof );
	  aold[j][i] = ( *mysolution->_SolOld[indVAR3[j]] ) ( idof );
        }
      }

      // pressure dof values
      for ( int i = 0; i < nP; i++ ) {
        dofsVAR[nBlocks * dim + i].resize ( nve1 );
        Soli[indexVAR[nBlocks * dim + i]].resize ( nve1 );
        aRhs[indexVAR[nBlocks * dim + i]].resize ( nve1 );
      }
      for ( unsigned i = 0; i < nve1; i++ ) {
        unsigned idof = mymsh->GetSolutionDof ( i, iel, SolType[nBlocks * dim ] );
        for ( int j = nBlocks * dim; j < nBlocks * dim + nP; j++ ) {
          dofsVAR[j][i] = myLinEqSolver->GetSystemDof ( indVAR[j], indexVAR[j], i, iel );
          Soli[indexVAR[j]][i] = ( *mysolution->_Sol[indVAR[j]] ) ( idof );
          aRhs[indexVAR[j]][i] = 0.;
        }
      }

      // stuck all dofs
      dofsAll.resize ( 0 );
      for ( int i = 0; i < nBlocks * dim + nP; i++ ) {
        dofsAll.insert ( dofsAll.end(), dofsVAR[i].begin(), dofsVAR[i].end() );
      }

      if ( assembleMatrix ) stack.new_recording();

      for ( int j = 0; j < nve; j++ ) {
        for ( unsigned idim = 0; idim < dim; idim++ ) {
          vxOld[idim][j] = vx_hat[idim][j] + Soli_old[indexVAR[idim]][j];
          vxNew[idim][j] = vx_hat[idim][j] + Soli[indexVAR[idim]][j];
        }
      }

      //Boundary integral
      {
        vector < adept::adouble> normal ( dim, 0 );

        // loop on faces
        for ( unsigned jface = 0; jface < mymsh->GetElementFaceNumber ( iel ); jface++ ) {
          std::vector< double > xx ( dim, 0. );

          // look for boundary faces
          if ( myel->GetFaceElementIndex ( iel, jface ) < 0 ) {

            unsigned int face = - ( mymsh->el->GetFaceElementIndex ( iel, jface ) + 1 );

            for ( unsigned tip = 0; tip < TimeIntPoints; tip++ ) {

              double tau = 0.;

              if ( !ml_sol->GetBdcFunction() ( xx, "PS", tau, face, time + ( -1. + s[tip] ) * dt )  && ( tau != 0. || Kslip != 0. ) ) {

                unsigned nve = mymsh->GetElementFaceDofNumber ( iel, jface, SolType2 );
                const unsigned felt = mymsh->GetElementFaceType ( iel, jface );

                vector <vector < adept::adouble > > sol_face ( dim );
                for ( unsigned idim = 0; idim < dim; idim++ ) {
                  vx_face[idim].resize ( nve );
                  sol_face[idim].resize ( nve );
                }

                std::vector < unsigned > ilocal ( nve );

                unsigned shift = ( flag_mat == 4 ) ? 0 : 1;

                for ( unsigned i = 0; i < nve; i++ ) {
                  ilocal[i] = mymsh->GetLocalFaceVertexIndex ( iel, jface, i );

                  for ( unsigned idim = 0; idim < dim; idim++ ) {
                    vx_face[idim][i]  = vx_hat[idim][ilocal[i]] +
                                        ( ( 1. - s[tip] ) * Soli_old[indexVAR[idim]][ilocal[i]] + s[tip] * Soli[indexVAR[idim]][ilocal[i]] );
                    sol_face[idim][i] = ( ( 1. - s[tip] ) * Soli_old[indexVAR[dim * shift + idim]][ilocal[i]] + s[tip] * Soli[indexVAR[dim * shift + idim]][ilocal[i]] );
                  }
                }

                for ( unsigned igs = 0; igs < mymsh->_finiteElement[felt][SolType2]->GetGaussPointNumber(); igs++ ) {
                  mymsh->_finiteElement[felt][SolType2]->JacobianSur ( vx_face, igs, Weight, phi, gradphi, normal );

                  vector < adept::adouble > solVAR_face ( dim );
                  for ( unsigned idim = 0; idim < dim; idim++ ) {
                    solVAR_face[idim] = 0.;
                    for ( unsigned i = 0; i < nve; i++ ) {
                      solVAR_face[idim] += phi[i] * sol_face[idim][i];
                    }
                  }

                  // *** phi_i loop ***
                  for ( unsigned i = 0; i < nve; i++ ) {
                    std::vector < adept::adouble> value ( dim );
                    for ( unsigned idim = 0; idim < dim; idim++ ) {
                      value[idim] =  - theta[tip] * phi[i] * ( tau * normal[idim] +  Kslip * solVAR_face[idim] ) / rhof * Weight;
                    }

                    for ( unsigned idim = 0; idim < dim; idim++ ) {
                      if ( ( !solidmark[ilocal[i]] ) ) {
                        aRhs[indexVAR[dim + idim]][ilocal[i]] += value[idim];
                      }
                      else {   //if interface node it goes to solid
                        aRhs[indexVAR[idim]][ilocal[i]] += value[idim];
                      }
                    }
                  }
                }
                Kslip = 0.;
              }
            }
          }
        }
      }

      double area = 1.;
      // *** Volume integral ***
      for ( unsigned ig = 0; ig < mymsh->_finiteElement[ielt][SolType2]->GetGaussPointNumber(); ig++ ) {
        mymsh->_finiteElement[ielt][SolType2]->Jacobian ( vx_hat, ig, Weight_hat, phi_hat, gradphi_hat);
        phi1 = mymsh->_finiteElement[ielt][SolType1]->GetPhi ( ig );

        // ---------------------------------------------------------------------------
        // store displacenet gradient
        for ( int i = 0; i < dim; i++ ) {
          for ( int j = 0; j < dim; j++ ) {
            GradSolhatVAR[i][j] = 0.;
          }

          for ( unsigned inode = 0; inode < nve; inode++ ) {
            for ( int j = 0; j < dim; j++ ) {
              GradSolhatVAR[i][j]     += gradphi_hat[inode * dim + j] * Soli[indexVAR[i]][inode];
            }
          }
        }

        // store displacenet and velocity solution old and New
        for ( int i = 0; i < dim * nBlocks; i++ ) {
          SolVAROld[i] = 0.;
          SolVARNew[i] = 0.;

          for ( unsigned inode = 0; inode < nve; inode++ ) {
            SolVAROld[i] += phi_hat[inode] * Soli_old[indexVAR[i]][inode];
            SolVARNew[i] += phi_hat[inode] * Soli[indexVAR[i]][inode];
          }
        }

        // store pressure (we use only one pressure in the time interval)

        for ( unsigned i = 0; i < nP; i++ ) {
          SolVAR[nBlocks * dim + i] = 0.;
          for ( unsigned inode = 0; inode < nve1; inode++ ) {
            SolVAR[nBlocks * dim + i] += phi1[inode] * Soli[indexVAR[nBlocks * dim + i]][inode];
          }
        }
        
        vector < double > accVAROld(dim,0.);
	vector < adept::adouble > accVARNew(dim,0.);
        
        //store acceleration old and new in the gauss point
        for ( int i = 0; i < dim; i++ ) {
	  for ( unsigned inode = 0; inode < nve; inode++ ) {
	    accVAROld[i] += phi_hat[inode] * aold[i][inode];
	  }
	  accVARNew[i] =  ( SolVARNew[i] - SolVAROld[i] )/( dt * dt * beta ) 
			- SolVAROld[i + dim]/( dt * beta ) 
			- (1. - 2. * beta)/(2. * beta) * accVAROld[i];
	}
          
        
        

        adept::adouble J_hat;
        adept::adouble F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
        double scaleWeightNoJacobian = 1.;

        //BEGIN SOLID ASSEMBLY ============
        if ( flag_mat == 4 ) {
	  
          //BEGIN redidual d_t - v = 0 in fixed domain
          for ( unsigned i = 0; i < nve; i++ ) {
            for ( int idim = 0; idim < dim; idim++ ) {
//               aRhs[indexVAR[dim + idim]][i] += ((SolVARNew[idim] - SolVAROld[idim]) / dt -
//                                                 (1. * SolVARNew[dim + idim] + 0 * SolVAROld[dim + idim])
//                                                ) * phi_hat[i] * Weight_hat;

              aRhs[indexVAR[dim + idim]][i] += ( ( Soli[idim][i] - Soli_old[idim][i] ) / dt -
                                                 Soli[idim + dim][i]
                                               ) * phi_hat[i] * Weight_hat;
	      
//               adept::adouble anew = ( Soli[idim][i] - Soli_old[idim][i] )/( dt* dt * beta ) - Soli_old[idim + dim][i]/( dt * beta ) - (1 - 2. * beta)/(2. * beta) * aold[idim][i];
//               adept::adouble vnew = Soli_old[idim + dim][i] + dt * ( (1 - Gamma) * aold[idim][i] + Gamma * anew );
// 	      
// 	      aRhs[indexVAR[dim + idim]][i] += (vnew - Soli[idim + dim][i]) * phi_hat[i] * Weight_hat;

	      
            }
          }
          //END redidual d_t - v = 0 in fixed domain

          //BEGIN Jacobian in the undeformed configuration
          //physical quantity
          for ( int i = 0; i < dim; i++ ) {
            for ( int j = 0; j < dim; j++ ) {
              F[i][j] += GradSolhatVAR[i][j];
            }
          }

          J_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                   - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];
          //END Jacobian in the undeformed configuration

	  double bkm = bulkModule;
	  bool incpr  = incompressible;	   	    
	  if( elementGroup == 15 ){
	    bkm = bulkModule1;	   
	    incpr = incompressible1;	   	    
	  }
	  	  
		   
          //BEGIN continuity block
          for ( unsigned i = 0; i < nve1; i++ ) {
            aRhs[indexVAR[nBlocks * dim]][i] += - phi1[i] * ( J_hat - 1.) * Weight_hat / rhof;
	    if(!incpr) aRhs[indexVAR[nBlocks * dim]][i] += phi1[i] * ( - SolVAR[nBlocks * dim] / bkm ) * Weight_hat / rhof;
            //if (twoPressure) aRhs[indexVAR[nBlocks * dim + jP]][i] +=  phi1[i] * SolVAR[nBlocks * dim + jP] * Weight_hat;   // fake continuity of the fluid in the solid
            if ( twoPressure ) aRhs[indexVAR[nBlocks * dim + jP]][i] +=  phi1[i] * Soli[indexVAR[nBlocks * dim + jP]][i] * Weight_hat; // fake continuity of the fluid in the solid
          }
          //END continuity block
        }
        //END SOLID ASSEMBLY ============
        else {
          // store gauss point coordinates old and new and mesh velocity old
          for ( int i = 0; i < dim; i++ ) {
            vxOld_ig[i] = 0.;
            vxNew_ig[i] = 0.;
            meshVelOld[i] = 0.;

            for ( unsigned inode = 0; inode < nve; inode++ ) {
              vxOld_ig[i] += phi_hat[inode] * vxOld[i][inode];
              vxNew_ig[i] += phi_hat[inode] * vxNew[i][inode];
              meshVelOld[i] += phi_hat[inode] * meshVelOldNode[i][inode];
            }
          }
          // get scaling factor for Weight_nojac
          std::vector<double> xc ( dim, 0. );

          if ( valve ) {
            if ( dim == 2 ) {
              xc[0] = -0.000193; // vein_valve_closed valve2_corta2.neu
              xc[1] = 0.067541 * 0;
            }
            else if ( dim == 3 ) {
              xc[0] = 0.0015; // vein_valve_closed
              xc[1] = 0.0;
              xc[2] = 0.01;
            }

            double distance = 0.;

            for ( unsigned k = 0; k < dim; k++ ) {
              distance += ( vx_hat[k][ nve - 1] - xc[k] ) * ( vx_hat[k][nve - 1] - xc[k] );
            }

            distance = sqrt ( distance );
            scaleWeightNoJacobian = 1. / ( 1 + 10000 * distance );

            if ( elementGroup == 16 ) scaleWeightNoJacobian *= 100;
          }
        }

        if ( !valve ) {
          if ( flag_mat == 2 || flag_mat == 3  || iel == mymsh->_elementOffset[iproc] ) {
            if ( ig == 0 ) {
              double GaussWeight = mymsh->_finiteElement[ielt][SolType2]->GetGaussWeight ( ig );
              area = Weight_hat / GaussWeight;

              if ( iel == mymsh->_elementOffset[iproc] ) {
                area_elem_first->set ( mymsh->processor_id(), area );
                area_elem_first->close();
                rapresentative_area = area_elem_first->l1_norm() / nprocs;
              }
            }
            scaleWeightNoJacobian = rapresentative_area / area;
          }
        }


        for ( unsigned tip = 0; tip <= TimeIntPoints; tip++ ) {

          bool end = ( tip == TimeIntPoints ) ? true : false;
          bool middle = ( !end ) ? true : false;

          // store displacenet and velocity at the current time
          for ( int i = 0; i < dim * nBlocks; i++ ) {
            SolVAR[i] = SolVAROld[i] * ( 1. - s[tip] ) +  SolVARNew[i] * s[tip];
          }
          
          // store acceleration at the current time
          std::vector <adept::adouble> accVAR(dim,0.);
          for ( int i = 0; i < dim; i++ ) {
            accVAR[i] = accVAROld[i] * ( 1. - s[tip] ) +  accVARNew[i] * s[tip];
          }
          
          // store nodes coordinates and Gauss point coordinates at the current time
          for ( unsigned i = 0; i < dim; i++ ) {
            for ( int j = 0; j < nve; j++ ) {
              vx[i][j] = vxOld[i][j] * ( 1. - s[tip] ) +  vxNew[i][j] * s[tip];
            }
          }
          if ( middle && flag_mat != 4 ) {
            mymsh->_finiteElement[ielt][SolType2]->Jacobian ( vx, ig, Weight, phi, gradphi, nablaphi );
          }
          else {
            mymsh->_finiteElement[ielt][SolType2]->Jacobian ( vx, ig, Weight, phi, gradphi);
          }

          // store solution gradient and laplace at the current time
          for ( int i = 0; i < nBlocks * dim; i++ ) {
            for ( int j = 0; j < dim; j++ ) {
              GradSolVAR[i][j] = 0.;
            }
            for ( unsigned inode = 0; inode < nve; inode++ ) {
              for ( int j = 0; j < dim; j++ ) {
                GradSolVAR[i][j]     += gradphi[inode * dim + j]  * ( s[tip] * Soli[indexVAR[i]][inode] + ( 1. - s[tip] ) * Soli_old[indexVAR[i]][inode] );
              }
            }
          }
      
          //BEGIN FLUID and POROUS MEDIA ASSEMBLY ============
          if ( flag_mat != 4 ) {
            if ( middle ) {
              vector < adept::adouble > gradPf ( dim, 0. );
              std::vector <std::vector <adept::adouble> > jacobianMatrix;
              mymsh->_finiteElement[ielt][SolType2]->GetJacobian ( vx, ig, Weight, jacobianMatrix );

              for ( unsigned i = 0; i < dim; i++ ) {
                for ( unsigned j = 0; j < dim; j++ ) {
                  gradPf[i] +=  Soli[indexVAR[nBlocks * dim + jP]][ 1 + j ] * jacobianMatrix[j][i];
                }
              }


              // store solution gradient and laplace at the current time
              for ( int i = dim; i < nBlocks * dim; i++ ) {

                for ( int j = 0; j < nabla_dim; j++ ) {
                  NablaSolVAR[i][j] = 0.;
                }
                for ( unsigned inode = 0; inode < nve; inode++ ) {

                  for ( int j = 0; j < nabla_dim; j++ ) {
                    NablaSolVAR[i][j]     += nablaphi[inode * nabla_dim + j] * ( s[tip] * Soli[indexVAR[i]][inode] + ( 1. - s[tip] ) * Soli_old[indexVAR[i]][inode] );
                  }
                }
              }

              // get Gauss point coordinates at the current time
              for ( unsigned i = 0; i < dim; i++ ) {
                vx_ig[i] = vxOld_ig[i] * ( 1. - s[tip] ) +  vxNew_ig[i] * s[tip];
              }

              adept::adouble springStiffness = 0.;
              if ( valve ) springStiffness = 2.e-4 * exp ( ( vx_ig[0] - ( - 1e-5 ) ) / 0.0001 ) * ( dim == 2 ) + // if 2D
                                             2.e-4 * exp ( ( vx_ig[0] - ( - 1e-5 ) ) / 0.0001 ) * ( dim == 3 ); //if 3D


              // get mesh velocity and gradient at current time
              for ( unsigned i = 0; i < dim; i++ ) {
                meshVel[i] = ( 1. - s[tip] ) * meshVelOld[i]  + s[tip] * ( SolVARNew[i] - SolVAROld[i] ) / dt ;
// 		meshVel[i] = ( 1. - s[tip] ) * meshVelOld[i]  
// 			    + s[tip] * ( meshVelOld[i] + dt * ( (1. - Gamma) * accVAROld[i] + Gamma * accVARNew[i] ) );
              }

              // speed
              adept::adouble speed2 = 0.;
              for ( int i = 0; i < dim; i++ ) {
                speed2 += SolVAR[i + dim] * SolVAR[i + dim];
              }
              adept::adouble speed = sqrt ( speed2 );

              // phiSupg evaluation
              double sqrtlambdak = ( *mysolution->_Sol[indLmbd] ) ( iel );
              adept::adouble tauSupg = 1. / ( sqrtlambdak * sqrtlambdak * 4.*IRe );
              adept::adouble Rek = speed / ( 4.*sqrtlambdak * IRe );

              if ( Rek > 1.0e-15 ) {
                adept::adouble xiRek = ( Rek >= 1. ) ? 1. : Rek;
                tauSupg   = xiRek / ( speed * sqrtlambdak );
              }

              //tauSupg=0;
              vector < adept::adouble > phiSupg ( nve, 0. );

              for ( unsigned i = 0; i < nve; i++ ) {
                for ( unsigned j = 0; j < dim; j++ ) {
                  phiSupg[i] += ( ( SolVAR[j + dim] - meshVel[j] ) * gradphi[i * dim + j] ) * tauSupg;
                }
              }

              //BEGIN ALE + Momentum (Navier-Stokes)
              for ( unsigned i = 0; i < nve; i++ ) {
                if ( flag_mat == 2 ) {
                  //BEGIN residual Navier-Stokes in moving domain
                  adept::adouble Lapdisp[3] = {0., 0., 0.};
                  adept::adouble LapvelVAR[3] = {0., 0., 0.};
                  adept::adouble LapStrong[3] = {0., 0., 0.};
                  adept::adouble AdvaleVAR[3] = {0., 0., 0.};


                  for ( int idim = 0.; idim < dim; idim++ ) {
                    for ( int jdim = 0.; jdim < dim; jdim++ ) {
                      unsigned kdim;

                      if ( idim == jdim ) kdim = jdim;
                      else if ( 1 == idim + jdim ) kdim = dim;    // xy
                      else if ( 2 == idim + jdim ) kdim = dim + 2; // xz
                      else if ( 3 == idim + jdim ) kdim = dim + 1; // yz

                      //laplaciano debole
                      LapvelVAR[idim]     += ( GradSolVAR[dim + idim][jdim] + GradSolVAR[dim + jdim][idim] ) * gradphi[i * dim + jdim];
                      Lapdisp[idim]     += ( idim == 0 ) * ( GradSolVAR[idim][jdim] + GradSolVAR[jdim][idim] ) * gradphi[i * dim + jdim];
                      //laplaciano strong
                      LapStrong[idim]     += ( NablaSolVAR[dim + idim][jdim] + NablaSolVAR[dim + jdim][kdim] ) * phiSupg[i];
                      AdvaleVAR[idim]	+= ( ( SolVAR[dim + jdim] - meshVel[jdim] ) * GradSolVAR[dim + idim][jdim]
                                           + ( dim == 2 ) * GradSolVAR[dim + jdim][jdim] * SolVAR[dim + idim] ) * ( phi[i] + phiSupg[i] );
                    }
                  }

                  for ( int idim = 0; idim < dim; idim++ ) {
                    adept::adouble timeDerivative = 0.;
                    adept::adouble value = 0.;


                    timeDerivative = theta[tip] * ( SolVARNew[dim + idim] - SolVAROld[dim + idim] ) * ( phi[i] + phiSupg[i] ) * Weight / dt;

                    value =  theta[tip] * (
                               + AdvaleVAR[idim]      	             // advection term
                               + IRe * LapvelVAR[idim]	             // viscous dissipation
                               + ( idim == 0 ) * springStiffness * Lapdisp[idim] // penalty spring
                               - IRe * LapStrong[idim]
                               - 1. / rhof * SolVAR[nBlocks * dim + jP] * gradphi[i * dim + idim] // pressure gradient weak
                               + 1. / rhof * gradPf[idim] * phiSupg[i] // pressure gradient supg
                             ) * Weight;

                    if ( !solidmark[i] ) {
                      aRhs[indexVAR[dim + idim]][i] += timeDerivative + value;
                    }
                    else {
                      aRhs[indexVAR[idim]][i] += timeDerivative + value;
                    }
                  }
                  //END redidual Navier-Stokes in moving domain
                }
                else if ( flag_mat == 3 ) {

                  double eps = 1.0e-20;
                  speed = sqrt ( speed2 + eps );

                  //BEGIN redidual Porous Media in moving domain
                  double DE = 0.;

                  if ( dim == 2 ) {
                    DE = 0.0002; // AAA_thrombus_2D
                    //DE = 0.00006; // turek2D
                  }
                  else if ( dim == 3 ) {
                    DE = 0.000112; // porous3D
                  }

                  double b = 4188;
                  double a = 1452;
                  double K = DE * IRe * rhof / b; // alpha = mu/b * De
                  double C2 = 2 * a / ( rhof * DE );

                  for ( int idim = 0; idim < dim; idim++ ) {
                    adept::adouble value = 0.;


                    value = theta[tip] * ( ( SolVAR[dim + idim] - meshVel[idim] ) * ( IRe / K + 0.5 * C2 * speed ) * ( phi[i] + phiSupg[i] )
                                           - 1. / rhof * SolVAR[nBlocks * dim + jP] * gradphi[i * dim + idim] // pressure gradient
                                           + 1. / rhof * gradPf[idim] * phiSupg[i] // pressure gradient supg
                                         ) * Weight;

                    if ( !solidmark[i] ) {
                      aRhs[indexVAR[dim + idim]][i] += value;
                    }
                    else {
                      aRhs[indexVAR[idim]][i] += value;
                    }
                  }
                  //END redidual Porous Media in moving domain
                }
              }
              //END Momentum (Navier-Stokes)
            }
            else {
              //BEGIN redidual Laplacian ALE map in the reference domain
              for ( unsigned i = 0; i < nve; i++ ) {
                if ( !solidmark[i] ) {
                  adept::adouble LapmapVAR[3] = {0., 0., 0.};
                  adept::adouble Weight_nojac = Weight * scaleWeightNoJacobian;

                  for ( int idim = 0; idim < dim; idim++ ) {
                    for ( int jdim = 0; jdim < dim; jdim++ ) {
                      LapmapVAR[idim] += ( GradSolVAR[idim][jdim] + 0.*GradSolVAR[jdim][idim] ) * gradphi[i * dim + jdim];
                    }
                  }

                  for ( int idim = 0; idim < dim; idim++ ) {
                    aRhs[indexVAR[idim]][i] += ( LapmapVAR[idim] * Weight_nojac );
                  }
                }
              }
              //END residual Laplacian ALE map in reference domain

              //BEGIN continuity block
              adept::adouble div_vel = 0.;

              for ( int i = 0; i < dim; i++ ) {
                div_vel += GradSolVAR[dim + i][i];
              }

              for ( unsigned i = 0; i < nve1; i++ ) {
                //if (twoPressure) aRhs[indexVAR[nBlocks * dim] ][i] +=  phi1[i] * SolVAR[nBlocks * dim ] * Weight_hat;   // fake continuity of the solid in the fluid
                if ( twoPressure ) aRhs[indexVAR[nBlocks * dim] ][i] +=  phi1[i] * Soli[indexVAR[nBlocks * dim]][i] * Weight_hat; // fake continuity of the solid in the fluid
                aRhs[indexVAR[nBlocks * dim + jP]][i] += ( -phi1[i] * div_vel ) * Weight / rhof;
              }
              //END continuity block ===========================
            }
          }
          //END FLUID ASSEMBLY ============

          //BEGIN SOLID ASSEMBLY ============
          else {
            if ( middle ) {
              //BEGIN build Chauchy Stress in moving domain
              adept::adouble Cauchy[3][3];
              //double Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};

              adept::adouble B[3][3];
              adept::adouble I1_B = 0.;

              for ( int I = 0; I < 3; ++I ) {
                for ( int J = 0; J < 3; ++J ) {
                  B[I][J] = 0.;

                  for ( int K = 0; K < 3; ++K ) {
                    //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
                    B[I][J] += F[I][K] * F[J][K];
                  }
                }
              }

              I1_B = B[0][0] + B[1][1] + B[2][2];
	      
	      double C1 = mus / 3.;
	      bool incpr  = incompressible;	   	    
	      if( elementGroup == 15 ){
		C1 = mus1 / 3.;
		incpr = incompressible1;	   	    
	      }
              double C2 = C1 / 2.;

	      if( incpr ){ //incompressible
		adept::adouble detB =   B[0][0] * ( B[1][1] * B[2][2] - B[2][1] * B[1][2] )
                                      - B[0][1] * ( B[2][2] * B[1][0] - B[1][2] * B[2][0] )
                                      + B[0][2] * ( B[1][0] * B[2][1] - B[2][0] * B[1][1] );

		adept::adouble invdetB = 1. / detB;
		adept::adouble invB[3][3];

		invB[0][0] = ( B[1][1] * B[2][2] - B[1][2] * B[2][1] ) * invdetB;
		invB[1][0] = - ( B[0][1] * B[2][2] - B[0][2] * B[2][1] ) * invdetB;
		invB[2][0] = ( B[0][1] * B[1][2] - B[0][2] * B[1][1] ) * invdetB;
		invB[0][1] = - ( B[1][0] * B[2][2] - B[1][2] * B[2][0] ) * invdetB;
		invB[1][1] = ( B[0][0] * B[2][2] - B[0][2] * B[2][0] ) * invdetB;
		invB[2][1] = - ( B[0][0] * B[1][2] - B[1][0] * B[0][2] ) * invdetB;
		invB[0][2] = ( B[1][0] * B[2][1] - B[2][0] * B[1][1] ) * invdetB;
		invB[1][2] = - ( B[0][0] * B[2][1] - B[2][0] * B[0][1] ) * invdetB;
		invB[2][2] = ( B[0][0] * B[1][1] - B[1][0] * B[0][1] ) * invdetB;
		
		for ( int I = 0; I < 3; ++I ) {
		  for ( int J = 0; J < 3; ++J ) {
		    Cauchy[I][J] =  2.* ( C1 * B[I][J] - C2 * invB[I][J] );
		  }
		  Cauchy[I][I] -= 1. / rhof * SolVAR[nBlocks * dim ];
		}
	      }
	      else { //slightly compressible
		adept::adouble B2[3][3];
		for (int I = 0; I < 3; ++I) {
		  for (int J = 0; J < 3; ++J) {
		    B2[I][J] = 0.;
		    for (int K = 0; K < 3; ++K) {
		      B2[I][J] += B[I][K] * B[K][J]; 
		    }
		  }
		}

		adept::adouble I2_B = B2[0][0] + B2[1][1] + B2[2][2];
	      
		adept::adouble J_hatm2over3 = pow(J_hat, -2. / 3.);
		adept::adouble J_hatm4over3 = J_hatm2over3 * J_hatm2over3;
		
// 		for (int I = 0; I < 3; ++I) {
// 		  for (int J = 0; J < 3; ++J) {
// 		    Cauchy[I][J] =  1. / J_hat * (   2. * J_hatm2over3 * (C1 + J_hatm2over3 * I1_B * C2) * B[I][J] - 2. * J_hatm4over3 * C2 * B2[I][J] 
// 						  - 2. / 3. * (C1 * J_hatm2over3 * I1_B + 2 * C2 * J_hatm4over3 * I2_B) * Id2th[I][J] 
// 						  - SolVAR[nBlocks * dim ] * Id2th[I][J] / rhof
// 						);
// 		  }
// 		}
		
		adept::adouble twoC1 = 2. * J_hatm2over3 * (C1 + J_hatm2over3 * I1_B * C2) / J_hat;
		adept::adouble twoC2 = 2. * J_hatm4over3 * C2 / J_hat;
		adept::adouble C3 = 2. / 3. * (C1 * J_hatm2over3 * I1_B + 2 * C2 * J_hatm4over3 * I2_B) / J_hat;
		adept::adouble C4 = 1. / ( J_hat * rhof); 
		
		for (int I = 0; I < 3; ++I) {
		  for (int J = 0; J < 3; ++J) {
		    Cauchy[I][J] =  twoC1 * B[I][J] - twoC2 * B2[I][J]; 
		  }
		  Cauchy[I][I] -= C3 + C4 * SolVAR[nBlocks * dim ];
		}
              }
              //END build Cauchy Stress in moving domain


              //BEGIN Momentum (Solid)
	      
              for ( unsigned i = 0; i < nve; i++ ) {
                adept::adouble CauchyDIR[3] = {0., 0., 0.};
                for ( int idim = 0.; idim < dim; idim++ ) {
                  for ( int jdim = 0.; jdim < dim; jdim++ ) {
                    CauchyDIR[idim] += gradphi[i * dim + jdim] * Cauchy[idim][jdim];
                  }
                }
                for ( int idim = 0; idim < dim; idim++ ) {
                  adept::adouble timeDerivative = 0.;
                  adept::adouble value = 0.;
                  timeDerivative = theta[tip] * rhos * ( SolVARNew[dim + idim] - SolVAROld[dim + idim] ) * phi[i] * Weight / dt;

		  //timeDerivative = theta[tip] * rhos * accVAR[idim] * phi[i] * Weight;
		  
                  value =  theta[tip] * ( - rhos * phi[i] * _gravity[idim]    // body force
                                          + CauchyDIR[idim]		     // stress
                                        ) * Weight;
                  aRhs[indexVAR[idim]][i] += timeDerivative + value;
                }
              }
              //END Momentum (Solid)
            }
          }
          //END SOLID ASSEMBLY ============
        }
      }

      //BEGIN local to global assembly
      //copy adouble -aRhs into double Rhs
      for ( unsigned i = 0; i < nBlocks * dim; i++ ) {
        Rhs[indexVAR[i]].resize ( nve );
        for ( int j = 0; j < nve; j++ ) {
          Rhs[indexVAR[i]][j] = -aRhs[indexVAR[i]][j].value();
        }
      }
      for ( unsigned i = 0; i < nP; i++ ) {
        Rhs[indexVAR[nBlocks * dim + i]].resize ( nve1 );
        for ( unsigned j = 0; j < nve1; j++ ) {
          Rhs[indexVAR[nBlocks * dim + i]][j] = -aRhs[indexVAR[nBlocks * dim + i]][j].value();
        }
      }
      for ( int i = 0; i < nBlocks * dim + nP; i++ ) {
        myRES->add_vector_blocked ( Rhs[indexVAR[i]], dofsVAR[i] );
      }

      if ( assembleMatrix ) {
        //Store equations
        for ( int i = 0; i < nBlocks * dim; i++ ) {
          stack.dependent ( &aRhs[indexVAR[i]][0], nve );
          stack.independent ( &Soli[indexVAR[i]][0], nve );
        }
        for ( int i = 0; i < nP; i++ ) {
          stack.dependent ( &aRhs[indexVAR[nBlocks * dim + i]][0], nve1 );
          stack.independent ( &Soli[indexVAR[nBlocks * dim + i]][0], nve1 );
        }
        Jac.resize ( ( nBlocks * dim * nve + nP * nve1 ) * ( nBlocks * dim * nve + nP * nve1 ) );

        start_time = clock();
        stack.jacobian_reverse ( &Jac[0], true );
        AdeptTime += clock() - start_time;

        myKK->add_matrix_blocked ( Jac, dofsAll, dofsAll );

        stack.clear_independents();
        stack.clear_dependents();

        //END local to global assembly
      }
    } //end list of elements loop

    if ( assembleMatrix ) myKK->close();

    myRES->close();

    delete area_elem_first;
    // *************************************
    std::cout << "   ********* Level Max " << igrid + 1 << " ADEPT TIME:\t" << \
              static_cast<double> ( AdeptTime ) / CLOCKS_PER_SEC << std::endl;
    // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************
  }

//****************************************************************************************


  void SetLambdaNew ( MultiLevelSolution & mlSol, const unsigned & level, const  FEOrder & order, Operator operatorType )
  {

    unsigned SolType;

    if ( order < FIRST || order > SECOND ) {
      std::cout << "Wong Solution Order" << std::endl;
      exit ( 0 );
    }
    else if ( order == FIRST ) SolType = 0;
    else if ( order == SERENDIPITY ) SolType = 1;
    else if ( order == SECOND ) SolType = 2;



    clock_t GetLambdaTime = 0;
    clock_t start_time, end_time;
    start_time = clock();

    adept::Stack& adeptStack = FemusInit::_adeptStack;

    Solution* mysolution = mlSol.GetSolutionLevel ( level );
    Mesh* mymsh	=  mlSol._mlMesh->GetLevel ( level );
    elem* myel	=  mymsh->el;

    unsigned indLmbd = mlSol.GetIndex ( "lmbd" );

    const unsigned geoDim = mymsh->GetDimension();
    const unsigned nablaGoeDim = ( 3 * ( geoDim - 1 ) + ! ( geoDim - 1 ) );
    const unsigned max_size = static_cast< unsigned > ( ceil ( pow ( 3, geoDim ) ) );


    const char varname[3][4] = {"DX", "DY", "DZ"};
    vector <unsigned> indVAR ( geoDim );

    for ( unsigned ivar = 0; ivar < geoDim; ivar++ ) {
      indVAR[ivar] = mlSol.GetIndex ( &varname[ivar][0] );
    }


    bool diffusion, elasticity;

    if ( operatorType == DIFFUSION ) {
      diffusion  = true;
      elasticity = false;
    }

    if ( operatorType == ELASTICITY ) {
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
    vector<vector<adept::adouble> > GradSolVAR ( varDim );
    vector<vector<adept::adouble> > NablaSolVAR ( varDim );

    for ( int ivar = 0; ivar < varDim; ivar++ ) {
      GradSolVAR[ivar].resize ( geoDim );
      NablaSolVAR[ivar].resize ( nablaGoeDim );
    }

    vector <double > phi;
    vector <adept::adouble> gradphi;
    vector <adept::adouble> nablaphi;
    adept::adouble Weight;

    phi.reserve ( max_size );
    gradphi.reserve ( max_size * geoDim );
    nablaphi.reserve ( max_size * nablaGoeDim );

    vector <vector < adept::adouble> > vx ( geoDim );

    for ( int ivar = 0; ivar < geoDim; ivar++ ) {
      vx[ivar].reserve ( max_size );
    }

    unsigned SolTypeVx = 2.;

    vector< vector< adept::adouble > > Soli ( varDim );
    vector< vector< adept::adouble > > aRhs ( varDim );
    vector< vector< adept::adouble > > aLhs ( varDim );

    for ( int ivar = 0; ivar < varDim; ivar++ ) {
      Soli[ivar].reserve ( max_size );
      aRhs[ivar].reserve ( max_size );
      aLhs[ivar].reserve ( max_size );
    }

    vector < double > K;
    K.reserve ( ( max_size * varDim ) * ( max_size * varDim ) );
    vector < double > M;
    M.reserve ( ( max_size * varDim ) * ( max_size * varDim ) );

    // mesh and procs
    unsigned nel    = mymsh->GetNumberOfElements();
    unsigned iproc  = mymsh->processor_id();

    // *** element loop ***
    for ( int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++ ) {

      unsigned kel        = iel;
      short unsigned kelt = mymsh->GetElementType ( kel );
      unsigned nve        = mymsh->GetElementDofNumber ( kel, SolType ) - 1;
      unsigned nveVx      = mymsh->GetElementDofNumber ( kel, SolTypeVx );

      // -------------- resize --------------
      for ( int ivar = 0; ivar < varDim; ivar++ ) {
        Soli[ivar].resize ( nve );
        aRhs[ivar].resize ( nve );
        aLhs[ivar].resize ( nve );
      }

      M.resize ( ( varDim * nve ) * ( varDim * nve ) );
      K.resize ( ( varDim * nve ) * ( varDim * nve ) );
      // ------------------------------------

      // ------------ get coordinates -------
      for ( int i = 0; i < geoDim; i++ ) {
        vx[i].resize ( nveVx );
        //vx1[i].resize(nveVx);
      }

      for ( unsigned i = 0; i < nveVx; i++ ) {
        unsigned inodeVx_Metis = mymsh->GetSolutionDof ( i, iel, SolTypeVx );

        for ( int j = 0; j < geoDim; j++ ) {
          //coordinates
          vx[j][i] = ( *mymsh->_topology->_Sol[j] ) ( inodeVx_Metis ) + ( *mysolution->_Sol[indVAR[j]] ) ( inodeVx_Metis );
        }
      }

      // ------------------------------------

      // ------------ init ------------------
      for ( unsigned i = 0; i < nve; i++ ) {
        for ( int ivar = 0; ivar < varDim; ivar++ ) {
          Soli[ivar][i] = 1.;
          aRhs[ivar][i] = 0.;
          aLhs[ivar][i] = 0.;
        }
      }

      // ------------------------------------

      adeptStack.new_recording();
      double hk = 1.;

      for ( unsigned ig = 0; ig < mymsh->_finiteElement[kelt][SolType]->GetGaussPointNumber(); ig++ ) {
        // *** get Jacobian and test function and test function derivatives in the moving frame***
        mymsh->_finiteElement[kelt][SolType]->Jacobian ( vx, ig, Weight, phi, gradphi, nablaphi );

        if ( ig == 0 ) {
          double referenceElementScale[6] = {8., 1. / 6., 1., 4., 1., 2.};
          double GaussWeight = mymsh->_finiteElement[kelt][SolType]->GetGaussWeight ( ig );
          double area = referenceElementScale[kelt] * Weight.value() / GaussWeight;
          hk = pow ( area, 1. / geoDim );

          //cout<<hk<<endl;
          if ( 0 == SolType ) break;
        }

        for ( int ivar = 0; ivar < varDim; ivar++ ) {
          for ( int jvar = 0; jvar < geoDim; jvar++ ) {
            GradSolVAR[ivar][jvar] = 0.;
          }

          for ( int jvar = 0; jvar < nablaGoeDim; jvar++ ) {
            NablaSolVAR[ivar][jvar] = 0.;
          }

          for ( unsigned inode = 0; inode < nve; inode++ ) {
            adept::adouble soli = Soli[ivar][inode];

            for ( int jvar = 0; jvar < geoDim; jvar++ ) {
              GradSolVAR[ivar][jvar] += gradphi[inode * geoDim + jvar] * soli;
            }

            for ( int jvar = 0; jvar < nablaGoeDim; jvar++ ) {
              NablaSolVAR[ivar][jvar] += nablaphi[inode * nablaGoeDim + jvar] * soli;
            }
          }
        }


        vector < adept::adouble > divGradSol ( varDim, 0. );

        for ( unsigned ivar = 0; ivar < varDim; ivar++ ) {
          for ( unsigned jvar = 0; jvar < geoDim; jvar++ ) {
            if ( diffusion ) {
              divGradSol[ivar] += NablaSolVAR[ivar][jvar];
            }
            else if ( elasticity ) {
              unsigned kvar;

              if ( ivar == jvar ) kvar = jvar;
              else if ( 1 == ivar + jvar ) kvar = geoDim;  // xy
              else if ( 2 == ivar + jvar ) kvar = geoDim + 2; // xz
              else if ( 3 == ivar + jvar ) kvar = geoDim + 1; // yz

              divGradSol[ivar]   += 0.5 * ( NablaSolVAR[ivar][jvar] + NablaSolVAR[jvar][kvar] );
            }
          }
        }

        //BEGIN local assembly
        for ( unsigned i = 0; i < nve; i++ ) {
          for ( unsigned ivar = 0; ivar < varDim; ivar++ ) {
            for ( unsigned jvar = 0; jvar < geoDim; jvar++ ) {
              aRhs[ivar][i] += gradphi[i * geoDim + jvar] * ( GradSolVAR[ivar][jvar] ) * Weight;

              if ( diffusion ) {
                aLhs[ivar][i] +=  divGradSol[ivar] * nablaphi[i * nablaGoeDim + jvar] * Weight;
                //aRhs[ivar][i] += gradphi[i*geoDim+jvar]*(GradSolVAR[ivar][jvar]) * Weight;
              }
              else if ( elasticity ) {
                unsigned kvar;

                if ( ivar == jvar ) kvar = jvar;
                else if ( 1 == ivar + jvar ) kvar = geoDim;  // xy
                else if ( 2 == ivar + jvar ) kvar = geoDim + 2; // xz
                else if ( 3 == ivar + jvar ) kvar = geoDim + 1; // yz

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
      double lambdak = 6. / ( hk * hk ); //if SolType is linear

      if ( SolType == 1 || SolType == 2 ) { // only if solType is quadratic or biquadratic
        for ( int ivar = 0; ivar < varDim; ivar++ ) {
          adeptStack.independent ( &Soli[ivar][0], nve );
        }

        //Store RHS in M
        for ( int ivar = 0; ivar < varDim; ivar++ ) {
          adeptStack.dependent ( &aRhs[ivar][0], nve );
        }

        adeptStack.jacobian ( &M[0] );
        adeptStack.clear_dependents();

        //Store LHS in K
        for ( int ivar = 0; ivar < varDim; ivar++ ) {
          adeptStack.dependent ( &aLhs[ivar][0], nve );
        }

        adeptStack.jacobian ( &K[0] );
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
        for ( int k = 0 ; k < indSize - 1 ; k++ ) {
          for ( int i = k + 1 ; i < indSize ; i++ ) {
            M[i * matSize + k] /= M[k * matSize + k];

            for ( int j = k + 1 ; j < indSize ; j++ ) {
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


        vector < double > x ( matSize, 1. );
        vector < double > y ( matSize );

        double phik = x[0] + x[1];
        lambdak = 1.;
        double error = 1.;

        while ( error > 1.0e-10 ) {
          double phikm1 = phik;
          double lambdakm1 = lambdak;

          // y = K x
          for ( int i = 0; i < indSize; i++ ) {
            y[i] = 0.;

            for ( int j = 0; j < indSize; j++ ) {
              y[i] += K[ i * matSize + j ] * x[j];
            }
          }

          // y = L^(-1) y
          for ( int i = 0; i < indSize; i++ ) {
            for ( int j = 0; j < i; j++ ) {
              y[i] -= M[i * matSize + j] * y[j];
            }
          }

          // x <--  y = U^(-1) y
          double l2norm = 0.;

          for ( int i = indSize - 1; i >= 0; i-- ) {
            x[i] = y[i];

            for ( int j = i + 1; j < indSize; j++ ) {
              x[i] -= M[ i * matSize + j] * x[j];
            }

            x[i] /= M[i * matSize + i];
            l2norm += x[i] * x[i];
          }

          l2norm = sqrt ( l2norm );

          phik = ( x[0] + x[1] );
          lambdak =  phik / phikm1;

          for ( int i = 0; i < indSize; i++ ) {
            x[i] /= l2norm;
          }

          phik /= l2norm;
          error = fabs ( ( lambdak - lambdakm1 ) / lambdak );
        }
      }

      //std::cout << lambdak*hk*hk << std::endl;
      mysolution->_Sol[indLmbd]->set ( iel, sqrt ( lambdak ) );
      //abort();
    } //end list of elements loop


    mysolution->_Sol[indLmbd]->close();
    // *************************************
    end_time = clock();
    GetLambdaTime += ( end_time - start_time );

    std::cout << "GetLambda Time = " << GetLambdaTime / CLOCKS_PER_SEC << std::endl;

  }

//***************************************************************************************************************

  void StoreMeshVelocity ( MultiLevelProblem & ml_prob )
  {

    TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ( "Fluid-Structure-Interaction" );
    const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();

    MultiLevelSolution* ml_sol = ml_prob._ml_sol;

    Solution* solution = ml_sol->GetSolutionLevel ( level );
    Mesh* msh = ml_prob._ml_msh->GetLevel ( level );

    //const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();
    const unsigned dim = msh->GetDimension();
    const char varname[9][4] = {"DX", "U", "Um", "DY", "V", "Vm", "DZ", "W", "Wm"};
    const char varname2[3][4] = {"AX", "AY", "AZ"};

    vector <unsigned> indVAR ( 3 * dim );
    vector <unsigned> indVAR2 ( dim );

    for ( unsigned ivar = 0; ivar < 3 * dim; ivar++ ) {
      indVAR[ivar] = ml_sol->GetIndex ( &varname[ivar][0] );
    }
    
    for ( unsigned ivar = 0; ivar < dim; ivar++ ) {
      indVAR2[ivar] = ml_sol->GetIndex ( &varname2[ivar][0] );
    }

    double dt =  my_nnlin_impl_sys.GetIntervalTime();

    
    int  iproc;
    MPI_Comm_rank ( MPI_COMM_WORLD, &iproc );

    for ( unsigned idim = 0; idim < dim; idim++ ) {
      unsigned ivar = idim * 3;

      for ( unsigned jdof = msh->_dofOffset[2][iproc]; jdof < msh->_dofOffset[2][iproc + 1]; jdof++ ) {
        bool solidmark = msh->GetSolidMark ( jdof );
        double vnew;
	double anew;

        if ( solidmark != solidmark ) {
          vnew = ( *solution->_Sol[indVAR[ivar + 1]] ) ( jdof ); // solid: mesh velocity equals solid velocity
	  anew = 0.;
        }
        else {
          double unew = ( *solution->_Sol[indVAR[ivar]] ) ( jdof );
          double uold = ( *solution->_SolOld[indVAR[ivar]] ) ( jdof );
          double vold = ( *solution->_Sol[indVAR[ivar + 2]] ) ( jdof );
	  double aold = ( *solution->_SolOld[indVAR2[idim]] ) ( jdof );
	  anew = (unew - uold)/(dt * dt * beta) - vold/(dt * beta) - (1. - 2. * beta)/(2. * beta) * aold;
          vnew = 1. / dt * ( unew - uold ) - 0. * vold;
	  //vnew = vold + dt * ( (1. - Gamma) * aold + Gamma * anew );
        }
        solution->_Sol[indVAR[ivar + 2]]->set ( jdof, vnew );
	solution->_Sol[indVAR2[idim]]->set ( jdof, anew );
      }

      solution->_Sol[indVAR[ivar + 2]]->close();
      solution->_Sol[indVAR2[idim]]->close();
    }

  }

}





