/** tutorial/Ex1
 * This example shows how to:
 * initialize a femus application;
 * define the multilevel-mesh object mlMsh;
 * read from the file ./input/square.neu the coarse-level mesh and associate it to mlMsh;
 * add in mlMsh uniform refined level-meshes;
 * define the multilevel-solution object mlSol associated to mlMsh;
 * add in mlSol different types of finite element solution variables;
 * initialize the solution varables;
 * define vtk and gmv writer objects associated to mlSol;
 * print vtk and gmv binary-format files in ./output directory.
 **/

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "PetscMatrix.hpp"

#include "TransientSystem.hpp"
#include "LinearImplicitSystem.hpp"

#include "slepceps.h"
#include <slepcmfn.h>

using namespace femus;

double dt = 1.; 

double k_v = 0.01;

double pi = acos ( -1. );
double k_h = 1 / ( 10 * pi );

const unsigned NumberOfLayers = 20;

clock_t start_time = clock();

bool wave = false;
bool twostage = false;
bool assembly = true; //assembly must be left always true

const double hRest[20] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

double InitalValueV ( const std::vector < double >& x ) {
  return 1. / 10.;
}

double InitalValueH ( const std::vector < double >& x ) {
  return hRest[0];
}

double InitalValueT ( const std::vector < double >& x ) {
  double pi = acos ( -1. );
  if ( x[0] < 5 ) return 5;
  else return 30;
//  return (- sin(pi*x[0]));
}


double InitalValueB ( const std::vector < double >& x ) {
  return 10.; //( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );
}


bool SetBoundaryCondition ( const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time ) {
  bool dirichlet = false;
  if ( !strcmp ( SolName, "HT" ) ) {
    if ( facename == 1 || facename == 2 ) {
      dirichlet = true;
      value = 0.;
    }
  }
  if ( !strcmp ( SolName, "T" ) ) {
    if ( facename == 1 || facename == 2 ) {
      dirichlet = true;
      value = 0.;
    }
  }
  return dirichlet;
}


void ETD ( MultiLevelProblem& ml_prob );

void RK4 ( MultiLevelProblem& ml_prob, const bool & implicitEuler );


int main ( int argc, char** args ) {

  SlepcInitialize ( &argc, &args, PETSC_NULL, PETSC_NULL );

  // init Petsc-MPI communicator
  FemusInit mpinit ( argc, args, MPI_COMM_WORLD );

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;

  unsigned nx = static_cast<unsigned> ( floor ( pow ( 2.,/*11*/4 ) + 0.5 ) ); //Grid cell size = 3.90625 m
  nx += 3;

  double length = 10.; //2 * 1465700.;

  //mlMsh.GenerateCoarseBoxMesh ( nx, 0, 0, -length / 2, length / 2, 0., 0., 0., 0., EDGE3, "seventh" );
  mlMsh.GenerateCoarseBoxMesh ( nx, 0, 0, 0, length, 0., 0., 0., 0., EDGE3, "seventh" );

  //mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol ( &mlMsh );

  for ( unsigned i = 0; i < NumberOfLayers; i++ ) {
    char name[10];
    sprintf ( name, "h%d", i );
    mlSol.AddSolution ( name, DISCONTINUOUS_POLYNOMIAL, ZERO, 2 );
    sprintf ( name, "v%d", i );
    mlSol.AddSolution ( name, LAGRANGE, FIRST, 2 );
    sprintf ( name, "T%d", i );
    mlSol.AddSolution ( name, DISCONTINUOUS_POLYNOMIAL, ZERO, 2 );
    sprintf ( name, "HT%d", i );
    mlSol.AddSolution ( name, DISCONTINUOUS_POLYNOMIAL, ZERO, 2 );
  }

  mlSol.AddSolution ( "b", DISCONTINUOUS_POLYNOMIAL, ZERO, 1, false );

  mlSol.AddSolution ( "eta", DISCONTINUOUS_POLYNOMIAL, ZERO, 1, false );

  mlSol.Initialize ( "All" );

  for ( unsigned i = 0; i < NumberOfLayers; i++ ) {
    char name[10];
    sprintf ( name, "h%d", i );
    mlSol.Initialize ( name, InitalValueH );
  }

  for ( unsigned i = 0; i < NumberOfLayers; i++ ) {
    char name[10];
    sprintf ( name, "T%d", i );
    mlSol.Initialize ( name, InitalValueT );
  }

  for ( unsigned i = 0; i < NumberOfLayers; i++ ) {
    char name[10];
    sprintf ( name, "v%d", i );
    mlSol.Initialize ( name, InitalValueV );
  }

  mlSol.Initialize ( "b", InitalValueB );

  mlSol.AttachSetBoundaryConditionFunction ( SetBoundaryCondition );
  mlSol.GenerateBdc ( "All" );

  MultiLevelProblem ml_prob ( &mlSol );

  // ******* Add FEM system to the MultiLevel problem *******
  TransientLinearImplicitSystem& system = ml_prob.add_system < TransientLinearImplicitSystem > ( "SWt" );
  for ( unsigned i = 0; i < NumberOfLayers; i++ ) {
    char name[10];
    sprintf ( name, "HT%d", i );
    system.AddSolutionToSystemPDE ( name );
  }
  system.init();

  mlSol.SetWriter ( VTK );
  std::vector<std::string> print_vars;
  print_vars.push_back ( "All" );
  //mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write ( DEFAULT_OUTPUTDIR, "linear", print_vars, 0 );

  unsigned numberOfTimeSteps = 3209; //17h=1020 with dt=60, 17h=10200 with dt=6
  dt = 1./10.;
  bool implicitEuler = true;
  for ( unsigned i = 0; i < numberOfTimeSteps; i++ ) {
    if ( wave == true ) assembly = ( i == 0 ) ? true : false;
    system.CopySolutionToOldSolution();
    //ETD ( ml_prob );
    RK4 ( ml_prob, implicitEuler );
    mlSol.GetWriter()->Write ( DEFAULT_OUTPUTDIR, "linear", print_vars, ( i + 1 ) / 1 );
  }
  std::cout << " TOTAL TIME:\t" << \
            static_cast<double> ( clock() - start_time ) / CLOCKS_PER_SEC << std::endl;
  return 0;
}


void ETD ( MultiLevelProblem& ml_prob ) {

  const unsigned& NLayers = NumberOfLayers;

  adept::Stack& s = FemusInit::_adeptStack;

  if ( assembly ) s.continue_recording();
  else s.pause_recording();

  TransientLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<TransientLinearImplicitSystem> ( "SWt" ); // pointer to the linear implicit system named "Poisson"

  unsigned level = ml_prob._ml_msh->GetNumberOfLevels() - 1u;

  Mesh* msh = ml_prob._ml_msh->GetLevel ( level ); // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel ( level ); // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  SparseMatrix* KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  NumericVector* EPS = pdeSys->_EPS; // pointer to the global residual vector object in pdeSys (level)

  NumericVector* RES2;
  RES2 = NumericVector::build().release();
  RES2->init ( *RES );

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  //solution variable
  std::vector < unsigned > solIndexh ( NLayers );

  std::vector < unsigned > solIndexv ( NLayers );

  std::vector < unsigned > solIndexHT ( NLayers );
  std::vector < unsigned > solPdeIndexHT ( NLayers );

  std::vector < unsigned > solIndexT ( NLayers );

  vector< int > l2GMapRow; // local to global mapping
  vector< int > l2GMapColumn; // local to global mapping

  for ( unsigned i = 0; i < NLayers; i++ ) {
    char name[10];
    sprintf ( name, "h%d", i );
    solIndexh[i] = mlSol->GetIndex ( name ); // get the position of "hi" in the sol object

    sprintf ( name, "v%d", i );
    solIndexv[i] = mlSol->GetIndex ( name ); // get the position of "vi" in the sol object

    sprintf ( name, "HT%d", i );
    solIndexHT[i] = mlSol->GetIndex ( name ); // get the position of "Ti" in the sol object
    solPdeIndexHT[i] = mlPdeSys->GetSolPdeIndex ( name ); // get the position of "Ti" in the pdeSys object

    sprintf ( name, "T%d", i );
    solIndexT[i] = mlSol->GetIndex ( name ); // get the position of "Ti" in the sol object

  }

  unsigned solTypeh = mlSol->GetSolutionType ( solIndexh[0] ); // get the finite element type for "hi"
  unsigned solTypev = mlSol->GetSolutionType ( solIndexv[0] ); // get the finite element type for "vi"
  unsigned solTypeHT = mlSol->GetSolutionType ( solIndexHT[0] ); // get the finite element type for "Ti"

  if ( assembly ) KK->zero();
  RES->zero();

  MatSetOption ( ( static_cast<PetscMatrix*> ( KK ) )->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE );

  for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
    for ( unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++ ) {
      double valueT = ( *sol->_SolOld[solIndexT[k]] ) ( i );
      double valueH = ( *sol->_Sol[solIndexh[k]] ) ( i );

      double valueHT = valueT * valueH;

      sol->_Sol[solIndexHT[k]]->set ( i, valueHT );
    }
    sol->_Sol[solIndexHT[k]]->close();
  }

  std::vector < double > maxW ( NLayers, -1.e6 );
  maxW[0] = 0.;

  unsigned start = msh->_dofOffset[solTypeHT][iproc];
  unsigned end = msh->_dofOffset[solTypeHT][iproc + 1];
  for ( unsigned i =  start; i <  end; i++ ) {

    vector < double > solhm ( NLayers );
    vector < double > solh ( NLayers ); // local coordinates
    vector < double > solhp ( NLayers );
    vector < double > solvm ( NLayers ); // local coordinates
    vector < double > solvp ( NLayers ); // local coordinates
    vector < adept::adouble > solHTm ( NLayers ); // local coordinates
    vector < adept::adouble > solHT ( NLayers ); // local coordinates
    vector < adept::adouble > solHTp ( NLayers ); // local coordinates

    vector < adept::adouble > solHTmm ( NLayers ); // local coordinates
    vector < adept::adouble > solHTpp ( NLayers ); // local coordinates

    vector< adept::adouble > aResHT ( NLayers );

    unsigned bc1 = ( i == start ) ? 0 : 1;
    unsigned bc2 = ( i == end - 1 ) ? 0 : 1;

    unsigned bc3 = ( i > start + 1 ) ? 1 : 0;
    unsigned bc4 = ( i < end - 2 ) ? 1 : 0;

    l2GMapRow.resize ( NLayers );
    l2GMapColumn.resize ( ( 1 + bc1 + bc2 ) * NLayers );

    std::fill ( aResHT.begin(), aResHT.end(), 0 ); //set aRes to zero

    for ( unsigned j = 0; j < NLayers; j++ ) {

      solh[j] = ( *sol->_Sol[solIndexh[j]] ) ( i );
      solHT[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i );
      l2GMapRow[/*NLayers +*/ j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i );

      l2GMapColumn[/*NLayers +*/ j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i );

      solvm[j] = ( *sol->_Sol[solIndexv[j]] ) ( i );
      solvp[j] = ( *sol->_Sol[solIndexv[j]] ) ( i + 1 );

      if ( i > start ) {
        solhm[j] = ( *sol->_Sol[solIndexh[j]] ) ( i - 1 );
        solHTm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 1 );

        l2GMapColumn[NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 1 );

      }

      if ( i < end - 1 ) {
        solhp[j] = ( *sol->_Sol[solIndexh[j]] ) ( i + 1 );
        solHTp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 1 );

        l2GMapColumn[ ( 1 + bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 1 );
      }

//       if ( i > start + 1 ) {
//         solHTmm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 2 );
//         if (i == end - 1) l2GMapColumn[( 1 + bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//         else l2GMapColumn[( (1 + bc1) + bc3 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//       }
//
//       if ( i < end - 2 ) {
//         solHTpp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 2 );
//         l2GMapColumn[( (1 + bc1) + bc3 + bc4 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 2 );
//       }

    }

    if ( assembly ) s.new_recording();

    vector < double > x ( 2 ); // local coordinates
    for ( unsigned j = 0; j < 2; j++ ) {
      unsigned xDof  = msh->GetSolutionDof ( j, i, 2 ); // global to global mapping between coordinates node and coordinate dof
      x[j] = ( *msh->_topology->_Sol[0] ) ( xDof ); // global extraction and local storage for the element coordinates
    }
    double dx = x[1] - x[0];

    double b = 10.; //( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );

    double hTot = 0.;
    for ( unsigned k = 0; k < NLayers; k++ ) {
      hTot += solh[k];
    }

    std::vector < double > hALE ( NLayers, 0. );

    hALE[0] = hRest[0] + ( hTot - b );
    for ( unsigned k = 1; k < NLayers; k++ ) {
      hALE[k] = hRest[k];
    }

    std::vector < double > w ( NLayers + 1, 0. );
    w[0] = 1.;

    std::vector < double > zMid ( NLayers );
    for ( unsigned k = 0; k < NLayers; k++ ) {
      zMid[k] = -b + solh[k] / 2.;
      for ( unsigned i = k + 1; i < NLayers; i++ ) {
        zMid[k] += solh[i];
      }
    }

    std::vector < double > psi2 ( NLayers );
    for ( unsigned k = 0; k < NLayers; k++ ) {
      psi2[k] = ( ( zMid[k] + 10 ) * zMid[k] ) * ( ( zMid[k] + 10 ) * zMid[k] ) / ( 25 * 25 );
      //psi2[k] = (zMid[k] + 10.)*zMid[k]/(5*5);
    }

    for ( unsigned k = NLayers; k > 1; k-- ) {
      //w[k-1] = - (10 - 2*x[0])/(5*5)*psi2[k];
      w[k - 1] = 1.;
      if ( maxW[k - 1] < w[k - 1] ) {
        maxW[k - 1] = w[k - 1];
      }
    }


    for ( unsigned k = 0; k < NLayers; k++ ) {

      //BEGIN FIRST ORDER
      if ( i > start ) {
        if ( solvm[k] > 0 ) {
          aResHT[k] += solHTm[k] * solvm[k] / dx;
        }
        else {
          aResHT[k] += solHT[k] * solvm[k] / dx;
        }
      }
      if ( i < end - 1 ) {
        if ( solvp[k] > 0 ) {
          aResHT[k] -= solHT[k] * solvp[k] / dx; //first order upwind
        }
        else {
          aResHT[k] -= solHTp[k] * solvp[k] / dx; //first order upwind
        }
      }
      //END

      //BEGIN THIRD ORDER
//       if ( i > start ) {
//         aResHT[k] += 0.5 * ( solHTm[k].value() + solHT[k].value() ) * solvm[k] / dx;
//         if ( solvm[k] > 0 ) {
//           if ( i > start + 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHT[k].value() - 2.*solHTm[k].value() + solHTmm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//       }
//       if ( i < end - 1 ) {
//         aResHT[k] -= 0.5 * ( solHTp[k].value() + solHT[k].value() ) * solvp[k] / dx;
//         if ( solvp[k] > 0 ) {
//           if (i > start) {
//             aResHT[k] -= - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvp[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 2 ) {
//             aResHT[k] -= - 1. / 6. * ( solHTpp[k].value() - 2.*solHTp[k].value() + solHT[k].value() ) * solvp[k]  / dx;
//           }
//         }
//       }
      //END

      if ( k < NLayers - 1 ) {
        //aResHT[k] += w[k + 1] * 0.5 * ( solHT[k] / solh[k] + solHT[k + 1]/ solh[k + 1] );
        aResHT[k] += w[k + 1] * ( solHT[k + 1] / solh[k + 1] );
      }
      if ( k >= 0 ) {
        //aResHT[k] -= w[k] * 0.5 * ( solHT[k - 1]/ solh[k - 1] + solHT[k] / solh[k] );
        aResHT[k] -= w[k] * ( solHT[k] / solh[k] );
      }

      adept::adouble deltaZt = 0.;
      adept::adouble deltaZb = 0.;
      adept::adouble ht = 0.;
      adept::adouble hb = 0.;
      if ( k > 0 ) {
        ht = ( solhm[k - 1] + solhm[k] + solhp[k - 1] + solhp[k] ) / 4.;
        deltaZt = ( solHT[k - 1] - solHT[k] ) / ht;
      }
      else {
        ht = 0.5 * ( solhm[k] + solhp[k] );
        deltaZt = 0.* ( 0. - solHT[k] ) / ht;
      }
      if ( k < NLayers - 1 ) {
        hb = ( solhm[k] + solhm[k + 1] + solhp[k] + solhp[k + 1] ) / 4.;
        deltaZb = ( solHT[k] - solHT[k + 1] ) / hb;
      }
      else {
        hb = 0.5 * ( solhm[k] + solhp[k] );
        deltaZb = 0.* ( solHT[k] - 0. ) / hb;
      }

      //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAA"<<deltaZt - deltaZb<<std::endl;

      aResHT[k] += solhm[k] * k_v * ( deltaZt - deltaZb ) / ( ( ht + hb ) / 2. ); // vertical diffusion

//       //aResHT[k] += ((solhp[k] - solhm[k]) * k_v * (solHTp[k] - solHTm[k])) / (dx*dx); // horizontal diffusion
//       aResHT[k] += k_h * solh[k] * (solHTm[k] - solHT[k])/(dx*dx); // horizontal diffusion
//       aResHT[k] += k_h * solh[k] * (solHTp[k] - solHT[k])/(dx*dx); // horizontal diffusion

    }

    vector< double > Res ( NLayers ); // local redidual vector
    vector< double > solht ( NLayers ); // local redidual vector
    for ( unsigned k = 0; k < NLayers; k++ ) {
      Res[k] =  aResHT[k].value();
      solht[k] = solHT[k].value();
    }

    RES->add_vector_blocked ( Res, l2GMapRow );

    if ( assembly ) {
      s.dependent ( &aResHT[0], NLayers );

      // define the independent variables
      s.independent ( &solHT[0], NLayers );
      if ( i > start ) {
        s.independent ( &solHTm[0], NLayers );
      }
      if ( i < end - 1 ) {
        s.independent ( &solHTp[0], NLayers );
      }
      /*    if ( i > start + 1) {
        s.independent ( &solHTmm[0], NLayers );
        }
        if ( i < end - 2 ) {
          s.independent ( &solHTpp[0], NLayers );
        } */

      // get the jacobian matrix (ordered by row major )
      vector < double > Jac ( NLayers * NLayers * ( 1 + bc1 + bc2 ) );
      s.jacobian ( &Jac[0], true );

      //store K in the global matrix KK
      KK->add_matrix_blocked ( Jac, l2GMapRow, l2GMapColumn );

      s.clear_independents();
      s.clear_dependents();
    }
  }

  RES->close();
  if ( assembly ) KK->close();

  // printing of the max value of w in every layer
//   for ( unsigned k = 0; k < NLayers; k++ ) {
//     std::cout << "layer " << k << " " << maxW[k] << std::endl;
//   }

//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,900,900,&viewer);
//   PetscObjectSetName((PetscObject)viewer,"FSI matrix");
//   PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(),viewer);
//   double a;
//   std::cin>>a;

//  abort();
  //std::cout << "dt = " << dt << std::endl;
  //std::cout << "first stage " << std::endl;

  MFN mfn;
  Mat A = ( static_cast<PetscMatrix*> ( KK ) )->mat();
  FN f;

  Vec v = ( static_cast< PetscVector* > ( RES ) )->vec();
  Vec y = ( static_cast< PetscVector* > ( EPS ) )->vec();

  MFNCreate ( PETSC_COMM_WORLD, &mfn );

  MFNSetOperator ( mfn, A );
  MFNGetFN ( mfn, &f );

  FNPhiSetIndex ( f, 1 );
  FNSetType ( f, FNPHI );
// FNView(f,PETSC_VIEWER_STDOUT_WORLD);

  FNSetScale ( f, dt, dt );
  MFNSetFromOptions ( mfn );

  MFNSolve ( mfn, v, y );
  MFNDestroy ( &mfn );

  //VecView(y,PETSC_VIEWER_STDOUT_WORLD);
  
  //std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
  
  sol->UpdateSol ( mlPdeSys->GetSolPdeIndex(), EPS, pdeSys->KKoffset );

  if ( twostage == true ) {

    std::cout << "second stage " << std::endl;

    RES2->zero();

    for ( unsigned i =  start; i <  end; i++ ) {

      vector < double > solhm ( NLayers );
      vector < double > solh ( NLayers ); // local coordinates
      vector < double > solhp ( NLayers );
      vector < double > solvm ( NLayers ); // local coordinates
      vector < double > solvp ( NLayers ); // local coordinates
      vector < double > solHTm ( NLayers ); // local coordinates
      vector < double > solHT ( NLayers ); // local coordinates
      vector < double > solHTp ( NLayers ); // local coordinates

      vector < double > solHTmm ( NLayers ); // local coordinates
      vector < double > solHTpp ( NLayers ); // local coordinates

      vector< double > aResHT ( NLayers, 0. );

      unsigned bc1 = ( i == start ) ? 0 : 1;
      unsigned bc2 = ( i == end - 1 ) ? 0 : 1;

      unsigned bc3 = ( i > start + 1 ) ? 1 : 0;
      unsigned bc4 = ( i < end - 2 ) ? 1 : 0;

      l2GMapRow.resize ( NLayers );
      l2GMapColumn.resize ( ( 1 + bc1 + bc2 ) * NLayers );

      //std::fill ( aResHT.begin(), aResHT.end(), 0 ); //set aRes to zero

      for ( unsigned j = 0; j < NLayers; j++ ) {

        solh[j] = ( *sol->_Sol[solIndexh[j]] ) ( i );
        solHT[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i );
        l2GMapRow[/*NLayers +*/ j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i );

        l2GMapColumn[/*NLayers +*/ j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i );

        solvm[j] = ( *sol->_Sol[solIndexv[j]] ) ( i );
        solvp[j] = ( *sol->_Sol[solIndexv[j]] ) ( i + 1 );


        if ( i > start ) {
          solhm[j] = ( *sol->_Sol[solIndexh[j]] ) ( i - 1 );
          solHTm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 1 );

          l2GMapColumn[NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 1 );

        }

        if ( i < end - 1 ) {
          solhp[j] = ( *sol->_Sol[solIndexh[j]] ) ( i + 1 );
          solHTp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 1 );

          l2GMapColumn[ ( 1 + bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 1 );
        }

//       if ( i > start + 1 ) {
//         solHTmm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 2 );
//         if (i == end - 1) l2GMapColumn[( 1 + bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//         else l2GMapColumn[( (1 + bc1) + bc3 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//       }
//
//       if ( i < end - 2 ) {
//         solHTpp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 2 );
//         l2GMapColumn[( (1 + bc1) + bc3 + bc4 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 2 );
//       }

      }

      // s.new_recording();

      vector < double > x ( 2 ); // local coordinates
      for ( unsigned j = 0; j < 2; j++ ) {
        unsigned xDof  = msh->GetSolutionDof ( j, i, 2 ); // global to global mapping between coordinates node and coordinate dof
        x[j] = ( *msh->_topology->_Sol[0] ) ( xDof ); // global extraction and local storage for the element coordinates
      }
      double dx = x[1] - x[0];

      double b = 10.; //( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );

      double hTot = 0.;
      for ( unsigned k = 0; k < NLayers; k++ ) {
        hTot += solh[k]/*.value()*/;
      }

      std::vector < double > hALE ( NLayers, 0. );

      hALE[0] = hRest[0] + ( hTot - b );
      for ( unsigned k = 1; k < NLayers; k++ ) {
        hALE[k] = hRest[k];
      }

      std::vector < double > w ( NLayers + 1, 0. );
      w[0] = 1.;

      std::vector < double > zMid ( NLayers );
      for ( unsigned k = 0; k < NLayers; k++ ) {
        zMid[k] = -b + solh[k] / 2.;
        for ( unsigned i = k + 1; i < NLayers; i++ ) {
          zMid[k] += solh[i];
        }
      }

      std::vector < double > psi2 ( NLayers );
      for ( unsigned k = 0; k < NLayers; k++ ) {
        psi2[k] = ( ( zMid[k] + 10 ) * zMid[k] ) * ( ( zMid[k] + 10 ) * zMid[k] ) / ( 25 * 25 );
        //psi2[k] = (zMid[k] + 10.)*zMid[k]/(5*5);
      }

      for ( unsigned k = NLayers; k > 1; k-- ) {
        //w[k-1] = - (10 - 2*x[0])/(5*5)*psi2[k];
        w[k - 1] = 1.;
        if ( maxW[k - 1] < w[k - 1] ) {
          maxW[k - 1] = w[k - 1];
        }
      }


      for ( unsigned k = 0; k < NLayers; k++ ) {

        //BEGIN FIRST ORDER
        if ( i > start ) {
          //aResHT[k] += 0.5 * (solHTm[k] + solHT[k]) * solvm[k]  / dx; //second order
          if ( solvm[k] > 0 ) {
            aResHT[k] += solHTm[k] * solvm[k] / dx;
          }
          else {
            aResHT[k] += solHT[k] * solvm[k] / dx;
          }
        }
        if ( i < end - 1 ) {
          //aResHT[k] -= 0.5 * (solHT[k] + solHTp[k]) * solvp[k]  / dx; //second order
          if ( solvp[k] > 0 ) {
            aResHT[k] -= solHT[k] * solvp[k] / dx; //first order upwind
          }
          else {
            aResHT[k] -= solHTp[k] * solvp[k] / dx; //first order upwind
          }
        }
//       else{
//         aResHT[k] -= solHT[k] /*.value()*/  / dx; //first order upwind
//       }
        //END

        //BEGIN THIRD ORDER
//       if ( i > start ) {
//         aResHT[k] += 0.5 * ( solHTm[k].value() + solHT[k].value() ) * solvm[k] / dx;
//         if ( solvm[k] > 0 ) {
//           if ( i > start + 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHT[k].value() - 2.*solHTm[k].value() + solHTmm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//       }
//       if ( i < end - 1 ) {
//         aResHT[k] -= 0.5 * ( solHTp[k].value() + solHT[k].value() ) * solvp[k] / dx;
//         if ( solvp[k] > 0 ) {
//           if (i > start) {
//             aResHT[k] -= - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvp[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 2 ) {
//             aResHT[k] -= - 1. / 6. * ( solHTpp[k].value() - 2.*solHTp[k].value() + solHT[k].value() ) * solvp[k]  / dx;
//           }
//         }
//       }
        //END

        if ( k < NLayers - 1 ) {
          //aResHT[k] += w[k + 1] * 0.5 * ( solHT[k] / solh[k] + solHT[k + 1]/ solh[k + 1] );
          aResHT[k] += w[k + 1] * ( solHT[k + 1] / solh[k + 1] );
        }
        if ( k >= 0 ) {
          //aResHT[k] -= w[k] * 0.5 * ( solHT[k - 1]/ solh[k - 1] + solHT[k] / solh[k] );
          aResHT[k] -= w[k] * ( solHT[k] / solh[k] );
        }

        double deltaZt = 0.;
        double deltaZb = 0.;
        double ht = 0.;
        double hb = 0.;
        if ( k > 0 ) {
          ht = ( solhm[k - 1] + solhm[k] + solhp[k - 1] + solhp[k] ) / 4.;
          deltaZt = ( solHT[k - 1] - solHT[k] ) / ht;
        }
        else {
          ht = 0.5 * ( solhm[k] + solhp[k] );
          deltaZt = 0.* ( 0. - solHT[k] ) / ht;
        }
        if ( k < NLayers - 1 ) {
          hb = ( solhm[k] + solhm[k + 1] + solhp[k] + solhp[k + 1] ) / 4.;
          deltaZb = ( solHT[k] - solHT[k + 1] ) / hb;
        }
        else {
          hb = 0.5 * ( solhm[k] + solhp[k] );
          deltaZb = 0.* ( solHT[k] - 0. ) / hb;
        }

        //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAA"<<deltaZt - deltaZb<<std::endl;

        aResHT[k] += solhm[k] * k_v * (deltaZt - deltaZb) / ( (ht + hb) / 2. ); // vertical diffusion
//
//       //aResHT[k] += ((solhp[k] - solhm[k]) * k_v * (solHTp[k] - solHTm[k])) / (dx*dx); // horizontal diffusion
//       aResHT[k] += k_h * solh[k] * (solHTm[k] - solHT[k])/(dx*dx); // horizontal diffusion
//       aResHT[k] += k_h * solh[k] * (solHTp[k] - solHT[k])/(dx*dx); // horizontal diffusion

      }

      RES2->add_vector_blocked ( aResHT, l2GMapRow );

    }
    RES2->close();

    //BEGIN ASSEMBLY: R2 = RES2 - RES - KK * EPS = RESnew - RESold - KK * (Vnew - Vold) = (ResNew - KK * Vnew) - (ResOld - KK * Vold) = 0 - 0 ;
    RES2->scale ( -1. );
    RES2->add ( *RES );
    RES2->add_vector ( *EPS, *KK );
    RES2->scale ( -1. );    
    //END ASSEMBLY R2

    EPS->zero();
    
    //SLEPC
    MFN mfn;
    Mat A2 = ( static_cast<PetscMatrix*> ( KK ) )->mat();
    FN f2;

    //std::cout << "dt = " << dt << std::endl;

    Vec v2 = ( static_cast< PetscVector* > ( RES2 ) )->vec();
    
    //VecView(v2,PETSC_VIEWER_STDOUT_WORLD);
    
    Vec y2 = ( static_cast< PetscVector* > ( EPS ) )->vec();

    MFNCreate ( PETSC_COMM_WORLD, &mfn );

    MFNSetOperator ( mfn, A2 );
    MFNGetFN ( mfn, &f2 );

    FNPhiSetIndex ( f2, 2 );
    FNSetType ( f2, FNPHI );
    //FNView(f,PETSC_VIEWER_STDOUT_WORLD);

    FNSetScale ( f2, dt, dt );
    MFNSetFromOptions ( mfn );

    MFNSolve ( mfn, v2, y2 );
    MFNDestroy ( &mfn );

    //VecView(y2,PETSC_VIEWER_STDOUT_WORLD);
    
    sol->UpdateSol ( mlPdeSys->GetSolPdeIndex(), EPS, pdeSys->KKoffset );

  }

  //PARAVIEW
  unsigned solIndexeta = mlSol->GetIndex ( "eta" );
  unsigned solIndexb = mlSol->GetIndex ( "b" );
  sol->_Sol[solIndexeta]->zero();
  for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
    sol->_Sol[solIndexeta]->add ( *sol->_Sol[solIndexh[k]] );
  }
  sol->_Sol[solIndexeta]->add ( -1, *sol->_Sol[solIndexb] );

  for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
    for ( unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++ ) {
      double valueHT = ( *sol->_Sol[solIndexHT[k]] ) ( i );
      double valueH = ( *sol->_Sol[solIndexh[k]] ) ( i );

      double valueT = valueHT / valueH;
      //if ( i == 10 ) std::cout << "temperature " << valueT << std::endl;
      //if (i == 0) valueT = 0.;
      //if (i == msh->_dofOffset[solTypeHT][iproc + 1] - 1 ) valueT = 0.;


      sol->_Sol[solIndexT[k]]->set ( i, valueT );
    }

    sol->_Sol[solIndexT[k]]->close();

  }

  delete RES2;
}


void RK4 ( MultiLevelProblem& ml_prob, const bool & implicitEuler ) {

  const unsigned& NLayers = NumberOfLayers;

  adept::Stack& s = FemusInit::_adeptStack;

  TransientLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<TransientLinearImplicitSystem> ( "SWt" ); // pointer to the linear implicit system named "Poisson"

  unsigned level = ml_prob._ml_msh->GetNumberOfLevels() - 1u;

  Mesh* msh = ml_prob._ml_msh->GetLevel ( level ); // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel ( level ); // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  //solution variable
  std::vector < unsigned > solIndexh ( NLayers );

  std::vector < unsigned > solIndexv ( NLayers );

  std::vector < unsigned > solIndexHT ( NLayers );
  std::vector < unsigned > solPdeIndexHT ( NLayers );

  std::vector < unsigned > solIndexT ( NLayers );

  vector< int > l2GMapRow; // local to global mapping
  vector< int > l2GMapColumn; // local to global mapping

  for ( unsigned i = 0; i < NLayers; i++ ) {
    char name[10];
    sprintf ( name, "h%d", i );
    solIndexh[i] = mlSol->GetIndex ( name ); // get the position of "hi" in the sol object

    sprintf ( name, "v%d", i );
    solIndexv[i] = mlSol->GetIndex ( name ); // get the position of "vi" in the sol object

    sprintf ( name, "HT%d", i );
    solIndexHT[i] = mlSol->GetIndex ( name ); // get the position of "Ti" in the sol object
    solPdeIndexHT[i] = mlPdeSys->GetSolPdeIndex ( name ); // get the position of "Ti" in the pdeSys object

    sprintf ( name, "T%d", i );
    solIndexT[i] = mlSol->GetIndex ( name ); // get the position of "Ti" in the sol object

  }

  unsigned solTypeh = mlSol->GetSolutionType ( solIndexh[0] ); // get the finite element type for "hi"
  unsigned solTypev = mlSol->GetSolutionType ( solIndexv[0] ); // get the finite element type for "vi"
  unsigned solTypeHT = mlSol->GetSolutionType ( solIndexHT[0] ); // get the finite element type for "Ti"

  for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
    for ( unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++ ) {
      double valueT = ( *sol->_SolOld[solIndexT[k]] ) ( i );
      double valueH = ( *sol->_Sol[solIndexh[k]] ) ( i );

      double valueHT = valueT * valueH;

      sol->_Sol[solIndexHT[k]]->set ( i, valueHT );
    }
    sol->_Sol[solIndexHT[k]]->close();
  }

  std::vector < double > maxW ( NLayers, -1.e6 );
  maxW[0] = 0.;

  unsigned start = msh->_dofOffset[solTypeHT][iproc];
  unsigned end = msh->_dofOffset[solTypeHT][iproc + 1];
  for ( unsigned i =  start; i <  end; i++ ) {

    vector < double > solhm ( NLayers );
    vector < double > solh ( NLayers ); // local coordinates
    vector < double > solhp ( NLayers );
    vector < double > solvm ( NLayers ); // local coordinates
    vector < double > solvp ( NLayers ); // local coordinates
    vector < double > solHTm ( NLayers ); // local coordinates
    vector < double > solHT ( NLayers ); // local coordinates
    vector < double > solHTp ( NLayers ); // local coordinates

    vector < adept::adouble > solHTmm ( NLayers ); // local coordinates
    vector < adept::adouble > solHTpp ( NLayers ); // local coordinates

    vector< adept::adouble > aResHT ( NLayers );

    unsigned bc1 = ( i == start ) ? 0 : 1;
    unsigned bc2 = ( i == end - 1 ) ? 0 : 1;

    unsigned bc3 = ( i > start + 1 ) ? 1 : 0;
    unsigned bc4 = ( i < end - 2 ) ? 1 : 0;

    l2GMapRow.resize ( NLayers );
    l2GMapColumn.resize ( ( 1 + bc1 + bc2 ) * NLayers );

    std::fill ( aResHT.begin(), aResHT.end(), 0 ); //set aRes to zero

    for ( unsigned j = 0; j < NLayers; j++ ) {

      solh[j] = ( *sol->_Sol[solIndexh[j]] ) ( i );
      solHT[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i );
      l2GMapRow[/*NLayers +*/ j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i );

      l2GMapColumn[/*NLayers +*/ j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i );

      solvm[j] = ( *sol->_Sol[solIndexv[j]] ) ( i );
      solvp[j] = ( *sol->_Sol[solIndexv[j]] ) ( i + 1 );

      if ( i > start ) {
        solhm[j] = ( *sol->_Sol[solIndexh[j]] ) ( i - 1 );
        solHTm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 1 );

        l2GMapColumn[NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 1 );

      }

      if ( i < end - 1 ) {
        solhp[j] = ( *sol->_Sol[solIndexh[j]] ) ( i + 1 );
        solHTp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 1 );

        l2GMapColumn[ ( 1 + bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 1 );
      }

//       if ( i > start + 1 ) {
//         solHTmm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 2 );
//         if (i == end - 1) l2GMapColumn[( 1 + bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//         else l2GMapColumn[( (1 + bc1) + bc3 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//       }
//
//       if ( i < end - 2 ) {
//         solHTpp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 2 );
//         l2GMapColumn[( (1 + bc1) + bc3 + bc4 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 2 );
//       }

    }

//   s.new_recording();

    vector < double > x ( 2 ); // local coordinates
    for ( unsigned j = 0; j < 2; j++ ) {
      unsigned xDof  = msh->GetSolutionDof ( j, i, 2 ); // global to global mapping between coordinates node and coordinate dof
      x[j] = ( *msh->_topology->_Sol[0] ) ( xDof ); // global extraction and local storage for the element coordinates
    }
    double dx = x[1] - x[0];

    double b = 10.; //( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );

    double hTot = 0.;
    for ( unsigned k = 0; k < NLayers; k++ ) {
      hTot += solh[k];
    }

    std::vector < double > hALE ( NLayers, 0. );

    hALE[0] = hRest[0] + ( hTot - b );
    for ( unsigned k = 1; k < NLayers; k++ ) {
      hALE[k] = hRest[k];
    }

    std::vector < double > w ( NLayers + 1, 0. );
    w[0] = 1.;

    std::vector < double > zMid ( NLayers );
    for ( unsigned k = 0; k < NLayers; k++ ) {
      zMid[k] = -b + solh[k] / 2.;
      for ( unsigned i = k + 1; i < NLayers; i++ ) {
        zMid[k] += solh[i];
      }
    }

    std::vector < double > psi2 ( NLayers );
    for ( unsigned k = 0; k < NLayers; k++ ) {
      psi2[k] = ( ( zMid[k] + 10 ) * zMid[k] ) * ( ( zMid[k] + 10 ) * zMid[k] ) / ( 25 * 25 );
      //psi2[k] = (zMid[k] + 10.)*zMid[k]/(5*5);
    }

    for ( unsigned k = NLayers; k > 1; k-- ) {
      //w[k-1] = - (10 - 2*x[0])/(5*5)*psi2[k];
      w[k - 1] = 1.;
      if ( maxW[k - 1] < w[k - 1] ) {
        maxW[k - 1] = w[k - 1];
      }
    }


    std::vector < double > k1_RK ( NLayers, 0. );
    std::vector < double > k2_RK ( NLayers, 0. );
    std::vector < double > k3_RK ( NLayers, 0. );
    std::vector < double > k4_RK ( NLayers, 0. );

    for ( unsigned RK_step = 0; RK_step < 4; RK_step++ ) {
      for ( unsigned k = 0; k < NLayers; k++ ) {
        double LHS = 0.;
        double addition = 0.;
        if ( RK_step == 1 ) {
          addition = k1_RK[k] * 0.5;
        }
        else if ( RK_step == 2 ) {
          addition = k2_RK[k] * 0.5;
        }

        else if ( RK_step == 3 ) {
          addition = k3_RK[k];
        }

        //BEGIN FIRST ORDER
        if ( i > start ) {
          if ( solvm[k] > 0 ) {
            LHS += ( solHTm[k] + addition ) * solvm[k] / dx;
          }
          else {
            LHS += ( solHT[k] + addition ) * solvm[k] / dx;
          }
        }
        if ( i < end - 1 ) {
          if ( solvp[k] > 0 ) {
            LHS -= ( solHT[k] + addition ) * solvp[k] / dx; //first order upwind
          }
          else {
            LHS -= ( solHTp[k] + addition ) * solvp[k] / dx; //first order upwind
          }
        }
        //END

        //BEGIN THIRD ORDER
//       if ( i > start ) {
//         aResHT[k] += 0.5 * ( solHTm[k].value() + solHT[k].value() ) * solvm[k] / dx;
//         if ( solvm[k] > 0 ) {
//           if ( i > start + 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHT[k].value() - 2.*solHTm[k].value() + solHTmm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//       }
//       if ( i < end - 1 ) {
//         aResHT[k] -= 0.5 * ( solHTp[k].value() + solHT[k].value() ) * solvp[k] / dx;
//         if ( solvp[k] > 0 ) {
//           if (i > start) {
//             aResHT[k] -= - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvp[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 2 ) {
//             aResHT[k] -= - 1. / 6. * ( solHTpp[k].value() - 2.*solHTp[k].value() + solHT[k].value() ) * solvp[k]  / dx;
//           }
//         }
//       }
        //END

        if ( k < NLayers - 1 ) {
          //aResHT[k] += w[k + 1] * 0.5 * ( solHT[k] / solh[k] + solHT[k + 1]/ solh[k + 1] );
          LHS += w[k + 1] * ( ( solHT[k + 1] + addition ) / solh[k + 1] );
        }
        if ( k >= 0 ) {
          //aResHT[k] -= w[k] * 0.5 * ( solHT[k - 1]/ solh[k - 1] + solHT[k] / solh[k] );
          LHS -= w[k] * ( ( solHT[k] + addition ) / solh[k] );
        }

        if ( RK_step == 0 ) {
          k1_RK[k] = LHS * dt;
        }
        else if ( RK_step == 1 ) {
          k2_RK[k] = LHS * dt;
        }
        else if ( RK_step == 2 ) {
          k3_RK[k] = LHS * dt;
        }
        else {
          k4_RK[k] = LHS * dt;
        }

      }

    }

    for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
      double valueHT = solHT[k] + 1. / 6. * ( k1_RK[k] + 2.*k2_RK[k] + 2.*k3_RK[k] + k4_RK[k] );
      sol->_Sol[solIndexHT[k]]->set ( i, valueHT );
      sol->_Sol[solIndexHT[k]]->close();
    }


    if ( implicitEuler == false ) {
      //BEGIN forward Euler for vertical diffusion
      std::vector < double > vert_diff ( NLayers, 0. );
      for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
        double deltaZt = 0.;
        double deltaZb = 0.;
        double ht = 0.;
        double hb = 0.;
        if ( k > 0 ) {
          ht = ( solhm[k - 1] + solhm[k] + solhp[k - 1] + solhp[k] ) / 4.;
          deltaZt = ( solHT[k - 1] - solHT[k] ) / ht;
        }
        else {
          ht = 0.5 * ( solhm[k] + solhp[k] );
          deltaZt = 0.* ( 0. - solHT[k] ) / ht;
        }
        if ( k < NLayers - 1 ) {
          hb = ( solhm[k] + solhm[k + 1] + solhp[k] + solhp[k + 1] ) / 4.;
          deltaZb = ( solHT[k] - solHT[k + 1] ) / hb;
        }
        else {
          hb = 0.5 * ( solhm[k] + solhp[k] );
          deltaZb = 0.* ( solHT[k] - 0. ) / hb;
        }
        //std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAA" << deltaZt - deltaZb << std::endl;
        vert_diff[k] = solhm[k] * k_v * ( deltaZt - deltaZb ) / ( ( ht + hb ) / 2. ); // vertical diffusion
      }

      for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
        double valueHT = ( *sol->_Sol[solIndexHT[k]] ) ( i );
        double valueH = ( *sol->_Sol[solIndexh[k]] ) ( i );

        double valueT = valueHT / valueH;
        valueT = valueT + dt * vert_diff[k];

        sol->_Sol[solIndexT[k]]->set ( i, valueT );
        sol->_Sol[solIndexT[k]]->close();
      }
      //END
    }

    else if ( implicitEuler == true ) {

      std::vector <double> Trhs (NLayers, 0.);
      std::vector < std::vector < double > > sysMatrix ( NLayers );
      
      for ( unsigned k = 0; k < NLayers; k++ ) {
        sysMatrix[k].assign( NLayers, 0. );

        double A = 0.;
	double C = 0.;
	double ht = 0.;
        double hb = 0.;

	if ( k > 0 ) {
          ht = ( solhm[k - 1] + solhm[k] + solhp[k - 1] + solhp[k] ) / 4.;
	  A = solhm[k] * k_v / ht ;
        }
        else {
          ht = 0.5 * ( solhm[k] + solhp[k] );
        }
        if ( k < NLayers - 1 ) {
          hb = ( solhm[k] + solhm[k + 1] + solhp[k] + solhp[k + 1] ) / 4.;
	  C = solhm[k] * k_v / hb;
	  C /= (ht + hb) * 0.5 ;
	  if(k > 0){
	    A /= (ht + hb) * 0.5 ;
	  }
        }
        else {
          hb = 0.5 * ( solhm[k] + solhp[k] );
	  A /= (ht + hb) * 0.5 ;
        }
	
	double B = 1. - A - C;
	
	sysMatrix[k][k] = B;
	if(k > 0) sysMatrix[k][k-1] = A;
	if(k < NLayers - 1) sysMatrix[k][k+1] = C;
	
	Trhs[k] = ( *sol->_Sol[solIndexHT[k]] ) ( i )/( *sol->_Sol[solIndexh[k]] ) ( i );
	
      }
           
      //risolvere il sistema Nlayer X Nlayer
      KSP                solver;
      Mat                triDiagA;
      Vec                b,x;
      PetscInt           k,j,nlayers;
      PetscErrorCode     ierr;
      
      nlayers = static_cast<PetscInt> (NLayers);
      ierr = VecCreate(PETSC_COMM_WORLD, &x);
      ierr = VecSetSizes(x, PETSC_DECIDE, nlayers);
      ierr = VecSetFromOptions(x);
      ierr = VecDuplicate(x,&b);
      ierr = MatCreate(PETSC_COMM_WORLD,&triDiagA);
      ierr = MatSetSizes(triDiagA,PETSC_DECIDE,PETSC_DECIDE, nlayers, nlayers);
      ierr = MatSetFromOptions(triDiagA);
      ierr = MatSetUp(triDiagA);
      
      for (k=0; k<nlayers; k++) {
        MatSetValues(triDiagA,1,&k,1,&k,&sysMatrix[k][k],INSERT_VALUES);
	VecSetValues(b, 1, &k, &Trhs[k], INSERT_VALUES );
	if(k>0){
	  j=k-1;
	  MatSetValues(triDiagA,1,&k,1,&j,& sysMatrix[k-1][k],INSERT_VALUES);
	}
	if(k<nlayers-1){
	  j=k+1;
	  MatSetValues(triDiagA,1,&k,1,&j,&sysMatrix[k][k+1],INSERT_VALUES);
	}
      }
      
      ierr = MatAssemblyBegin(triDiagA,MAT_FINAL_ASSEMBLY);
      ierr = MatAssemblyEnd(triDiagA,MAT_FINAL_ASSEMBLY);
      
      
//       PetscViewer    viewer;
//       PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,900,900,&viewer);
//       PetscObjectSetName((PetscObject)viewer,"implicit Euler matrix");
//       PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
//       MatView(triDiagA,viewer);
//       double a;
//       std::cin>>a;
      
      
      ierr = KSPCreate(PETSC_COMM_WORLD, &solver);
      ierr = KSPSetOperators(solver,triDiagA,triDiagA);
      ierr = KSPSetType(solver, KSPRICHARDSON);
      ierr = KSPSolve(solver, b, x);
      
      //1. aggiornare solT con x 
      for (k = 0; k < nlayers; k++ ) {
	PetscScalar valueT = 0.;
	ierr = VecGetValues(x, 1, &k, &valueT);
        sol->_Sol[solIndexT[k]]->set ( i, valueT );
        sol->_Sol[solIndexT[k]]->close();
	
        //2. aggiornare solHT
        double valueH = ( *sol->_Sol[solIndexh[k]] ) ( i );
        double valueHT = valueH * valueT;
        sol->_Sol[solIndexHT[k]]->set ( i, valueHT );
        sol->_Sol[solIndexHT[k]]->close();
      }
      //3. checkare che runni, fare pulizia e scegliere LU solver
    }

  }

  //BEGIN no vertical diffusion
//   for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
//     for ( unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++ ) {
//       double valueHT = ( *sol->_Sol[solIndexHT[k]] ) ( i );
//       double valueH = ( *sol->_Sol[solIndexh[k]] ) ( i );
//
//       double valueT = valueHT / valueH;
//
//       sol->_Sol[solIndexT[k]]->set ( i, valueT );
//     }
//
//     sol->_Sol[solIndexT[k]]->close();
//
//   }
  //END

//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,900,900,&viewer);
//   PetscObjectSetName((PetscObject)viewer,"FSI matrix");
//   PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(),viewer);
//   double a;
//   std::cin>>a;

//  abort();

//   unsigned solIndexeta = mlSol->GetIndex ( "eta" );
//   unsigned solIndexb = mlSol->GetIndex ( "b" );
//   sol->_Sol[solIndexeta]->zero();
//   for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
//     sol->_Sol[solIndexeta]->add ( *sol->_Sol[solIndexh[k]] );
//   }
//   sol->_Sol[solIndexeta]->add ( -1, *sol->_Sol[solIndexb] );


}




