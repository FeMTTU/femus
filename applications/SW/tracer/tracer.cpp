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

double dt = 0.01 /*1./10.*/; //= dx / maxWaveSpeed * 0.85;
unsigned counter = 0;

double k_v = 0.0001;

double pi = acos ( -1. );
double k_h = 1/(10*pi);

const unsigned NumberOfLayers = 5;

const double hRest[10]={1,1,1,1,1};
//const double hRest[10]={2,2,2,2,2,2,2,2,2,2};
//const double hRest[20] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

double InitalValueV ( const std::vector < double >& x ) {
  return 1.;
}

double InitalValueH ( const std::vector < double >& x ) {
  return hRest[0];
}


double InitalValueT ( const std::vector < double >& x ) {
  double pi = acos ( -1. );
//   return 17.5 + 25/pi * atan(x[0]/100.);
//   if ( x[0] < 0 ) return 5;
//   else return 30;
  return (- sin(pi*x[0]));
}


double InitalValueB ( const std::vector < double >& x ) {
  return 5; //( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );
}


bool SetBoundaryCondition ( const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time ) {
  bool dirichlet = false;
  if(!strcmp(SolName, "HT")) {
    if ( facename == 1 || facename == 2 ) {
      dirichlet = true;
      value = 0.;
    }
  }
  if(!strcmp(SolName, "T")) {
    if ( facename == 1 || facename == 2 ) {
      dirichlet = true;
      value = 0.;
    }
  }
   return dirichlet;
}


void ETD2 ( MultiLevelProblem& ml_prob );

void RK ( MultiLevelProblem& ml_prob, const unsigned & numberOfTimeSteps );


int main ( int argc, char** args ) {

  SlepcInitialize ( &argc, &args, PETSC_NULL, PETSC_NULL );

  // init Petsc-MPI communicator
  FemusInit mpinit ( argc, args, MPI_COMM_WORLD );

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;

  unsigned nx = static_cast<unsigned> ( floor ( pow ( 2.,/*11*/6 ) + 0.5 ) ); //Grid cell size = 0.03125 m

  double length = 2; //2 * 1465700.;

  mlMsh.GenerateCoarseBoxMesh ( nx, 0, 0, -length / 2, length / 2, 0., 0., 0., 0., EDGE3, "seventh" );

  //mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol ( &mlMsh );

  for ( unsigned i = 0; i < NumberOfLayers; i++ ) {
    char name[10];
    sprintf ( name, "h%d", i );
    mlSol.AddSolution ( name, DISCONTINUOUS_POLYNOMIAL, ZERO, 2);
    sprintf ( name, "v%d", i );
    mlSol.AddSolution ( name, LAGRANGE, FIRST, 2 );
    sprintf ( name, "T%d", i );
    mlSol.AddSolution ( name, DISCONTINUOUS_POLYNOMIAL, ZERO, 2);
    sprintf ( name, "HT%d", i );
    mlSol.AddSolution ( name, DISCONTINUOUS_POLYNOMIAL, ZERO, 2);
  }

  mlSol.AddSolution ( "b", DISCONTINUOUS_POLYNOMIAL, ZERO, 1, false );

  mlSol.AddSolution ( "eta", DISCONTINUOUS_POLYNOMIAL, ZERO, 1, false );

  mlSol.Initialize ( "All" );

  for ( unsigned i = 0; i < NumberOfLayers; i++ ) {
    char name[10];
    sprintf ( name, "v%d", i );
    mlSol.Initialize ( name, InitalValueV );
  }

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

  mlSol.Initialize ( "b", InitalValueB );

  mlSol.AttachSetBoundaryConditionFunction ( SetBoundaryCondition );
  mlSol.GenerateBdc ( "All" );

  MultiLevelProblem ml_prob ( &mlSol );

  // ******* Add FEM system to the MultiLevel problem *******
  //TransientLinearImplicitSystem& system = ml_prob.add_system < TransientLinearImplicitSystem > ( "SWhv" );
  TransientLinearImplicitSystem& system2 = ml_prob.add_system < TransientLinearImplicitSystem > ( "SWt" );
  for ( unsigned i = 0; i < NumberOfLayers; i++ ) {
    char name[10];
//     sprintf ( name, "h%d", i );
//     system.AddSolutionToSystemPDE ( name );
//     sprintf ( name, "v%d", i );
//     system.AddSolutionToSystemPDE ( name );
    sprintf ( name, "HT%d", i );
    system2.AddSolutionToSystemPDE ( name );
  }
  //system.init();
  system2.init();
  
  mlSol.SetWriter ( VTK );
  std::vector<std::string> print_vars;
  print_vars.push_back ( "All" );
  //mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write ( DEFAULT_OUTPUTDIR, "linear", print_vars, 0 );

  unsigned numberOfTimeSteps = 5000; //17h=1020 with dt=60, 17h=10200 with dt=6
  dt = 0.02;
  for ( unsigned i = 0; i < numberOfTimeSteps; i++ ) {
     system2.CopySolutionToOldSolution();
//     dt = 60.;
//     ETD ( ml_prob );
//     dt = 60.;
    //ETD2 ( ml_prob );
    RK ( ml_prob, numberOfTimeSteps );
    mlSol.GetWriter()->Write ( DEFAULT_OUTPUTDIR, "linear", print_vars, ( i + 1 ) / 1 );
    counter++;
  }
  return 0;
}


void ETD2 ( MultiLevelProblem& ml_prob ) {

  const unsigned& NLayers = NumberOfLayers;

  adept::Stack& s = FemusInit::_adeptStack;

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

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  //solution variable
  std::vector < unsigned > solIndexh ( NLayers );
  //std::vector < unsigned > solPdeIndexh ( NLayers );

  std::vector < unsigned > solIndexv ( NLayers );
  //std::vector < unsigned > solPdeIndexv ( NLayers );

  std::vector < unsigned > solIndexHT ( NLayers );
  std::vector < unsigned > solPdeIndexHT ( NLayers );

  std::vector < unsigned > solIndexT ( NLayers );

  vector< int > l2GMapRow; // local to global mapping
  vector< int > l2GMapColumn; // local to global mapping

  for ( unsigned i = 0; i < NLayers; i++ ) {
    char name[10];
    sprintf ( name, "h%d", i );
    solIndexh[i] = mlSol->GetIndex ( name ); // get the position of "hi" in the sol object
    //solPdeIndexh[i] = mlPdeSys->GetSolPdeIndex ( name ); // get the position of "hi" in the pdeSys object

    sprintf ( name, "v%d", i );
    solIndexv[i] = mlSol->GetIndex ( name ); // get the position of "vi" in the sol object
    //solPdeIndexv[i] = mlPdeSys->GetSolPdeIndex ( name ); // get the position of "vi" in the pdeSys object

    sprintf ( name, "HT%d", i );
    solIndexHT[i] = mlSol->GetIndex ( name ); // get the position of "Ti" in the sol object
    solPdeIndexHT[i] = mlPdeSys->GetSolPdeIndex ( name ); // get the position of "Ti" in the pdeSys object

    sprintf ( name, "T%d", i );
    solIndexT[i] = mlSol->GetIndex ( name ); // get the position of "Ti" in the sol object

  }

  unsigned solTypeh = mlSol->GetSolutionType ( solIndexh[0] ); // get the finite element type for "hi"
  unsigned solTypev = mlSol->GetSolutionType ( solIndexv[0] ); // get the finite element type for "vi"
  unsigned solTypeHT = mlSol->GetSolutionType ( solIndexHT[0] ); // get the finite element type for "Ti"

  KK->zero();
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
    
    vector < adept::adouble > solHTmm(NLayers);    // local coordinates
    vector < adept::adouble > solHTpp(NLayers);    // local coordinates

    //vector< adept::adouble > aResh ( NLayers );
    //vector< adept::adouble > aResv ( NLayers );
    vector< adept::adouble > aResHT ( NLayers );

    unsigned bc1 = ( i == start ) ? 0 : 1;
    unsigned bc2 = ( i == end - 1 ) ? 0 : 1;
    
    unsigned bc3 = ( i > start + 1) ? 1 : 0;
    unsigned bc4 = ( i < end - 2 ) ? 1 : 0;

    l2GMapRow.resize ( NLayers );
    l2GMapColumn.resize ( ( 1 + bc1 + bc2 ) * NLayers );
    //l2GMapColumn.resize ( ( 1 + bc1 + bc2 + bc3 + bc4) * NLayers );

    //std::fill ( aResh.begin(), aResh.end(), 0 ); //set aRes to zero
    std::fill ( aResHT.begin(), aResHT.end(), 0 ); //set aRes to zero

    for ( unsigned j = 0; j < NLayers; j++ ) {

      solh[j] = ( *sol->_Sol[solIndexh[j]] ) ( i );
      solHT[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i );
      //l2GMapRow[j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i );
      l2GMapRow[/*NLayers +*/ j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i );

      //l2GMapColumn[j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i );
      l2GMapColumn[/*NLayers +*/ j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i );

      solvm[j] = ( *sol->_Sol[solIndexv[j]] ) ( i );
      solvp[j] = ( *sol->_Sol[solIndexv[j]] ) ( i + 1 );

      //l2GMapColumn[2 * NLayers + j] = pdeSys->GetSystemDof ( solIndexv[j], solPdeIndexv[j], 0, i );
      //l2GMapColumn[3 * NLayers + j] = pdeSys->GetSystemDof ( solIndexv[j], solPdeIndexv[j], 1, i );

      if ( i > start ) {
        solhm[j] = ( *sol->_Sol[solIndexh[j]] ) ( i - 1 );
        solHTm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 1 );

        //l2GMapColumn[4 * NLayers + j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i - 1 );
        l2GMapColumn[NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 1 );

      }

      if ( i < end - 1 ) {
        solhp[j] = ( *sol->_Sol[solIndexh[j]] ) ( i + 1 );
        solHTp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 1 );

        //l2GMapColumn[ ( 4 + 2 * bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexh[j], solPdeIndexh[j], 0, i + 1 );
        l2GMapColumn[ ( 1 + bc1) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 1 );
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
    
    s.new_recording();

    vector < double > x ( 2 ); // local coordinates
    for ( unsigned j = 0; j < 2; j++ ) {
      unsigned xDof  = msh->GetSolutionDof ( j, i, 2 ); // global to global mapping between coordinates node and coordinate dof
      x[j] = ( *msh->_topology->_Sol[0] ) ( xDof ); // global extraction and local storage for the element coordinates
    }
    double dx = x[1] - x[0];

    double b = 5; //( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );

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
    
//     NEW w
//     for ( unsigned k = NLayers; k > 1; k-- ) {
//       w[k - 1] = w[k] - ( hALE[k - 1] - solh[k - 1]) / dt; 
//       if(bc2){ 
//         w[k - 1] -=   0.5 * ( solh[k - 1] + solhp[k - 1] ) * solvp[k - 1] /dx;
//       }
//       else{
//         w[k - 1] -=   solh[k - 1] * 1 /dx;  
//       }
//       if(bc1){
//         w[k - 1] +=   0.5 * ( solh[k - 1] + solhm[k - 1] ) * solvm[k - 1] /dx;
//       }
//       else{
//         w[k - 1] +=   solh[k - 1] * 1 /dx;   
//       }
//       //std::cout<< w[k-1] << " ";
//     }
//     //std::cout<<std::endl;
    
//     OLD w     
//     for(unsigned k = 1; k < NLayers; k++){
//       w[k] = w[k-1] + (  0.5 * ( solh[k-1].value() + solhp[k-1].value() ) * solvp[k-1].value()
// 		       - 0.5 * ( solh[k-1].value() + solhm[k-1].value() ) * solvm[k-1].value() )/dx
// 		    + ( hALE[k-1] - solh[k-1].value()) / dt;//TODO
// 		    //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"<<w[k-1]<<std::endl;
//     }
    
    std::vector < double > zMid ( NLayers );
    for ( unsigned k = 0; k < NLayers; k++ ) {
      zMid[k] = -b + solh[k] / 2.;
      for ( unsigned i = k + 1; i < NLayers; i++ ) {
        zMid[k] += solh[i];
      }
    }
    
    for ( unsigned k = NLayers; k > 1; k-- ) {
      w[k-1] = - 0.1 * zMid[k-1];
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
        aResHT[k] += w[k + 1] * 0.5 * ( solHT[k] / solh[k] + solHT[k + 1] / solh[k + 1] );
      }
      if ( k > 0 ) {
        aResHT[k] -= w[k] * 0.5 * ( solHT[k - 1] / solh[k - 1] + solHT[k] / solh[k] );
      }
      
      adept::adouble deltaZt = 0.;
      adept::adouble deltaZb = 0.;
      adept::adouble ht = 0.;
      adept::adouble hb = 0.;
      if ( k > 0 ){
	ht = (solhm[k-1] + solhm[k] + solhp[k-1] + solhp[k]) / 4.;
	deltaZt = ( solHT[k-1] - solHT[k] ) / ht;
	//aResv[k] -= 0.5 * w[k] * deltaZt;
      }
      else{
	ht = 0.5 * (solhm[k] + solhp[k]);
	deltaZt = 0.*( 0. - solHT[k] ) / ht;
      }
      if (k < NLayers - 1) {
	hb = (solhm[k] + solhm[k+1] + solhp[k] + solhp[k+1] ) / 4.;
	deltaZb = (solHT[k] - solHT[k+1]) / hb;
	//aResv[k] -= 0.5 * w[k+1] * deltaZb; 
      }
      else{
	hb = 0.5 * (solhm[k] + solhp[k]);
	deltaZb = 0.*(solHT[k] - 0.) / hb;
      }
      
      //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAA"<<deltaZt - deltaZb<<std::endl;
      
      aResHT[k] += solhm[k] * k_v * (deltaZt - deltaZb) / ( (ht + hb) / 2. ); // vertical diffusion

      //aResHT[k] += ((solhp[k] - solhm[k]) * k_v * (solHTp[k] - solHTm[k])) / (dx*dx); // horizontal diffusion
      aResHT[k] += k_h * solh[k] * (solHTm[k] - solHT[k])/(dx*dx); // horizontal diffusion
      aResHT[k] += k_h * solh[k] * (solHTp[k] - solHT[k])/(dx*dx); // horizontal diffusion

    }

    vector< double > Res ( NLayers ); // local redidual vector
    for ( unsigned k = 0; k < NLayers; k++ ) {
      //Res[k] =  aResh[k].value();
      Res[k] =  aResHT[k].value();
      //std::cout<< "Res["<<k<<"] = " << Res[k] <<std::endl;
      //std::cout<< "Res["<<NLayers+k<<"] = " << Res[NLayers+k] <<std::endl;
    }

    RES->add_vector_blocked ( Res, l2GMapRow );

    //s.dependent ( &aResh[0], NLayers );
    s.dependent ( &aResHT[0], NLayers );

    // define the independent variables
    //s.independent ( &solh[0], NLayers );
    s.independent ( &solHT[0], NLayers );
    //s.independent ( &solvm[0], NLayers );
    //s.independent ( &solvp[0], NLayers );
    if ( i > start ) {
      //s.independent ( &solhm[0], NLayers );
      s.independent ( &solHTm[0], NLayers );
    }
    if ( i < end - 1 ) {
      //s.independent ( &solhp[0], NLayers );
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
    //vector < double > Jac ( NLayers * NLayers * ( 1 + bc1 + bc2 + bc3 +bc4 ) );
    s.jacobian ( &Jac[0], true );

    //store K in the global matrix KK
    KK->add_matrix_blocked ( Jac, l2GMapRow, l2GMapColumn );

    s.clear_independents();
    s.clear_dependents();

  }

  RES->close();
  KK->close();


//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,900,900,&viewer);
//   PetscObjectSetName((PetscObject)viewer,"FSI matrix");
//   PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(),viewer);
//   double a;
//   std::cin>>a;
//

//  abort();
  MFN mfn;
  Mat A = ( static_cast<PetscMatrix*> ( KK ) )->mat();
  FN f, f1, f2, f3 , f4;

  //std::cout << "dt = " << dt << " dx = "<< dx << " maxWaveSpeed = "<<maxWaveSpeed << std::endl;
  std::cout << "dt = " << dt << std::endl;

  //dt = 100.;

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

  sol->UpdateSol ( mlPdeSys->GetSolPdeIndex(), EPS, pdeSys->KKoffset );

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
      if (i == 0) valueT = 0.;
//       if (i == msh->_dofOffset[solTypeHT][iproc + 1] - 1 ) valueT = 0.;
      
      
      sol->_Sol[solIndexT[k]]->set ( i, valueT );
    }

    sol->_Sol[solIndexT[k]]->close();

  }


}

void RK ( MultiLevelProblem& ml_prob, const unsigned & numberOfTimeSteps )
{

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

//     for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
//         for ( unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++ ) {
//             double valueT = ( *sol->_SolOld[solIndexT[k]] ) ( i );
//             //double valueT = ( *sol->_Sol[solIndexT[k]] ) ( i );
//             double valueH = ( *sol->_Sol[solIndexh[k]] ) ( i );
// 
//             double valueHT = valueT * valueH;
// 
//             sol->_Sol[solIndexHT[k]]->set ( i, valueHT );
//         }
//         sol->_Sol[solIndexHT[k]]->close();
//     }

//     std::vector < double > maxW ( NLayers, -1.e6 );
//     maxW[0] = 0.;

    unsigned start = msh->_dofOffset[solTypeHT][iproc];
    unsigned end = msh->_dofOffset[solTypeHT][iproc + 1];
    for ( unsigned i =  start; i <  end; i++ ) {

        vector < double > solhm ( NLayers );
        vector < double > solh ( NLayers ); // local coordinates
        vector < double > solhp ( NLayers );
        vector < double > solvm ( NLayers ); // local coordinates
        vector < double > solvp ( NLayers ); // local coordinates
        
//         vector < double > solHT ( NLayers ); // local coordinates
//         vector < double > solHTp ( NLayers ); // local coordinates
//         vector < double > solHTm ( NLayers ); // local coordinates
        
        vector < double > solTm ( NLayers ); // local coordinates
        vector < double > solT ( NLayers ); // local coordinates
        vector < double > solTp ( NLayers ); // local coordinates
        
//         vector < double > solHTmm ( NLayers ); // local coordinates
//         vector < double > solHTpp ( NLayers ); // local coordinates
        
//         vector < double > solTmm ( NLayers ); // local coordinates
//         vector < double > solTpp ( NLayers ); // local coordinates
//         vector < double > solhmm ( NLayers ); // local coordinates
//         vector < double > solhpp ( NLayers ); // local coordinates

        unsigned bc1 = ( i == start ) ? 0 : 1;
        unsigned bc2 = ( i == end - 1 ) ? 0 : 1;
        unsigned bc3 = ( i > start + 1 ) ? 1 : 0;
        unsigned bc4 = ( i < end - 2 ) ? 1 : 0;

        for ( unsigned j = 0; j < NLayers; j++ ) {

            solh[j] = ( *sol->_Sol[solIndexh[j]] ) ( i );
            solT[j] = ( *sol->_Sol[solIndexT[j]] ) ( i );
            //solHT[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i );

            solvm[j] = ( *sol->_Sol[solIndexv[j]] ) ( i );
            solvp[j] = ( *sol->_Sol[solIndexv[j]] ) ( i + 1 );

            if ( i > start ) {
                solhm[j] = ( *sol->_Sol[solIndexh[j]] ) ( i - 1 );
                solTm[j] = ( *sol->_Sol[solIndexT[j]] ) ( i - 1 );
                //solHTm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 1 );

           }

            if ( i < end - 1 ) {
                solhp[j] = ( *sol->_Sol[solIndexh[j]] ) ( i + 1 );
                solTp[j] = ( *sol->_Sol[solIndexT[j]] ) ( i + 1 );
                if ( (i+1) ==  end-1 ) solTp[j] = 0.;
                //solHTp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 1 );

            }

//             if ( i > start + 1 ) {
//                 solTmm[j] = ( *sol->_Sol[solIndexT[j]] ) ( i - 2 );
//                 solhmm[j] = ( *sol->_Sol[solIndexh[j]] ) ( i - 2 );
//                 //solHTmm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 2 );
//             }
// 
//             if ( i < end - 2 ) {
//                 solTpp[j] = ( *sol->_Sol[solIndexT[j]] ) ( i + 2 );
//                 if ( (i+2) ==  end-1 ) solTpp[j] = 0.;
//                 solhpp[j] = ( *sol->_Sol[solIndexh[j]] ) ( i + 2 );
//                 //solHTpp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 2 );
//             }

        }


        vector < double > x ( 2 ); // local coordinates
        for ( unsigned j = 0; j < 2; j++ ) {
            unsigned xDof  = msh->GetSolutionDof ( j, i, 2 ); // global to global mapping between coordinates node and coordinate dof
            x[j] = ( *msh->_topology->_Sol[0] ) ( xDof ); // global extraction and local storage for the element coordinates
        }
        double dx = x[1] - x[0];

        double b = 5.;

        std::vector < double > w ( NLayers + 1, 0. );

 
//         double xmid = 0.5 * ( x[1] + x[0] );
//         for ( unsigned k = NLayers; k > 1; k-- ) {
//             w[k - 1] = - 0.1 * zTop[k-1];
// //             if ( maxW[k - 1] < w[k - 1] ) {
// //                 maxW[k - 1] = w[k - 1];
// //             }
//         }
        
        std::vector < double > zMid ( NLayers );
        for ( unsigned k = 0; k < NLayers; k++ ) {
          zMid[k] = -b + solh[k] / 2.;
          for ( unsigned i = k + 1; i < NLayers; i++ ) {
            zMid[k] += solh[i];
          }
        }
        
        std::vector < double > zTop ( NLayers );
        zTop[0] = 0;
        for ( unsigned k = 1; k < NLayers; k++ ) {
            zTop[k] = zTop[k-1] - solh[k];
        }
        
        for ( unsigned k = NLayers; k > 1; k-- ) {
          w[k-1] = - 0.1 * zTop[k-1];
        }        

        std::vector < double > k1_RK ( NLayers, 0. );
        std::vector < double > k2_RK ( NLayers, 0. );
        std::vector < double > k3_RK ( NLayers, 0. );
        std::vector < double > k4_RK ( NLayers, 0. );

        unsigned RK_order = 2;

        for ( unsigned RK_step = 0; RK_step < RK_order; RK_step++ ) {
            for ( unsigned k = 0; k < NLayers; k++ ) {
                double LHS = 0.;
                double addition = 0.;
                if ( RK_step == 1 ) {
                    if ( RK_order == 2 ) {
                        addition = k1_RK[k] ;
                    } else if ( RK_order == 3 || RK_order == 4 ) {
                        addition = k1_RK[k] * 0.5;
                    }
                } else if ( RK_step == 2 ) {
                    if ( RK_order == 2 || RK_order == 4 ) {
                        addition = k2_RK[k] * 0.5;
                    } else if ( RK_order == 3 ) {
                        addition = - k1_RK[k] + k2_RK[k] * 2.;
                    }
                } else if ( RK_step == 3 ) {
                    addition = k3_RK[k];
                }


//                 if(addition > 1.e-2 || addition < - 1.e-2) std::cout<< "additon = " << addition <<  std::endl;

                //BEGIN HORIZONTAL ADVECTION
                
                
                
                //BEGIN CENTRAL DIFF h and T separate
                if ( i > start ) {
                     LHS +=  ( 0.5 * (solhm[k] + addition + solh[k] + addition ) ) * solvm[k] * ( 0.5 * (solTm[k] + addition + solT[k] + addition ) ) / dx;                      
                }
                if ( i < end - 1 ) {
                      LHS -= ( 0.5 * (solh[k] + addition + solhp[k] + addition) ) * solvp[k] * ( 0.5 * (solT[k] + addition + solTp[k] + addition ) ) / dx; 
                }
                //END
    
                
                //BEGIN CENTRAL DIFF hT multiplied
/*                if ( i > start ) {
                    LHS += 0.5 * (solHTm[k] + solHT[k]) * solvm[k]  / dx; 
                }
                if ( i < end - 1 ) {
                    LHS -= 0.5 * (solHT[k] + solHTp[k]) * solvp[k]  / dx; 
                }       */         
                //END
                
                
                //BEGIN FIRST ORDER UPWINDING h and T separate
//                 if ( i > start ) {
//                     if ( solvm[k] > 0 ) {
//                         LHS += ( ( solTm[k] + addition ) * ( solhm[k] + addition ) ) * solvm[k] / dx;
//                     } else {
//                         LHS += ( ( solT[k] + addition ) * ( solh[k] + addition ) ) * solvm[k] / dx;
//                     }
//                 }
//                 if ( i < end - 1 ) {
//                     if ( solvp[k] > 0 ) {
//                         LHS -= ( ( solT[k] + addition ) * ( solh[k] + addition ) ) * solvp[k] / dx; 
//                     } else {
//                         LHS -= ( ( solTp[k] + addition ) * ( solhp[k] + addition ) ) * solvp[k] / dx; 
//                     }
//                 }
                //END
                
                
                //BEGIN FIRST ORDER UPWINDING hT multiplied
//                 if ( i > start ) {
//                     if ( solvm[k] > 0 ) {
//                         LHS += ( solHTm[k] + addition ) * solvm[k] / dx;
//                     } else {
//                         LHS += ( solHT[k] + addition ) * solvm[k] / dx;
//                     }
//                 }
//                 if ( i < end - 1 ) {
//                     if ( solvp[k] > 0 ) {
//                         LHS -= ( solHT[k] + addition ) * solvp[k] / dx; 
//                     } else {
//                         LHS -= ( solHTp[k] + addition ) * solvp[k] / dx; 
//                     }
//                 }
                //END
                
                
                //BEGIN THIRD ORDER UPWINDING h and T separate
                
//                 if ( i > start ) {
//                     LHS += 0.5 * ( (solTm[k] + addition) * (solhm[k] + addition) + (solT[k] + addition) * (solh[k] + addition) ) * solvm[k] / dx;
//                     if ( solvm[k] > 0 ) {
//                         if ( i > start + 1 ) {
//                             LHS += - 1. / 6. * ( (solT[k] + addition) * (solh[k] + addition) - 2.* (solTm[k] + addition) * (solhm[k] + addition) + (solTmm[k] + addition) * (solhmm[k] + addition) ) * solvm[k]  / dx;
//                         }
//                     } else {
//                         if ( i < end - 1 ) {
//                             LHS += - 1. / 6. * ( (solTp[k] + addition) * (solhp[k] + addition) - 2.* (solT[k] + addition) * (solh[k] + addition) + (solTm[k] + addition) * (solhm[k] + addition) ) * solvm[k]  / dx;
//                         }
//                     }
//                 }
//                 if ( i < end - 1 ) {
//                     LHS -= 0.5 * ( (solTp[k] + addition) * (solhp[k] + addition) + (solT[k] + addition) * (solh[k] + addition) ) * solvp[k] / dx;
//                     if ( solvp[k] > 0 ) {
//                         if ( i > start ) {
//                             LHS -= - 1. / 6. * ( (solTp[k] + addition) * (solhp[k] + addition) - 2.* (solT[k] + addition) * (solh[k] + addition) + (solTm[k] + addition) * (solhm[k] + addition) ) * solvp[k]  / dx;
//                         }
//                     } else {
//                         if ( i < end - 2 ) {
//                             LHS -= - 1. / 6. * ( (solTpp[k] + addition) * (solhpp[k] + addition) - 2.* (solTp[k] + addition) * (solhp[k] + addition) + (solT[k] + addition) * (solh[k] + addition) ) * solvp[k]  / dx;
//                         }
//                     }
//                 }
                
                //END
                
                
                
                //END HORIZONTAL ADVECTION

                
                
                
                //BEGIN VERTICAL ADVECTION
                
                
                //BEGIN CENTRAL DIFF h and T separate
                if ( k < NLayers - 1 ) {
                     LHS += w[k + 1] * 0.5 * ( solT[k] + addition + solT[k + 1] + addition ); 
                }
                if ( k > 0 ) {
                    LHS -= w[k] * 0.5 * ( solT[k - 1] + addition + solT[k] + addition ) ; 
                }
                //END
                
                //BEGIN CENTRAL DIFF hT multiplied
//                 if ( k < NLayers - 1 ) {
//                     LHS += w[k + 1] * 0.5 * ( ( solHT[k] + addition ) / solh[k] + ( solHT[k + 1] + addition ) / solh[k + 1] ); 
//                 }
//                 if ( k > 0 ) {
//                     LHS -= w[k] * 0.5 * ( ( solHT[k - 1] + addition ) / solh[k - 1] + ( solHT[k] + addition ) / solh[k] ); 
//                 }
                //END
                
                //BEGIN FIRST ORDER UPWIND h and T separate
//                 if ( k < NLayers - 1 ) {
//                     if ( w[k + 1] > 0 ) {
//                         LHS += w[k + 1] * ( solT[k + 1] + addition ) ;
//                     } else {
//                         LHS += w[k + 1] * ( solT[k] + addition ) ;
//                     }
//                 }
//                 if ( k > 0 ) {
//                     if ( w[k] > 0 ) {
//                         LHS -= w[k] * ( solT[k] + addition );
//                     } else {
//                         LHS -= w[k] * ( solT[k - 1] + addition );
//                     }
//                 }
               //END
               
               //BEGIN FIRST ORDER UPWIND hT multiplied
//                 if ( k < NLayers - 1 ) {
//                     if ( w[k + 1] > 0 ) {
//                         LHS += w[k + 1] * ( ( solHT[k + 1] + addition ) / solh[k + 1]  );
//                     } else {
//                         LHS += w[k + 1] * ( ( solHT[k] + addition ) / solh[k] );
//                     }
//                 }
//                 if ( k > 0 ) {
//                     if ( w[k] > 0 ) {
//                         LHS -= w[k] * ( ( solHT[k] + addition ) / solh[k] );
//                     } else {
//                         LHS -= w[k] * ( ( solHT[k - 1] + addition ) / solh[k - 1] );
//                     }
//                 }
              //END
                
                  
                  
               //END VERTICAL ADVECTION
               
                // horizontal diffusion
                if ( i > start ) {
                  LHS -= k_h * ( ( 0.5 * ( solhm[k] + addition + solh[k] + addition ) ) * ( ( solT[k] + addition )  - ( solTm[k] + addition ) ) / dx ) / dx; 
                  //LHS += k_h * (0.5 * (solhm[k] + solh[k])) * ((solHTm[k] / solhm[k] + addition)  - (solHT[k] / solh[k] + addition))/(dx*dx); 
                }
                if ( i < end-1 ) {
                  LHS += k_h * ( ( 0.5 * ( solh[k] + addition + solhp[k] + addition ) )  * ( ( solTp[k] + addition )  - ( solT[k] + addition ) ) / dx ) / dx; 
                  //LHS += k_h * (0.5 * (solhp[k] + solh[k])) * ((solHTp[k] / solhp[k] + addition)  - (solHT[k] / solh[k] + addition))/(dx*dx); 
                }


                if ( RK_step == 0 ) {
                    k1_RK[k] = LHS * dt;
                } else if ( RK_step == 1 ) {
                    k2_RK[k] = LHS * dt;
                } else if ( RK_step == 2 ) {
                    k3_RK[k] = LHS * dt;
                } else {
                    k4_RK[k] = LHS * dt;
                }

            }

        }

        for ( unsigned k = 0; k < NumberOfLayers; k++ ) {

            double valueHT;
            double valueT;

            if ( RK_order == 1 ) {                
//                 valueHT = solHT[k] + k1_RK[k] ;
//                 valueT = solHT[k]/solh[k] + k1_RK[k];
                valueHT = solT[k] * solh[k] + k1_RK[k] ;
                valueT = solT[k] + k1_RK[k] / solh[k];
            }

            else if ( RK_order == 2 ) {                
//                 valueHT = solHT[k] + 0.5 * ( k1_RK[k] + k2_RK[k] );
//                 valueT = solHT[k]/solh[k] + ( 0.5 * ( k1_RK[k] + k2_RK[k] ) );
                valueHT = solT[k] * solh[k] + 0.5 * ( k1_RK[k] + k2_RK[k] );
                valueT = solT[k]  + ( 0.5 * ( k1_RK[k] + k2_RK[k] ) ) / solh[k];
            }

            else if ( RK_order == 3 ) {                
//                 valueHT = solHT[k] + 1. / 6. * ( k1_RK[k] + 4. * k2_RK[k] + k3_RK[k] );
//                 valueT = solHT[k]/solh[k] + ( 1. / 6. * ( k1_RK[k] + 4. * k2_RK[k] + k3_RK[k] ) );
                valueHT = solT[k] * solh[k] + 1. / 6. * ( k1_RK[k] + 4. * k2_RK[k] + k3_RK[k] );
                valueT = solT[k] + ( 1. / 6. * ( k1_RK[k] + 4. * k2_RK[k] + k3_RK[k] ) ) /  solh[k];
            }

            else if ( RK_order == 4 ) {                
//                 valueHT = solHT[k] + 1. / 6. * ( k1_RK[k] + 2.*k2_RK[k] + 2.*k3_RK[k] + k4_RK[k] );
//                 valueT = solHT[k]/solh[k] + ( 1. / 6. * ( k1_RK[k] + 2.*k2_RK[k] + 2.*k3_RK[k] + k4_RK[k] ) );
                valueHT = solT[k] * solh[k] + 1. / 6. * ( k1_RK[k] + 2.*k2_RK[k] + 2.*k3_RK[k] + k4_RK[k] );
                valueT = solT[k] + ( 1. / 6. * ( k1_RK[k] + 2.*k2_RK[k] + 2.*k3_RK[k] + k4_RK[k] ) ) /  solh[k];
            }

            sol->_Sol[solIndexHT[k]]->set ( i, valueHT );
            if (i == start || i == end - 1 ) valueT = 0.;
            sol->_Sol[solIndexT[k]]->set ( i, valueT );

            sol->_Sol[solIndexHT[k]]->close();
            sol->_Sol[solIndexT[k]]->close();
        }

    }

    //BEGIN no vertical diffusion
    for ( unsigned k = 0; k < NumberOfLayers; k++) {
        for ( unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++) {
            double valueT = ( *sol->_Sol[solIndexT[k]] ) ( i );

            std::cout.precision ( 14 );
            if ( counter==numberOfTimeSteps - 1) {
              std::cout << valueT << std::endl;
            }
        }
    }
    //END


}

