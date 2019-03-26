#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelMesh.hpp"
#include "TransientSystem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Parameter.hpp"
#include "FemusInit.hpp"
#include "SparseMatrix.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "SolvertypeEnum.hpp"
#include "FElemTypeEnum.hpp"
#include "OprtrTypeEnum.hpp"
#include "Files.hpp"

using std::cout;
using std::endl;

using namespace femus;

void AssembleMatrixResNS ( MultiLevelProblem &ml_prob );
void AssembleMatrixResT ( MultiLevelProblem &ml_prob );

void SetLambda ( MultiLevelSolution &mlSol, const unsigned &level, const  FEOrder &order, Operator operatorType );


double InitVariableU ( const std::vector < double >& x );

bool SetBoundaryConditionTurek ( const std::vector < double >& x, const char name[],
                                 double &value, const int FaceName, const double time );

bool SetBoundaryConditionCavityFlow ( const std::vector < double >& x, const char name[],
                                      double &value, const int FaceName, const double time );

int main ( int argc, char **args ) {

  /// Init Petsc-MPI communicator
  FemusInit mpinit ( argc, args, MPI_COMM_WORLD );

  Files files;
  files.CheckIODirectories();
  //files.RedirectCout();

  bool Gmres = 0, Asm = 0;
  if ( argc >= 2 ) {
    if ( !strcmp ( "gmres", args[1] ) ) 	Gmres = 1;
    else if ( !strcmp ( "asm", args[1] ) ) 	Asm = 1;

    if ( Gmres + Asm == 0 ) {
      cout << "wrong input arguments!" << endl;
      exit ( 0 );
    }
  }
  else {
    cout << "No input argument set default smoother = Gmres" << endl;
    Gmres = 1;
  }

  /// INIT MESH =================================

  unsigned short nm, nr;
  nm = 4;
  std::cout << "MULTIGRID levels: " << nm << endl;

  nr = 0;
  std::cout << "MAX_REFINEMENT levels: " << nr << endl << endl;

  int tmp = nm;
  nm += nr;
  nr = tmp;

  char *infile = new char [50];

  sprintf ( infile, "./input/box10x10.neu" );

  //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;

  //Steadystate NonLinearMultiLevelProblem
  //MultiLevelMesh ml_msh(nm,nr,infile,"seventh",Lref,SetRefinementFlag);
  MultiLevelMesh ml_msh;
  ml_msh.ReadCoarseMesh ( infile, "seventh", Lref );
  ml_msh.RefineMesh ( nm, nr, NULL );

  ml_msh.EraseCoarseLevels ( nm - 1 );

  MultiLevelSolution ml_sol ( &ml_msh );

  // generate solution vector

  FEOrder orderPre = FIRST;
  FEOrder orderVel = FIRST;
  FEOrder orderTemp = FIRST;

  ml_sol.AddSolution ( "U", LAGRANGE, orderVel );
  ml_sol.AddSolution ( "V", LAGRANGE, orderVel );
  ml_sol.AddSolution ( "lmbd", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false );

  // the pressure variable should be the last for the Schur decomposition
  // ml_sol.AddSolution("P",DISCONTINUOUS_POLYNOMIAL,FIRST);
  ml_sol.AddSolution ( "P", LAGRANGE, orderPre );
  ml_sol.AssociatePropertyToSolution ( "P", "Pressure" );

  ml_sol.AddSolution ( "T", LAGRANGE, orderTemp );

  //Initialize (update Init(...) function)
  ml_sol.Initialize ( "U" );
  ml_sol.Initialize ( "V" );
  ml_sol.Initialize ( "P" );
  ml_sol.Initialize ( "T" );
  ml_sol.Initialize ( "lmbd" );
  //Set Boundary (update Dirichlet(...) function)
  ml_sol.AttachSetBoundaryConditionFunction ( SetBoundaryConditionCavityFlow );
  ml_sol.GenerateBdc ( "U" );
  ml_sol.GenerateBdc ( "V" );
  ml_sol.GenerateBdc ( "P" );
  ml_sol.GenerateBdc ( "T" );

  SetLambda ( ml_sol, 0, orderVel, ELASTICITY );
  //SetLambda(ml_sol, 0, 2,"diffusion");

  MultiLevelProblem ml_prob ( &ml_sol );

  // add fluid material
  Parameter parameter ( Lref, Uref );

  // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  Fluid fluid ( parameter, 0.001, 1, "Newtonian", 0.001, 1. );
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;

  ml_prob.parameters.set<Fluid> ( "Fluid" ) = fluid;


  //BEGIN Navier-Stokes Multilevel Problem
  std::cout << std::endl;
  std::cout << " *********** Navier-Stokes ************  " << std::endl;

  NonLinearImplicitSystem & system1 = ml_prob.add_system<NonLinearImplicitSystem> ( "Navier-Stokes" );
  system1.AddSolutionToSystemPDE ( "U" );
  system1.AddSolutionToSystemPDE ( "V" );
  system1.AddSolutionToSystemPDE ( "P" );

  // Set MG Options
  system1.SetAssembleFunction ( AssembleMatrixResNS );
  system1.SetMaxNumberOfNonLinearIterations ( 90 );
  system1.SetMaxNumberOfLinearIterations ( 2 );
  system1.SetAbsoluteLinearConvergenceTolerance ( 1.e-10 );
  system1.SetNonLinearConvergenceTolerance ( 1.e-10 );
  system1.SetMgType ( F_CYCLE );
  system1.SetNumberPreSmoothingStep ( 1 );
  system1.SetNumberPostSmoothingStep ( 1 );

  //Set Smoother Options
  if ( Gmres ) 		system1.SetLinearEquationSolverType ( FEMuS_DEFAULT );
  else if ( Asm ) 		system1.SetLinearEquationSolverType ( FEMuS_ASM );
  //else if(Vanka)	system1.SetLinearEquationSolverType(VANKA_SMOOTHER);

  system1.init();
  //common smoother options
//   system1.AddStabilization(true);
  system1.SetSolverFineGrids ( GMRES );
  system1.SetPreconditionerFineGrids ( ILU_PRECOND );
  system1.SetTolerances ( 1.e-12, 1.e-20, 1.e+50, 4 );

  system1.ClearVariablesToBeSolved();
  //system1.AddVariableToBeSolved("All");
  system1.AddVariableToBeSolved ( "U" );
  system1.AddVariableToBeSolved ( "V" );
  system1.AddVariableToBeSolved ( "P" );
  //for Vanka and ASM smoothers
  system1.SetNumberOfSchurVariables ( 0 );
  system1.SetElementBlockNumber ( 4 );
  //system1.SetElementBlockNumber("All",1);
  //for Gmres smoother
  system1.SetDirichletBCsHandling ( PENALTY );
  //system1.SetDirichletBCsHandling(ELIMINATION);

  // Solve Navier-Stokes system
  ml_prob.get_system("Navier-Stokes").SetOuterSolver(PREONLY);
  ml_prob.get_system ( "Navier-Stokes" ).MGsolve();
  //END Navier-Stokes Multilevel Problem


  //BEGIN Temperature MultiLevel Problem
  std::cout << std::endl;
  std::cout << " *********** Temperature ************* " << std::endl;

  LinearImplicitSystem & system2 = ml_prob.add_system<LinearImplicitSystem> ( "Temperature" );
  system2.AddSolutionToSystemPDE ( "T" );


  // Set MG Options
  system2.SetAssembleFunction ( AssembleMatrixResT );
  system2.SetMaxNumberOfLinearIterations ( 6 );
  system2.SetAbsoluteLinearConvergenceTolerance ( 1.e-9 );
  system2.SetMgType ( V_CYCLE );
  system2.SetNumberPreSmoothingStep ( 1 );
  system2.SetNumberPostSmoothingStep ( 1 );

  //Set Smoother Options
  if ( Gmres ) 		system2.SetLinearEquationSolverType ( FEMuS_DEFAULT );
  else if ( Asm ) 		system2.SetLinearEquationSolverType ( FEMuS_ASM );

  system2.init();
  //common smoother option
  system2.SetSolverFineGrids ( GMRES );
  system2.SetTolerances ( 1.e-12, 1.e-20, 1.e+50, 4 );
  system2.SetPreconditionerFineGrids ( ILU_PRECOND );
  //for Vanka and ASM smoothers
  system2.ClearVariablesToBeSolved();
  system2.AddVariableToBeSolved ( "All" );
  system2.SetNumberOfSchurVariables ( 0 );
  system2.SetElementBlockNumber ( 4 );
  //for Gmres smoother
  system2.SetDirichletBCsHandling ( PENALTY );
  //system2.SetDirichletBCsHandling(ELIMINATION);


  // Solve Temperature system
  ml_prob.get_system ( "Temperature" ).SetOuterSolver(PREONLY);
  ml_prob.get_system ( "Temperature" ).MGsolve();
  //END Temperature Multilevel Problem

  /// Print all solutions
  std::vector<std::string> print_vars;
  print_vars.push_back ( "U" );
  print_vars.push_back ( "V" );
  print_vars.push_back ( "P" );
  print_vars.push_back ( "T" );

  VTKWriter vtkio ( &ml_sol );
  vtkio.Write ( files.GetOutputPath(), "biquadratic", print_vars );
  //vtkio.write(DEFAULT_OUTPUTDIR,"biquadratic",print_vars);

  GMVWriter gmvio ( &ml_sol );
  gmvio.Write ( DEFAULT_OUTPUTDIR, "biquadratic", print_vars );
  // gmvio.write(files.GetOutputPath(),"biquadratic",print_vars);

  //   XDMFWriter xdmfio(ml_sol);
  //   xdmfio.write(files.GetOutputPath(),"biquadratic",print_vars);

  //Destroy all the new systems
  ml_prob.clear();

  delete [] infile;
  return 0;
}

//--------------------------------------------------------------------------------------------------------------

double InitVariableU ( const std::vector < double >& x ) {
  double um = 0.2;
  double  value = 1.5 * um * ( 4.0 / ( 0.1681 ) ) * x[1] * ( 0.41 - x[1] );
  return value;
}

//-------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionTurek ( const std::vector < double >& x, const char name[],
                                 double &value, const int FaceName, const double time ) {
  bool test = 1; //Dirichlet
  value = 0.;
  //   cout << "Time bdc : " <<  time << endl;
  if ( !strcmp ( name, "U" ) ) {
    if ( 1 == FaceName ) { //inflow
      test = 1;
      double um = 0.2; // U/Uref
      value = 1.5 * 0.2 * ( 4.0 / ( 0.1681 ) ) * x[1] * ( 0.41 - x[1] );
    }
    else if ( 2 == FaceName ) { //outflow
      test = 0;
      //    test=1;
      value = 0.;
    }
    else if ( 3 == FaceName ) { // no-slip fluid wall
      test = 1;
      value = 0.;
    }
    else if ( 4 == FaceName ) { // no-slip solid wall
      test = 1;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "V" ) ) {
    if ( 1 == FaceName ) {      //inflow
      test = 1;
      value = 0.;
    }
    else if ( 2 == FaceName ) { //outflow
      test = 0;
      //    test=1;
      value = 0.;
    }
    else if ( 3 == FaceName ) { // no-slip fluid wall
      test = 1;
      value = 0;
    }
    else if ( 4 == FaceName ) { // no-slip solid wall
      test = 1;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "W" ) ) {
    if ( 1 == FaceName ) {
      test = 1;
      value = 0.;
    }
    else if ( 2 == FaceName ) {
      test = 1;
      value = 0.;
    }
    else if ( 3 == FaceName ) {
      test = 1;
      value = 0.;
    }
    else if ( 4 == FaceName ) {
      test = 1;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "P" ) ) {
    if ( 1 == FaceName ) {
      test = 0;
      value = 0.;
    }
    else if ( 2 == FaceName ) {
      test = 0;
      value = 0.;
    }
    else if ( 3 == FaceName ) {
      test = 0;
      value = 0.;
    }
    else if ( 4 == FaceName ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp ( name, "T" ) ) {
    if ( 1 == FaceName ) { //inflow
      test = 1;
      value = 1;
    }
    else if ( 2 == FaceName ) { //outflow
      test = 0;
      value = 0.;
    }
    else if ( 3 == FaceName ) { // no-slip fluid wall
      test = 0;
      value = 0.;
    }
    else if ( 4 == FaceName ) { // no-slip solid wall
      test = 1;
      value = 5.;
    }
  }

  return test;
}


bool SetBoundaryConditionCavityFlow ( const std::vector < double >& x, const char name[],
                                      double& value, const int FaceName, const double time ) {
  bool test = 1; //Dirichlet
  value = 0.;
  if ( !strcmp ( name, "V" ) ) {
    if ( 1 == FaceName ) {      //inflow
      test = 1;
      if ( x[1] < 0.5 && x[1] > -0.5 ) value = 1.; //4*(0.5-y)*(y+0.5);
    }
  }

  if ( !strcmp ( name, "P" ) ) {
    test = 0.;
    value = 0.;
    if ( x[0] < -.5 + 1.e-08 && x[1] < -.5 + 1.e-08 ) {
      test = 1;
    }
  }

  return test;
}


static unsigned counter = 0;

void AssembleMatrixResNS ( MultiLevelProblem &ml_prob ) {

  adept::Stack & adeptStack = FemusInit::_adeptStack;

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  //pointers and references
  NonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<NonLinearImplicitSystem> ( "Navier-Stokes" );
  const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();

  MultiLevelSolution*	 ml_sol	               = ml_prob._ml_sol;
  Solution*	 mysolution  	               = ml_sol->GetSolutionLevel ( level );

  LinearEquationSolver*  myLinEqSolver       = my_nnlin_impl_sys._LinSolver[level];

  Mesh		*mymsh		=  ml_prob._ml_msh->GetLevel ( level );
  elem		*myel		=  mymsh->el;
  SparseMatrix	*myKK		=  myLinEqSolver->_KK;
  NumericVector 	*myRES		=  myLinEqSolver->_RES;

  const unsigned dim = mymsh->GetDimension();
  const unsigned nabla_dim = 3 * ( dim - 1 );
  const unsigned max_size = static_cast< unsigned > ( ceil ( pow ( 3, dim ) ) );

  // local objects
  vector<adept::adouble> SolVAR ( dim + 1 );
  vector<vector<adept::adouble> > GradSolVAR ( dim + 1 );
  vector<vector<adept::adouble> > NablaSolVAR ( dim + 1 );

  for ( int i = 0; i < dim + 1; i++ ) {
    GradSolVAR[i].resize ( dim );
    NablaSolVAR[i].resize ( nabla_dim );
  }

  vector <double > phi;
  vector <adept::adouble> gradphi;
  vector <adept::adouble> nablaphi;
  adept::adouble Weight;

  phi.reserve ( max_size );
  gradphi.reserve ( max_size * dim );
  nablaphi.reserve ( max_size * nabla_dim );

  vector <double > phi1;
  vector <adept::adouble> gradphi1;
  vector <adept::adouble> nablaphi1;
  adept::adouble Weight1;

  phi1.reserve ( max_size );
  gradphi1.reserve ( max_size * dim );
  nablaphi1.reserve ( max_size * nabla_dim );

  vector <vector < adept::adouble> > vx ( dim );
  vector <vector < adept::adouble> > vx_face ( dim );

  for ( int i = 0; i < dim; i++ ) {
    vx[i].reserve ( max_size );
    vx_face[i].resize ( 9 );
  }

  vector< vector< adept::adouble > > Soli ( dim + 1 );
  vector< vector< int > > dofsVAR ( dim + 1 );
  for ( int i = 0; i < dim + 1; i++ ) {
    Soli[i].reserve ( max_size );
    dofsVAR[i].reserve ( max_size );
  }

  vector< vector< double > > Rhs ( dim + 1 );
  vector< vector< adept::adouble > > aRhs ( dim + 1 );
  for ( int i = 0; i < dim + 1; i++ ) {
    aRhs[i].reserve ( max_size );
    Rhs[i].reserve ( max_size );
  }

  vector < int > dofsAll;
  dofsAll.reserve ( max_size * ( dim + 1 ) );

  vector < double > KKloc;
  KKloc.reserve ( dim * max_size * ( dim + 1 ) *dim * max_size * ( dim + 1 ) );

  vector < double > Jac;
  Jac.reserve ( dim * max_size * ( dim + 1 ) *dim * max_size * ( dim + 1 ) );

  // ------------------------------------------------------------------------
  // Physical parameters
  double rhof	 	= ml_prob.parameters.get<Fluid> ( "Fluid" ).get_density();
  double IRe 		= ml_prob.parameters.get<Fluid> ( "Fluid" ).get_IReynolds_number();
  double betans	= 1.;

  // gravity
  double _gravity[3] = {0., 0., 0.};

  double DRe = 1 + ( counter * counter ) * 5;
  //double DRe=150;
  IRe = ( DRe * ( counter + 1 ) < 10000 ) ? 1. / ( DRe * ( counter + 1 ) ) : 1. / 10000.;


  cout << "iteration=" << counter << " Reynolds Number = " << 1. / IRe << endl;
  counter++;
  // -----------------------------------------------------------------
  // space discretization parameters
  unsigned SolType2 = ml_sol->GetSolutionType ( ml_sol->GetIndex ( "U" ) );

  unsigned SolType1 = ml_sol->GetSolutionType ( ml_sol->GetIndex ( "P" ) );

  unsigned SolTypeVx = 2;

  // mesh and procs
  unsigned nel    = mymsh->GetNumberOfElements();
  unsigned igrid  = mymsh->GetLevel();
  unsigned iproc  = mymsh->processor_id();

  //----------------------------------------------------------------------------------
  //variable-name handling
  const char varname[4][3] = {"U", "V", "W", "P"};
  vector <unsigned> indexVAR ( dim + 1 );
  vector <unsigned> indVAR ( dim + 1 );
  vector <unsigned> SolType ( dim + 1 );

  for ( unsigned ivar = 0; ivar < dim; ivar++ ) {
    indVAR[ivar] = ml_sol->GetIndex ( &varname[ivar][0] );
    SolType[ivar] = ml_sol->GetSolutionType ( &varname[ivar][0] );
    indexVAR[ivar] = my_nnlin_impl_sys.GetSolPdeIndex ( &varname[ivar][0] );
  }
  indVAR[dim] = ml_sol->GetIndex ( &varname[3][0] );
  SolType[dim] = ml_sol->GetSolutionType ( &varname[3][0] );
  indexVAR[dim] = my_nnlin_impl_sys.GetSolPdeIndex ( &varname[3][0] );

  unsigned indLmbd = ml_sol->GetIndex ( "lmbd" );

  //----------------------------------------------------------------------------------

  start_time = clock();

  myKK->zero();

  // *** element loop ***
  for ( int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++ ) {

    unsigned kel        = iel;
    short unsigned kelt = mymsh->GetElementType ( kel );
    unsigned nve2       = mymsh->GetElementDofNumber ( kel, SolType2 );
    unsigned nve1       = mymsh->GetElementDofNumber ( kel, SolType1 );
    unsigned nveVx      = mymsh->GetElementDofNumber ( kel, SolTypeVx );

    // *******************************************************************************************************
    //cout<<SolType1<<" "<<SolType2<<" "<<SolTypeVx<<" "<<nve1<<" "<<nve2<<" "<<nveVx<<endl;
    //Rhs
    for ( int i = 0; i < dim; i++ ) {
      dofsVAR[i].resize ( nve2 );
      Soli[indexVAR[i]].resize ( nve2 );
      aRhs[indexVAR[i]].resize ( nve2 );
      Rhs[indexVAR[i]].resize ( nve2 );
    }
    dofsVAR[dim].resize ( nve1 );
    Soli[indexVAR[dim]].resize ( nve1 );
    aRhs[indexVAR[dim]].resize ( nve1 );
    Rhs[indexVAR[dim]].resize ( nve1 );

    dofsAll.resize ( 0 );

    KKloc.resize ( ( dim * nve2 + nve1 ) * ( dim * nve2 + nve1 ) );
    Jac.resize ( ( dim * nve2 + nve1 ) * ( dim * nve2 + nve1 ) );

    // ----------------------------------------------------------------------------------------

    // coordinates
    for ( int i = 0; i < dim; i++ ) {
      vx[i].resize ( nveVx );
    }
    for ( unsigned i = 0; i < nveVx; i++ ) {
      unsigned inode_Metis = mymsh->GetSolutionDof ( i, iel, SolTypeVx );
      for ( int j = 0; j < dim; j++ ) {
        //coordinates
        vx[j][i] = ( *mymsh->_topology->_Sol[j] ) ( inode_Metis );
      }
    }

    // Velocity
    for ( unsigned i = 0; i < nve2; i++ ) {
      unsigned inode_Metis = mymsh->GetSolutionDof ( i, iel, SolType2 );
      for ( int j = 0; j < dim; j++ ) {
        // velocity dofs
        Soli[indexVAR[j]][i] = ( *mysolution->_Sol[indVAR[j]] ) ( inode_Metis );
        dofsVAR[j][i] = myLinEqSolver->GetSystemDof ( indVAR[j], indexVAR[j], i, iel );
        aRhs[indexVAR[j]][i] = 0.;
      }
    }

    // Pressure dofs
    for ( unsigned i = 0; i < nve1; i++ ) {
      unsigned inode_Metis = mymsh->GetSolutionDof ( i, iel, SolType1 );
      Soli[indexVAR[dim]][i] = ( *mysolution->_Sol[indVAR[dim]] ) ( inode_Metis );
      dofsVAR[dim][i] = myLinEqSolver->GetSystemDof ( indVAR[dim], indexVAR[dim], i, iel );
      aRhs[indexVAR[dim]][i] = 0.;
    }

    // build dof composition
    for ( int idim = 0; idim < dim; idim++ ) {
      dofsAll.insert ( dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end() );
    }
    dofsAll.insert ( dofsAll.end(), dofsVAR[dim].begin(), dofsVAR[dim].end() );

    {

      adeptStack.new_recording();

      // *** Gauss point loop ***

      adept::adouble hk;
      for ( unsigned ig = 0; ig < mymsh->_finiteElement[kelt][SolType2]->GetGaussPointNumber(); ig++ ) {
        // *** get Jacobian and test function and test function derivatives in the moving frame***
        mymsh->_finiteElement[kelt][SolType2]->Jacobian ( vx, ig, Weight, phi, gradphi, nablaphi );
        mymsh->_finiteElement[kelt][SolType1]->Jacobian ( vx, ig, Weight1, phi1, gradphi1, nablaphi1 );

        for ( int i = 0; i < nve1; i++ ) {
          for ( int j = 0; j < nabla_dim; j++ ) {
            //cout<<nablaphi[i*nabla_dim+j]<<" ";
          }
          //cout<<endl;
        }
        //cout<<endl;


        if ( ig == 0 ) {
          double referenceElementArea[6] = {8, 1. / 6., 1., 4., 1., 2.};
          double GaussWeight = mymsh->_finiteElement[kelt][SolType2]->GetGaussWeight ( ig );
          double area = referenceElementArea[kelt] * Weight.value() / GaussWeight;
          hk = pow ( area, 1. / dim );
          //cout<<hk<<" ";
        }

        // velocity: solution, gradient and laplace
        for ( int i = 0; i < dim; i++ ) {
          SolVAR[i] = 0.;
          for ( int j = 0; j < dim; j++ ) {
            GradSolVAR[i][j] = 0.;
          }
          for ( int j = 0; j < nabla_dim; j++ ) {
            NablaSolVAR[i][j] = 0.;
          }
          for ( unsigned inode = 0; inode < nve2; inode++ ) {
            adept::adouble soli = Soli[indexVAR[i]][inode];
            SolVAR[i] += phi[inode] * soli;
            for ( int j = 0; j < dim; j++ ) {
              GradSolVAR[i][j] += gradphi[inode * dim + j] * soli;
            }
            for ( int j = 0; j < nabla_dim; j++ ) {
              NablaSolVAR[i][j] += nablaphi[inode * nabla_dim + j] * soli;
            }
          }
        }

        // pressure, solution and gradient
        SolVAR[dim] = 0.;
        for ( int j = 0; j < dim; j++ ) {
          GradSolVAR[dim][j] = 0.;
        }
        for ( unsigned inode = 0; inode < nve1; inode++ ) {
          adept::adouble soli = Soli[indexVAR[dim]][inode];
          SolVAR[dim] += phi1[inode] * soli;
          for ( int j = 0; j < dim; j++ ) {
            GradSolVAR[dim][j] += gradphi1[inode * dim + j] * soli;
          }
        }

        bool Tezduyare = 0;
        bool FrancaAndFrey = !Tezduyare;

        if ( Tezduyare ) {
          //BEGIN TAU_SUPG, TAU_PSPG EVALUATION
          // ******** T.E. Tezduyar, S. Mittal, S.E. Ray and R. Shih *****************************
          // Computer Methods in Applied Mechanics and Engineering 95 (1992) 221-242 North-Holland
          // *************************************************************************************
          // velocity
          vector < adept::adouble > u ( dim );
          for ( int ivar = 0; ivar < dim; ivar++ ) {
            u[ivar] = SolVAR[ivar];
          }

          // speed
          adept::adouble uL2Norm = 0.;
          for ( int ivar = 0; ivar < dim; ivar++ ) {
            uL2Norm += u[ivar] * u[ivar];
          }
          uL2Norm = sqrt ( uL2Norm );
          adept::adouble tauSupg = 0.;
          adept::adouble tauPspg = 0.;
          if ( uL2Norm / ( 2.*IRe ) > 1.0e-10 ) {
            // velocity direction s = u/|u|
            vector < adept::adouble > s ( dim );
            for ( int ivar = 0; ivar < dim; ivar++ )
              s[ivar] = u[ivar] / uL2Norm;

            // element lenght h(s) = 2. ( \sum_i |s . gradphi_i | )^(-1)
            adept::adouble h = 0.;
            for ( unsigned i = 0; i < nve2; i++ ) {
              adept::adouble sDotGradphi = 0.;
              for ( int ivar = 0; ivar < dim; ivar++ )
                sDotGradphi += s[ivar] * gradphi[i * dim + ivar];

              //Start WARNING!!! do not change the following if into h += fabs(sDotGradphi).
              //Adept does not like it, and convergence is bad
              if ( sDotGradphi.value() < 0. )
                h -= sDotGradphi;
              else
                h += sDotGradphi;
              //End WARNING
            }
            h = 2. / h;
            //tauSupg
            adept::adouble Reu   = 1. / 3.* ( uL2Norm * h ) / ( 2.*IRe );
            adept::adouble zReu  = ( Reu.value() < 1. ) ? Reu : 1.;
            tauSupg = h / ( 2.*uL2Norm ) * zReu;

            uL2Norm = 0.01;
            h = 2.*hk / sqrt ( acos ( -1. ) );
            Reu   = 1. / 3.* ( uL2Norm * h ) / ( 2.*IRe );
            zReu  = ( Reu.value() < 1. ) ? Reu : 1.;
            tauPspg = h / ( 2.*uL2Norm ) * zReu;
          }
          //END TAU_SUPG, TAU_PSPG EVALUATION

          //BEGIN FLUID ASSEMBLY
          {
            vector < adept::adouble > Res ( dim, 0. );
            for ( unsigned ivar = 0; ivar < dim; ivar++ ) {
              Res[ivar] += 0. - GradSolVAR[dim][ivar];
              for ( unsigned jvar = 0; jvar < dim; jvar++ ) {
                unsigned kvar;
                if ( ivar == jvar ) kvar = jvar;
                else if ( 1 == ivar + jvar ) kvar = dim;  // xy
                else if ( 2 == ivar + jvar ) kvar = dim + 2; // xz
                else if ( 3 == ivar + jvar ) kvar = dim + 1; // yz
                Res[ivar] += - SolVAR[jvar] * GradSolVAR[ivar][jvar]
                             + IRe * ( NablaSolVAR[ivar][jvar] + NablaSolVAR[jvar][kvar] );
              }
            }

            adept::adouble div_vel = 0.;
            for ( int ivar = 0; ivar < dim; ivar++ ) {
              div_vel += GradSolVAR[ivar][ivar];
            }

            //BEGIN redidual momentum block
            for ( unsigned i = 0; i < nve2; i++ ) {
              for ( unsigned ivar = 0; ivar < dim; ivar++ ) {
                adept::adouble Advection = 0.;
                adept::adouble Laplacian = 0.;
                adept::adouble phiSupg = 0.;
                for ( unsigned jvar = 0; jvar < dim; jvar++ ) {
                  unsigned kvar;
                  if ( ivar == jvar ) kvar = jvar;
                  else if ( 1 == ivar + jvar ) kvar = dim;  // xy
                  else if ( 2 == ivar + jvar ) kvar = dim + 2; // xz
                  else if ( 3 == ivar + jvar ) kvar = dim + 1; // yz
                  Advection += SolVAR[jvar] * GradSolVAR[ivar][jvar] * phi[i];
                  Laplacian += IRe * gradphi[i * dim + jvar] * ( GradSolVAR[ivar][jvar] + GradSolVAR[jvar][ivar] );
                  phiSupg   += ( SolVAR[jvar] * gradphi[i * dim + jvar] ) * tauSupg;
                }
                aRhs[indexVAR[ivar]][i] += ( - Advection - Laplacian
                                             + SolVAR[dim] * gradphi[i * dim + ivar]
                                             + Res[ivar] * phiSupg ) * Weight;
              }
            }
            //END redidual momentum block

            //BEGIN continuity block
            for ( unsigned i = 0; i < nve1; i++ ) {
              adept::adouble MinusGradPhi1DotRes = 0.;
              for ( int ivar = 0; ivar < dim; ivar++ ) {
                MinusGradPhi1DotRes += -gradphi1[i * dim + ivar] * Res[ivar] * tauPspg;
              }
              aRhs[indexVAR[dim]][i] += ( - ( -div_vel ) * phi1[i] + MinusGradPhi1DotRes ) * Weight;
            }
            //END continuity block
          }
          //END FLUID ASSEMBLY
        }
        else if ( FrancaAndFrey ) {
          //BEGIN TAU_SUPG EVALUATION
          // ********************************* Tau_Supg FRANCA and FREY **************************
          // Computer Methods in Applied Mechanics and Engineering 99 (1992) 209-233 North-Holland
          // *************************************************************************************

          // Element diameter free stability parameters for stabilized methods applied to fluids
          // by Leopoldo P. Franca and Alexandre L. Madureira,
          // Computer Methods in Applied Mechanics and Engineering, Vol. 105 (1993) 395-403

          // velocity
          // double Ck[6][3]={{1.,1.,1.},{1.,1.,1.},{1.,1.,1.},{0.5, 11./270., 11./270.},{0.,1./42.,1./42.},{1.,1.,1.}};

          vector < adept::adouble > a ( dim );
          for ( int ivar = 0; ivar < dim; ivar++ ) {
            a[ivar] = SolVAR[ivar];
          }

          // speed
          adept::adouble aL2Norm = 0.;
          for ( int ivar = 0; ivar < dim; ivar++ ) {
            aL2Norm += a[ivar] * a[ivar];
          }
          aL2Norm = sqrt ( aL2Norm );

          adept::adouble deltaSupg = 0.;

          double sqrtlambdak = ( *mysolution->_Sol[indLmbd] ) ( iel );
          adept::adouble tauSupg = 1. / ( sqrtlambdak * sqrtlambdak * 4.*IRe );
          adept::adouble Rek   = aL2Norm / ( 4.*sqrtlambdak * IRe );

          /*adept::adouble mk=std::min(1./3., 2.*Ck[kelt][SolType2]);
          adept::adouble Rek   = mk * aL2Norm * hk / ( 4.*IRe);
          adept::adouble tauSupg = (mk*hk*hk/2.) / ( 4.*IRe);
          */
          if ( Rek > 1.0e-15 ) {
            adept::adouble xiRek = ( Rek >= 1. ) ? 1. : Rek;

            tauSupg   = xiRek / ( aL2Norm * sqrtlambdak );
            deltaSupg = ( xiRek * aL2Norm ) / sqrtlambdak;

// 	      tauSupg = hk / (2.*aL2Norm)*xiRek;
// 	      double lambda=1.;
// 	      deltaSupg = lambda*aL2Norm*hk*xiRek;
          }

          //END TAU_SUPG EVALUATION ============


          //BEGIN FLUID ASSEMBLY ============
          {
            vector < adept::adouble > Res ( dim, 0. );
            for ( unsigned ivar = 0; ivar < dim; ivar++ ) {
              Res[ivar] += 0. - GradSolVAR[dim][ivar];
              for ( unsigned jvar = 0; jvar < dim; jvar++ ) {
                unsigned kvar;
                if ( ivar == jvar ) kvar = jvar;
                else if ( 1 == ivar + jvar ) kvar = dim;  // xy
                else if ( 2 == ivar + jvar ) kvar = dim + 2; // xz
                else if ( 3 == ivar + jvar ) kvar = dim + 1; // yz
                Res[ivar] += - SolVAR[jvar] * GradSolVAR[ivar][jvar] // inconsistent
                             + IRe * ( NablaSolVAR[ivar][jvar] + NablaSolVAR[jvar][kvar] ); //consistent
              }

            }

            adept::adouble div_vel = 0.;
            for ( int ivar = 0; ivar < dim; ivar++ ) {
              div_vel += GradSolVAR[ivar][ivar];
            }

            //BEGIN redidual momentum block
            for ( unsigned i = 0; i < nve2; i++ ) {
              for ( unsigned ivar = 0; ivar < dim; ivar++ ) {

                adept::adouble Advection = 0.;
                adept::adouble Laplacian = 0.;
                adept::adouble phiSupg = 0.;
                for ( unsigned jvar = 0; jvar < dim; jvar++ ) {
                  unsigned kvar;
                  if ( ivar == jvar ) kvar = jvar;
                  else if ( 1 == ivar + jvar ) kvar = dim;  // xy
                  else if ( 2 == ivar + jvar ) kvar = dim + 2; // xz
                  else if ( 3 == ivar + jvar ) kvar = dim + 1; // yz

                  Advection += SolVAR[jvar] * GradSolVAR[ivar][jvar] * phi[i];
                  Laplacian += IRe * gradphi[i * dim + jvar] * ( GradSolVAR[ivar][jvar] + GradSolVAR[jvar][ivar] );
                  phiSupg   += ( SolVAR[jvar] * gradphi[i * dim + jvar] ) * tauSupg;

                  aRhs[indexVAR[ivar]][i] += 	 Res[ivar] * ( -IRe * nablaphi[i * nabla_dim + jvar] ) * tauSupg * Weight; //only in least square
                  aRhs[indexVAR[jvar]][i] += 	 Res[ivar] * ( -IRe * nablaphi[i * nabla_dim + kvar] ) * tauSupg * Weight; //only in least square
                }
                aRhs[indexVAR[ivar]][i] += ( - Advection - Laplacian
                                             + ( SolVAR[dim] - deltaSupg * div_vel ) * gradphi[i * dim + ivar]
                                             + Res[ivar] * phiSupg ) * Weight;
              }
            }
            //END redidual momentum block

            //BEGIN continuity block
            for ( unsigned i = 0; i < nve1; i++ ) {
              adept::adouble MinusGradPhi1DotRes = 0.;
              for ( int ivar = 0; ivar < dim; ivar++ ) {
                MinusGradPhi1DotRes += -gradphi1[i * dim + ivar] * Res[ivar] * tauSupg;
              }
              aRhs[indexVAR[dim]][i] += ( - ( -div_vel ) * phi1[i] + MinusGradPhi1DotRes ) * Weight;
            }
            //END continuity block ===========================
          }
          //END FLUID ASSEMBLY ============
        }
      }
    }

    //BEGIN local to global assembly
    //copy adouble aRhs into double Rhs
    for ( unsigned i = 0; i < dim; i++ ) {
      for ( int j = 0; j < nve2; j++ ) {
        Rhs[indexVAR[i]][j] = aRhs[indexVAR[i]][j].value();
      }
    }
    for ( unsigned j = 0; j < nve1; j++ ) {
      Rhs[indexVAR[dim]][j] = aRhs[indexVAR[dim]][j].value();
    }
    for ( int i = 0; i < dim + 1; i++ ) {
      myRES->add_vector_blocked ( Rhs[indexVAR[i]], dofsVAR[i] );
    }

    //Store equations
    for ( int i = 0; i < dim; i++ ) {
      adeptStack.dependent ( &aRhs[indexVAR[i]][0], nve2 );
      adeptStack.independent ( &Soli[indexVAR[i]][0], nve2 );
    }
    adeptStack.dependent ( &aRhs[indexVAR[dim]][0], nve1 );
    adeptStack.independent ( &Soli[indexVAR[dim]][0], nve1 );
    adeptStack.jacobian ( &Jac[0] );
    unsigned nveAll = ( dim * nve2 + nve1 );
    for ( int inode = 0; inode < nveAll; inode++ ) {
      for ( int jnode = 0; jnode < nveAll; jnode++ ) {
        KKloc[inode * nveAll + jnode] = -Jac[jnode * nveAll + inode];
      }
    }
    myKK->add_matrix_blocked ( KKloc, dofsAll, dofsAll );
    adeptStack.clear_independents();
    adeptStack.clear_dependents();

    //END local to global assembly

  } //end list of elements loop

  myKK->close();
  myRES->close();

  // *************************************
  end_time = clock();
  AssemblyTime += ( end_time - start_time );
  // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}


//   void LU( double *a , const int n ){
//
//     for (int k = 0 ; k < n-1 ; k++) {
//       for ( int i = k+1 ; i < n ; i++) {
// 	a[i*n+k] /= a[k*n+k];
// 	for ( int j= k+1 ; j < n ; j++) {
// 	  a[i*n+j] -=  a[i*n+k]*a[k*n+j] ;
// 	}
//       }
//     }
//   }

void SetLambda ( MultiLevelSolution &mlSol, const unsigned &level, const  FEOrder &order, Operator operatorType ) {

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

  adept::Stack & adeptStack = FemusInit::_adeptStack;

  Solution *mysolution = mlSol.GetSolutionLevel ( level );
  Mesh *mymsh	=  mlSol._mlMesh->GetLevel ( level );
  elem *myel	=  mymsh->el;


  unsigned indLmbd = mlSol.GetIndex ( "lmbd" );

  const unsigned geoDim = mymsh->GetDimension();
  const unsigned nablaGoeDim = ( 3 * ( geoDim - 1 ) + ! ( geoDim - 1 ) );
  const unsigned max_size = static_cast< unsigned > ( ceil ( pow ( 3, geoDim ) ) );

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
    cout << "wrong operator name in SetLambda\n"
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
    }
    for ( unsigned i = 0; i < nveVx; i++ ) {
      unsigned inodeVx_Metis = mymsh->GetSolutionDof ( i, iel, SolTypeVx );
      for ( int j = 0; j < geoDim; j++ ) {
        //coordinates
        vx[j][i] = ( *mymsh->_topology->_Sol[j] ) ( inodeVx_Metis );
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
    cout << hk << endl;
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
      while ( error > 1.0e-10 && counter <= 500 ) {
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

    cout << lambdak*hk*hk << std::endl;
    mysolution->_Sol[indLmbd]->set ( iel, sqrt ( lambdak ) );
    //abort();
  } //end list of elements loop


  mysolution->_Sol[indLmbd]->close();
  // *************************************
  end_time = clock();
  GetLambdaTime += ( end_time - start_time );

  std::cout << "GetLambda Time = " << GetLambdaTime / CLOCKS_PER_SEC << std::endl;
  //abort();
}


//------------------------------------------------------------------------------------------------------------
void AssembleMatrixResT ( MultiLevelProblem &ml_prob ) {

  //pointers and references

  LinearImplicitSystem& mylin_impl_sys = ml_prob.get_system<LinearImplicitSystem> ( "Temperature" );
  const unsigned level = mylin_impl_sys.GetLevelToAssemble();

  Solution*      mysolution	       = ml_prob._ml_sol->GetSolutionLevel ( level );

  LinearEquationSolver*  mylsyspde     = mylin_impl_sys._LinSolver[level];
  Mesh*          mymsh		       = ml_prob._ml_msh->GetLevel ( level );
  elem*          myel		       = mymsh->el;
  SparseMatrix*  myKK		       = mylsyspde->_KK;
  NumericVector* myRES		       = mylsyspde->_RES;
  MultiLevelSolution* ml_sol           = ml_prob._ml_sol;


  //data
  const unsigned	dim	= mymsh->GetDimension();
  unsigned 		nel	= mymsh->GetNumberOfElements();
  unsigned 		igrid	= mymsh->GetLevel();
  unsigned 		iproc	= mymsh->processor_id();
  double		IPe	= 1. / ( ml_prob.parameters.get<Fluid> ( "Fluid" ).get_Peclet_number() );

  //solution variable
  unsigned SolIndex;
  unsigned SolPdeIndex;
  SolIndex = ml_sol->GetIndex ( "T" );
  SolPdeIndex = mylin_impl_sys.GetSolPdeIndex ( "T" );
  //solution order
  unsigned order_ind = ml_sol->GetSolutionType ( SolIndex );
  //unsigned end_ind   = order_ind;

  //coordinates
  vector< vector < double> > coordinates ( dim );
  //const char coordinate_name[3][2] = {"X","Y","Z"};
  //vector < unsigned > coordinate_Index(dim);
//   for(unsigned ivar=0; ivar<dim; ivar++) {
//     coordinate_Index[ivar]=ivar;//ml_prob.GetIndex(coordinate_name[ivar]);
//   }

  // declare
  vector< int > metis_node;
  vector< int > KK_dof;
  vector <double> phi;
  vector <double> gradphi;
  vector <double> nablaphi;

  double weight;
  vector< double > F;
  vector< double > B;

  // reserve
  const unsigned max_size = static_cast< unsigned > ( ceil ( pow ( 3, dim ) ) );
  metis_node.reserve ( max_size );
  KK_dof.reserve ( max_size );
  for ( int i = 0; i < dim; i++ )
    coordinates[i].reserve ( max_size );
  phi.reserve ( max_size );
  gradphi.reserve ( max_size * dim );
  nablaphi.reserve ( max_size * ( 3 * ( dim - 1 ) ) );
  F.reserve ( max_size );
  B.reserve ( max_size * max_size );

  // Set to zeto all the entries of the Global Matrix
  myKK->zero();

  // *** element loop ***

  for ( int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++ ) {

    unsigned kel = iel;
    short unsigned kelt = mymsh->GetElementType ( kel );
    unsigned nve = mymsh->GetElementDofNumber ( kel, order_ind );

    // resize
    metis_node.resize ( nve );
    KK_dof.resize ( nve );
    phi.resize ( nve );
    gradphi.resize ( nve * dim );
    nablaphi.resize ( nve * ( 3 * ( dim - 1 ) ) );
    for ( int i = 0; i < dim; i++ ) {
      coordinates[i].resize ( nve );
    }

    // set to zero all the entries of the FE matrices
    F.resize ( nve );
    memset ( &F[0], 0, nve * sizeof ( double ) );
    {
      B.resize ( nve * nve );
      memset ( &B[0], 0, nve * nve * sizeof ( double ) );
    }

    // get local to global mappings
    for ( unsigned i = 0; i < nve; i++ ) {
      unsigned inode_coord_metis = mymsh->GetSolutionDof ( i, iel, 2 );
      metis_node[i] = mymsh->GetSolutionDof ( i, iel, order_ind );
      for ( unsigned ivar = 0; ivar < dim; ivar++ ) {
        coordinates[ivar][i] = ( *mymsh->_topology->_Sol[ivar] ) ( inode_coord_metis );
      }
      KK_dof[i] = mylsyspde->GetSystemDof ( SolIndex, SolPdeIndex, i, iel );
    }

    {
      // *** Gauss poit loop ***
      for ( unsigned ig = 0; ig < ml_prob._ml_msh->_finiteElement[kelt][order_ind]->GetGaussPointNumber(); ig++ ) {
        // *** get Jacobian and test function and test function derivatives ***
        ml_prob._ml_msh->_finiteElement[kelt][order_ind]->Jacobian ( coordinates, ig, weight, phi, gradphi, nablaphi );
        //Temperature and velocity current solution
        double SolT = 0;
        vector < double > gradSolT ( dim, 0. );
        for ( unsigned ivar = 0; ivar < dim; ivar++ ) {
          gradSolT[ivar] = 0;
        }
        vector < double > SolU ( dim, 0. );
        vector < unsigned > SolIndexU ( dim );
        SolIndexU[0] = ml_sol->GetIndex ( "U" );
        SolIndexU[1] = ml_sol->GetIndex ( "V" );
        if ( dim == 3 ) SolIndexU[2] = ml_sol->GetIndex ( "W" );

        unsigned SolType = ml_sol->GetSolutionType ( "T" );
        for ( unsigned i = 0; i < nve; i++ ) {
          double soli = ( *mysolution->_Sol[SolIndex] ) ( metis_node[i] );
          SolT += phi[i] * soli;
          for ( unsigned ivar = 0; ivar < dim; ivar++ ) gradSolT[ivar] += gradphi[i * dim + ivar] * soli;
          for ( int j = 0; j < dim; j++ )  {
            SolU[j] += phi[i] * ( *mysolution->_Sol[SolIndexU[j]] ) ( metis_node[i] );
          }
        }
        // *** phi_i loop ***
        for ( unsigned i = 0; i < nve; i++ ) {
          //BEGIN RESIDUALS A block ===========================
          double Adv_rhs = 0;
          double Lap_rhs = 0;
          for ( unsigned ivar = 0; ivar < dim; ivar++ ) {
            Lap_rhs += gradphi[i * dim + ivar] * gradSolT[ivar];
            Adv_rhs += SolU[ivar] * gradSolT[ivar];
          }
          F[i] += ( -IPe * Lap_rhs - Adv_rhs * phi[i] ) * weight;
          //END RESIDUALS A block ===========================
          {
            // *** phi_j loop ***
            for ( unsigned j = 0; j < nve; j++ ) {
              double Lap = 0;
              double Adv1 = 0;
              for ( unsigned ivar = 0; ivar < dim; ivar++ ) {
                // Laplacian
                Lap  += gradphi[i * dim + ivar] * gradphi[j * dim + ivar] * weight;
                // advection term I
                Adv1 += SolU[ivar] * gradphi[j * dim + ivar] * phi[i] * weight;
              }
              B[i * nve + j] += IPe * Lap + Adv1;
            } // end phij loop
          } // end phii loop
        } // endif assemble_matrix
      } // end gauss point loop
    } // endif single element not refined or fine grid loop
    //--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the global Matrix/vector

    myRES->add_vector_blocked ( F, KK_dof );
    myKK->add_matrix_blocked ( B, KK_dof, KK_dof );
  } //end list of elements loop for each subdomain

  myRES->close();
  myKK->close();

  // ***************** END ASSEMBLY *******************

}


