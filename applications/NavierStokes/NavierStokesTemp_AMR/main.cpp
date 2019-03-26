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
#include "Files.hpp"

using std::cout;
using std::endl;
using namespace femus;

void AssembleMatrixResNS(MultiLevelProblem &ml_prob);
void AssembleMatrixResT(MultiLevelProblem &ml_prob);


double InitVariableU(const std::vector < double >& x);


bool SetBoundaryCondition(const std::vector < double >& x, const char name[],
			  double &value, const int FaceName, const double time);

int main(int argc,char **args) {

  bool Gmres=0, Asm=0;
  if(argc >= 3) {
    if( !strcmp("gmres",args[2]))       Gmres=1;
    else if( !strcmp("asm",args[2]))    Asm=1;

    if(Gmres+Asm==0) {
      cout << "wrong input arguments!" << endl;
      abort();
    }
  }
  else {
    cout << "No input argument set default smoother = Gmres" << endl;
    Gmres=1;
  }

  /// Init Petsc-MPI communicator
  FemusInit mpinit(argc,args,MPI_COMM_WORLD);

  Files files;
        files.CheckIODirectories();
        files.RedirectCout();

  /// INIT MESH =================================

  unsigned short nm,nr;
  nm=2;
  std::cout<<"MULTIGRID levels: "<< nm << endl;

  nr=0;
  std::cout<<"MAX_REFINEMENT levels: " << nr << endl<< endl;

  int tmp=nm;  nm+=nr;  nr=tmp;

  char *infile = new char [50];

  sprintf(infile,"./input/nsbenc.neu");

  //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;

  //Steadystate NonLinearMultiLevelProblem
  //MultiLevelMesh ml_msh(nm,nr,infile,"seventh",Lref,SetRefinementFlag);
  MultiLevelMesh ml_msh;
  ml_msh.ReadCoarseMesh(infile,"seventh",Lref);
  ml_msh.RefineMesh(nm,nr,NULL);

  // ml_msh.EraseCoarseLevels(2);

  MultiLevelSolution ml_sol(&ml_msh);

  // generate solution vector
  ml_sol.AddSolution("T",LAGRANGE,SECOND);
  ml_sol.AddSolution("U",LAGRANGE,SECOND);
  ml_sol.AddSolution("V",LAGRANGE,SECOND);
  // the pressure variable should be the last for the Schur decomposition
  ml_sol.AddSolution("P",DISCONTINUOUS_POLYNOMIAL,FIRST);
  ml_sol.AssociatePropertyToSolution("P","Pressure");

  //Initialize (update Init(...) function)
  ml_sol.Initialize("U",InitVariableU);
  ml_sol.Initialize("V");
  ml_sol.Initialize("P");
  ml_sol.Initialize("T");

  //Set Boundary (update Dirichlet(...) function)
  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  ml_sol.GenerateBdc("U");
  ml_sol.GenerateBdc("V");
  ml_sol.GenerateBdc("P");
  ml_sol.GenerateBdc("T");

  MultiLevelProblem ml_prob(&ml_sol);

  // add fluid material
  Parameter parameter(Lref,Uref);

  // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  Fluid fluid(parameter,0.001,1,"Newtonian",0.001,1.);
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;

  ml_prob.parameters.set<Fluid>("Fluid") = fluid;


  //BEGIN Navier-Stokes Multilevel Problem
  std::cout << std::endl;
  std::cout << " *********** Navier-Stokes ************  " << std::endl;

  NonLinearImplicitSystem & system1 = ml_prob.add_system<NonLinearImplicitSystem> ("Navier-Stokes");
  system1.AddSolutionToSystemPDE("U");
  system1.AddSolutionToSystemPDE("V");
  system1.AddSolutionToSystemPDE("P");

  // Set MG Options
  system1.SetAssembleFunction(AssembleMatrixResNS);
  system1.SetMaxNumberOfNonLinearIterations(3);
  system1.SetMaxNumberOfLinearIterations(2);
  system1.SetAbsoluteLinearConvergenceTolerance(1.e-10);
  system1.SetNonLinearConvergenceTolerance(1.e-04);
  system1.SetMgType(F_CYCLE);
  system1.SetNumberPreSmoothingStep(1);
  system1.SetNumberPostSmoothingStep(1);

  //Set Smoother Options
  if(Gmres) 		system1.SetLinearEquationSolverType(FEMuS_DEFAULT);
  else if(Asm) 		system1.SetLinearEquationSolverType(FEMuS_ASM);

  system1.init();

  std::string  AMR  = "yes";
  unsigned int maxAMRlevels 	= 6;
  std::string  AMRnorm="l2";
  double       AMRthreshold 	=0.001;

  system1.SetAMRSetOptions(AMR,maxAMRlevels,AMRnorm,AMRthreshold);

  //common smoother options
  system1.SetSolverFineGrids(GMRES);
  system1.SetPreconditionerFineGrids(ILU_PRECOND);
  system1.SetTolerances(1.e-12,1.e-20,1.e+50,4);
  //for Vanka and ASM smoothers
  system1.ClearVariablesToBeSolved();
  system1.AddVariableToBeSolved("All");
  //system1.AddVariableToBeSolved("U");
  //system1.AddVariableToBeSolved("V");
  //system1.AddVariableToBeSolved("P");
  system1.SetNumberOfSchurVariables(1);
  //system1.SetElementBlockNumber(4);
  system1.SetElementBlockNumber("All",1);
  //for Gmres smoother
  system1.SetDirichletBCsHandling(PENALTY);

  // Solve Navier-Stokes system
  ml_prob.get_system("Navier-Stokes").SetOuterSolver(PREONLY);
  ml_prob.get_system("Navier-Stokes").MGsolve();
  //END Navier-Stokes Multilevel Problem


//   //BEGIN Temperature MultiLevel Problem
//   std::cout << std::endl;
//   std::cout << " *********** Temperature ************* " << std::endl;
//
//   LinearImplicitSystem & system2 = ml_prob.add_system<LinearImplicitSystem> ("Temperature");
//   system2.AddSolutionToSystemPDE("T");
//
//
//   // Set MG Options
//   system2.SetAssembleFunction(AssembleMatrixResT);
//   system2.SetMaxNumberOfLinearIterations(6);
//   system2.SetAbsoluteLinearConvergenceTolerance(1.e-9);
//   system2.SetMgType(V_CYCLE);
//   system2.SetNumberPreSmoothingStep(1);
//   system2.SetNumberPostSmoothingStep(1);
//
//   //Set Smoother Options
//   if(Gmres) 		system2.SetLinearEquationSolverType(FEMuS_DEFAULT);
//   else if(Asm) 		system2.SetLinearEquationSolverType(FEMuS_ASM);
//   else if(Vanka)	system2.SetLinearEquationSolverType(VANKA_SMOOTHER);
//
//   system2.init();
//   //common smoother option
//   system2.SetSolverFineGrids(GMRES);
//   system2.SetTolerances(1.e-12,1.e-20,1.e+50,4);
//   system2.SetPreconditionerFineGrids(ILU_PRECOND);
//   //for Vanka and ASM smoothers
//   system2.ClearVariablesToBeSolved();
//   system2.AddVariableToBeSolved("All");
//   system2.SetNumberOfSchurVariables(0);
//   system2.SetElementBlockNumber(4);
//   //for Gmres smoother
//   system2.SetDirichletBCsHandling(PENALTY);
//
//   // Solve Temperature system
//   ml_prob.get_system("Temperature").solve();
//   //END Temperature Multilevel Problem
//
//   double l2normvarU = ml_sol.GetSolutionLevel(3)->GetSolutionName("U")->l2_norm();
//
//   double l2normvarUStored = 16.313927822836003;
//
//   std::cout << "Solution U l2norm: " << l2normvarU << std::endl;
//
//   if( fabs((l2normvarU - l2normvarUStored )/l2normvarUStored) > 1.e-6)
//   {
//     exit(1);
//   }
//
//   double l2normvarV = ml_sol.GetSolutionLevel(3)->GetSolutionName("V")->l2_norm();
//
//   double l2normvarVStored = 6.0644257018060355;
//
//   std::cout << "Solution V l2norm: " << l2normvarV << std::endl;
//
//   if( fabs((l2normvarV - l2normvarVStored )/l2normvarVStored )> 1.e-6)
//   {
//     exit(1);
//   }
//
//   double l2normvarP = ml_sol.GetSolutionLevel(3)->GetSolutionName("P")->l2_norm();
//
//   double l2normvarPStored = 1.8202105018866834;
//
//   std::cout << "Solution P l2norm: " << l2normvarP << std::endl;
//
//   if( fabs((l2normvarP - l2normvarPStored )/l2normvarPStored) > 1.e-6)
//   {
//     exit(1);
//   }
//
//   double l2normvarT = ml_sol.GetSolutionLevel(3)->GetSolutionName("T")->l2_norm();
//
//   double l2normvarTStored = 219.68194612060503;
//
//   std::cout << "Solution T l2norm: " << l2normvarT <<std::endl;
//
//   if( fabs((l2normvarT - l2normvarTStored )/l2normvarTStored) > 1.e-6)
//   {
//     exit(1);
//   }

  std::vector<std::string> print_vars;
  print_vars.push_back("U");
  print_vars.push_back("V");
  print_vars.push_back("P");
  print_vars.push_back("T");

  GMVWriter gmvio(&ml_sol);
  gmvio.Write(files.GetOutputPath(),"biquadratic",print_vars);



  //Destroy all the new systems
  ml_prob.clear();

  delete [] infile;
  return 0;
}

//--------------------------------------------------------------------------------------------------------------

double InitVariableU(const std::vector < double >& x) {
   double um = 0.2;
   double  value=1.5*um*(4.0/(0.1681))*x[1]*(0.41-x[1]);
   return value;
}

//-------------------------------------------------------------------------------------------------------------------

bool SetBoundaryCondition(const std::vector < double >& x,const char name[],
			  double &value, const int FaceName, const double time){
  bool test=1; //Dirichlet
  value=0.;
  //   cout << "Time bdc : " <<  time << endl;
  if(!strcmp(name,"U")) {
    if(1==FaceName){   //inflow
      test=1;
      double um = 0.2; // U/Uref
      value=1.5*0.2*(4.0/(0.1681))*x[1]*(0.41-x[1]);
    }
    else if(2==FaceName ){  //outflow
      test=0;
      //    test=1;
      value=0.;
    }
    else if(3==FaceName ){  // no-slip fluid wall
      test=1;
      value=0.;
    }
    else if(4==FaceName ){  // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if(!strcmp(name,"V")){
    if(1==FaceName){            //inflow
      test=1;
      value=0.;
    }
    else if(2==FaceName ){      //outflow
      test=0;
      //    test=1;
      value=0.;
    }
    else if(3==FaceName ){      // no-slip fluid wall
      test=1;
      value=0;
    }
    else if(4==FaceName ){      // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if(!strcmp(name,"W")){
    if(1==FaceName){
      test=1;
      value=0.;
    }
    else if(2==FaceName ){
      test=1;
      value=0.;
    }
    else if(3==FaceName ){
      test=1;
      value=0.;
    }
    else if(4==FaceName ){
      test=1;
      value=0.;
    }
  }
  else if(!strcmp(name,"P")){
    if(1==FaceName){
      test=0;
      value=0.;
    }
    else if(2==FaceName ){
      test=0;
      value=0.;
    }
    else if(3==FaceName ){
      test=0;
      value=0.;
    }
    else if(4==FaceName ){
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"T")) {
    if(1==FaceName){   //inflow
      test=1;
      value=1;
    }
    else if(2==FaceName ){  //outflow
      test=0;
      value=0.;
    }
    else if(3==FaceName ){  // no-slip fluid wall
      test=0;
      value=0.;
    }
    else if(4==FaceName ){  // no-slip solid wall
      test=1;
      value=5.;
    }
  }

  return test;
}

// //------------------------------------------------------------------------------------------------------------


void AssembleMatrixResNS(MultiLevelProblem &ml_prob){

  //pointers
  NonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<NonLinearImplicitSystem>("Navier-Stokes");
  const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();

  Solution*	 mysolution  	             = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearEquationSolver*  mylsyspde	     = my_nnlin_impl_sys._LinSolver[level];
  const char* pdename                        = my_nnlin_impl_sys.name().c_str();

  MultiLevelSolution* ml_sol=ml_prob._ml_sol;


  Mesh*		 mymsh    	= ml_prob._ml_msh->GetLevel(level);
  elem*		 myel		= mymsh->el;
  SparseMatrix*	 myKK		= mylsyspde->_KK;
  NumericVector* myRES 		= mylsyspde->_RES;

  //data
  const unsigned dim = mymsh->GetDimension();
  unsigned nel= mymsh->GetNumberOfElements();
  unsigned igrid= mymsh->GetLevel();
  unsigned iproc = mymsh->processor_id();
  double ILambda= 0;
  double IRe = ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
  bool penalty = true; //mylsyspde->GetStabilization();
  const bool symm_mat = false;//mylsyspde->GetMatrixProperties();
  const bool NavierStokes = true;
  unsigned nwtn_alg = 2;
  bool newton = (nwtn_alg==0) ? 0:1;

  // solution and coordinate variables
  const char Solname[4][2] = {"U","V","W","P"};
  vector < unsigned > SolPdeIndex(dim+1);
  vector < unsigned > SolIndex(dim+1);

  //const char coordinate_name[3][2] = {"X","Y","Z"};
  //vector < unsigned > coordinate_Index(dim);
  vector< vector < double> > coordinates(dim);

  for(unsigned ivar=0; ivar<dim; ivar++) {
    SolPdeIndex[ivar]=my_nnlin_impl_sys.GetSolPdeIndex(&Solname[ivar][0]);
    SolIndex[ivar]=ml_sol->GetIndex(&Solname[ivar][0]);
    //coordinate_Index[ivar]=ivar;//ml_prob.GetIndex(&coordinate_name[ivar][0]);
  }
  SolPdeIndex[dim]=my_nnlin_impl_sys.GetSolPdeIndex(&Solname[3][0]);
  SolIndex[dim]=ml_sol->GetIndex(&Solname[3][0]);
  //solution order
  unsigned order_ind2 = ml_sol->GetSolutionType(SolIndex[0]);
  unsigned order_ind1 = ml_sol->GetSolutionType(SolIndex[dim]);

  // declare
  vector < int > metis_node2;
  vector < int > node1;
  vector< vector< int > > KK_dof(dim+1);
  vector <double> phi2;
  vector <double> gradphi2;
  vector <double> nablaphi2;
  const double *phi1;
  double Weight2;
  double normal[3];
  vector< vector< double > > F(dim+1);
  vector< vector< vector< double > > > B(dim+1);

  // reserve
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  metis_node2.reserve(max_size);
  node1.reserve( static_cast< unsigned > (ceil(pow(2,dim))));
  for(int i=0;i<dim;i++) {
    coordinates[i].reserve(max_size);
  }
  phi2.reserve(max_size);
  gradphi2.reserve(max_size*dim);
  nablaphi2.reserve(max_size*(3*(dim-1)));
  for(int i=0;i<dim;i++) {
    KK_dof[i].reserve(max_size);
  }

  for(int i=0;i<dim+1;i++) F[i].reserve(max_size);

  {
    for(int i=0;i<dim+1;i++){
      B[i].resize(dim+1);
      for(int j=0;j<dim+1;j++){
	B[i][j].reserve(max_size*max_size);
      }
    }
  }

  vector < double > SolVAR(dim+1);
  vector < vector < double > > gradSolVAR(dim);
  for(int i=0;i<dim;i++) {
    gradSolVAR[i].resize(dim);
  }

  // Set to zeto all the entries of the matrix
  myKK->zero();

  // *** element loop ***

  for (int iel=mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc+1]; iel++) {

    unsigned kel = iel;
    short unsigned kelt=mymsh->GetElementType(kel);
    unsigned nve2=mymsh->GetElementDofNumber(kel,order_ind2);
    unsigned nve1=mymsh->GetElementDofNumber(kel,order_ind1);

    //set to zero all the entries of the FE matrices
    metis_node2.resize(nve2);
    node1.resize(nve1);
    phi2.resize(nve2);
    nablaphi2.resize(nve2*(3*(dim-1)));
    gradphi2.resize(nve2*dim);
    for(int ivar=0; ivar<dim; ivar++) {
      coordinates[ivar].resize(nve2);
      KK_dof[ivar].resize(nve2);

      F[SolPdeIndex[ivar]].resize(nve2);
      memset(&F[SolPdeIndex[ivar]][0],0,nve2*sizeof(double));

      {
	B[SolPdeIndex[ivar]][SolPdeIndex[ivar]].resize(nve2*nve2);
	B[SolPdeIndex[ivar]][SolPdeIndex[dim]].resize(nve2*nve1);
	B[SolPdeIndex[dim]][SolPdeIndex[ivar]].resize(nve1*nve2);
	memset(&B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][0],0,nve2*nve2*sizeof(double));
	memset(&B[SolPdeIndex[ivar]][SolPdeIndex[dim]][0],0,nve2*nve1*sizeof(double));
	memset(&B[SolPdeIndex[dim]][SolPdeIndex[ivar]][0],0,nve1*nve2*sizeof(double));
      }
    }
    KK_dof[dim].resize(nve1);
    F[SolPdeIndex[dim]].resize(nve1);
    memset(&F[SolPdeIndex[dim]][0],0,nve1*sizeof(double));


    if(nwtn_alg==2){
      for(int ivar=0; ivar<dim; ivar++) {
	for(int ivar2=1; ivar2<dim; ivar2++) {
	  B[SolPdeIndex[ivar]][SolPdeIndex[(ivar+ivar2)%dim]].resize(nve2*nve2);
	  memset(&B[SolPdeIndex[ivar]][SolPdeIndex[(ivar+ivar2)%dim]][0],0,nve2*nve2*sizeof(double));
	}
      }
    }

    if(penalty){
      B[SolPdeIndex[dim]][SolPdeIndex[dim]].resize(nve1*nve1,0.);
      memset(&B[SolPdeIndex[dim]][SolPdeIndex[dim]][0],0,nve1*nve1*sizeof(double));
    }

    for( unsigned i=0;i<nve2;i++){
      unsigned inode_metis=mymsh->GetSolutionDof(i, iel, 2);
      metis_node2[i]=inode_metis;
      for(unsigned ivar=0; ivar<dim; ivar++) {
	coordinates[ivar][i]=(*mymsh->_topology->_Sol[ivar])(inode_metis);
	KK_dof[ivar][i]=mylsyspde->GetSystemDof(SolIndex[ivar],SolPdeIndex[ivar],i, iel);
      }
    }

    for(unsigned i=0;i<nve1;i++) {
      unsigned inode=mymsh->GetSolutionDof(i, iel, order_ind1);
      node1[i]=inode;
      KK_dof[dim][i]=mylsyspde->GetSystemDof(SolIndex[dim],SolPdeIndex[dim],i, iel);
    }

    {
      // *** Gauss poit loop ***
      for(unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->Jacobian(coordinates,ig,Weight2,phi2,gradphi2,nablaphi2);
	phi1=ml_prob._ml_msh->_finiteElement[kelt][order_ind1]->GetPhi(ig);

	//velocity variable
	for(unsigned ivar=0; ivar<dim; ivar++) {
	  SolVAR[ivar]=0;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){
	    gradSolVAR[ivar][ivar2]=0;
	  }
	  unsigned SolIndex=ml_sol->GetIndex(&Solname[ivar][0]);
	  unsigned SolType=ml_sol->GetSolutionType(&Solname[ivar][0]);
	  for(unsigned i=0; i<nve2; i++) {
	    double soli = (*mysolution->_Sol[SolIndex])(metis_node2[i]);
	    SolVAR[ivar]+=phi2[i]*soli;
	    for(unsigned ivar2=0; ivar2<dim; ivar2++){
	      gradSolVAR[ivar][ivar2] += gradphi2[i*dim+ivar2]*soli;
	    }
	  }
	}
	//pressure variable
	SolVAR[dim]=0;
	unsigned SolIndex=ml_sol->GetIndex(&Solname[3][0]);
	unsigned SolType=ml_sol->GetSolutionType(&Solname[3][0]);
	for(unsigned i=0; i<nve1; i++){
	  unsigned sol_dof = node1[i];
	  double soli = (*mysolution->_Sol[SolIndex])(sol_dof);
	  SolVAR[dim]+=phi1[i]*soli;
	}

	// *** phi_i loop ***
	for(unsigned i=0; i<nve2; i++){

	  //BEGIN RESIDUALS A block ===========================
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    double Adv_rhs=0;
	    double Lap_rhs=0;
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      Lap_rhs += gradphi2[i*dim+ivar2]*gradSolVAR[ivar][ivar2];
	      Adv_rhs += SolVAR[ivar2]*gradSolVAR[ivar][ivar2];
	    }
	    F[SolPdeIndex[ivar]][i]+= (-IRe*Lap_rhs-NavierStokes*Adv_rhs*phi2[i]+SolVAR[dim]*gradphi2[i*dim+ivar])*Weight2;
	  }
	  //END RESIDUALS A block ===========================

	  {
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve2; j++) {
	      double Lap=0;
	      double Adv1=0;
	      double Adv2 = phi2[i]*phi2[j]*Weight2;
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		// Laplacian
		Lap  += gradphi2[i*dim+ivar]*gradphi2[j*dim+ivar]*Weight2;
		// advection term I
		Adv1 += SolVAR[ivar]*gradphi2[j*dim+ivar]*phi2[i]*Weight2;
	      }

	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][i*nve2+j] += IRe*Lap+ NavierStokes*newton*Adv1;
		if(nwtn_alg==2){
		  // Advection term II
		  B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][i*nve2+j]       += Adv2*gradSolVAR[ivar][ivar];
		  for(unsigned ivar2=1; ivar2<dim; ivar2++) {
		    B[SolPdeIndex[ivar]][SolPdeIndex[(ivar+ivar2)%dim]][i*nve2+j] += Adv2*gradSolVAR[ivar][(ivar+ivar2)%dim];
		  }
		}
	      }
  	    } //end phij loop

	    // *** phi1_j loop ***
	    for(unsigned j=0; j<nve1; j++){
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndex[ivar]][SolPdeIndex[dim]][i*nve1+j] -= gradphi2[i*dim+ivar]*phi1[j]*Weight2;
	      }
	    } //end phi1_j loop
	  } // endif assemble_matrix
	} //end phii loop


	// *** phi1_i loop ***
	for(unsigned i=0; i<nve1; i++){
	  //BEGIN RESIDUALS B block ===========================
	  double div = 0;
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    div += gradSolVAR[ivar][ivar];
	  }
	  F[SolPdeIndex[dim]][i]+= (phi1[i]*div +penalty*ILambda*phi1[i]*SolVAR[dim])*Weight2;
	  //END RESIDUALS  B block ===========================

	  {
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve2; j++) {
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndex[dim]][SolPdeIndex[ivar]][i*nve2+j]-= phi1[i]*gradphi2[j*dim+ivar]*Weight2;
	      }
	    }  //end phij loop
	  } // endif assemble_matrix
	}  //end phi1_i loop

	if( penalty){  //block nve1 nve1
	  // *** phi_i loop ***
	  for(unsigned i=0; i<nve1; i++){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve1; j++){
	      B[SolPdeIndex[dim]][SolPdeIndex[dim]][i*nve1+j]-= ILambda*phi1[i]*phi1[j]*Weight2;
	    }
	  }
	}   //end if penalty
      }  // end gauss point loop

      //--------------------------------------------------------------------------------------------------------
      // Boundary Integral --> to be addded
      //number of faces for each type of element
//       if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {
//
// 	unsigned nfaces = myel->GetElementFaceNumber(kel);
//
// 	// loop on faces
// 	for(unsigned jface=0;jface<nfaces;jface++){
//
// 	  // look for boundary faces
// 	  if(myel->GetFaceElementIndex(kel,jface)<0){
// 	    for(unsigned ivar=0; ivar<dim; ivar++) {
// 	      ml_prob.ComputeBdIntegral(pdename, &Solname[ivar][0], kel, jface, level, ivar);
// 	    }
// 	  }
// 	}
//       }
      //--------------------------------------------------------------------------------------------------------
    } // endif single element not refined or fine grid loop
    //--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the Global Matrix/Vector
    for(unsigned ivar=0; ivar<dim; ivar++) {
      myRES->add_vector_blocked(F[SolPdeIndex[ivar]],KK_dof[ivar]);
      {
	myKK->add_matrix_blocked(B[SolPdeIndex[ivar]][SolPdeIndex[ivar]],KK_dof[ivar],KK_dof[ivar]);
	myKK->add_matrix_blocked(B[SolPdeIndex[ivar]][SolPdeIndex[dim]],KK_dof[ivar],KK_dof[dim]);
	myKK->add_matrix_blocked(B[SolPdeIndex[dim]][SolPdeIndex[ivar]],KK_dof[dim],KK_dof[ivar]);
	if(nwtn_alg==2){
	  for(unsigned ivar2=1; ivar2<dim; ivar2++) {
	    myKK->add_matrix_blocked(B[SolPdeIndex[ivar]][SolPdeIndex[(ivar+ivar2)%dim]],KK_dof[ivar],KK_dof[(ivar+ivar2)%dim]);
	  }
	}
      }
    }
    //Penalty
    if(penalty) myKK->add_matrix_blocked(B[SolPdeIndex[dim]][SolPdeIndex[dim]],KK_dof[dim],KK_dof[dim]);
    myRES->add_vector_blocked(F[SolPdeIndex[dim]],KK_dof[dim]);
    //--------------------------------------------------------------------------------------------------------
  } //end list of elements loop for each subdomain


  myKK->close();
  myRES->close();
  // ***************** END ASSEMBLY *******************
}

//------------------------------------------------------------------------------------------------------------
void AssembleMatrixResT(MultiLevelProblem &ml_prob){

  //pointers and references
  LinearImplicitSystem& mylin_impl_sys = ml_prob.get_system<LinearImplicitSystem>("Temperature");
  const unsigned level = mylin_impl_sys.GetLevelToAssemble();

  Solution*      mysolution	       = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver*  mylsyspde     = mylin_impl_sys._LinSolver[level];
  Mesh*          mymsh		       = ml_prob._ml_msh->GetLevel(level);
  elem*          myel		       = mymsh->el;
  SparseMatrix*  myKK		       = mylsyspde->_KK;
  NumericVector* myRES		       = mylsyspde->_RES;
  MultiLevelSolution* ml_sol           = ml_prob._ml_sol;


  //data
  const unsigned	dim	= mymsh->GetDimension();
  unsigned 		nel	= mymsh->GetNumberOfElements();
  unsigned 		igrid	= mymsh->GetLevel();
  unsigned 		iproc	= mymsh->processor_id();
  double		IPe	= 1./(ml_prob.parameters.get<Fluid>("Fluid").get_Peclet_number());

  //solution variable
  unsigned SolIndex;
  unsigned SolPdeIndex;
  SolIndex=ml_sol->GetIndex("T");
  SolPdeIndex=mylin_impl_sys.GetSolPdeIndex("T");
  //solution order
  unsigned order_ind = ml_sol->GetSolutionType(SolIndex);
  //unsigned end_ind   = mymsh->GetEndIndex(order_ind);

  //coordinates
  vector< vector < double> > coordinates(dim);
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
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  metis_node.reserve(max_size);
  KK_dof.reserve(max_size);
  for(int i=0;i<dim;i++)
    coordinates[i].reserve(max_size);
  phi.reserve(max_size);
  gradphi.reserve(max_size*dim);
  nablaphi.reserve(max_size*(3*(dim-1)));
  F.reserve(max_size);
  B.reserve(max_size*max_size);

  // Set to zeto all the entries of the Global Matrix
  myKK->zero();

  // *** element loop ***

  for (int iel=mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc+1]; iel++) {

    unsigned kel = iel;
    short unsigned kelt = mymsh->GetElementType(kel);
    unsigned nve=mymsh->GetElementDofNumber(kel,order_ind);

    // resize
    metis_node.resize(nve);
    KK_dof.resize(nve);
    phi.resize(nve);
    gradphi.resize(nve*dim);
    nablaphi.resize(nve*(3*(dim-1)));
    for(int i=0;i<dim;i++){
      coordinates[i].resize(nve);
    }

    // set to zero all the entries of the FE matrices
    F.resize(nve);
    memset(&F[0],0,nve*sizeof(double));
    {
      B.resize(nve*nve);
      memset(&B[0],0,nve*nve*sizeof(double));
    }

    // get local to global mappings
    for( unsigned i=0;i<nve;i++){
      unsigned inode_metis=mymsh->GetSolutionDof(i, iel, 2);
      metis_node[i]=inode_metis;
      for(unsigned ivar=0; ivar<dim; ivar++) {
	coordinates[ivar][i]=(*mymsh->_topology->_Sol[ivar])(inode_metis);
      }
      KK_dof[i]=mylsyspde->GetSystemDof(SolIndex,SolPdeIndex,i, iel);
    }

    {
      // *** Gauss poit loop ***
      for(unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[kelt][order_ind]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	ml_prob._ml_msh->_finiteElement[kelt][order_ind]->Jacobian(coordinates,ig,weight,phi,gradphi,nablaphi);
	//Temperature and velocity current solution
	double SolT=0;
	vector < double > gradSolT(dim,0.);
	for(unsigned ivar=0; ivar<dim; ivar++){
	  gradSolT[ivar]=0;
	}
	vector < double > SolU(dim,0.);
	vector < unsigned > SolIndexU(dim);
	SolIndexU[0]=ml_sol->GetIndex("U");
	SolIndexU[1]=ml_sol->GetIndex("V");
	if(dim==3) SolIndexU[2]=ml_sol->GetIndex("W");

	unsigned SolType=ml_sol->GetSolutionType("T");
	for(unsigned i=0; i<nve; i++) {
	  double soli = (*mysolution->_Sol[SolIndex])(metis_node[i]);
	  SolT+=phi[i]*soli;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++) gradSolT[ivar2] += gradphi[i*dim+ivar2]*soli;
	  for(int j=0;j<dim;j++)  {
	    SolU[j]+=phi[i]*(*mysolution->_Sol[SolIndexU[j]])(metis_node[i]);
	  }
	}
	// *** phi_i loop ***
	for(unsigned i=0; i<nve; i++){
	  //BEGIN RESIDUALS A block ===========================
	  double Adv_rhs=0;
	  double Lap_rhs=0;
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    Lap_rhs += gradphi[i*dim+ivar]*gradSolT[ivar];
	    Adv_rhs += SolU[ivar]*gradSolT[ivar];
	  }
	  F[i]+= (-IPe*Lap_rhs-Adv_rhs*phi[i])*weight;
	  //END RESIDUALS A block ===========================
	  {
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve; j++) {
	      double Lap=0;
	      double Adv1=0;
	      for(unsigned ivar=0; ivar<dim; ivar++) {
	      // Laplacian
		Lap  += gradphi[i*dim+ivar]*gradphi[j*dim+ivar]*weight;
		// advection term I
		Adv1 += SolU[ivar]*gradphi[j*dim+ivar]*phi[i]*weight;
	      }
	      B[i*nve+j] += IPe*Lap + Adv1;
	    } // end phij loop
	  } // end phii loop
	} // endif assemble_matrix
      } // end gauss point loop
    } // endif single element not refined or fine grid loop
    //--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the global Matrix/vector

    myRES->add_vector_blocked(F,KK_dof);
    myKK->add_matrix_blocked(B,KK_dof,KK_dof);
  } //end list of elements loop for each subdomain

  myRES->close();
  myKK->close();

   // ***************** END ASSEMBLY *******************

}


