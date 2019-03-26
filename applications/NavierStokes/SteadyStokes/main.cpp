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
#include "LinearImplicitSystem.hpp"
#include "SolvertypeEnum.hpp"
#include "FElemTypeEnum.hpp"
#include "Files.hpp"

using std::cout;
using std::endl;

using namespace femus;

void AssembleMatrixResSteadyStokes(MultiLevelProblem &ml_prob);


double InitVariableU(const std::vector < double >& x);


bool SetBoundaryCondition(const std::vector < double >& x,const char name[],
			  double &value, const int FaceName, const double time);


int main(int argc,char **args) {

  bool Gmres=0, Asm=0;
  if(argc >= 2) {
    if( !strcmp("gmres",args[1])) 	Gmres=1;
    else if( !strcmp("asm",args[1])) 	Asm=1;

    if(Gmres+Asm==0) {
      cout << "wrong input arguments!" << endl;
      exit(0);
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
  nm=6;
  std::cout<<"MULTIGRID levels: "<< nm << endl;

  nr=0;
  std::cout<<"MAX_REFINEMENT levels: " << nr << endl<< endl;

  int tmp=nm;  nm+=nr;  nr=tmp;

  char *infile = new char [50];

  sprintf(infile,"./input/nsbenc.neu");

  //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;

  MultiLevelMesh ml_msh;
  ml_msh.ReadCoarseMesh(infile,"seventh",Lref);
  ml_msh.RefineMesh(nm,nr,NULL);

  ml_msh.PrintInfo();

  MultiLevelSolution ml_sol(&ml_msh);

  // generate solution vector
//   ml_sol.AddSolution("T",LAGRANGE,SECOND);
//   ml_sol.AddSolution("U",LAGRANGE,SECOND);
//   ml_sol.AddSolution("V",LAGRANGE,SECOND);
//   // the pressure variable should be the last for the Schur decomposition
//   ml_sol.AddSolution("P",DISCONTINUOUS_POLYNOMIAL,FIRST);

  ml_sol.AddSolution("T",LAGRANGE,FIRST);
  ml_sol.AddSolution("U",LAGRANGE,FIRST);
  ml_sol.AddSolution("V",LAGRANGE,FIRST);
  // the pressure variable should be the last for the Schur decomposition
  ml_sol.AddSolution("P",LAGRANGE,FIRST);
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


  //BEGIN Stokes Multilevel Problem
  std::cout << std::endl;
  std::cout << " *********** Stokes ************  " << std::endl;

  LinearImplicitSystem & system1 = ml_prob.add_system<LinearImplicitSystem> ("Stokes");
  system1.AddSolutionToSystemPDE("U");
  system1.AddSolutionToSystemPDE("V");
  system1.AddSolutionToSystemPDE("P");

  // Set MG Options
  system1.SetAssembleFunction(AssembleMatrixResSteadyStokes);
  system1.SetMaxNumberOfLinearIterations(2);
  system1.SetAbsoluteLinearConvergenceTolerance(1.e-10);
  system1.SetMgType(F_CYCLE);
  system1.SetNumberPreSmoothingStep(1);
  system1.SetNumberPostSmoothingStep(1);

  //Set Smoother Options
  if(Gmres) 		system1.SetLinearEquationSolverType(FEMuS_DEFAULT);
  else if(Asm) 		system1.SetLinearEquationSolverType(FEMuS_ASM);

  system1.init();
  //common smoother options
//   system1.AddStabilization(true);
  system1.SetSolverFineGrids(GMRES);
  system1.SetPreconditionerFineGrids(ILU_PRECOND);
  system1.SetTolerances(1.e-12,1.e-20,1.e+50,4);

  system1.ClearVariablesToBeSolved();
  //system1.AddVariableToBeSolved("All");
  system1.AddVariableToBeSolved("U");
  system1.AddVariableToBeSolved("V");
  system1.AddVariableToBeSolved("P");
  //for Vanka and ASM smoothers
  system1.SetNumberOfSchurVariables(0);
  system1.SetElementBlockNumber(4);
  //system1.SetElementBlockNumber("All",1);
  //for Gmres smoother
  system1.SetDirichletBCsHandling(PENALTY);
  //system1.SetDirichletBCsHandling(ELIMINATION);

  // Solve Navier-Stokes system
  ml_prob.get_system("Navier-Stokes").SetOuterSolver(PREONLY);
  ml_prob.get_system("Stokes").MGsolve();
  //END Stokes Multilevel Problem

  /// Print all solutions
  std::vector<std::string> print_vars;
  print_vars.push_back("U");
  print_vars.push_back("V");
  print_vars.push_back("P");


  VTKWriter vtkio(&ml_sol);
  vtkio.Write(files.GetOutputPath(),"biquadratic",print_vars);

//   XDMFWriter xdmfio(ml_prob);
//   xdmfio.write("biquadratic",print_vars);

//   GMVWriter gmvio(ml_sol);
//   gmvio.write("biquadratic",print_vars);

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

  return test;
}

// //------------------------------------------------------------------------------------------------------------


void AssembleMatrixResSteadyStokes(MultiLevelProblem &ml_prob){

  //pointers
  LinearImplicitSystem& my_lin_impl_sys    = ml_prob.get_system<LinearImplicitSystem>("Stokes");
  const unsigned level = my_lin_impl_sys.GetLevelToAssemble();

  Solution*	 mysolution  	             = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearEquationSolver*  mylsyspde	     = my_lin_impl_sys._LinSolver[level];
  const char* pdename                        = my_lin_impl_sys.name().c_str();

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
  bool penalty = true;
  const bool symm_mat = false;

  // solution and coordinate variables
  const char Solname[4][2] = {"U","V","W","P"};
  vector < unsigned > SolPdeIndex(dim+1);
  vector < unsigned > SolIndex(dim+1);

  vector< vector < double> > coordinates(dim);

  for(unsigned ivar=0; ivar<dim; ivar++) {
    SolPdeIndex[ivar]=my_lin_impl_sys.GetSolPdeIndex(&Solname[ivar][0]);
    SolIndex[ivar]=ml_sol->GetIndex(&Solname[ivar][0]);
  }
  SolPdeIndex[dim]=my_lin_impl_sys.GetSolPdeIndex(&Solname[3][0]);
  SolIndex[dim]=ml_sol->GetIndex(&Solname[3][0]);
  //solution order
  unsigned order_ind_vel = ml_sol->GetSolutionType(SolIndex[0]);
  unsigned order_ind_p = ml_sol->GetSolutionType(SolIndex[dim]);

  double alpha = 0.;
  if(order_ind_p == order_ind_vel && order_ind_vel == 0) // if pressure and velocity are both linear, we need stabilization
  {
    alpha = 0.013333;
  }

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
    unsigned nve2=mymsh->GetElementDofNumber(kel,order_ind_vel);
    unsigned nve1=mymsh->GetElementDofNumber(kel,order_ind_p);

    //set to zero all the entries of the FE matrices
    metis_node2.resize(nve2);
    node1.resize(nve1);
    phi2.resize(nve2);
    gradphi2.resize(nve2*dim);
    nablaphi2.resize(nve2*(3*(dim-1)));
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


    if(penalty){
      B[SolPdeIndex[dim]][SolPdeIndex[dim]].resize(nve1*nve1,0.);
      memset(&B[SolPdeIndex[dim]][SolPdeIndex[dim]][0],0,nve1*nve1*sizeof(double));
    }

    for( unsigned i=0;i<nve2;i++){
      unsigned inode_coord_metis=mymsh->GetSolutionDof(i, iel, 2);
      metis_node2[i]=mymsh->GetSolutionDof(i, iel, order_ind_vel);
      for(unsigned ivar=0; ivar<dim; ivar++) {
	coordinates[ivar][i]=(*mymsh->_topology->_Sol[ivar])(inode_coord_metis);
	KK_dof[ivar][i]=mylsyspde->GetSystemDof(SolIndex[ivar],SolPdeIndex[ivar],i, iel);
      }
    }

    double hk = sqrt( (coordinates[0][2] - coordinates[0][0])*(coordinates[0][2] - coordinates[0][0]) +
      (coordinates[1][2] - coordinates[1][0])*(coordinates[1][2] - coordinates[1][0]) );

    for(unsigned i=0;i<nve1;i++) {
      unsigned inode = mymsh->GetSolutionDof(i, kel, order_ind_p);
      node1[i]=inode;
      KK_dof[dim][i]=mylsyspde->GetSystemDof(SolIndex[dim],SolPdeIndex[dim],i, iel);
    }

    {
      // *** Gauss poit loop ***
      for(unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[kelt][order_ind_vel]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	ml_prob._ml_msh->_finiteElement[kelt][order_ind_vel]->Jacobian(coordinates,ig,Weight2,phi2,gradphi2,nablaphi2);
	phi1=ml_prob._ml_msh->_finiteElement[kelt][order_ind_p]->GetPhi(ig);

	double GradSolP[3] = {0.,0.,0.};
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
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){
	    GradSolP[ivar2] += gradphi2[i*dim+ivar2]*soli;
	  }
	}

// 	std::cout << "   0: " << GradSolP[0] << "    1: " << GradSolP[1] << std::endl;

	// *** phi_i loop ***
	for(unsigned i=0; i<nve2; i++){

	  //BEGIN RESIDUALS A block ===========================
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    double Lap_rhs=0;
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      Lap_rhs += gradphi2[i*dim+ivar2]*gradSolVAR[ivar][ivar2];
	    }
	    F[SolPdeIndex[ivar]][i]+= (-IRe*Lap_rhs + SolVAR[dim]*gradphi2[i*dim+ivar])*Weight2;
	  }
	  //END RESIDUALS A block ===========================

	  {
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve2; j++) {
	      double Lap=0;
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		// Laplacian
		Lap  += gradphi2[i*dim+ivar]*gradphi2[j*dim+ivar]*Weight2;
	      }

	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][i*nve2+j] += IRe*Lap;
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
	  F[SolPdeIndex[dim]][i]+= (phi1[i]*div + /*penalty*ILambda*phi1[i]*SolVAR[dim]*/
	                             + 0.*((hk*hk)/(4.*IRe))*alpha*(GradSolP[0]*gradphi2[i*dim + 0] + GradSolP[1]*gradphi2[i*dim + 1]) )*Weight2;


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

	if(penalty){  //block nve1 nve1
	  // *** phi_i loop ***
	  for(unsigned i=0; i<nve1; i++){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve1; j++){
	      //B[SolPdeIndex[dim]][SolPdeIndex[dim]][i*nve1+j]-= ILambda*phi1[i]*phi1[j]*Weight2;
	      for(unsigned ivar=0; ivar<dim; ivar++) {
	        B[SolPdeIndex[dim]][SolPdeIndex[dim]][i*nve1+j] -= ((hk*hk)/(4.*IRe))*alpha*(gradphi2[i*dim + ivar]*gradphi2[j*dim + ivar])*Weight2;
	      }
	    }
	  }
	}   //end if penalty
      }  // end gauss point loop

      //--------------------------------------------------------------------------------------------------------
      // Boundary Integral --> to be added
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
