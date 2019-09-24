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
#include "NonLinearImplicitSystem.hpp"
#include "FElemTypeEnum.hpp"
#include "Files.hpp"

using std::cout;
using std::endl;
using namespace femus;

void AssembleMatrixResNS(MultiLevelProblem &ml_prob);

bool SetBoundaryCondition(const std::vector < double >& x, const char name[],
                          double &value, const int FaceName, const double time);

int main(int argc,char **args) {

  /// Init Petsc-MPI communicator
  FemusInit mpinit(argc,args,MPI_COMM_WORLD);

  Files files;
        files.CheckIODirectories();
        files.RedirectCout();

  /// INIT MESH =================================

  unsigned short nm,nr;
  nm=4;
  std::cout<<"MULTIGRID levels: "<< nm << endl;

  nr=0;
  std::cout<<"MAX_REFINEMENT levels: " << nr << endl<< endl;

  int tmp=nm;  nm+=nr;  nr=tmp;

  char *infile = new char [50];

  sprintf(infile,"./input/nsbenchreg.neu");

  //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;

  MultiLevelMesh ml_msh(nm,nr,infile,"seventh",Lref,NULL);

  MultiLevelSolution ml_sol(&ml_msh);

  // generate solution vector
  ml_sol.AddSolution("U",LAGRANGE,SECOND,2);
  ml_sol.AddSolution("V",LAGRANGE,SECOND,2);
  // the pressure variable should be the last for the Schur decomposition
  ml_sol.AddSolution("P",DISCONTINUOUS_POLYNOMIAL,FIRST,1);
  ml_sol.AssociatePropertyToSolution("P","Pressure");

  //Initialize (update Init(...) function)
  ml_sol.Initialize("All");

  //Set Boundary (update Dirichlet(...) function)
  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  ml_sol.GenerateBdc("U","Time_dependent");
  ml_sol.GenerateBdc("V");
  ml_sol.GenerateBdc("P");


  MultiLevelProblem ml_prob(&ml_sol);


  Parameter parameter(Lref,Uref);

  // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  Fluid fluid(parameter,0.001,1.,"Newtonian");
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;

  // add fluid material
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;



  //create systems
  // add the system Navier-Stokes to the MultiLevel problem
  TransientNonlinearImplicitSystem & system = ml_prob.add_system<TransientNonlinearImplicitSystem> ("Navier-Stokes");
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  system.AddSolutionToSystemPDE("P");

  // init all the systems
  system.init();

  // System Navier-Stokes
  system.SetAssembleFunction(AssembleMatrixResNS);
  system.SetMaxNumberOfLinearIterations(1);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-8);
  system.SetMgType(V_CYCLE);
  system.SetMaxNumberOfNonLinearIterations(15);

  // time loop parameter
  system.SetIntervalTime(0.1);
  const unsigned int n_timesteps = 20;
  const unsigned int write_interval = 1;

  for (unsigned time_step = 0; time_step < n_timesteps; time_step++) {

    // Solving Navier-Stokes system
    std::cout << std::endl;
    std::cout << " *********** Navier-Stokes ************  " << std::endl;
    ml_prob.get_system("Navier-Stokes").SetOuterSolver(PREONLY);
    ml_prob.get_system("Navier-Stokes").MGsolve();

    //update Solution
    ml_prob.get_system<TransientNonlinearImplicitSystem>("Navier-Stokes").CopySolutionToOldSolution();

    // print solution
    if ( !(time_step%write_interval) ) {

      //print solution
      std::vector<std::string> print_vars;
      print_vars.push_back("U");
      print_vars.push_back("V");
      print_vars.push_back("P");

//       ml_prob.printsol_vtu_inline("biquadratic",print_vars,time_step);
      VTKWriter vtkio(&ml_sol);
      vtkio.Write(files.GetOutputPath(),"biquadratic",print_vars,time_step);
    }

  } //end loop timestep


  // Destroy all the new systems
  ml_prob.clear();


  delete[] infile;
  return 0;
}

//-------------------------------------------------------------------------------------------------------------------

bool SetBoundaryCondition(const std::vector < double >& x,const char name[],
			  double &value, const int FaceName, const double time){
  bool test=1; //Dirichlet
  value=0.;

  if(!strcmp(name,"U")) {
    if (1==FaceName) { //inflow
      test=1;
      double um = 0.2; // U/Uref
      if(time < 2.0) {
        value=1.5*um*(4.0/(0.1681))*x[1]*(0.41-x[1])*0.5*(1. - cos(0.5*3.141592653589793*time) );
      } else {
        value=1.5*um*(4.0/(0.1681))*x[1]*(0.41-x[1]);
      }
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


//------------------------------------------------------------------------------------------------------------
void AssembleMatrixResNS(MultiLevelProblem &ml_prob){

  //pointers

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem>("Navier-Stokes");
  const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();

  MultiLevelSolution *ml_sol			      = ml_prob._ml_sol;
  Solution*	 mysolution  	                      = ml_sol->GetSolutionLevel(level);
  LinearEquationSolver*  mylsyspde	              = my_nnlin_impl_sys._LinSolver[level];
  const char* pdename                                 = my_nnlin_impl_sys.name().c_str();

  Mesh*		 mymsh    	   = ml_prob._ml_msh->GetLevel(level);
  elem*		 myel		   = mymsh->el;
  SparseMatrix*	 myKK	 	   = mylsyspde->_KK;
  NumericVector* myRES 		   = mylsyspde->_RES;

  //data
  double dt = my_nnlin_impl_sys.GetIntervalTime();
  double theta = 0.5;
  const unsigned dim = mymsh->GetDimension();
  unsigned nel= mymsh->GetNumberOfElements();
  unsigned igrid= mymsh->GetLevel();
  unsigned iproc = mymsh->processor_id();
  double ILambda = 0.;
  double IRe = ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
  bool penalty = true;         // mylsyspde->GetStabilization();
  const bool symm_mat = false; // mylsyspde->GetMatrixProperties();
  const bool NavierStokes = true;
  unsigned nwtn_alg = 1;
  bool newton = (nwtn_alg==0) ? 0:1;

  // solution and coordinate variables
  const char Solname[4][2] = {"U","V","W","P"};
  vector < unsigned > SolPdeIndex(dim+1);
  vector < unsigned > SolIndex(dim+1);

  vector< vector < double> > coordinates(dim);

  for(unsigned ivar=0; ivar<dim; ivar++) {
    SolPdeIndex[ivar]=my_nnlin_impl_sys.GetSolPdeIndex(&Solname[ivar][0]);
    SolIndex[ivar]=ml_sol->GetIndex(&Solname[ivar][0]);
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


  for(int i=0;i<dim+1;i++){
    B[i].resize(dim+1);
    for(int j=0;j<dim+1;j++){
      B[i][j].reserve(max_size*max_size);
    }
  }


  vector < double > Sol(dim+1);
  vector < vector < double > > gradSol(dim);
  for(int i=0;i<dim;i++) {
    gradSol[i].resize(dim);
  }

  vector < double > SolOld(dim+1);
  vector < vector < double > > gradSolOld(dim);
  for(int i=0;i<dim;i++) {
    gradSolOld[i].resize(dim);
  }
  //vector < double > AccSol(dim);

  // Set to zeto all the entries of the matrix
  myKK->zero();

  // *** element loop ***

  for (int iel=mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc+1]; iel++) {

    unsigned kel = iel;
    short unsigned kelt= mymsh->GetElementType(kel);
    unsigned nve2=mymsh->GetElementDofNumber(kel,order_ind2);
    unsigned nve1=mymsh->GetElementDofNumber(kel,order_ind1);

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


      B[SolPdeIndex[ivar]][SolPdeIndex[ivar]].resize(nve2*nve2);
      B[SolPdeIndex[ivar]][SolPdeIndex[dim]].resize(nve2*nve1);
      B[SolPdeIndex[dim]][SolPdeIndex[ivar]].resize(nve1*nve2);
      memset(&B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][0],0,nve2*nve2*sizeof(double));
      memset(&B[SolPdeIndex[ivar]][SolPdeIndex[dim]][0],0,nve2*nve1*sizeof(double));
      memset(&B[SolPdeIndex[dim]][SolPdeIndex[ivar]][0],0,nve1*nve2*sizeof(double));

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
	  Sol[ivar]=0.;
	  SolOld[ivar]=0.;
	  //AccSol[ivar]=0.;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){
	    gradSol[ivar][ivar2]=0;
	    gradSolOld[ivar][ivar2]=0;
	  }
	  unsigned SolIndex=ml_sol->GetIndex(&Solname[ivar][0]);
	  unsigned SolType =ml_sol->GetSolutionType(&Solname[ivar][0]);
	  for(unsigned i=0; i<nve2; i++) {
	    double soli    = (*mysolution->_Sol[SolIndex])(metis_node2[i]);
	    double sololdi = (*mysolution->_SolOld[SolIndex])(metis_node2[i]);
	    Sol[ivar]+=phi2[i]*soli;
	    SolOld[ivar]+=phi2[i]*sololdi;
	    for(unsigned ivar2=0; ivar2<dim; ivar2++){
	      gradSol[ivar][ivar2]    += gradphi2[i*dim+ivar2]*soli;
	      gradSolOld[ivar][ivar2] += gradphi2[i*dim+ivar2]*sololdi;
	    }
	  }
	}
	//pressure variable
	Sol[dim]=0;
	unsigned SolIndex=ml_sol->GetIndex(&Solname[3][0]);
	unsigned SolType=ml_sol->GetSolutionType(&Solname[3][0]);
	for(unsigned i=0; i<nve1; i++){
	  unsigned sol_dof = node1[i];
	  double soli = (*mysolution->_Sol[SolIndex])(sol_dof);
	  Sol[dim]+=phi1[i]*soli;
	}

	// *** phi_i loop ***
	for(unsigned i=0; i<nve2; i++){

	  //BEGIN RESIDUALS A block ===========================
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    double Adv_rhs=0;
	    double Adv_old_rhs=0;
	    double Lap_rhs=0;
	    double Lap_old_rhs=0;
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      Lap_rhs 	  += gradphi2[i*dim+ivar2]*gradSol[ivar][ivar2];
	      Lap_old_rhs += gradphi2[i*dim+ivar2]*gradSolOld[ivar][ivar2];
	      Adv_rhs 	  += Sol[ivar2]*gradSol[ivar][ivar2];
	      Adv_old_rhs += SolOld[ivar2]*gradSolOld[ivar][ivar2];
	    }
	    F[SolPdeIndex[ivar]][i]+= ( -theta*dt*IRe*Lap_rhs                           // Laplacian
					-(1.-theta)*dt*IRe*Lap_old_rhs                  // Laplacian
					-theta*dt*NavierStokes*Adv_rhs*phi2[i]          // advection
					-(1.-theta)*dt*NavierStokes*Adv_old_rhs*phi2[i] // advection
					+dt*Sol[dim]*gradphi2[i*dim+ivar]                  // pressure
					-(Sol[ivar] - SolOld[ivar])*phi2[i]        // acceleration
	    )*Weight2;
	  }
	  //END RESIDUALS A block ===========================

	  {
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve2; j++) {
	      double Lap=0;
	      double Adv1=0;
	      double Adv2 = phi2[i]*phi2[j]*Weight2;
	      double Mass = phi2[i]*phi2[j]*Weight2;
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		// Laplacian
		Lap  += gradphi2[i*dim+ivar]*gradphi2[j*dim+ivar]*Weight2;
		// advection term I
		Adv1 += Sol[ivar]*gradphi2[j*dim+ivar]*phi2[i]*Weight2;
	      }

	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][i*nve2+j] += theta*dt*IRe*Lap + theta*dt*NavierStokes*newton*Adv1 + Mass;
		if(nwtn_alg==2){
		  // Advection term II
		  B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][i*nve2+j]       += theta*dt*Adv2*gradSol[ivar][ivar];
		  for(unsigned ivar2=1; ivar2<dim; ivar2++) {
		    B[SolPdeIndex[ivar]][SolPdeIndex[(ivar+ivar2)%dim]][i*nve2+j] += theta*dt*Adv2*gradSol[ivar][(ivar+ivar2)%dim];
		  }
		}
	      }
  	    } //end phij loop

	    // *** phi1_j loop ***
	    for(unsigned j=0; j<nve1; j++){
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndex[ivar]][SolPdeIndex[dim]][i*nve1+j] -= dt*gradphi2[i*dim+ivar]*phi1[j]*Weight2;
	      }
	    } //end phi1_j loop
	  } // endif assemble_matrix
	} //end phii loop


	// *** phi1_i loop ***
	for(unsigned i=0; i<nve1; i++){
	  //BEGIN RESIDUALS B block ===========================
	  double div = 0;
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    div += gradSol[ivar][ivar];
	  }
	  F[SolPdeIndex[dim]][i]+= dt*(phi1[i]*div +penalty*ILambda*phi1[i]*Sol[dim])*Weight2;
	  //END RESIDUALS  B block ===========================

	  {
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve2; j++) {
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndex[dim]][SolPdeIndex[ivar]][i*nve2+j]-= dt*phi1[i]*gradphi2[j*dim+ivar]*Weight2;
	      }
	    }  //end phij loop
	  } // endif assemble_matrix
	}  //end phi1_i loop

	if(penalty){  //block nve1 nve1
	  // *** phi_i loop ***
	  for(unsigned i=0; i<nve1; i++){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve1; j++){
	      B[SolPdeIndex[dim]][SolPdeIndex[dim]][i*nve1+j]-=  dt*ILambda*phi1[i]*phi1[j]*Weight2;
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


/*

static char help[] = "Illustrates use of the preconditioner ASM.\n\
The Additive Schwarz Method for solving a linear system in parallel with KSP.  The\n\
code indicates the procedure for setting user-defined subdomains.  Input\n\
parameters include:\n\
  -user_set_subdomain_solvers:  User explicitly sets subdomain solvers\n\
  -user_set_subdomains:  Activate user-defined subdomains\n\n";

/*
   Note:  This example focuses on setting the subdomains for the ASM
   preconditioner for a problem on a 2D rectangular grid.  See ex1.c
   and ex2.c for more detailed comments on the basic usage of KSP
   (including working with matrices and vectors).

   The ASM preconditioner is fully parallel, but currently the routine
   PCASMCreateSubdomains2D(), which is used in this example to demonstrate
   user-defined subdomains (activated via -user_set_subdomains), is
   uniprocessor only.

   This matrix in this linear system arises from the discretized Laplacian,
   and thus is not very interesting in terms of experimenting with variants
   of the ASM preconditioner.
*/

/*T
   Concepts: KSP^Additive Schwarz Method (ASM) with user-defined subdomains
   Processors: n
T*/

/*
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/

// #include <petscksp.h>
//
// #undef __FUNCT__
// #define __FUNCT__ "main"
// int main(int argc,char **args)
// {
//   Vec            x,b,u;                 /* approx solution, RHS, exact solution */
//   Mat            A;                       /* linear system matrix */
//   KSP            ksp;                    /* linear solver context */
//   PC             pc;                      /* PC context */
//   IS             *is,*is_local;           /* array of index sets that define the subdomains */
//   PetscInt       overlap = 1;             /* width of subdomain overlap */
//   PetscInt       Nsub;                    /* number of subdomains */
//   PetscInt       m = 15,n = 17;          /* mesh dimensions in x- and y- directions */
//   PetscInt       M = 2,N = 1;            /* number of subdomains in x- and y- directions */
//   PetscInt       i,j,Ii,J,Istart,Iend;
//   PetscErrorCode ierr;
//   PetscMPIInt    size;
//   PetscBool      flg;
//   PetscBool      user_subdomains = PETSC_FALSE;
//   PetscScalar    v, one = 1.0;
//   PetscReal      e;
//
//   PetscInitialize(&argc,&args,(char*)0,help);
//   ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
//   ierr = PetscOptionsGetInt(NULL,"-m",&m,NULL);CHKERRQ(ierr);
//   ierr = PetscOptionsGetInt(NULL,"-n",&n,NULL);CHKERRQ(ierr);
//   ierr = PetscOptionsGetInt(NULL,"-Mdomains",&M,NULL);CHKERRQ(ierr);
//   ierr = PetscOptionsGetInt(NULL,"-Ndomains",&N,NULL);CHKERRQ(ierr);
//   ierr = PetscOptionsGetInt(NULL,"-overlap",&overlap,NULL);CHKERRQ(ierr);
//   ierr = PetscOptionsGetBool(NULL,"-user_set_subdomains",&user_subdomains,NULL);CHKERRQ(ierr);
//
//   /* -------------------------------------------------------------------
//          Compute the matrix and right-hand-side vector that define
//          the linear system, Ax = b.
//      ------------------------------------------------------------------- */
//
//   /*
//      Assemble the matrix for the five point stencil, YET AGAIN
//   */
//   ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
//   ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m*n,m*n);CHKERRQ(ierr);
//   ierr = MatSetFromOptions(A);CHKERRQ(ierr);
//   ierr = MatSetUp(A);CHKERRQ(ierr);
//   ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
//   for (Ii=Istart; Ii<Iend; Ii++) {
//     v = -1.0; i = Ii/n; j = Ii - i*n;
//     if (i>0)   {J = Ii - n; ierr = MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);CHKERRQ(ierr);}
//     if (i<m-1) {J = Ii + n; ierr = MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);CHKERRQ(ierr);}
//     if (j>0)   {J = Ii - 1; ierr = MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);CHKERRQ(ierr);}
//     if (j<n-1) {J = Ii + 1; ierr = MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);CHKERRQ(ierr);}
//     v = 4.0; ierr = MatSetValues(A,1,&Ii,1,&Ii,&v,INSERT_VALUES);CHKERRQ(ierr);
//   }
//   ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//   ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//
//   /*
//      Create and set vectors
//   */
//   ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
//   ierr = VecSetSizes(b,PETSC_DECIDE,m*n);CHKERRQ(ierr);
//   ierr = VecSetFromOptions(b);CHKERRQ(ierr);
//   ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
//   ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
//   ierr = VecSet(u,one);CHKERRQ(ierr);
//   ierr = MatMult(A,u,b);CHKERRQ(ierr);
//
//   /*
//      Create linear solver context
//   */
//   ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
//
//   /*
//      Set operators. Here the matrix that defines the linear system
//      also serves as the preconditioning matrix.
//   */
//   ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
//
//   /*
//      Set the default preconditioner for this program to be ASM
//   */
//   ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
//   ierr = PCSetType(pc,PCASM);CHKERRQ(ierr);
//
//   /* -------------------------------------------------------------------
//                   Define the problem decomposition
//      ------------------------------------------------------------------- */
//
//   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        Basic method, should be sufficient for the needs of many users.
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//      Set the overlap, using the default PETSc decomposition via
//          PCASMSetOverlap(pc,overlap);
//      Could instead use the option -pc_asm_overlap <ovl>
//
//      Set the total number of blocks via -pc_asm_blocks <blks>
//      Note:  The ASM default is to use 1 block per processor.  To
//      experiment on a single processor with various overlaps, you
//      must specify use of multiple blocks!
//   */
//
//   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        More advanced method, setting user-defined subdomains
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//      Firstly, create index sets that define the subdomains.  The utility
//      routine PCASMCreateSubdomains2D() is a simple example (that currently
//      supports 1 processor only!).  More generally, the user should write
//      a custom routine for a particular problem geometry.
//
//      Then call either PCASMSetLocalSubdomains() or PCASMSetTotalSubdomains()
//      to set the subdomains for the ASM preconditioner.
//   */

//   if (!user_subdomains) { /* basic version */
//     ierr = PCASMSetOverlap(pc,overlap);CHKERRQ(ierr);
//   }
//   else { /* advanced version */
//     if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"PCASMCreateSubdomains2D() is currently a uniprocessor routine only!");
//     ierr = PCASMCreateSubdomains2D(m,n,M,N,1,overlap,&Nsub,&is,&is_local);CHKERRQ(ierr);
//     ierr = PCASMSetLocalSubdomains(pc,Nsub,is,is_local);CHKERRQ(ierr);
//     flg  = PETSC_FALSE;
//     ierr = PetscOptionsGetBool(NULL,"-subdomain_view",&flg,NULL);CHKERRQ(ierr);
//     if (flg) {
//       ierr = PetscPrintf(PETSC_COMM_SELF,"Nmesh points: %D x %D; subdomain partition: %D x %D; overlap: %D; Nsub: %D\n",m,n,M,N,overlap,Nsub);CHKERRQ(ierr);
//       ierr = PetscPrintf(PETSC_COMM_SELF,"IS:\n");CHKERRQ(ierr);
//       for (i=0; i<Nsub; i++) {
//         ierr = PetscPrintf(PETSC_COMM_SELF,"  IS[%D]\n",i);CHKERRQ(ierr);
//         ierr = ISView(is[i],PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
//       }
//       ierr = PetscPrintf(PETSC_COMM_SELF,"IS_local:\n");CHKERRQ(ierr);
//       for (i=0; i<Nsub; i++) {
//         ierr = PetscPrintf(PETSC_COMM_SELF,"  IS_local[%D]\n",i);CHKERRQ(ierr);
//         ierr = ISView(is_local[i],PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
//       }
//     }
//   }
//
//   /* -------------------------------------------------------------------
//                 Set the linear solvers for the subblocks
//      ------------------------------------------------------------------- */
//
//   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        Basic method, should be sufficient for the needs of most users.
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//      By default, the ASM preconditioner uses the same solver on each
//      block of the problem.  To set the same solver options on all blocks,
//      use the prefix -sub before the usual PC and KSP options, e.g.,
//           -sub_pc_type <pc> -sub_ksp_type <ksp> -sub_ksp_rtol 1.e-4
//
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         Advanced method, setting different solvers for various blocks.
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//      Note that each block's KSP context is completely independent of
//      the others, and the full range of uniprocessor KSP options is
//      available for each block.
//
//      - Use PCASMGetSubKSP() to extract the array of KSP contexts for
//        the local blocks.
//      - See ex7.c for a simple example of setting different linear solvers
//        for the individual blocks for the block Jacobi method (which is
//        equivalent to the ASM method with zero overlap).
//   */
//
//   flg  = PETSC_FALSE;
//   ierr = PetscOptionsGetBool(NULL,"-user_set_subdomain_solvers",&flg,NULL);CHKERRQ(ierr);
//   if (flg) {
//     KSP       *subksp;        /* array of KSP contexts for local subblocks */
//     PetscInt  nlocal,first;   /* number of local subblocks, first local subblock */
//     PC        subpc;          /* PC context for subblock */
//     PetscBool isasm;
//
//     ierr = PetscPrintf(PETSC_COMM_WORLD,"User explicitly sets subdomain solvers.\n");CHKERRQ(ierr);
//
//     /*
//        Set runtime options
//     */
//     ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
//
//     /*
//        Flag an error if PCTYPE is changed from the runtime options
//      */
//     ierr = PetscObjectTypeCompare((PetscObject)pc,PCASM,&isasm);CHKERRQ(ierr);
//     if (!isasm) SETERRQ(PETSC_COMM_WORLD,1,"Cannot Change the PCTYPE when manually changing the subdomain solver settings");
//
//     /*
//        Call KSPSetUp() to set the block Jacobi data structures (including
//        creation of an internal KSP context for each block).
//
//        Note: KSPSetUp() MUST be called before PCASMGetSubKSP().
//     */
//     ierr = KSPSetUp(ksp);CHKERRQ(ierr);
//
//     /*
//        Extract the array of KSP contexts for the local blocks
//     */
//     ierr = PCASMGetSubKSP(pc,&nlocal,&first,&subksp);CHKERRQ(ierr);
//
//     /*
//        Loop over the local blocks, setting various KSP options
//        for each block.
//     */
//     for (i=0; i<nlocal; i++) {
//       ierr = KSPGetPC(subksp[i],&subpc);CHKERRQ(ierr);
//       ierr = PCSetType(subpc,PCILU);CHKERRQ(ierr);
//       ierr = KSPSetType(subksp[i],KSPGMRES);CHKERRQ(ierr);
//       ierr = KSPSetTolerances(subksp[i],1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
//     }
//   } else {
//     /*
//        Set runtime options
//     */
//     ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
//   }
//
//   /* -------------------------------------------------------------------
//                       Solve the linear system
//      ------------------------------------------------------------------- */
//
//   ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
//
//   /* -------------------------------------------------------------------
//                       Compare result to the exact solution
//      ------------------------------------------------------------------- */
//   ierr = VecAXPY(x,-1.0,u);CHKERRQ(ierr);
//   ierr = VecNorm(x,NORM_INFINITY, &e);CHKERRQ(ierr);
//
//   flg  = PETSC_FALSE;
//   ierr = PetscOptionsGetBool(NULL,"-print_error",&flg,NULL);CHKERRQ(ierr);
//   if (flg) {
//     ierr = PetscPrintf(PETSC_COMM_WORLD, "Infinity norm of the error: %G\n", e);CHKERRQ(ierr);
//   }
//
//   /*
//      Free work space.  All PETSc objects should be destroyed when they
//      are no longer needed.
//   */
//
//   if (user_subdomains) {
//     for (i=0; i<Nsub; i++) {
//       ierr = ISDestroy(&is[i]);CHKERRQ(ierr);
//       ierr = ISDestroy(&is_local[i]);CHKERRQ(ierr);
//     }
//     ierr = PetscFree(is);CHKERRQ(ierr);
//     ierr = PetscFree(is_local);CHKERRQ(ierr);
//   }
//   ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
//   ierr = VecDestroy(&u);CHKERRQ(ierr);
//   ierr = VecDestroy(&x);CHKERRQ(ierr);
//   ierr = VecDestroy(&b);CHKERRQ(ierr);
//   ierr = MatDestroy(&A);CHKERRQ(ierr);
//   ierr = PetscFinalize();
//   return 0;
// }









