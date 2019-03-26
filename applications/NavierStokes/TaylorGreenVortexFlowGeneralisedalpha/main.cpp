#include "ElemTypeEnum.hpp"
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

double InitVariableU(const std::vector < double >& x);
double InitVariableV(const std::vector < double >& x);
double InitVariableP(const std::vector < double >& x);

bool SetBoundaryCondition(const std::vector < double >& x,const char name[],
			  double &value, const int FaceName, const double time);

int main(int argc,char **args) {

  /// Init Petsc-MPI communicator
  FemusInit mpinit(argc,args,MPI_COMM_WORLD);

  Files files;
        files.CheckIODirectories();
        files.RedirectCout();

  /// INIT MESH =================================

  unsigned short nm,nr;
  nm=5;
  std::cout<<"MULTIGRID levels: "<< nm << endl;

  nr=0;
  std::cout<<"MAX_REFINEMENT levels: " << nr << endl<< endl;

  int tmp=nm;  nm+=nr;  nr=tmp;

  char *infile = new char [50];

  sprintf(infile,"./input/square.neu");

  //Adimensional quantity (Lref,Uref)
  double Lref = 1.0;
  double Ladimref = 1./(2.*3.1415926535897932);
  double Uref = 1.0;
  //MultiLevelMesh ml_msh(nm,nr,infile,"seventh",Ladimref,SetRefinementFlag);

  MultiLevelMesh ml_msh;
  //ml_msh.ReadCoarseMesh(infile,"seventh",Lref);
  ml_msh.GenerateCoarseBoxMesh(4,4,0,0.,2.*3.1415926535897932,0.,2.*3.1415926535897932,0.,0.,QUAD9,"seventh");
  ml_msh.RefineMesh(nm,nr,NULL);

  MultiLevelSolution ml_sol(&ml_msh);

  MultiLevelProblem ml_prob(&ml_sol);

  // generate solution vector
  ml_sol.AddSolution("U",LAGRANGE,SECOND,2);
  ml_sol.AddSolution("V",LAGRANGE,SECOND,2);
  ml_sol.AddSolution("AX",LAGRANGE,SECOND,1,0);
  ml_sol.AddSolution("AY",LAGRANGE,SECOND,1,0);
  // the pressure variable should be the last for the Schur decomposition
  ml_sol.AddSolution("P",LAGRANGE,FIRST,1);

  //Initialize (update Init(...) function)
  ml_sol.Initialize("U",InitVariableU);
  ml_sol.Initialize("V",InitVariableV);
  ml_sol.Initialize("AX");
  ml_sol.Initialize("AY");
  ml_sol.Initialize("P",InitVariableP);

  //Set Boundary (update Dirichlet(...) function)
  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  ml_sol.GenerateBdc("U","Time_dependent");
  ml_sol.GenerateBdc("V","Time_dependent");
  ml_sol.GenerateBdc("P","Time_dependent");
  ml_sol.GenerateBdc("AX","Steady");
  ml_sol.GenerateBdc("AY","Steady");


  Parameter parameter(Lref,Uref);

  // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  Fluid fluid(parameter,0.01,1.,"Newtonian");
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
  const unsigned int n_timesteps = 5;
  const unsigned int write_interval = 1;

  VTKWriter vtkio(&ml_sol);

  for (unsigned time_step = 0; time_step < n_timesteps; time_step++) {

    // Solving Navier-Stokes system
    std::cout << std::endl;
    std::cout << " *********** Navier-Stokes ************  " << std::endl;
    ml_prob.get_system("Navier-Stokes").SetOuterSolver(PREONLY);
    ml_prob.get_system("Navier-Stokes").MGsolve();

    //The update of the acceleration must be done before the update of the other variables
    ml_prob.get_system<TransientNonlinearImplicitSystem>("Navier-Stokes").NewmarkAccUpdate();

    //update Solution
    ml_prob.get_system<TransientNonlinearImplicitSystem>("Navier-Stokes").CopySolutionToOldSolution();

    // print solution
    if ( !(time_step%write_interval) ) {

      //print solution
      std::vector<std::string> print_vars;
      print_vars.push_back("U");
      print_vars.push_back("V");
      print_vars.push_back("P");
      print_vars.push_back("AX");
      print_vars.push_back("AY");

//       ml_prob.printsol_vtu_inline("biquadratic",print_vars,time_step);
      vtkio.Write(files.GetOutputPath(),"biquadratic",print_vars,time_step);
    }

  } //end loop timestep

  // Destroy all the new systems
  ml_prob.clear();

  delete[] infile;
  return 0;
}



//--------------------------------------------------------------------------------------------------------------

double InitVariableU(const std::vector < double >& x){
  double value = sin(x[0])*cos(x[1]);
  return value;
}

double InitVariableV(const std::vector < double >& x){
  double value = -cos(x[0])*sin(x[1]);
  return value;
}

double InitVariableP(const std::vector < double >& x){
  double value = 1.*0.25*(cos(2.*x[0])+cos(2.*x[1]));
  return value;
}

//-------------------------------------------------------------------------------------------------------------------

bool SetBoundaryCondition(const std::vector < double >& x,const char name[],
			  double &value, const int FaceName, const double time){
  bool test=1; //Dirichlet
  value=0.;

  if(!strcmp(name,"U")) {
    test=1;
    value=sin(x[0])*cos(x[1])*exp(-2.*0.01*time);
  }
  else if(!strcmp(name,"V")) {
    test=1;
    value=-cos(x[0])*sin(x[1])*exp(-2.*0.01*time);
  }
  else if(!strcmp(name,"P")) {
    test=0;
    value=0.;
    if(x[0] < 1.e-08 && x[1] < 1.e-08) {
      test=1;
      value = 1.*0.25*(cos(2.*x[0])+cos(2.*x[1]))*exp(-4.*0.01*time);
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

  Mesh*		 mymsh    	= ml_prob._ml_msh->GetLevel(level);
  elem*		 myel		= mymsh->el;
  SparseMatrix*	 myKK		= mylsyspde->_KK;
  NumericVector* myRES 		= mylsyspde->_RES;

  //data
  double dt = my_nnlin_impl_sys.GetIntervalTime();
  double theta = 0.5;
  const unsigned dim = mymsh->GetDimension();
  unsigned nel= mymsh->GetNumberOfElements();
  unsigned igrid= mymsh->GetLevel();
  unsigned iproc = mymsh->processor_id();
  double ILambda = 0.;
  double IRe = ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
  bool penalty = true; // mylsyspde->GetStabilization();
  const bool symm_mat = false; //mylsyspde->GetMatrixProperties();
  const bool NavierStokes = true;
  unsigned nwtn_alg = 1;
  bool newton = (nwtn_alg==0) ? 0:1;

  // solution and coordinate variables
  const char Solname[4][2] = {"U","V","W","P"};
  vector < unsigned > SolPdeIndex(dim+1);
  vector < unsigned > SolIndex(dim+1);

  const char AccSolname[3][3] = {"AX","AY","AZ"};


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
  vector < double > AccSol(dim);

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
      unsigned inode=mymsh->GetSolutionDof(i, iel,order_ind1);
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
	  AccSol[ivar]=0.;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){
	    gradSol[ivar][ivar2]=0;
	    gradSolOld[ivar][ivar2]=0;
	  }
	  unsigned SolIndex=ml_sol->GetIndex(&Solname[ivar][0]);
	  unsigned SolType =ml_sol->GetSolutionType(&Solname[ivar][0]);
	  unsigned AccSolIndex=ml_sol->GetIndex(&AccSolname[ivar][0]);
	  for(unsigned i=0; i<nve2; i++) {
	    double soli    = (*mysolution->_Sol[SolIndex])(metis_node2[i]);
	    double sololdi = (*mysolution->_SolOld[SolIndex])(metis_node2[i]);
	    double accsoli = (*mysolution->_Sol[AccSolIndex])(metis_node2[i]);
	    Sol[ivar]+=phi2[i]*soli;
	    SolOld[ivar]+=phi2[i]*sololdi;
	    AccSol[ivar]+=phi2[i]*accsoli;
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

	if( penalty){  //block nve1 nve1
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

