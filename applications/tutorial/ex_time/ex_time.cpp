#include "MultiLevelProblem.hpp"
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


#define NSUB_X  16
#define NSUB_Y  16

using std::cout;
using std::endl;
using namespace femus;

void AssembleMatrixResNS(MultiLevelProblem &ml_prob);

bool SetBoundaryCondition(const std::vector < double >& x, const char name[],
                          double &value, const int FaceName, const double time);

double SetInitialCondition (const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {
         
           double value = 0.;

             if(!strcmp(name,"U")) {
                 value = 0.;
             }

           
      return value;   
}


int main(int argc,char **args) {

    // ======= Init ========================
  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
    // ======= Files ========================
  Files files; 
        files.CheckIODirectories();
        files.RedirectCout();

  // ======= Quad Rule ===================
  std::string fe_quad_rule("seventh");

  // ======= Mesh ========================
  /// INIT MESH =================================
  //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;

  MultiLevelMesh ml_msh;
  ml_msh.GenerateCoarseBoxMesh(NSUB_X,NSUB_Y,0,0.,1.,0.,1.,0.,0.,QUAD9,fe_quad_rule.c_str());
  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  ml_msh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  ml_msh.PrintInfo();


    // ======= Solution ========================
  MultiLevelSolution mlSol(&ml_msh);  // define the multilevel solution and attach the mlMsh object to it
  mlSol.AddSolution("U",   LAGRANGE, FIRST,2);

    // ======= Problem ========================
  MultiLevelProblem ml_prob(&mlSol);  // define the multilevel problem attach the mlSol object to it

  ml_prob.SetFilesHandler(&files);
  
    // ======= Initial values ========================
  mlSol.Initialize("All");    // initialize all variables to zero
  mlSol.Initialize("U",   SetInitialCondition, &ml_prob);

    // ======= Boundary Conditions ========================
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);  // attach the boundary condition function and generate boundary data

//   mlSol.GenerateBdc("All");  //this would do it also for the non-equation-related variables
  mlSol.GenerateBdc("U"); //"Time_dependent");

    // ======= Parameters ========================
  Parameter parameter(Lref,Uref);

  // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  Fluid fluid(parameter,0.001,1.,"Newtonian");
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;

  // add fluid material
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;



  //create systems
  // add the system Navier-Stokes to the MultiLevel problem
  TransientNonlinearImplicitSystem & system = ml_prob.add_system<TransientNonlinearImplicitSystem> ("Timedep");
  system.AddSolutionToSystemPDE("U");

  // init all the systems
  system.init();

  // System Navier-Stokes
  system.SetAssembleFunction(AssembleMatrixResNS);
  system.SetMaxNumberOfLinearIterations(1);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-8);
  system.SetMgType(V_CYCLE);
  system.SetMaxNumberOfNonLinearIterations(15);

  //**************
   mlSol.SetWriter(VTK);   //need to move this here for the DebugNonlinear function
//    mlSol.GetWriter()->SetDebugOutput(true);
  //**************

  
  
  // time loop parameter
  system.SetIntervalTime(0.1);
  const unsigned int n_timesteps = 20;
  const unsigned int write_interval = 1;

  for (unsigned time_step = 0; time_step < n_timesteps; time_step++) {

    // Solving Navier-Stokes system
    std::cout << std::endl;
    std::cout << " *********** Timedep ************  " << std::endl;
    ml_prob.get_system("Timedep").MLsolve();

    //update Solution
    ml_prob.get_system<TransientNonlinearImplicitSystem>("Timedep").CopySolutionToOldSolution();

    // print solution
    if ( !(time_step%write_interval) ) {

    // ======= Final Print ========================
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
  mlSol.GetWriter()->Write(files.GetOutputPath()/*DEFAULT_OUTPUTDIR*/, "biquadratic", variablesToBePrinted,time_step);    // print solutions

    }

  } //end loop timestep


  // Destroy all the new systems
  ml_prob.clear();


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

  TransientNonlinearImplicitSystem* mlPdeSys = & ml_prob.get_system<TransientNonlinearImplicitSystem>("Timedep");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  MultiLevelSolution *ml_sol			      = ml_prob._ml_sol;
  Solution*	 mysolution  	                      = ml_sol->GetSolutionLevel(level);
  LinearEquationSolver*  mylsyspde	              = mlPdeSys->_LinSolver[level];
  const char* pdename                                 = mlPdeSys->name().c_str();

  Mesh*		 mymsh    	   = ml_prob._ml_msh->GetLevel(level);
  elem*		 myel		   = mymsh->el;
  SparseMatrix*	 myKK	 	   = mylsyspde->_KK;
  NumericVector* myRES 		   = mylsyspde->_RES;

  //data
  double dt = mlPdeSys->GetIntervalTime();
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
  const unsigned int n_unknowns = mlPdeSys->GetSolPdeIndex().size();

  const char Solname[1][2] = {"U"};
  vector < unsigned > SolPdeIndex(n_unknowns);
  vector < unsigned > SolIndex(n_unknowns);

  vector< vector < double> > coordinates(dim);

  for(unsigned ivar=0; ivar< SolPdeIndex.size(); ivar++) {
    SolPdeIndex[ivar]=mlPdeSys->GetSolPdeIndex(&Solname[ivar][0]);
    SolIndex[ivar]=ml_sol->GetIndex(&Solname[ivar][0]);
  }

  //solution order
  unsigned order_ind2 = ml_sol->GetSolutionType(SolIndex[0]);

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

    //set to zero all the entries of the FE matrices
    metis_node2.resize(nve2);
    phi2.resize(nve2);
    gradphi2.resize(nve2*dim);
    nablaphi2.resize(nve2*(3*(dim-1)));

    for(int ivar=0; ivar< dim; ivar++) {
       coordinates[ivar].resize(nve2);
     }
        

    
    for(int ivar=0; ivar<n_unknowns; ivar++) {
      KK_dof[ivar].resize(nve2);

      F[SolPdeIndex[ivar]].resize(nve2);
      memset(&F[SolPdeIndex[ivar]][0],0,nve2*sizeof(double));


      B[SolPdeIndex[ivar]][SolPdeIndex[ivar]].resize(nve2*nve2);
      memset(&B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][0],0,nve2*nve2*sizeof(double));

    }


    if(nwtn_alg==2){
      for(int ivar=0; ivar<dim; ivar++) {
	for(int ivar2=1; ivar2<dim; ivar2++) {
	  B[SolPdeIndex[ivar]][SolPdeIndex[(ivar+ivar2)%dim]].resize(nve2*nve2);
	  memset(&B[SolPdeIndex[ivar]][SolPdeIndex[(ivar+ivar2)%dim]][0],0,nve2*nve2*sizeof(double));
	}
      }
    }



    for( unsigned i=0;i<nve2;i++) {
      unsigned inode_metis=mymsh->GetSolutionDof(i, iel, 2);
      metis_node2[i]=inode_metis;
      for(unsigned ivar=0; ivar<dim; ivar++) {
	coordinates[ivar][i]=(*mymsh->_topology->_Sol[ivar])(inode_metis);
      }	
      for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
           KK_dof[ivar][i]=mylsyspde->GetSystemDof(SolIndex[ivar],SolPdeIndex[ivar],i, iel);
        }
    }



      // *** Gauss poit loop ***
      for(unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->Jacobian(coordinates,ig,Weight2,phi2,gradphi2,nablaphi2);

	//velocity variable
	for(unsigned ivar=0; ivar<n_unknowns; ivar++) {
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


	// *** phi_i loop ***
	for(unsigned i=0; i<nve2; i++){

	  //BEGIN RESIDUALS A block ===========================
	  for(unsigned ivar=0; ivar<n_unknowns; ivar++) {
	    double Adv_rhs=0;
	    double Adv_old_rhs=0;
	    double Lap_rhs=0;
	    double Lap_old_rhs=0;
	    for(unsigned ivar2=0; ivar2< dim; ivar2++) {
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

	      for(unsigned ivar=0; ivar< n_unknowns; ivar++) {
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

	  } // endif assemble_matrix
	} //end phii loop


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
    //Sum the local matrices/vectors into the Global Matrix/Vector
    for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
      myRES->add_vector_blocked(F[SolPdeIndex[ivar]],KK_dof[ivar]);
       myKK->add_matrix_blocked(B[SolPdeIndex[ivar]][SolPdeIndex[ivar]],KK_dof[ivar],KK_dof[ivar]);
    }
    //--------------------------------------------------------------------------------------------------------
  } //end list of elements loop for each subdomain


  myKK->close();
  myRES->close();
  // ***************** END ASSEMBLY *******************
}


