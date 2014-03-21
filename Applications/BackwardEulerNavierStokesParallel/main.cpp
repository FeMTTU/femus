#include "ElemType.hpp"
#include "NonLinearTimeDependentMultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "LinearEquationSolver.hpp"
#include "Fluid.hpp"
#include "Parameter.hpp"
#include "FemTTUInit.hpp"
#include "SparseMatrix.hpp"
using std::cout;
using std::endl;

void AssembleMatrixResNS(MultiLevelProblem &nl_td_ml_prob, unsigned level, const unsigned &gridn, const unsigned &ipde, const bool &assembe_matrix);

bool SetBoundaryCondition(const double &x, const double &y, const double &z,const char name[], 
			  double &value, const int FaceName, const double time);

bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber,const int &level);

int main(int argc,char **args) {
  
  /// Init Petsc-MPI communicator
  FemTTUInit mpinit(argc,args,MPI_COMM_WORLD);
  
  /// INIT MESH =================================  
  
  unsigned short nm,nr;
  nm=2;
  std::cout<<"MULTIGRID levels: "<< nm << endl;

  nr=2;
  std::cout<<"MAX_REFINEMENT levels: " << nr << endl<< endl;
  
  int tmp=nm;  nm+=nr;  nr=tmp;
  
  char *infile = new char [50];
 
  sprintf(infile,"./input/nsbenchreg.neu");
  
  //Adimensional quantity (Lref,Uref)
  double Lref = 0.1;
  double Uref = 0.2;
  Parameter parameter(Lref,Uref);
  
  // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  Fluid fluid(parameter,0.001,1.,"Newtonian");
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;
  
  //Steadystate NonLinearMultiLevelProblem  
  NonLinearTimeDependentMultiLevelProblem nl_td_ml_prob(nm,nr,infile,"seventh",Lref,SetRefinementFlag);
  
  /// END MESH =================================  
  
  ///Start System Variables;===========================
  //Focus here is on VARIABLES first, rather than on Equations
 
  // add fluid material
  nl_td_ml_prob.parameters.set<Fluid>("Fluid") = fluid;
  
  // generate solution vector
  nl_td_ml_prob.AddSolution("U","biquadratic");
  nl_td_ml_prob.AddSolution("V","biquadratic");
  // the pressure variable should be the last for the Schur decomposition
  nl_td_ml_prob.AddSolution("P","disc_linear",1);
  nl_td_ml_prob.AssociatePropertyToSolution("P","Pressure");
  
  //Initialize (update Init(...) function)
  nl_td_ml_prob.Initialize("All");
  
  
  //Set Boundary (update Dirichlet(...) function)
  nl_td_ml_prob.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  nl_td_ml_prob.GenerateBdc("U","Time_dependent");
//   nl_td_ml_prob.GenerateBdc("U");
  nl_td_ml_prob.GenerateBdc("V");
  nl_td_ml_prob.GenerateBdc("P");
  
  //Set Time step information
  nl_td_ml_prob.SetTimeStep(0.1);
  nl_td_ml_prob.SetPrintTimeStep(1);
  nl_td_ml_prob.SetSaveTimeStep(33300);
  nl_td_ml_prob.SetNumTimeSteps(100);  
  
  ///End System Variables; ==============================

  // START EQUATIONS =================================  
  
  /// Start Navier-Stokes Muligrid Block
    
  //start Multigrid for UVWP
  nl_td_ml_prob.ClearSolPdeIndex();
  nl_td_ml_prob.AddPde("NS");
  nl_td_ml_prob.AddSolutionToSolPdeIndex("NS","U"); 
  nl_td_ml_prob.AddSolutionToSolPdeIndex("NS","V");
  nl_td_ml_prob.AddSolutionToSolPdeIndex("NS","P");
  
  // create Multigrid (PRLO, REST, MAT, VECs) based on SolPdeIndex
  nl_td_ml_prob.CreatePdeStructure();
  
  // create index of solutions to be to used in the Vanka Smoother  
  nl_td_ml_prob.ClearVankaIndex();
  nl_td_ml_prob.AddToVankaIndex("NS","U"); 
  nl_td_ml_prob.AddToVankaIndex("NS","V"); 
  nl_td_ml_prob.AddToVankaIndex("NS","P"); 
    
  //Equation
  nl_td_ml_prob.AttachAssembleFunction(AssembleMatrixResNS);
  nl_td_ml_prob.SetNonLinearAlgorithm(true,"Newton",1.e-06);  //Navier-Stokes (Quasi-Newton - Newton)
  nl_td_ml_prob.SetMatrixProperties("NS","Symmetric");
  nl_td_ml_prob.AddStabilization("NS",true);
  
  nl_td_ml_prob.SetDirichletBCsHandling("NS","Penalty");
  
  //Solver Configuration 
  //Solver I (Gmres)
  nl_td_ml_prob.SetSmoother("Gmres");
  nl_td_ml_prob.SetTolerances("NS",1.e-12,1.e-20,1.e+50,15);
  
 // Solver II (Vanka - MPSC)
//   nl_td_ml_prob.SetSmoother("Vanka");
//   nl_td_ml_prob.SetVankaSchurOptions(true,true,1);
//   nl_td_ml_prob.SetSolverFineGrids("GMRES");
//   nl_td_ml_prob.SetPreconditionerFineGrids("ILU");
//   nl_td_ml_prob.SetTolerances(1.e-12,1.e-20,1.e+50,15);
//   nl_td_ml_prob.SetSchurTolerances(1.e-12,1.e-20,1.e+50,5);
//   nl_td_ml_prob.SetDimVankaBlock(3);                             //2^lev 1D 4^lev 2D 8^lev 3D
  
    for (unsigned time_step = nl_td_ml_prob.GetInitTimeStep(); time_step < nl_td_ml_prob.GetInitTimeStep() + nl_td_ml_prob.GetNumTimeSteps(); 
       time_step++) {
   
    //Solve with V-cycle or F-cycle
    if(time_step==0) nl_td_ml_prob.Solve("NS",15,1,1,"F-Cycle");
    else nl_td_ml_prob.Solve("NS",15,1,1,"V-Cycle");
  
//     //The update of the acceleration must be done before the update of the other variables
//     //update time step
//     nl_td_ml_prob._NewmarkAccUpdate();

    //update Solution
    nl_td_ml_prob._UpdateSolution();

    //print solution for restart
    if ( !(time_step%nl_td_ml_prob.GetSaveTimeStep()) ) {
      nl_td_ml_prob.SaveData();
    }

    // print solution
    if ( !(time_step%nl_td_ml_prob.GetPrintTimeStep()) ) {
       
      ///print solution 
      std::vector<std::string> print_vars;
      print_vars.resize(3);
      print_vars[0] = "U";
      print_vars[1] = "V";
      print_vars[2] = "P";
      
      nl_td_ml_prob.printsol_vtu_inline("biquadratic",print_vars);
//       nl_td_ml_prob.printsol_xdmf_hdf5("biquadratic",print_vars);
       
      bool DebugOn=1;
      //nl_ml_prob.printsol_gmv_binary("quadratic",nm-1,DebugOn);
      nl_td_ml_prob.printsol_gmv_binary("quadratic",nm,DebugOn);
  
      
      
    }
  
  } //end loop timestep
  
 
  // Delete Multigrid (PRLO, REST, MAT, VECs) based on SolPdeIndex
  nl_td_ml_prob.DeletePdeStructure();
  
  /// End Navier-Stokes Muligrid Block
   
  /// Destroy the last PETSC objects
  nl_td_ml_prob.FreeMultigrid(); 
   
  delete [] infile;
  return(0);
}

//-----------------------------------------------------------------------------------------------------------------

bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber, const int &level) {
   bool refine=0;
   // refinemenet based on Elemen Group Number
   if(ElemGroupNumber==5) refine=1;
   if(ElemGroupNumber==6 && level<2) refine=1;
    if(ElemGroupNumber==7) refine=0;
//    if(x>0 && x<5) refine=1;
   return refine;
}

//-------------------------------------------------------------------------------------------------------------------

bool SetBoundaryCondition(const double &x, const double &y, const double &z,const char name[], 
			  double &value, const int FaceName, const double time){
  bool test=1; //Dirichlet
  value=0.;
// //   cout << "Time bdc : " <<  time << endl;
  if(!strcmp(name,"U")) {
    if (1==FaceName) { //inflow
      test=1;
      double um = 1.; // U/Uref
      if(time < 2.0) {
        value=1.5*um*(4.0/(0.1681*100))*y*(0.41/0.1-y)*0.5*(1. - cos(0.5*3.141592653589793*time) );
      } else {
        value=1.5*um*(4.0/(0.1681*100))*y*(0.41/0.1-y);
      }
//     value = 1.e-04; //non-adimensional formula
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
void AssembleMatrixResNS(MultiLevelProblem &nl_ml_prob, unsigned level, const unsigned &gridn, const unsigned &ipde, const bool &assembe_matrix){
     
  NonLinearTimeDependentMultiLevelProblem& nl_td_ml_prob = static_cast<NonLinearTimeDependentMultiLevelProblem&>(nl_ml_prob);
  
  const char* pdename=nl_td_ml_prob.GetThisPdeName(ipde);
  
  //pointers 
  Solution*	 mysolution  	= nl_td_ml_prob._solution[level];
  LinearEquationSolver*  mylsyspde 	= nl_td_ml_prob._LinSolver[ipde][level];
  mesh*		 mymsh    	= nl_td_ml_prob._msh[level];
  elem*		 myel		= mymsh->el;
  SparseMatrix*	 myKK		= mylsyspde->_KK;
  NumericVector* myRES 		= mylsyspde->_RES;
    
  //data
  double dt = nl_td_ml_prob.GetTimeStep();
  double theta = 0.5;
  const unsigned dim = mymsh->GetDimension();
  unsigned nel= mymsh->GetElementNumber();
  unsigned igrid= mymsh->GetGridNumber();
  unsigned iproc = mymsh->GetProcID();
  unsigned nprocs = mymsh->GetNumProcs();
  double ILambda= mylsyspde->GetCompressibility();
  double IRe = nl_td_ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
  bool penalty = mylsyspde->GetStabilization();
  const bool symm_mat = mylsyspde->GetMatrixProperties();
  const bool NavierStokes = nl_td_ml_prob.GetNonLinearCase();
  unsigned nwtn_alg = nl_td_ml_prob.GetNonLinearAlgorithm();
  bool newton = (nwtn_alg==0) ? 0:1;
  
  // solution and coordinate variables
  const char Solname[4][2] = {"U","V","W","P"};
  vector < unsigned > SolPdeIndex(dim+1);
  vector < unsigned > SolIndex(dim+1);  

  const char coordinate_name[3][2] = {"X","Y","Z"};
  vector < unsigned > coordinate_Index(dim);
  vector< vector < double> > coordinates(dim);
  
  for(unsigned ivar=0; ivar<dim; ivar++) {
    SolPdeIndex[ivar]=nl_td_ml_prob.GetSolPdeIndex(pdename,&Solname[ivar][0]);
    SolIndex[ivar]=nl_td_ml_prob.GetIndex(&Solname[ivar][0]);
    coordinate_Index[ivar]=nl_td_ml_prob.GetIndex(&coordinate_name[ivar][0]);
  }
  SolPdeIndex[dim]=nl_td_ml_prob.GetSolPdeIndex(pdename, &Solname[3][0]);
  SolIndex[dim]=nl_td_ml_prob.GetIndex(&Solname[3][0]);       
  //solution order
  unsigned order_ind2 = nl_td_ml_prob.SolType[SolIndex[0]];
  unsigned end_ind2   = mymsh->GetEndIndex(order_ind2);
  unsigned order_ind1 = nl_td_ml_prob.SolType[SolIndex[dim]];
  unsigned end_ind1   = mymsh->GetEndIndex(order_ind1);
  
  // declare 
  vector < int > metis_node2; 
  vector < int > node1;
  vector< vector< int > > KK_dof(dim+1); 
  vector <double> phi2;
  vector <double> gradphi2;
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
  for(int i=0;i<dim;i++) {
    KK_dof[i].reserve(max_size);
  }
   
  for(int i=0;i<dim+1;i++) F[i].reserve(max_size);
    
  if(assembe_matrix){
    for(int i=0;i<dim+1;i++){
      B[i].resize(dim+1);
      for(int j=0;j<dim+1;j++){
	B[i][j].reserve(max_size*max_size);
      }
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
  if(assembe_matrix) myKK->zero();
  
  // *** element loop ***
 
  for (int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

    unsigned kel = mymsh->IS_Mts2Gmt_elem[iel];
    short unsigned kelt=myel->GetElementType(kel);
    unsigned nve2=myel->GetElementDofNumber(kel,end_ind2);
    unsigned nve1=myel->GetElementDofNumber(kel,end_ind1);
    
    //set to zero all the entries of the FE matrices
    metis_node2.resize(nve2);
    node1.resize(nve1);
    phi2.resize(nve2);
    gradphi2.resize(nve2*dim);
    for(int ivar=0; ivar<dim; ivar++) {
      coordinates[ivar].resize(nve2);
      KK_dof[ivar].resize(nve2);
      
      F[SolPdeIndex[ivar]].resize(nve2);
      memset(&F[SolPdeIndex[ivar]][0],0,nve2*sizeof(double));
      
      if(assembe_matrix){
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
      
      
    if(assembe_matrix*nwtn_alg==2){
      for(int ivar=0; ivar<dim; ivar++) {
	for(int ivar2=1; ivar2<dim; ivar2++) {
	  B[SolPdeIndex[ivar]][SolPdeIndex[(ivar+ivar2)%dim]].resize(nve2*nve2);
	  memset(&B[SolPdeIndex[ivar]][SolPdeIndex[(ivar+ivar2)%dim]][0],0,nve2*nve2*sizeof(double));
	}
      }
    }
  
    if(assembe_matrix*penalty){
      B[SolPdeIndex[dim]][SolPdeIndex[dim]].resize(nve1*nve1,0.);
      memset(&B[SolPdeIndex[dim]][SolPdeIndex[dim]][0],0,nve1*nve1*sizeof(double));
    }
    
    for( unsigned i=0;i<nve2;i++){
      unsigned inode=myel->GetElementVertexIndex(kel,i)-1u;
      unsigned inode_metis=mymsh->GetMetisDof(inode,2);
      metis_node2[i]=inode_metis;
      for(unsigned ivar=0; ivar<dim; ivar++) {
	coordinates[ivar][i]=(*mysolution->_Sol[coordinate_Index[ivar]])(inode_metis);
	KK_dof[ivar][i]=mylsyspde->GetKKDof(SolIndex[ivar],SolPdeIndex[ivar],inode);
      }
    }
    
    for(unsigned i=0;i<nve1;i++) {
      unsigned inode=(order_ind1<dim)?(myel->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
      node1[i]=inode;
      KK_dof[dim][i]=mylsyspde->GetKKDof(SolIndex[dim],SolPdeIndex[dim],inode);
    }
   
    if(igrid==gridn || !myel->GetRefinedElementIndex(kel)) {
      // *** Gauss poit loop ***
      for(unsigned ig=0;ig < nl_td_ml_prob.type_elem[kelt][order_ind2]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	(nl_td_ml_prob.type_elem[kelt][order_ind2]->*(nl_td_ml_prob.type_elem[kelt][order_ind2])->Jacobian_ptr)(coordinates,ig,Weight2,phi2,gradphi2);
	phi1=nl_td_ml_prob.type_elem[kelt][order_ind1]->GetPhi(ig);

	//velocity variable
	for(unsigned ivar=0; ivar<dim; ivar++) {
	  Sol[ivar]=0.;
	  SolOld[ivar]=0.;
	  //AccSol[ivar]=0.;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSol[ivar][ivar2]=0;
	    gradSolOld[ivar][ivar2]=0;
	  }
	  unsigned SolIndex=nl_td_ml_prob.GetIndex(&Solname[ivar][0]);
	  unsigned SolType =nl_td_ml_prob.GetSolType(&Solname[ivar][0]);
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
	unsigned SolIndex=nl_td_ml_prob.GetIndex(&Solname[3][0]);
	unsigned SolType=nl_td_ml_prob.GetSolType(&Solname[3][0]);
	for(unsigned i=0; i<nve1; i++){
	  unsigned sol_dof = mymsh->GetMetisDof(node1[i],SolType);
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
	  
	  if(assembe_matrix){
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
	  } // endif assembe_matrix
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
	  
	  if(assembe_matrix){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve2; j++) {
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndex[dim]][SolPdeIndex[ivar]][i*nve2+j]-= dt*phi1[i]*gradphi2[j*dim+ivar]*Weight2;
	      }
	    }  //end phij loop
	  } // endif assembe_matrix
	}  //end phi1_i loop
	
	if(assembe_matrix * penalty){  //block nve1 nve1
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
      // Boundary Integral
      //number of faces for each type of element
      if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {
     
	unsigned nfaces = myel->GetElementFaceNumber(kel);

	// loop on faces
	for(unsigned jface=0;jface<nfaces;jface++){ 
	  
	  // look for boundary faces
	  if(myel->GetFaceElementIndex(kel,jface)<0){
	    for(unsigned ivar=0; ivar<dim; ivar++) {
	      nl_td_ml_prob.ComputeBdIntegral(pdename, &Solname[ivar][0], kel, jface, level, ivar);
	    }
	  }
	}	
      }
      //--------------------------------------------------------------------------------------------------------
    } // endif single element not refined or fine grid loop
    //--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the Global Matrix/Vector
    for(unsigned ivar=0; ivar<dim; ivar++) {
      myRES->add_vector_blocked(F[SolPdeIndex[ivar]],KK_dof[ivar]);
      if(assembe_matrix){
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
    if(assembe_matrix*penalty) myKK->add_matrix_blocked(B[SolPdeIndex[dim]][SolPdeIndex[dim]],KK_dof[dim],KK_dof[dim]);
    myRES->add_vector_blocked(F[SolPdeIndex[dim]],KK_dof[dim]);
    //--------------------------------------------------------------------------------------------------------  
  } //end list of elements loop for each subdomain
  
  
  if(assembe_matrix) myKK->close();
  myRES->close();
  // ***************** END ASSEMBLY *******************
}
