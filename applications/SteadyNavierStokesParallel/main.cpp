#include "ElemType.hpp"
#include "MultiLevelProblem.hpp"
#include "Fluid.hpp"
#include "Parameter.hpp"
#include "FemTTUInit.hpp"
#include "NumericVector.hpp"
#include "LinearEquationSolver.hpp"
#include "SparseMatrix.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "LinearImplicitSystem.hpp"
#include "VTKOutput.hpp"

using std::cout;
using std::endl;
   
void AssembleMatrixResNS(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix);
void AssembleMatrixResT(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix);


double InitVariableU(const double &x, const double &y, const double &z);


bool SetBoundaryCondition(const double &x, const double &y, const double &z,const char name[], 
			  double &value, const int FaceName, const double time);

bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber,const int &level);

int main(int argc,char **args) {

  bool linear=1;
  bool vanka=1;
  if(argc == 2) {
    if( strcmp("vanka",args[1])) vanka=0;
  }
  else {
    cout << "No input arguments!" << endl;
    exit(0);
  }
  
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
 
  sprintf(infile,"./input/nsbenc.neu");
  
  //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;
  
  //Steadystate NonLinearMultiLevelProblem  
  MultiLevelProblem ml_prob(nm,nr,infile,"seventh",Lref,SetRefinementFlag);
  
  // add fluid material
  Parameter parameter(Lref,Uref);
  
  // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  //Fluid fluid(parameter,0.001,1.,"Newtonian",0.001,1.);
  
  Fluid fluid(parameter,0.001,1,"Newtonian",0.001,1.);
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;
  
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;
  
  // generate solution vector
  ml_prob.AddSolution("T","biquadratic");
  ml_prob.AddSolution("U","biquadratic");
  ml_prob.AddSolution("V","biquadratic");
  // the pressure variable should be the last for the Schur decomposition
  ml_prob.AddSolution("P","disc_linear");
  ml_prob.AssociatePropertyToSolution("P","Pressure");

 
  //Initialize (update Init(...) function)
  ml_prob.Initialize("U",InitVariableU);
  ml_prob.Initialize("V");
  ml_prob.Initialize("P");
  ml_prob.Initialize("T");
  
  //Set Boundary (update Dirichlet(...) function)
  ml_prob.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  ml_prob.GenerateBdc("U");
  ml_prob.GenerateBdc("V");
  ml_prob.GenerateBdc("P");
  ml_prob.GenerateBdc("T");
  
  //create systems
  // add the system Navier-Stokes to the MultiLevel problem
  NonLinearImplicitSystem & system1 = ml_prob.add_system<NonLinearImplicitSystem> ("Navier-Stokes");
  system1.AddSolutionToSytemPDE("U");
  system1.AddSolutionToSytemPDE("V");
  system1.AddSolutionToSytemPDE("P");
  
  // add the system Temperature to the MultiLevel problem
  LinearImplicitSystem & system2 = ml_prob.add_system<LinearImplicitSystem> ("Temperature");
  system2.AddSolutionToSytemPDE("T");

  // init all the systems
  ml_prob.init();
   
  // System Navier-Stokes
  system1.AttachAssembleFunction(AssembleMatrixResNS);  
  system1.SetMaxNumberOfNonLinearIterations(3);
  system1.SetMaxNumberOfLinearIterations(1);
  system1.SetAbsoluteConvergenceTolerance(1.e-10);  
  system1.SetMgType(F_CYCLE);
  system1.SetNumberPreSmoothingStep(1);
  system1.SetNumberPostSmoothingStep(1);
 
  if(!vanka){
    system1.SetMgSmoother(GMRES_SMOOTHER);
    system1.SetTolerances(1.e-12,1.e-20,1.e+50,10);
  }
  else{
    system1.SetMgSmoother(VANKA_SMOOTHER);
    system1.AddStabilization(true);
    system1.ClearVankaIndex();
    system1.AddVariableToVankaIndex("U");
    system1.AddVariableToVankaIndex("V");
    system1.AddVariableToVankaIndex("P");
    system1.SetSolverFineGrids(GMRES);
    system1.SetPreconditionerFineGrids(ILU_PRECOND); 
    system1.SetVankaSchurOptions(false,1);
    system1.SetTolerances(1.e-12,1.e-20,1.e+50,10);
    //system1.SetSchurTolerances(1.e-12,1.e-20,1.e+50,1);
    system1.SetDimVankaBlock(3);                
  }
  
  // Solving Navier-Stokes system
  std::cout << std::endl;
  std::cout << " *********** Navier-Stokes ************  " << std::endl;
  ml_prob.get_system("Navier-Stokes").solve();
  
  
  // System Temperature
  system2.AttachAssembleFunction(AssembleMatrixResT);
  system2.SetMaxNumberOfLinearIterations(6);
  system2.SetAbsoluteConvergenceTolerance(1.e-10);  
  system2.SetMgType(F_CYCLE);
   
  // Solving Temperature system
  std::cout << std::endl;
  std::cout << " *********** Temperature ************* " << std::endl;
  ml_prob.get_system("Temperature").solve();
   
  /// Print all solutions
  std::vector<std::string> print_vars;
  print_vars.push_back("U");
  print_vars.push_back("V");
  print_vars.push_back("P");
  print_vars.push_back("T");
  
//   ml_prob.printsol_vtu_inline("biquadratic",print_vars);
//     
  ml_prob.printsol_gmv_binary("biquadratic",0,1);
  
  VTKOutput vtkio(ml_prob);
  vtkio.write_system_solutions("biquadratic",print_vars);
  
  VTKOutput vtkio2(ml_prob);
  vtkio2.write_system_solutions("biquadratic",print_vars);
  
  // Destroy all the new systems
  ml_prob.clear();
  
  delete [] infile;
  return 0;
}

//-----------------------------------------------------------------------------------------------------------------

bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber, const int &level) {
  bool refine=0;
  // refinemenet based on Elemen Group Number
  if(ElemGroupNumber==5 ) {
    refine=1;
  }
  if(ElemGroupNumber==6 && level<2) {
    refine=1;
  }
  if(ElemGroupNumber==7 ) {
    refine=0;
  }

  return refine;
}

//--------------------------------------------------------------------------------------------------------------

double InitVariableU(const double &x, const double &y, const double &z) { 
   double um = 0.2;
   double  value=1.5*um*(4.0/(0.1681))*y*(0.41-y); 
   return value;
}

//-------------------------------------------------------------------------------------------------------------------

bool SetBoundaryCondition(const double &x, const double &y, const double &z,const char name[], 
			  double &value, const int FaceName, const double time){
  bool test=1; //Dirichlet
  value=0.;
  //   cout << "Time bdc : " <<  time << endl;
  if(!strcmp(name,"U")) {
    if(1==FaceName){   //inflow
      test=1;
      double um = 0.2; // U/Uref
      value=1.5*0.2*(4.0/(0.1681))*y*(0.41-y);
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


void AssembleMatrixResNS(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix){
     
  //pointers 
  Solution*	 mysolution  	             = ml_prob._solution[level];
  NonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<NonLinearImplicitSystem>("Navier-Stokes");
  LinearEquationSolver*  mylsyspde	     = my_nnlin_impl_sys._LinSolver[level];   
  const char* pdename                        = my_nnlin_impl_sys.name().c_str();
  
  mesh*		 mymsh    	= ml_prob._msh[level];
  elem*		 myel		= mymsh->el;
  SparseMatrix*	 myKK		= mylsyspde->_KK;
  NumericVector* myRES 		= mylsyspde->_RES;
    
  //data
  const unsigned dim = mymsh->GetDimension();
  unsigned nel= mymsh->GetElementNumber();
  unsigned igrid= mymsh->GetGridNumber();
  unsigned iproc = mymsh->GetProcID();
  unsigned nprocs = mymsh->GetNumProcs();
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
    SolIndex[ivar]=ml_prob.GetIndex(&Solname[ivar][0]);
    //coordinate_Index[ivar]=ivar;//ml_prob.GetIndex(&coordinate_name[ivar][0]);
  }
  SolPdeIndex[dim]=my_nnlin_impl_sys.GetSolPdeIndex(&Solname[3][0]);
  SolIndex[dim]=ml_prob.GetIndex(&Solname[3][0]);       
  //solution order
  unsigned order_ind2 = ml_prob.SolType[SolIndex[0]];
  unsigned end_ind2   = mymsh->GetEndIndex(order_ind2);
  unsigned order_ind1 = ml_prob.SolType[SolIndex[dim]];
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
    
  vector < double > SolVAR(dim+1);
  vector < vector < double > > gradSolVAR(dim);
  for(int i=0;i<dim;i++) {
    gradSolVAR[i].resize(dim);  
  }
  
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
	coordinates[ivar][i]=(*mymsh->_coordinate->_Sol[ivar])(inode_metis);
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
      for(unsigned ig=0;ig < ml_prob.type_elem[kelt][order_ind2]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	(ml_prob.type_elem[kelt][order_ind2]->*(ml_prob.type_elem[kelt][order_ind2])->Jacobian_ptr)(coordinates,ig,Weight2,phi2,gradphi2);
	phi1=ml_prob.type_elem[kelt][order_ind1]->GetPhi(ig);

	//velocity variable
	for(unsigned ivar=0; ivar<dim; ivar++) {
	  SolVAR[ivar]=0;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSolVAR[ivar][ivar2]=0; 
	  }
	  unsigned SolIndex=ml_prob.GetIndex(&Solname[ivar][0]);
	  unsigned SolType=ml_prob.GetSolType(&Solname[ivar][0]);
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
	unsigned SolIndex=ml_prob.GetIndex(&Solname[3][0]);
	unsigned SolType=ml_prob.GetSolType(&Solname[3][0]);
	for(unsigned i=0; i<nve1; i++){
	  unsigned sol_dof = mymsh->GetMetisDof(node1[i],SolType);
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
	  
	  if(assembe_matrix){
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
	  } // endif assembe_matrix
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
	  
	  if(assembe_matrix){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve2; j++) {
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndex[dim]][SolPdeIndex[ivar]][i*nve2+j]-= phi1[i]*gradphi2[j*dim+ivar]*Weight2;
	      }
	    }  //end phij loop
	  } // endif assembe_matrix
	}  //end phi1_i loop
	
	if(assembe_matrix * penalty){  //block nve1 nve1
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
      // Boundary Integral
      //number of faces for each type of element
      if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {
     
	unsigned nfaces = myel->GetElementFaceNumber(kel);

	// loop on faces
	for(unsigned jface=0;jface<nfaces;jface++){ 
	  
	  // look for boundary faces
	  if(myel->GetFaceElementIndex(kel,jface)<0){
	    for(unsigned ivar=0; ivar<dim; ivar++) {
	      ml_prob.ComputeBdIntegral(pdename, &Solname[ivar][0], kel, jface, level, ivar);
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

//------------------------------------------------------------------------------------------------------------
void AssembleMatrixResT(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix){
  
  //pointers and references
  Solution*      mysolution	       = ml_prob._solution[level];
  LinearImplicitSystem& mylin_impl_sys = ml_prob.get_system<LinearImplicitSystem>("Temperature");
  LinearEquationSolver*  mylsyspde     = mylin_impl_sys._LinSolver[level];   
  mesh*          mymsh		       = ml_prob._msh[level];
  elem*          myel		       = mymsh->el;
  SparseMatrix*  myKK		       = mylsyspde->_KK;
  NumericVector* myRES		       = mylsyspde->_RES;
   
  //data
  const unsigned	dim	= mymsh->GetDimension();
  unsigned 		nel	= mymsh->GetElementNumber();
  unsigned 		igrid	= mymsh->GetGridNumber();
  unsigned 		iproc	= mymsh->GetProcID();
  double		IPe	= 1./(ml_prob.parameters.get<Fluid>("Fluid").get_Peclet_number());  
  
  //solution variable
  unsigned SolIndex;  
  unsigned SolPdeIndex;
  SolIndex=ml_prob.GetIndex("T");
  SolPdeIndex=mylin_impl_sys.GetSolPdeIndex("T");
  //solution order
  unsigned order_ind = ml_prob.SolType[SolIndex];
  unsigned end_ind   = mymsh->GetEndIndex(order_ind);
  
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
  F.reserve(max_size);
  B.reserve(max_size*max_size);
  
  // Set to zeto all the entries of the Global Matrix
  if(assembe_matrix) myKK->zero();
  
  // *** element loop ***
 
  for (int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

    unsigned kel = mymsh->IS_Mts2Gmt_elem[iel];
    short unsigned kelt=myel->GetElementType(kel);
    unsigned nve=myel->GetElementDofNumber(kel,end_ind);
    
    // resize
    metis_node.resize(nve);
    KK_dof.resize(nve);
    phi.resize(nve);
    gradphi.resize(nve*dim);	
    for(int i=0;i<dim;i++){
      coordinates[i].resize(nve);
    }
    
    // set to zero all the entries of the FE matrices
    F.resize(nve);
    memset(&F[0],0,nve*sizeof(double));
    if(assembe_matrix){
      B.resize(nve*nve);
      memset(&B[0],0,nve*nve*sizeof(double));
    }
    
    // get local to global mappings
    for( unsigned i=0;i<nve;i++){
      unsigned inode=myel->GetElementVertexIndex(kel,i)-1u;	
      unsigned inode_metis=mymsh->GetMetisDof(inode,2);
      metis_node[i]=inode_metis;
      for(unsigned ivar=0; ivar<dim; ivar++) {
	coordinates[ivar][i]=(*mymsh->_coordinate->_Sol[ivar])(inode_metis);
      }
      KK_dof[i]=mylsyspde->GetKKDof(SolIndex,SolPdeIndex,inode);
    }
        
    if(igrid==gridn || !myel->GetRefinedElementIndex(kel)) {
      // *** Gauss poit loop ***
      for(unsigned ig=0;ig < ml_prob.type_elem[kelt][order_ind]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	(ml_prob.type_elem[kelt][order_ind]->*(ml_prob.type_elem[kelt][order_ind])->Jacobian_ptr)(coordinates,ig,weight,phi,gradphi);
	//Temperature and velocity current solution
	double SolT=0;
	vector < double > gradSolT(dim,0.);
	for(unsigned ivar=0; ivar<dim; ivar++){
	  gradSolT[ivar]=0; 
	}
	vector < double > SolU(dim,0.);
	vector < unsigned > SolIndexU(dim);
	SolIndexU[0]=ml_prob.GetIndex("U");
	SolIndexU[1]=ml_prob.GetIndex("V");
	if(dim==3) SolIndexU[2]=ml_prob.GetIndex("W");
	  	  
	unsigned SolType=ml_prob.GetSolType("T");
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
	  if(assembe_matrix){
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
	} // endif assembe_matrix
      } // end gauss point loop
    } // endif single element not refined or fine grid loop
    //--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the global Matrix/vector
      
    myRES->add_vector_blocked(F,KK_dof);
    if(assembe_matrix) myKK->add_matrix_blocked(B,KK_dof,KK_dof);  
  } //end list of elements loop for each subdomain
    
  myRES->close();
  if(assembe_matrix) myKK->close();
  
   // ***************** END ASSEMBLY *******************
  
}


