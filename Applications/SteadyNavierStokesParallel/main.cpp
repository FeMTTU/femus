#include "ElemType.hpp"
#include "NonLinearMultiLevelProblem.hpp"
#include "Fluid.hpp"
#include "Parameter.hpp"
#include "FemTTUInit.hpp"
#include "NumericVector.hpp"
#include "LinearSolver.hpp"
#include "SparseMatrix.hpp"

using std::cout;
using std::endl;
   
int AssembleMatrixResNS(NonLinearMultiLevelProblem &nl_ml_prob, unsigned level, const unsigned &gridn, const unsigned &ipde, const bool &assembe_matrix);
int AssembleMatrixResT(NonLinearMultiLevelProblem &nl_ml_prob, unsigned level, const unsigned &gridn, const unsigned &ipde, const bool &assembe_matrix);


double InitVariables(const double &x, const double &y, const double &z,const char name[]);

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
  Parameter parameter(Lref,Uref);
  
  // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  Fluid fluid(parameter,0.001,1.,"Newtonian",0.001,1.);
  
  //Steadystate NonLinearMultiLevelProblem  
  NonLinearMultiLevelProblem nl_ml_prob(nm,nr,infile,"seventh",Lref,SetRefinementFlag);
  
  /// END MESH =================================  
  
  ///Start System Variables;===========================
  //Focus here is on VARIABLES first, rather than on Equations
 
  // add fluid material
  nl_ml_prob.Add_Fluid(&fluid);
  
  // generate solution vector
  nl_ml_prob.AddSolution("T","biquadratic");
  
  nl_ml_prob.AddSolution("U","biquadratic");
  nl_ml_prob.AddSolution("V","biquadratic");
  // the pressure variable should be the last for the Schur decomposition
  nl_ml_prob.AddSolution("P","disc_linear");
  // nl_ml_prob.AddSolution("P","linear");
  nl_ml_prob.AssociatePropertyToSolution("P","Pressure");
  //nl_ml_prob.AssociatePropertyToSolution("P","Default");
 
  //Initialize (update Init(...) function)
  nl_ml_prob.AttachInitVariableFunction(InitVariables);
  nl_ml_prob.Initialize("All");
  
  //Set Boundary (update Dirichlet(...) function)
  nl_ml_prob.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  nl_ml_prob.GenerateBdc("U");
  nl_ml_prob.GenerateBdc("V");
  nl_ml_prob.GenerateBdc("P");
  nl_ml_prob.GenerateBdc("T");
  
  ///End System Variables; ==============================

  // START PdeS =================================  
  
  nl_ml_prob.AddPde("NS1");
  nl_ml_prob.AddPde("NS2");
  nl_ml_prob.AddPde("Temp");
    
  /// Start Navier-Stokes Muligrid Block
  //start Multigrid for UVWP
  nl_ml_prob.ClearSolPdeIndex();
  nl_ml_prob.AddSolutionToSolPdeIndex("NS1","U"); 
  nl_ml_prob.AddSolutionToSolPdeIndex("NS1","V");
  nl_ml_prob.AddSolutionToSolPdeIndex("NS1","P");
  
  nl_ml_prob.AddSolutionToSolPdeIndex("NS2","U"); 
  nl_ml_prob.AddSolutionToSolPdeIndex("NS2","V");
  nl_ml_prob.AddSolutionToSolPdeIndex("NS2","P");
  
  nl_ml_prob.AddSolutionToSolPdeIndex("Temp","T");
  
  // create Multigrid (PRLO, REST, MAT, VECs) based on SolPdeIndex
  nl_ml_prob.CreatePdeStructure();
   
  
  nl_ml_prob.SetDirichletBCsHandling("NS1","Penalty");
  nl_ml_prob.SetDirichletBCsHandling("NS2","Penalty");
  nl_ml_prob.SetDirichletBCsHandling("Temp","Penalty");
 
//   nl_ml_prob.SetDirichletBCsHandling("NS1","Elimination");
//   nl_ml_prob.SetDirichletBCsHandling("NS2","Elimination");
//   nl_ml_prob.SetDirichletBCsHandling("Temp","Elimination");
  
  //Equation 1
  nl_ml_prob.AttachAssembleFunction(AssembleMatrixResNS);
  nl_ml_prob.SetNonLinearAlgorithm(true,"Newton",1.e-07);  //Navier-Stokes (Quasi-Newton - Newton)
  nl_ml_prob.SetMatrixProperties("NS1","Symmetric");
  nl_ml_prob.AddStabilization("NS1",true);
  
  //Solver Configuration 
  //Solver I (Gmres)
//   nl_ml_prob.SetSmoother("Gmres");
//   nl_ml_prob.SetTolerances("NS1",1.e-12,1.e-20,1.e+50,10);
//   // Solving
//   nl_ml_prob.FullMultiGrid("NS1",2,1,1,"F-Cycle");
   
  //Equation 2
  nl_ml_prob.AttachAssembleFunction(AssembleMatrixResNS);
  nl_ml_prob.SetNonLinearAlgorithm(true,"Newton",1.e-07);  //Navier-Stokes (Quasi-Newton - Newton)
  nl_ml_prob.SetMatrixProperties("NS2","Symmetric");
  nl_ml_prob.AddStabilization("NS2",true);
 
  //Solver Configuration 
  
  if(!vanka){
    nl_ml_prob.SetSmoother("Gmres"); 
    nl_ml_prob.SetTolerances("NS2",1.e-12,1.e-20,1.e+50,10);
    // Solving
    nl_ml_prob.Solve("NS2",3,1,1,"F-Cycle",!linear);
  }  
  else{ 
    // Solver II (Vanka - MPSC)
    // create index of solutions to be to used in the Vanka Smoother  
    nl_ml_prob.ClearVankaIndex();
    nl_ml_prob.AddToVankaIndex("NS2","U"); 
    nl_ml_prob.AddToVankaIndex("NS2","V"); 
    nl_ml_prob.AddToVankaIndex("NS2","P"); 
  
    nl_ml_prob.SetSmoother("Vanka");
    nl_ml_prob.SetVankaSchurOptions(false,1);
    nl_ml_prob.SetSolverFineGrids("NS2","GMRES");
    nl_ml_prob.SetPreconditionerFineGrids("NS2","ILU");
    nl_ml_prob.SetTolerances("NS2",1.e-12,1.e-20,1.e+50,10);
    nl_ml_prob.SetSchurTolerances("NS2",1.e-12,1.e-20,1.e+50,1);
    nl_ml_prob.SetDimVankaBlock("NS2",3);                             //2^lev 1D 4^lev 2D 8^lev 3D
    // Solving
    nl_ml_prob.Solve("NS2",3,1,1,"F-Cycle", !linear);
  }
  
  //Equation 3
  nl_ml_prob.AttachAssembleFunction(AssembleMatrixResT);
  nl_ml_prob.SetNonLinearAlgorithm(true,"Newton",1.e-07);  //Navier-Stokes (Quasi-Newton - Newton)
  nl_ml_prob.SetMatrixProperties("Temp","Symmetric");
  nl_ml_prob.AddStabilization("Temp",true);
  
  if(!vanka){
    //Solver Configuration 
    // Solver I (Gmres)
    nl_ml_prob.SetSmoother("Gmres");
    nl_ml_prob.SetTolerances("Temp",1.e-12,1.e-20,1.e+50,10);
    // Solving
    nl_ml_prob.Solve("Temp",3,1,1,"F-Cycle",linear);
  }
  else{ 
    nl_ml_prob.ClearVankaIndex();
    nl_ml_prob.AddToVankaIndex("Temp","T");

    nl_ml_prob.SetSmoother("Vanka");
    nl_ml_prob.SetVankaSchurOptions(false,0);
    nl_ml_prob.SetSolverFineGrids("Temp","GMRES");
    nl_ml_prob.SetPreconditionerFineGrids("Temp","ILU");
    nl_ml_prob.SetTolerances("Temp",1.e-12,1.e-20,1.e+50,10);
    nl_ml_prob.SetSchurTolerances("Temp",1.e-12,1.e-20,1.e+50,1);
    nl_ml_prob.SetDimVankaBlock("Temp",3);                             //2^lev 1D 4^lev 2D 8^lev 3D
    // Solving
    nl_ml_prob.Solve("Temp",3,1,1,"F-Cycle",linear);
  }
   
  // Delete Multigrid (PRLO, REST, MAT, VECs) based on SolPdeIndex
  nl_ml_prob.DeletePdeStructure();
  /// End Navier-Stokes Muligrid Block
  
  //Post processing
  //   double cforce[3];
  //   nl_ml_prob.ComputeBdStress(4,cforce);
  //   std::cout << "forcex:  " << cforce[0] << "   forcey: " << cforce[1] << std::endl;
  //   
  //   double cl[2];
  //   cl[0] = (2.*cforce[0])/(1.*0.2*0.2*0.1);
  //   cl[1] = (2.*cforce[1])/(1.*0.2*0.2*0.1);
  //   std::cout << "drag_force:  " << cl[0] << "   lift_force: " << cl[1] << std::endl;
  //   
  /// Print all solutions
  std::vector<std::string> print_vars;
  print_vars.resize(3);
  print_vars[0] = "U";
  print_vars[1] = "V";
  print_vars[2] = "P";
  
  bool DebugOn=1;
  //nl_ml_prob.printsol_gmv_binary("quadratic",nm-1,DebugOn);
  nl_ml_prob.printsol_gmv_binary("quadratic",nm,DebugOn);
  
  //nl_ml_prob.printsol_gmv_binary("linear",nm,DebugOn);
  //nl_ml_prob.printsol_gmv_binary("linear",nm-1,DebugOn);
  //nl_ml_prob.printsol_gmv_binary("linear",nm-2,DebugOn);
  
  nl_ml_prob.printsol_vtu_inline("biquadratic",print_vars);
  nl_ml_prob.printsol_vtu_inline("quadratic",print_vars);
  nl_ml_prob.printsol_vtu_inline("linear",print_vars);
  //nl_ml_prob.printsol_xdmf_hdf5("biquadratic",print_vars);
  
  /// Destroy the last PETSC objects
  nl_ml_prob.FreeMultigrid(); 
   
  delete [] infile;
  return(0);
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
  // if(ElemGroupNumber==7 && level<=3) refine=1;
  //    if(x>0 && x<5) refine=1;
  return refine;
}

//--------------------------------------------------------------------------------------------------------------

double InitVariables(const double &x, const double &y, const double &z,const char name[]){
  double value=0.;
  if(!strcmp(name,"U")){
    double um = 0.2;
    value=1.5*um*(4.0/(0.1681))*y*(0.41-y);
  }
  else if(!strcmp(name,"V")){
    value=0;
  } 
  else if(!strcmp(name,"W")){
    value=0;
  }
  else if(!strcmp(name,"P")){
    value=0;
  }
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


int AssembleMatrixResNS(NonLinearMultiLevelProblem &nl_ml_prob, unsigned level, const unsigned &gridn, const unsigned &ipde, const bool &assembe_matrix){
     
  const char* pdename=nl_ml_prob.GetThisPdeName(ipde);
  
  //pointers 
  Solution*	 mysolution  	= nl_ml_prob._solution[level];
  LinearSolver*  mylsyspde 	= nl_ml_prob._LinSolver[ipde][level];
  mesh*		 mymsh    	= nl_ml_prob._msh[level];
  elem*		 myel		= mymsh->el;
  SparseMatrix*	 myKK		= mylsyspde->_KK;
  NumericVector* myRES 		= mylsyspde->_RES;
    
  //parameters
  const unsigned dim = mymsh->GetDimension();
  unsigned nel= mymsh->GetElementNumber();
  unsigned igrid= mymsh->GetGridNumber();
  unsigned iproc = mymsh->GetProcID();
  unsigned nprocs = mymsh->GetNumProcs();
  double ILambda= mylsyspde->GetCompressibility();
  double IRe = nl_ml_prob._fluid->get_IReynolds_number();
  bool penalty = mylsyspde->GetStabilization();
  const bool symm_mat = mylsyspde->GetMatrixProperties();
  const bool NavierStokes = nl_ml_prob.GetNonLinearCase();
  unsigned nwtn_alg = nl_ml_prob.GetNonLinearAlgorithm();
  bool newton = (nwtn_alg==0) ? 0:1;
     
  // declare and reseve
  vector < int > node2; 
  vector < int > node1;
  vector< vector < double> > coord(dim);
  vector <double> phi2;
  vector <double> gradphi2;
  const double *phi1;
  vector< vector< int > > nodeVAR(dim+1); 
     
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  
  node2.reserve(max_size);
  node1.reserve( static_cast< unsigned > (ceil(pow(2,dim))));
  for(int i=0;i<dim;i++) {
    coord[i].reserve(max_size);
  }
  phi2.reserve(max_size);
  gradphi2.reserve(max_size*dim);	
  for(int i=0;i<dim;i++) {
    nodeVAR[i].reserve(max_size);
  }
  
  double Weight2;
  double normal[3];
   
  vector< vector< double > > F(dim+1);
  for(int i=0;i<dim+1;i++) F[i].reserve(max_size);
  
  vector< vector< vector< double > > > B(dim+1); 
  if(assembe_matrix){
    for(int i=0;i<dim+1;i++){
      B[i].resize(dim+1);
      for(int j=0;j<dim+1;j++){
	B[i][j].reserve(max_size*max_size);
      }
    }
  }
    
  //variable-name handling
  const char varname[4][2] = {"U","V","W","P"};
  const char coordname[3][2] = {"X","Y","Z"};
  vector < unsigned > SolPdeIndexVAR(dim+1);
  vector < unsigned > SolIndexCOORD(dim);
  vector < unsigned > SolIndexVAR(dim+1);  
    
  for(unsigned ivar=0; ivar<dim; ivar++) {
    SolPdeIndexVAR[ivar]=nl_ml_prob.GetSolPdeIndex(pdename,&varname[ivar][0]);
    SolIndexCOORD[ivar]=nl_ml_prob.GetIndex(&coordname[ivar][0]);
    SolIndexVAR[ivar]=nl_ml_prob.GetIndex(&varname[ivar][0]);
  }
  SolPdeIndexVAR[dim]=nl_ml_prob.GetSolPdeIndex(pdename, &varname[3][0]);
  SolIndexVAR[dim]=nl_ml_prob.GetIndex(&varname[3][0]);
        
  //unknown order
  unsigned order_ind2 = nl_ml_prob.SolType[SolIndexVAR[0]];
  unsigned end_ind2   = mymsh->GetEndIndex(order_ind2);
  unsigned order_ind1 = nl_ml_prob.SolType[SolIndexVAR[dim]];
  unsigned end_ind1   = mymsh->GetEndIndex(order_ind1);
  
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
    node2.resize(nve2);
    node1.resize(nve1);
    phi2.resize(nve2);
    gradphi2.resize(nve2*dim);
    for(int ivar=0; ivar<dim; ivar++) {
      coord[ivar].resize(nve2);
      nodeVAR[ivar].resize(nve2);
      
      F[SolPdeIndexVAR[ivar]].resize(nve2);
      memset(&F[SolPdeIndexVAR[ivar]][0],0,nve2*sizeof(double));
      
      if(assembe_matrix){
	B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[ivar]].resize(nve2*nve2);
	B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[dim]].resize(nve2*nve1);
	B[SolPdeIndexVAR[dim]][SolPdeIndexVAR[ivar]].resize(nve1*nve2);
	memset(&B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[ivar]][0],0,nve2*nve2*sizeof(double));
	memset(&B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[dim]][0],0,nve2*nve1*sizeof(double));
	memset(&B[SolPdeIndexVAR[dim]][SolPdeIndexVAR[ivar]][0],0,nve1*nve2*sizeof(double));
      }
    }
    nodeVAR[dim].resize(nve1);
    F[SolPdeIndexVAR[dim]].resize(nve1);
    memset(&F[SolPdeIndexVAR[dim]][0],0,nve1*sizeof(double));
      
      
    if(assembe_matrix*nwtn_alg==2){
      for(int ivar=0; ivar<dim; ivar++) {
	for(int ivar2=1; ivar2<dim; ivar2++) {
	  B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[(ivar+ivar2)%dim]].resize(nve2*nve2);
	  memset(&B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[(ivar+ivar2)%dim]][0],0,nve2*nve2*sizeof(double));
	}
      }
    }
  
    if(assembe_matrix*penalty){
      B[SolPdeIndexVAR[dim]][SolPdeIndexVAR[dim]].resize(nve1*nve1,0.);
      memset(&B[SolPdeIndexVAR[dim]][SolPdeIndexVAR[dim]][0],0,nve1*nve1*sizeof(double));
    }
    
    for( unsigned i=0;i<nve2;i++){
      unsigned inode=myel->GetElementVertexIndex(kel,i)-1u;	
      node2[i]=inode;
      unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
      for(unsigned ivar=0; ivar<dim; ivar++) {
	coord[ivar][i]=(*mysolution->_Sol[SolIndexCOORD[ivar]])(inode_Metis);
	nodeVAR[ivar][i]=mylsyspde->GetKKDof(SolIndexVAR[ivar],SolPdeIndexVAR[ivar],inode);
      }
    }
    
    for(unsigned i=0;i<nve1;i++) {
      unsigned inode=(order_ind1<dim)?(myel->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
      node1[i]=inode;
      nodeVAR[dim][i]=mylsyspde->GetKKDof(SolIndexVAR[dim],SolPdeIndexVAR[dim],inode);
    }
   
    if(igrid==gridn || !myel->GetRefinedElementIndex(kel)) {
      // *** Gauss poit loop ***
      for(unsigned ig=0;ig < nl_ml_prob.type_elem[kelt][order_ind2]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	(nl_ml_prob.type_elem[kelt][order_ind2]->*(nl_ml_prob.type_elem[kelt][order_ind2])->Jacobian_ptr)(coord,ig,Weight2,phi2,gradphi2);
	phi1=nl_ml_prob.type_elem[kelt][order_ind1]->GetPhi(ig);

	//velocity variable
	for(unsigned ivar=0; ivar<dim; ivar++) {
	  SolVAR[ivar]=0;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSolVAR[ivar][ivar2]=0; 
	  }
	  unsigned SolIndex=nl_ml_prob.GetIndex(&varname[ivar][0]);
	  unsigned SolType=nl_ml_prob.GetSolType(&varname[ivar][0]);
	  for(unsigned i=0; i<nve2; i++) {
	    unsigned sol_dof = mymsh->GetMetisDof(node2[i],SolType);
	    double soli = (*mysolution->_Sol[SolIndex])(sol_dof);
	    SolVAR[ivar]+=phi2[i]*soli;
	    for(unsigned ivar2=0; ivar2<dim; ivar2++){
	      gradSolVAR[ivar][ivar2] += gradphi2[i*dim+ivar2]*soli; 
	    }
	  }
	}
	//pressure variable
	SolVAR[dim]=0;
	unsigned SolIndex=nl_ml_prob.GetIndex(&varname[3][0]);
	unsigned SolType=nl_ml_prob.GetSolType(&varname[3][0]);
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
	    F[SolPdeIndexVAR[ivar]][i]+= (-IRe*Lap_rhs-NavierStokes*Adv_rhs*phi2[i]+SolVAR[dim]*gradphi2[i*dim+ivar])*Weight2; 
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
		B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[ivar]][i*nve2+j] += IRe*Lap + NavierStokes*newton*Adv1;
		if(nwtn_alg==2){
		  // Advection term II
		  B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[ivar]][i*nve2+j]       += Adv2*gradSolVAR[ivar][ivar];
		  for(unsigned ivar2=1; ivar2<dim; ivar2++) {
		    B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[(ivar+ivar2)%dim]][i*nve2+j] += Adv2*gradSolVAR[ivar][(ivar+ivar2)%dim];
		  }
		}
	      }
  	    } //end phij loop
	    
	    // *** phi1_j loop ***
	    for(unsigned j=0; j<nve1; j++){
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[dim]][i*nve1+j] -= gradphi2[i*dim+ivar]*phi1[j]*Weight2;
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
	  F[SolPdeIndexVAR[dim]][i]+= (phi1[i]*div +penalty*ILambda*phi1[i]*SolVAR[dim])*Weight2;
	  //END RESIDUALS  B block ===========================
	  
	  if(assembe_matrix){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve2; j++) {
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndexVAR[dim]][SolPdeIndexVAR[ivar]][i*nve2+j]-= phi1[i]*gradphi2[j*dim+ivar]*Weight2;
	      }
	    }  //end phij loop
	  } // endif assembe_matrix
	}  //end phi1_i loop
	
	if(assembe_matrix * penalty){  //block nve1 nve1
	  // *** phi_i loop ***
	  for(unsigned i=0; i<nve1; i++){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve1; j++){
	      B[SolPdeIndexVAR[dim]][SolPdeIndexVAR[dim]][i*nve1+j]-= ILambda*phi1[i]*phi1[j]*Weight2;
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
	      nl_ml_prob.ComputeBdIntegral(pdename, &varname[ivar][0], kel, jface, level, ivar);
	    }
	  }
	}	
      }
      //--------------------------------------------------------------------------------------------------------
    } // endif single element not refined or fine grid loop
    //--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the Global Matrix/Vector
    for(unsigned ivar=0; ivar<dim; ivar++) {
      myRES->add_vector_blocked(F[SolPdeIndexVAR[ivar]],nodeVAR[ivar]);
      if(assembe_matrix){
	myKK->add_matrix_blocked(B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[ivar]],nodeVAR[ivar],nodeVAR[ivar]);  
	myKK->add_matrix_blocked(B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[dim]],nodeVAR[ivar],nodeVAR[dim]);
	myKK->add_matrix_blocked(B[SolPdeIndexVAR[dim]][SolPdeIndexVAR[ivar]],nodeVAR[dim],nodeVAR[ivar]);
	if(nwtn_alg==2){
	  for(unsigned ivar2=1; ivar2<dim; ivar2++) {
	    myKK->add_matrix_blocked(B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[(ivar+ivar2)%dim]],nodeVAR[ivar],nodeVAR[(ivar+ivar2)%dim]);  
	  }
	}
      }
    }
    //Penalty
    if(assembe_matrix*penalty) myKK->add_matrix_blocked(B[SolPdeIndexVAR[dim]][SolPdeIndexVAR[dim]],nodeVAR[dim],nodeVAR[dim]);
    myRES->add_vector_blocked(F[SolPdeIndexVAR[dim]],nodeVAR[dim]);
    //--------------------------------------------------------------------------------------------------------  
  } //end list of elements loop for each subdomain
  
  
  if(assembe_matrix) myKK->close();
  myRES->close();
  // ***************** END ASSEMBLY *******************
  return 0;
}

//------------------------------------------------------------------------------------------------------------
int AssembleMatrixResT(NonLinearMultiLevelProblem &nl_ml_prob, unsigned level, const unsigned &gridn, const unsigned &ipde, const bool &assembe_matrix){
  
  
    
  //pointers and references
  Solution*      mysolution	= nl_ml_prob._solution[level];
  LinearSolver*  mylsyspde	= nl_ml_prob._LinSolver[ipde][level];
  mesh*          mymsh		= nl_ml_prob._msh[level];
  elem*          myel		= mymsh->el;
  SparseMatrix*  myKK		= mylsyspde->_KK;
  NumericVector* myRES		= mylsyspde->_RES;
   
  //data
  const unsigned	dim	= mymsh->GetDimension();
  unsigned 		nel	= mymsh->GetElementNumber();
  unsigned 		igrid	= mymsh->GetGridNumber();
  unsigned 		iproc	= mymsh->GetProcID();
  double		IPe	= 1./nl_ml_prob._fluid->get_Peclet_number();  
  
  //solution variable
  unsigned SolIndex;  
  unsigned SolPdeIndex;
  const char* pdename=nl_ml_prob.GetThisPdeName(ipde);
  SolIndex=nl_ml_prob.GetIndex("T");
  SolPdeIndex=nl_ml_prob.GetSolPdeIndex(pdename,"T");
  //solution order
  unsigned order_ind = nl_ml_prob.SolType[SolIndex];
  unsigned end_ind   = mymsh->GetEndIndex(order_ind);
  
  //coordinates
  vector< vector < double> > coordinates(dim); 
  const char coordinate_name[3][2] = {"X","Y","Z"};
  vector < unsigned > coordinate_Index(dim);
  for(unsigned ivar=0; ivar<dim; ivar++) {
    coordinate_Index[ivar]=nl_ml_prob.GetIndex(coordinate_name[ivar]);
  }
  
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
	coordinates[ivar][i]=(*mysolution->_Sol[coordinate_Index[ivar]])(inode_metis);
      }
      KK_dof[i]=mylsyspde->GetKKDof(SolIndex,SolPdeIndex,inode);
    }
        
    if(igrid==gridn || !myel->GetRefinedElementIndex(kel)) {
      // *** Gauss poit loop ***
      for(unsigned ig=0;ig < nl_ml_prob.type_elem[kelt][order_ind]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	(nl_ml_prob.type_elem[kelt][order_ind]->*(nl_ml_prob.type_elem[kelt][order_ind])->Jacobian_ptr)(coordinates,ig,weight,phi,gradphi);
	//Temperature and velocity current solution
	double SolT=0;
	vector < double > gradSolT(dim,0.);
	for(unsigned ivar=0; ivar<dim; ivar++){
	  gradSolT[ivar]=0; 
	}
	vector < double > SolU(dim,0.);
	vector < unsigned > SolIndexU(dim);
	SolIndexU[0]=nl_ml_prob.GetIndex("U");
	SolIndexU[1]=nl_ml_prob.GetIndex("V");
	if(dim==3) SolIndexU[2]=nl_ml_prob.GetIndex("W");
	  	  
	unsigned SolType=nl_ml_prob.GetSolType("T");
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
  
  return 0;
}


