#include "ElemType.hpp"
#include "NonLinearTimeDependentMultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "PetscVector.hpp"
#include "LinearSolver.hpp"
#include "Fluid.hpp"
#include "Parameter.hpp"
#include "FemTTUInit.hpp"
using std::cout;
using std::endl;
   

int AssembleMatrixResNS(NonLinearMultiLevelProblem &nl_td_ml_prob, unsigned level, const unsigned &gridn);

double InitVariables(const double &x, const double &y, const double &z,const char name[]);

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
  
  //Steadystate NonLinearMultiLevelProblem  
  NonLinearTimeDependentMultiLevelProblem nl_td_ml_prob(nm,nr,infile,"seventh",Lref,SetRefinementFlag);
  
  /// END MESH =================================  
  
  ///Start System Variables;===========================
  //Focus here is on VARIABLES first, rather than on Equations
 
  // add fluid material
  nl_td_ml_prob.Add_Fluid(&fluid);
  
  // generate solution vector
  nl_td_ml_prob.AddSolution("U","biquadratic");
  nl_td_ml_prob.AddSolution("V","biquadratic");
  // the pressure variable should be the last for the Schur decomposition
  nl_td_ml_prob.AddSolution("P","disc_linear",1);
  nl_td_ml_prob.AssociatePropertyToSolution("P","Pressure");
  
  //Initialize (update Init(...) function)
  nl_td_ml_prob.AttachInitVariableFunction(InitVariables);
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
    nl_td_ml_prob.FullMultiGrid("NS",15,1,1,"V-Cycle");
  
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

//--------------------------------------------------------------------------------------------------------------

double InitVariables(const double &x, const double &y, const double &z,const char name[]){
  double value=0.;
  if(!strcmp(name,"U")){
   value=0;
   //value=1.5*1*(4.0/(0.1681*100))*y*(0.41/0.1-y);
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

int AssembleMatrixResNS(NonLinearMultiLevelProblem &nl_td_ml_prob2, unsigned level, const unsigned &gridn){
    
  clock_t AssemblyTime=0; 
  clock_t start_time, end_time;
  start_time=clock();
  
  NonLinearTimeDependentMultiLevelProblem& nl_td_ml_prob = static_cast<NonLinearTimeDependentMultiLevelProblem&>(nl_td_ml_prob2);
  
  const char pdename[]="NS";
  unsigned ipde=nl_td_ml_prob.GetPdeIndex(pdename);
  
  
  //pointers and references
  Solution*       mysolution = nl_td_ml_prob2._solution[level];
  LinearSolverM*  mylsyspde  = nl_td_ml_prob2._LinSolver[ipde][level];
  mesh*           mymsh      = nl_td_ml_prob2._msh[level];
  elem*           myel       = mymsh->el;
  Mat&            myKK       = mylsyspde->KK;
  Vec&            myRES      = mylsyspde->RES;
  vector <int>&   myKKIndex  = mylsyspde->KKIndex; 
    
  // Allocation
  PetscInt node2[27];
  PetscInt node1[27];
  PetscInt nodeVAR[4][27];
  double B[4][4][27*27];
  double F[4][27];
  double coord[3][27];
  double phi2[27],gradphi2[27][3],Weight2;
  double normal[3];
  const double * phi1;
  PetscErrorCode ierr;  

  //data
  double dt = nl_td_ml_prob.GetTimeStep();
  double theta = 0.5;
  const unsigned dim = mymsh->GetDimension();
  unsigned nel=mymsh->GetElementNumber();
  unsigned igrid=mymsh->GetGridNumber();
  unsigned iproc = mymsh->GetProcID();
  unsigned nprocs = mymsh->GetNumProcs();
  double _ILambda= mylsyspde->GetCompressibility();
  double _IRe = nl_td_ml_prob._fluid->get_IReynolds_number();
  bool _penalty = mylsyspde->GetStabilization();
  const bool symm_mat = mylsyspde->GetMatrixProperties();
  const bool _NavierStokes = nl_td_ml_prob.GetNonLinearCase();
  unsigned _nwtn_alg = nl_td_ml_prob.GetNonLinearAlgorithm();
  bool _newton;
  if(_nwtn_alg==0) {
   _newton=0;
  } else {
    _newton=1;
  }
  
  //variable-name handling
  const char varname[4][2] = {"U","V","W","P"};
  const char coordname[3][2] = {"X","Y","Z"};
  unsigned SolPdeIndexVAR[4];
  unsigned SolIndexCOORD[3];
  unsigned SolIndexVAR[4];  
  double SolVAR[4];
  double SolOldVAR[4];
  double gradSolVAR[3][3];
  double gradSolOldVAR[3][3];
  for(unsigned ivar=0; ivar<dim; ivar++) {
    SolPdeIndexVAR[ivar]=nl_td_ml_prob.GetSolPdeIndex("NS",&varname[ivar][0]);
    SolIndexCOORD[ivar]=nl_td_ml_prob.GetIndex(&coordname[ivar][0]);
    SolIndexVAR[ivar]=nl_td_ml_prob.GetIndex(&varname[ivar][0]);
  }
  SolPdeIndexVAR[3]=nl_td_ml_prob.GetSolPdeIndex("NS",&varname[3][0]);
  SolIndexVAR[3]=nl_td_ml_prob.GetIndex(&varname[3][0]);
  
  //unknown order
  unsigned order_ind2 = nl_td_ml_prob2.SolType[nl_td_ml_prob2.GetIndex(&varname[0][0])];
  unsigned end_ind2   = mymsh->GetEndIndex(order_ind2);
  unsigned order_ind1 = nl_td_ml_prob2.SolType[nl_td_ml_prob2.GetIndex(&varname[3][0])];
  unsigned end_ind1   = mymsh->GetEndIndex(order_ind1);
  
  // Set to zeto all the entries of the matrix
  ierr = MatZeroEntries(myKK);  CHKERRQ(ierr);
  
  /// *** element loop ***
  for(int isdom=iproc; isdom<iproc+1; isdom++) {
    for (int iel=mymsh->IS_Mts2Gmt_elem_offset[isdom]; iel < mymsh->IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {

    unsigned kel = mymsh->IS_Mts2Gmt_elem[iel];
    short unsigned kelt=myel->GetElementType(kel);
    unsigned nve2=myel->GetElementDofNumber(kel,end_ind2);
    unsigned nve1=myel->GetElementDofNumber(kel,end_ind1);
    
    //set to zero all the entries of the FE matrices
    for(int ivar=0; ivar<dim; ivar++) {
      memset(F[SolPdeIndexVAR[ivar]],0,nve2*sizeof(double));
      memset(B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[ivar]],0,nve2*nve2*sizeof(double));
      memset(B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[3]],0,nve2*nve1*sizeof(double));
      memset(B[SolPdeIndexVAR[3]][SolPdeIndexVAR[ivar]],0,nve1*nve2*sizeof(double));
    }
    memset(F[SolPdeIndexVAR[3]],0,nve1*sizeof(double));
    
    if(_nwtn_alg==2){
      for(int ivar=0; ivar<dim; ivar++) {
	for(int ivar2=1; ivar2<dim; ivar2++) {
	  memset(B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[(ivar+ivar2)%dim]],0,nve2*nve2*sizeof(double));
	}
      }
    }
  
    if(_penalty){
      memset(B[SolPdeIndexVAR[3]][SolPdeIndexVAR[3]],0,nve1*nve1*sizeof(double));
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
      unsigned inode=(order_ind1<3)?(myel->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
      node1[i]=inode;
      nodeVAR[3][i]=mylsyspde->GetKKDof(SolIndexVAR[3],SolPdeIndexVAR[3],inode);
    }
   
    if(igrid==gridn || !myel->GetRefinedElementIndex(kel)) {
      // *** Gauss poit loop ***
      for(unsigned ig=0;ig < nl_td_ml_prob.type_elem[kelt][order_ind2]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	(nl_td_ml_prob.type_elem[kelt][order_ind2]->*(nl_td_ml_prob.type_elem[kelt][order_ind2])->Jacobian_ptr)(coord,ig,Weight2,phi2,gradphi2);
	phi1=nl_td_ml_prob.type_elem[kelt][order_ind1]->GetPhi(ig);

	//velocity variable
	for(unsigned ivar=0; ivar<dim; ivar++) {
	  SolVAR[ivar]=0.;
	  SolOldVAR[ivar]=0.;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	    gradSolVAR[ivar][ivar2]=0; 
	    gradSolOldVAR[ivar][ivar2]=0;
	  }
	  unsigned SolIndex=nl_td_ml_prob.GetIndex(&varname[ivar][0]);
	  unsigned SolType=nl_td_ml_prob.GetSolType(&varname[ivar][0]);
	  for(unsigned i=0; i<nve2; i++) {
	    unsigned sol_dof = mymsh->GetMetisDof(node2[i],SolType);
	    double soli = (*mysolution->_Sol[SolIndex])(sol_dof);
	    double sololdi = (*mysolution->_SolOld[SolIndex])(sol_dof);
	    SolVAR[ivar]+=phi2[i]*soli;
	    SolOldVAR[ivar]+=phi2[i]*sololdi; 
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      gradSolVAR[ivar][ivar2]    += gradphi2[i][ivar2]*soli; 
	      gradSolOldVAR[ivar][ivar2] += gradphi2[i][ivar2]*sololdi;
	    }
	  }
	}
	//pressure variable
	SolVAR[3]=0.;
	unsigned SolIndex=nl_td_ml_prob.GetIndex(&varname[3][0]);
	unsigned SolType=nl_td_ml_prob.GetSolType(&varname[3][0]);
	for(unsigned i=0; i<nve1; i++){
	  unsigned sol_dof = mymsh->GetMetisDof(node1[i],SolType);
	  double soli = (*mysolution->_Sol[SolIndex])(sol_dof);
	  SolVAR[3]+=phi1[i]*soli;
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
	      Lap_rhs     += gradphi2[i][ivar2]*gradSolVAR[ivar][ivar2];
	      Lap_old_rhs += gradphi2[i][ivar2]*gradSolOldVAR[ivar][ivar2];
	      Adv_rhs     += SolVAR[ivar2]*gradSolVAR[ivar][ivar2];
	      Adv_old_rhs += SolOldVAR[ivar2]*gradSolOldVAR[ivar][ivar2];
	    }
	    F[SolPdeIndexVAR[ivar]][i] += ( -theta*dt*_IRe*Lap_rhs                           // Laplacian
	                              -(1.-theta)*dt*_IRe*Lap_old_rhs                  // Laplacian
	                              -theta*dt*_NavierStokes*Adv_rhs*phi2[i]          // advection 
	                              -(1.-theta)*dt*_NavierStokes*Adv_old_rhs*phi2[i] // advection 
	                              +dt*SolVAR[3]*gradphi2[i][ivar]                  // pressure
	                              -(SolVAR[ivar] - SolOldVAR[ivar])*phi2[i]        // acceleration
	    )*Weight2;
	  }
	  
	  //END RESIDUALS A block ===========================

	  // *** phi_j loop ***
	  for(unsigned j=0; j<nve2; j++) {
		
	    double Lap=0;
	    double Adv1=0;
	    double Adv2 = phi2[i]*phi2[j]*Weight2;
	    double Mass = phi2[i]*phi2[j]*Weight2;
	    for(unsigned ivar=0; ivar<dim; ivar++) {
	      // Laplacian
	      Lap  += gradphi2[i][ivar]*gradphi2[j][ivar]*Weight2;
	      // advection term I
	      Adv1 += SolVAR[ivar]*gradphi2[j][ivar]*phi2[i]*Weight2;
	    }

	    for(unsigned ivar=0; ivar<dim; ivar++) {    
	      B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[ivar]][i*nve2+j] += theta*dt*_IRe*Lap + theta*dt*_NavierStokes*_newton*Adv1 + Mass;
	      if(_nwtn_alg==2){
		// Advection term II
		B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[ivar]][i*nve2+j] += theta*dt*Adv2*gradSolVAR[ivar][ivar];
		for(unsigned ivar2=1; ivar2<dim; ivar2++) {
		  B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[(ivar+ivar2)%dim]][i*nve2+j] += theta*dt*Adv2*gradSolVAR[ivar][(ivar+ivar2)%dim];
		}
	      }
	    }
  	    
	  }   //end phij loop
	} //end phii loop
  

	// *** phi_i loop ***
	for(unsigned i=0; i<nve2; i++){
	  // *** phi_j loop ***
	  for(unsigned j=0; j<nve1; j++){
	    for(unsigned ivar=0; ivar<dim; ivar++) {
	      B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[3]][i*nve1+j] -= dt*gradphi2[i][ivar]*phi1[j]*Weight2;
	    }
	  }
	} 
	
	// *** phi_i loop ***
	for(unsigned i=0; i<nve1; i++){
	  //BEGIN RESIDUALS B block ===========================
	  
	  double div = 0;
	  for(unsigned ivar=0; ivar<dim; ivar++) div += gradSolVAR[ivar][ivar];
	  F[SolPdeIndexVAR[3]][i]+= dt*(phi1[i]*div +_penalty*_ILambda*phi1[i]*SolVAR[3])*Weight2;
	  	  	    
	  //END RESIDUALS  B block ===========================

	  // *** phi_j loop ***
	  for(unsigned j=0; j<nve2; j++) {
	    for(unsigned ivar=0; ivar<dim; ivar++) {
	      B[SolPdeIndexVAR[3]][SolPdeIndexVAR[ivar]][i*nve2+j]-= dt*phi1[i]*gradphi2[j][ivar]*Weight2;
	    }
	  }  //end phij loop
	}  //end nve1
	
	if(_penalty){  //block nve1 nve1
	  // *** phi_i loop ***
	  for(unsigned i=0; i<nve1; i++){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve1; j++){
	      B[SolPdeIndexVAR[3]][SolPdeIndexVAR[3]][i*nve1+j]-= dt*_ILambda*phi1[i]*phi1[j]*Weight2;
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
	      nl_td_ml_prob.ComputeBdIntegral(pdename, &varname[ivar][0], kel, jface, level, ivar);
	    }
          }
        }
    }
//--------------------------------------------------------------------------------------------------------
    
    }  // endif single element not refined or fine grid loop
   
   
   //--------------------------------------------------------------------------------------------------------
   //Sum the small matrices into the Big Matrix
    for(unsigned ivar=0; ivar<dim; ivar++) {
      ierr = VecSetValues(myRES,nve2,nodeVAR[ivar],F[SolPdeIndexVAR[ivar]],ADD_VALUES);CHKERRQ(ierr);
      ierr = MatSetValuesBlocked(myKK,nve2,nodeVAR[ivar],nve2,nodeVAR[ivar],
				 B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[ivar]],ADD_VALUES);CHKERRQ(ierr);
      ierr = MatSetValuesBlocked(myKK,nve2,nodeVAR[ivar],nve1,nodeVAR[3],B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[3]],
				 ADD_VALUES);CHKERRQ(ierr);
      ierr = MatSetValuesBlocked(myKK,nve1,nodeVAR[3],nve2,nodeVAR[ivar],B[SolPdeIndexVAR[3]][SolPdeIndexVAR[ivar]],
				 ADD_VALUES);CHKERRQ(ierr);
      if(_nwtn_alg==2){
	for(unsigned ivar2=1; ivar2<dim; ivar2++) {
	  ierr = MatSetValuesBlocked(myKK,nve2,nodeVAR[ivar],nve2,nodeVAR[(ivar+ivar2)%dim],B[SolPdeIndexVAR[ivar]][SolPdeIndexVAR[(ivar+ivar2)%dim]],
				   ADD_VALUES);CHKERRQ(ierr);
	}
      }
    }
     
    //Penalty
    if(_penalty){
      ierr = MatSetValuesBlocked(myKK,nve1,nodeVAR[3],nve1,nodeVAR[3],B[SolPdeIndexVAR[3]][SolPdeIndexVAR[3]],ADD_VALUES);CHKERRQ(ierr);
    }
    ierr = VecSetValues(myRES,nve1,nodeVAR[3],F[SolPdeIndexVAR[3]],ADD_VALUES);CHKERRQ(ierr);
    //--------------------------------------------------------------------------------------------------------  
    
    } //end list of elements loop for each subdomain
  } //end list of subdomain
  

  //BEGIN MATRIX ASSEMBLY ============ 
  ierr = MatAssemblyBegin(myKK,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(myKK,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  //END MATRIX ASSEMBLY   ============ 
  
  
//    PetscViewer viewer;
// // PetscViewerSetFormat(viewer,PETSC_VIEWER_DEFAULT);
//    PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,700,700,&viewer);
//    MatView(myKK,viewer);
//    double ff;
//    std::cin>>ff;

  
   //BEGIN RESIDUAL ASSEMBLY ============ 
  ierr = VecAssemblyBegin(myRES);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(myRES);CHKERRQ(ierr);
  
  if(nprocs!=1) {
//     ierr = VecGhostUpdateBegin(myRES,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
//     ierr = VecGhostUpdateEnd(myRES,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(myRES,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(myRES,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  }
  //END RESIDUAL ASSEMBLY   ============ 
 
  // *************************************
  end_time=clock();
  AssemblyTime+=(end_time-start_time);
  // ***************** END ASSEMBLY *******************
  
  return ierr;
  
}

