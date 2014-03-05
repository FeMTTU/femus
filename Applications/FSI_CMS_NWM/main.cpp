#include "ElemType.hpp"
#include "NonLinearMultiLevelProblem.hpp"
// Timedependent MultiGrid Header
#include "NonLinearTimeDependentMultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "PetscVector.hpp"
#include "LinearSolver.hpp"
#include "Solid.hpp"
#include "Fluid.hpp"
#include "Parameter.hpp"
#include <iostream>
#include "FemTTUInit.hpp"
using std::cout;
using std::endl;

// User defined functions
int AssembleMatrixResFSI(NonLinearMultiLevelProblem &nl_td_ml_prob, unsigned level, const unsigned &gridn, const unsigned &ipde, const bool &assembe_matrix);

double SetVariableTimeStep(const double time);

bool SetBoundaryCondition(const double &x, const double &y, const double &z,const char name[], 
			  double &value, const int FaceName, const double = 0.);

double InitVariables(const double &x, const double &y, const double &z,const char name[]);

bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber,const int &level);

//------------------------------------------------------------------------------------------------------------------

int main(int argc,char **args) {
  
  /// Init Petsc-MPI communicator
  FemTTUInit mpinit(argc,args,MPI_COMM_WORLD);

  unsigned short nm,nr;
  std::cout<<"#MULTIGRID levels? (>=1) \n";
  //std::cin>>nm;
  nm=3;

  std::cout<<"#MAX_REFINEMENT levels? (>=0) \n";
  //std::cin>>nr;
  nr=1;
  int tmp=nm;
  nm+=nr;
  nr=tmp;

  char *infile = new char [50];

  //input file
  sprintf(infile,"./input/mesh.comsolbenchmark");

  const double Lref = 1.;
  
  // Generate Parameter Object
  Parameter parameter(1.,1.);
  // Generate Solid Object
  Solid solid(parameter,200000,0.3,7850.,"Neo-Hookean");
  // Generate Fluid Object
  Fluid fluid(parameter,0.001,1000.,"Newtonian");

  NonLinearTimeDependentMultiLevelProblem nl_td_ml_prob(nm,nr,infile,"fifth",Lref,SetRefinementFlag);

  // END MESH =================================

  // Add fluid object
  nl_td_ml_prob.Add_Fluid(&fluid);
  // Add Solid Object
  nl_td_ml_prob.Add_Solid(&solid);

  //Start System Variables;===========================
  //Focus here is on VARIABLES first, rather than on Equations
  // generate solution vector
  nl_td_ml_prob.AddSolution("DX","biquadratic");
  nl_td_ml_prob.AddSolution("DY","biquadratic");
  //nl_td_ml_prob.AddSolution("DZ","biquadratic");
  nl_td_ml_prob.AddSolution("U","biquadratic");
  nl_td_ml_prob.AddSolution("V","biquadratic");
  //    nl_td_ml_prob.AddSolutionVector("W","biquadratic");
  nl_td_ml_prob.AddSolution("AX","biquadratic",1,0);
  nl_td_ml_prob.AddSolution("AY","biquadratic",1,0);
  //    nl_td_ml_prob.AddSolutionVector("AZ","biquadratic",1,0);
  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  nl_td_ml_prob.AddSolution("P","disc_linear",1);
  nl_td_ml_prob.AssociatePropertyToSolution("P","Pressure"); // Add this line

  //Initialize (update Init(...) function)
  nl_td_ml_prob.AttachInitVariableFunction(InitVariables);
  nl_td_ml_prob.Initialize("All");

  //Set Boundary (update Dirichlet(...) function)
  nl_td_ml_prob.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  nl_td_ml_prob.GenerateBdc("DX","Steady");
  nl_td_ml_prob.GenerateBdc("DY","Steady");
  //nl_td_ml_prob.GenerateBdc("DZ","Steady");
  nl_td_ml_prob.GenerateBdc("U","Time_dependent");
  nl_td_ml_prob.GenerateBdc("V","Steady");
  //    nl_td_ml_prob.GenerateBdc("W","Steady");
  nl_td_ml_prob.GenerateBdc("AX","Steady");
  nl_td_ml_prob.GenerateBdc("AY","Steady");
  //    nl_td_ml_prob.GenerateBdc("AZ","Steady");
  nl_td_ml_prob.GenerateBdc("P","Steady");

  //Set Time step information
  nl_td_ml_prob.SetTimeStep(0.005);
  nl_td_ml_prob.SetPrintTimeStep(1);
  nl_td_ml_prob.SetSaveTimeStep(33300);
  nl_td_ml_prob.SetNumTimeSteps(2);  //165   
  // nl_td_ml_prob.InitializeFromRestart(5);
  nl_td_ml_prob.AttachSetTimeStepFunction(SetVariableTimeStep);
  

  std::vector<std::string> mov_vars;
  mov_vars.resize(2);
  mov_vars[0] = "DX";
  mov_vars[1] = "DY";
  //mov_vars[2] = "DZ";
  nl_td_ml_prob.SetMovingMesh(mov_vars);
  nl_td_ml_prob.MarkStructureNode();
 

  //Solver Configuration 
 
  
  
  nl_td_ml_prob.AddPde("FSI");
  nl_td_ml_prob.AddSolutionToSolPdeIndex("FSI","DX"); 
  nl_td_ml_prob.AddSolutionToSolPdeIndex("FSI","DY");
  //nl_td_ml_prob.AddSolutionToSolPdeIndex("FSI","DZ");
  nl_td_ml_prob.AddSolutionToSolPdeIndex("FSI","U");
  nl_td_ml_prob.AddSolutionToSolPdeIndex("FSI","V");
  //nl_td_ml_prob.AddSolutionToSolPdeIndex("FSI","W");
  nl_td_ml_prob.AddSolutionToSolPdeIndex("FSI","P");
  
  // create Multigrid (PRLO, REST, MAT, VECs) based on MGIndex
  nl_td_ml_prob.CreatePdeStructure();
  
  nl_td_ml_prob.SetDirichletBCsHandling("FSI","Elimination");
  
  //Solver I (Gmres)
  //   nl_td_ml_prob.SetSmoother("Gmres");
  //   nl_td_ml_prob.SetTolerances("FSI",1.e-12,1.e-20,1.e+50,1);
  
  
  //Solver II (Vanka-smoother-MPSC)
  nl_td_ml_prob.AddStabilization("FSI",true);
  nl_td_ml_prob.SetSolverFineGrids("FSI","GMRES");
  nl_td_ml_prob.SetPreconditionerFineGrids("FSI","LU");
  nl_td_ml_prob.SetVankaSchurOptions(false,1);
  nl_td_ml_prob.SetTolerances("FSI",1.e-12,1.e-20,1.e+50,1);
  nl_td_ml_prob.SetSchurTolerances("FSI",1.e-12,1.e-20,1.e+50,4);
  nl_td_ml_prob.SetDimVankaBlock("FSI",4);                //2^lev 1D 4^lev 2D 8^lev 3D

  //End System Variables; ==============================

  // START EQUATIONS =================================

  // Start FSI Muligrid Block
  nl_td_ml_prob.AttachAssembleFunction(AssembleMatrixResFSI);


  
  // create index of solutions to be to used in the Vanka Smoother
  nl_td_ml_prob.ClearVankaIndex();
  nl_td_ml_prob.AddToVankaIndex("FSI","DX");
  nl_td_ml_prob.AddToVankaIndex("FSI","DY");
  nl_td_ml_prob.AddToVankaIndex("FSI","U");
  nl_td_ml_prob.AddToVankaIndex("FSI","V");
  nl_td_ml_prob.AddToVankaIndex("FSI","P");

  
  for (unsigned time_step = nl_td_ml_prob.GetInitTimeStep(); time_step < nl_td_ml_prob.GetInitTimeStep() + nl_td_ml_prob.GetNumTimeSteps(); 
       time_step++) {
   
    //Solve with V-cycle or F-cycle
    nl_td_ml_prob.Solve("FSI",15,0,3,"V-Cycle");
  
    //The update of the acceleration must be done before the update of the other variables
    //update time step
    nl_td_ml_prob._NewmarkAccUpdate();

    //update Solution
    nl_td_ml_prob._UpdateSolution();

    //print solution for restart
    if ( !(time_step%nl_td_ml_prob.GetSaveTimeStep()) ) {
      nl_td_ml_prob.SaveData();
    }

    // print solution
    if ( !(time_step%nl_td_ml_prob.GetPrintTimeStep()) ) {
       
      std::vector<std::string> print_vars;
      print_vars.resize(5);
      print_vars[0] = "DX";
      print_vars[1] = "DY";
      print_vars[2] = "U";
      print_vars[3] = "V";
      print_vars[4] = "P";

      nl_td_ml_prob.printsol_vtu_inline("biquadratic",print_vars);
      nl_td_ml_prob.printsol_vtu_inline("linear",print_vars);
      //nl_td_ml_prob.printsol_xdmf_hdf5("biquadratic",print_vars);
       
    }
  
  } //end loop timestep
  
  //print the XDMF time archive
  //nl_td_ml_prob.printsol_xdmf_archive("biquadratic");

  // Delete Multigrid (PRLO, REST, MAT, VECs) based on MGIndex
  nl_td_ml_prob.DeletePdeStructure();

  // End FSI Muligrid Block
  // Destroy the last PETSC objects
  nl_td_ml_prob.FreeMultigrid();

  delete [] infile;
  return(0);
}

//---------------------------------------------------------------------------------------------------------------------

bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber,const int &level) {
  bool refine=0;

  //refinemenet based on Elemen Group Number
  if (ElemGroupNumber==5) refine=1;
  if (ElemGroupNumber==6) refine=1;

  return refine;

}

//---------------------------------------------------------------------------------------------------------------------

double SetVariableTimeStep(const double time) {
  if (time < 0.75) {
    return 0.005;
  }
  else {
    return 0.25;
  }
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryCondition(const double &x, const double &y, const double &z,const char name[], double &value, const int FaceName, const double time) {
  bool test=1; //Dirichlet
  value=0.;
  //   cout << "Time bdc : " <<  time << endl;
  if (!strcmp(name,"U")) {
    if (1==FaceName) { //inflow
      test=1;
      //comsol Benchmark
      value = (0.05*time*time)/(sqrt( (0.04 - time*time)*(0.04 - time*time) + (0.1*time)*(0.1*time) ))*y*(0.0001-y)*4.*100000000;
    }
    else if (2==FaceName ) {  //outflow
      test=0;
      //    test=1;
      value=0.;
    }
    else if (3==FaceName ) {  // no-slip fluid wall
      test=1;
      value=0.;
    }
    else if (4==FaceName ) {  // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if (!strcmp(name,"V")) {
    if (1==FaceName) {          //inflow
      test=1;
      value=0.;
    }
    else if (2==FaceName ) {    //outflow
      test=0;
      //    test=1;
      value=0.;
    }
    else if (3==FaceName ) {    // no-slip fluid wall
      test=1;
      value=0;
    }
    else if (4==FaceName ) {    // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if (!strcmp(name,"W")) {
    if (1==FaceName) {
      test=1;
      value=0.;
    }
    else if (2==FaceName ) {
      test=1;
      value=0.;
    }
    else if (3==FaceName ) {
      test=1;
      value=0.;
    }
    else if (4==FaceName ) {
      test=1;
      value=0.;
    }
  }
  else if (!strcmp(name,"P")) {
    if (1==FaceName) {
      test=0;
      value=0.;
    }
    else if (2==FaceName ) {
      test=0;
      value=0.;
    }
    else if (3==FaceName ) {
      test=0;
      value=0.;
    }
    else if (4==FaceName ) {
      test=0;
      value=0.;
    }
  }
  else if (!strcmp(name,"DX")) {
    if (1==FaceName) {       //inflow
      test=1;
      value=0.;
    }
    else if (2==FaceName ) { //outflow
      test=1;
      value=0.;
    }
    else if (3==FaceName ) { // no-slip fluid wall
      test=0;
      value=0.;
    }
    else if (4==FaceName ) { // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if (!strcmp(name,"DY")) {
    if (1==FaceName) {       //inflow
      test=0;
      value=0.;
    }
    else if (2==FaceName ) { //outflow
      test=0;
      value=0.;
    }
    else if (3==FaceName ) { // no-slip fluid wall
      test=1;
      value=0.;
    }
    else if (4==FaceName ) { // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if (!strcmp(name,"DZ")) {
    if (1==FaceName) {       //inflow
      test=1;
      value=0.;
    }
    else if (2==FaceName ) { //outflow
      test=1;
      value=0.;
    }
    else if (3==FaceName ) { // no-slip fluid wall
      test=1;
      value=0.;
    }
    else if (4==FaceName ) { // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if (!strcmp(name,"AX")) {
    test=0;
    value=0;
  }
  else if (!strcmp(name,"AY")) {
    test=0;
    value=0;
  }
  else if (!strcmp(name,"AZ")) {
    test=0;
    value=0;
  }
  return test;
}

//---------------------------------------------------------------------------------------------------------------------

double InitVariables(const double &x, const double &y, const double &z,const char name[]) {
  double value=0.;
  if (!strcmp(name,"U")) {
    value=0;
  }
  else if (!strcmp(name,"V")) {
    value=0;
  }
  else if (!strcmp(name,"W")) {
    value=0;
  }
  else if (!strcmp(name,"P")) {
    value=0;
  }
  else if (!strcmp(name,"DX")) {
    value=0;
  }
  else if (!strcmp(name,"DY")) {
    value=0;
  }
  else if (!strcmp(name,"DZ")) {
    value=0;
  }
  else if (!strcmp(name,"AX")) {
    value=0;
  }
  else if (!strcmp(name,"AY")) {
    value=0;
  }
  else if (!strcmp(name,"AZ")) {
    value=0;
  }
  return value;
}

//--------------------------------------------------------------------------------------------------------------------

int AssembleMatrixResFSI(NonLinearMultiLevelProblem &nl_td_ml_prob2, unsigned level, const unsigned &gridn, const unsigned &ipde, const bool &assembe_matrix) {

  clock_t AssemblyTime=0;
  clock_t start_time, end_time;
  PetscErrorCode ierr;
  
  //Static conversion from the Base class (NonLinearMultiLevelProblem) to the Derived class (NonLinearMultiLevelProblemTimeLoop)
  NonLinearTimeDependentMultiLevelProblem& nl_td_ml_prob = static_cast<NonLinearTimeDependentMultiLevelProblem&>(nl_td_ml_prob2);
 
  const char* pdename= nl_td_ml_prob.GetThisPdeName(ipde);
  
  
  //pointers and references
  Solution*		mysolution	=	nl_td_ml_prob._solution[level];
  LinearSolver*		mylsyspde	=	nl_td_ml_prob._LinSolver[ipde][level];
  mesh*			mymsh		=	nl_td_ml_prob._msh[level];
  elem*			myel		=	mymsh->el;
  SparseMatrix*		myKK		=	mylsyspde->_KK;
  NumericVector*	myRES		=	mylsyspde->_RES;
  vector <int>&		myKKIndex	=	mylsyspde->KKIndex;
  
  
  const unsigned dim = mymsh->GetDimension();
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  
  // local objects
  vector<double> SolVAR(3*dim+1);
  vector<double> SolOldVAR(3*dim+1);
  double GradSolVAR[10][3];
  double GradSolhatVAR[10][3];
  double GradSolOldVAR[10][3];
  double GradSolOldhatVAR[10][3];
  
  vector <int> metis_node2(max_size);
  vector <int> metis_node1(max_size);
  vector <bool> solidmark(max_size);
  
  vector <double > phi(max_size);
  vector <double > phi_hat(max_size);
  vector <double > phi_old(max_size);
  
  const double * phi1 = NULL;
    
  double Weight=0.;
  double Weight_nojac=0.;
  double Weight_hat=0.;
  double Weight_old=0.;
  
  vector <vector < double> > vx(dim);
  vector <vector < double> > vx_hat(dim);
  vector <vector < double> > vx_old(dim);
  
  for(int i=0;i<dim;i++){
    vx[i].resize(max_size);
    vx_hat[i].resize(max_size);
    vx_old[i].resize(max_size);
  }
  
  vector <double> gradphi(max_size*dim);
  vector <double> gradphi_hat(max_size*dim);
  vector <double> gradphi_old(max_size*dim);
  
  vector< vector< double > > Rhs(7);
  vector< vector< vector< double > > > B(7); 
  vector< vector< int > > dofsVAR(7); 
  
  for(int i=0;i<7;i++){
    B[i].resize(7);
  }
   
  // algorithm parameters
  double eps_pen = 1.e40;
  bool   _newton = 0;
  bool   penalty = mylsyspde->GetStabilization();
  
  // ------------------------------------------------------------------------
  // Physical parameters
  double _rhof = nl_td_ml_prob._fluid->get_density();
  double _rhos = nl_td_ml_prob._solid->get_density();
  double _mu_lame = nl_td_ml_prob._solid->get_lame_shear_modulus();
  double _lambda_lame = nl_td_ml_prob._solid->get_lame_lambda();
  double _mus= _mu_lame/_rhof;
  double _IRe = nl_td_ml_prob._fluid->get_IReynolds_number();
  double _lambda= _lambda_lame / _rhof;
  double _betafsi= _rhos / _rhof;
  double _betans=1.;
  int    _solid_model= nl_td_ml_prob._solid->get_physical_model();

  //physical quantity
  double Jnp1_hat=0.;
  double Jn_hat=0.;
  double I_bleft=0.;
  double I_e=0.;
  double Cauchy[3][3] = {
    { 1.,0.,0.},
    { 0.,1.,0.},
    { 0.,0.,1.}
  };
  double Cauchy_old[3][3] = {
    { 1.,0.,0.},
    { 0.,1.,0.},
    { 0.,0.,1.}
  };
  double F[3][3] = {
    { 1.,0.,0.},
    { 0.,1.,0.},
    { 0.,0.,1.}
  };
  double b_left[3][3] = {
    { 1.,0.,0.},
    { 0.,1.,0.},
    { 0.,0.,1.}
  };
  double e[3][3] = {
    { 1.,0.,0.},
    { 0.,1.,0.},
    { 0.,0.,1.}
  };
  double e_old[3][3] = {
    { 1.,0.,0.},
    { 0.,1.,0.},
    { 0.,0.,1.}
  };
  double tg_stiff_matrix[3][3] = {
    { 0.,0.,0.},
    { 0.,0.,0.},
    { 0.,0.,0.}
  };
  const double Id2th[3][3] = {
    { 1.,0.,0.},
    { 0.,1.,0.},
    { 0.,0.,1.}
  };

  //initialization C tensor: Saint-Venaint Kirchoff model
  double C_mat[3][3][3][3];
  for (int I=0; I<3; ++I) {
    for (int J=0; J<3; ++J) {
      for (int K=0; K<3; ++K) {
        for (int L=0; L<3; ++L) {
          C_mat[I][J][K][L] = 2.*_mus*Id2th[I][K]*Id2th[J][L];
        }
      }
    }
  }

  // ale map
  double _lambda_map=0.;
  double _mu_ale[3] = {1.,1.,1.};
 
  // gravity
  double _gravity[3]={0.,0.,0.};

  // -----------------------------------------------------------------


  // time discretization algorithm paramters
  double dt = nl_td_ml_prob.GetTimeStep();
  const double gamma = 0.5;
  const double gammaratio = (1.-gamma)/gamma;
  double _theta_ns=1.0;

  // space discretization parameters
  unsigned order_ind2 = nl_td_ml_prob.SolType[nl_td_ml_prob.GetIndex("U")];  
  unsigned end_ind2   = mymsh->GetEndIndex(order_ind2);

  unsigned order_ind1 = nl_td_ml_prob.SolType[nl_td_ml_prob.GetIndex("P")];  
  unsigned end_ind1   = mymsh->GetEndIndex(order_ind1);

  // mesh and procs
  unsigned nel    = mymsh->GetElementNumber();
  unsigned igrid  = mymsh->GetGridNumber();
  unsigned iproc  = mymsh->GetProcID();
  unsigned nprocs = mymsh->GetNumProcs();

  //----------------------------------------------------------------------------------
  //variable-name handling
  const char varname[10][3] = {"DX","DY","DZ","U","V","W","P","AX","AY","AZ"};
  const char coordname[3][2] = {"X","Y","Z"};
  unsigned indexVAR[7];
  unsigned indCOORD[3];
  unsigned indVAR[10];  
  unsigned SolType[10];  
  
  for(unsigned ivar=0; ivar<dim; ivar++) {
    indCOORD[ivar]=nl_td_ml_prob.GetIndex(&coordname[ivar][0]);
    indVAR[ivar]=nl_td_ml_prob.GetIndex(&varname[ivar][0]);
    indVAR[ivar+dim]=nl_td_ml_prob.GetIndex(&varname[ivar+3][0]);
    indVAR[ivar+2*dim+1]=nl_td_ml_prob.GetIndex(&varname[ivar+7][0]);
    
    SolType[ivar]=nl_td_ml_prob.GetSolType(&varname[ivar][0]);
    SolType[ivar+dim]=nl_td_ml_prob.GetSolType(&varname[ivar+3][0]);
    SolType[ivar+2*dim+1]=nl_td_ml_prob.GetSolType(&varname[ivar+7][0]);
    
    indexVAR[ivar]=nl_td_ml_prob.GetSolPdeIndex("FSI",&varname[ivar][0]);
    indexVAR[ivar+dim]=nl_td_ml_prob.GetSolPdeIndex("FSI",&varname[ivar+3][0]);
  }
  indexVAR[2*dim]=nl_td_ml_prob.GetSolPdeIndex("FSI",&varname[6][0]);
  indVAR[2*dim]=nl_td_ml_prob.GetIndex(&varname[6][0]);
  SolType[2*dim]=nl_td_ml_prob.GetSolType(&varname[6][0]);
  //----------------------------------------------------------------------------------
  
  start_time=clock();
    
  myKK->zero();
  
  /// *** element loop ***
  for(int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

    unsigned kel        = mymsh->IS_Mts2Gmt_elem[iel]; 
    short unsigned kelt = myel->GetElementType(kel);
    unsigned nve        = myel->GetElementDofNumber(kel,end_ind2);
    unsigned nve1       = myel->GetElementDofNumber(kel,end_ind1);
    int flag_mat        = myel->GetElementMaterial(kel);

    //*******************************************************************************************************
    
    //initialization of everything is in common fluid and solid
    
    //Rhs
    for(int i=0; i<2*dim; i++) {
      dofsVAR[i].resize(nve);	
      Rhs[indexVAR[i]].resize(nve);
      memset(&Rhs[indexVAR[i]][0],0,nve*sizeof(double));
    }
    dofsVAR[2*dim].resize(nve1);
    Rhs[indexVAR[2*dim]].resize(nve1);
    memset(&Rhs[indexVAR[2*dim]][0],0,nve1*sizeof(double));
      
    //Kinematic relation (solid) and ALE Map (fluid)
    for(int i=0; i<dim; i++) {
      B[indexVAR[i]][indexVAR[i]].resize(nve*nve);
      memset(&B[indexVAR[i]][indexVAR[i]][0],0,nve*nve*sizeof(double));
    }

    //Stiffness Matrix (solid)
    for(int i=0; i<dim; i++) {
      for(int j=0; j<dim; j++) {
	B[indexVAR[dim+i]][indexVAR[j]].resize(nve*nve);
	memset(&B[indexVAR[dim+i]][indexVAR[j]][0],0,nve*nve*sizeof(double));
      }
    }
      
    //Mass Matrix (solid and fluid) and Diffusion Matrix (fluid only) 
    for(int i=0; i<dim; i++) {
      B[indexVAR[dim+i]][indexVAR[dim+i]].resize(nve*nve);
      memset(&B[indexVAR[dim+i]][indexVAR[dim+i]][0],0,nve*nve*sizeof(double));
    }

    //Pressure gradient (fluid and solid)
    for(int i=0; i<dim; i++) {
      B[indexVAR[dim+i]][indexVAR[2*dim]].resize(nve*nve1);
      memset(&B[indexVAR[dim+i]][indexVAR[2*dim]][0],0,nve*nve1*sizeof(double));
    }

    // Pressure Mass Matrix
    B[indexVAR[2*dim]][indexVAR[2*dim]].resize(nve1*nve1);
    memset(&B[indexVAR[2*dim]][indexVAR[2*dim]][0],0,nve1*nve1*sizeof(double));
    
    if (flag_mat==2) { //initialization for fluid only
      // Fluid Continuity Matrix: divergence of the velocity
      for(int i=0; i<dim; i++) {
	B[indexVAR[2*dim]][indexVAR[dim+i]].resize(nve1*nve);
	memset(&B[indexVAR[2*dim]][indexVAR[dim+i]][0],0,nve1*nve*sizeof(double));
      }
    }
    else{ // initialization for solid only
      // Kinematic relation
      for(int i=0; i<dim; i++) {
	B[indexVAR[i]][indexVAR[dim+i]].resize(nve*nve);
	memset(&B[indexVAR[i]][indexVAR[dim+i]][0],0,nve*nve*sizeof(double));
      }
      // Solid Continuity Matrix: divergence of the displacemnet
      for(int i=0; i<dim; i++) {
	B[indexVAR[2*dim]][indexVAR[i]].resize(nve1*nve);
	memset(&B[indexVAR[2*dim]][indexVAR[i]][0],0,nve1*nve*sizeof(double));
      }
    }
    
    // ----------------------------------------------------------------------------------------
    // coordinates, displacement, velocity dofs
        
    for (unsigned i=0;i<nve;i++) {
      // gambit nodes
      unsigned inode=myel->GetElementVertexIndex(kel,i)-1u;
      // dof metis
      unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
      metis_node2[i]=inode_Metis;
      
      //unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
      // flag to know if the node "inode" lays on the fluid-solid interface
      solidmark[i]=myel->GetNodeRegion(inode); // to check
      for(int j=0; j<dim; j++) {
	//Updated coordinates (Moving frame)
        vx[j][i]= (*mysolution->_Sol[indCOORD[j]])(inode_Metis) + (*mysolution->_Sol[indVAR[j]])(inode_Metis);
	//Old coordinates (Moving frame)
        vx_old[j][i]= (*mysolution->_Sol[indCOORD[j]])(inode_Metis) + (*mysolution->_SolOld[indVAR[j]])(inode_Metis);
	//Fixed coordinates (Reference frame)
	vx_hat[j][i]= (*mysolution->_Sol[indCOORD[j]])(inode_Metis);  
	// displacement dofs
	dofsVAR[j][i]= mylsyspde->GetKKDof(indVAR[j],indexVAR[j],inode); 
	// velocity dofs
	dofsVAR[j+dim][i]= mylsyspde->GetKKDof(indVAR[j+dim],indexVAR[j+dim],inode);   
      }
    }

    // pressure dofs
    for (unsigned i=0;i<nve1;i++) {
      unsigned inode=(order_ind1<3)?(myel->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
      metis_node1[i]=mymsh->GetMetisDof(inode,SolType[2*dim]);
      dofsVAR[2*dim][i]=mylsyspde->GetKKDof(indVAR[2*dim],indexVAR[2*dim],inode);
    }
    // ----------------------------------------------------------------------------------------
       
    if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {
      /// *** Gauss point loop ***
      for (unsigned ig=0;ig < nl_td_ml_prob.type_elem[kelt][order_ind2]->GetGaussPointNumber(); ig++) {

	// *** get Jacobian and test function and test function derivatives in the moving frame***
	(nl_td_ml_prob.type_elem[kelt][order_ind2]->*(nl_td_ml_prob.type_elem[kelt][order_ind2])->Jacobian_ptr)(vx,ig,Weight,phi,gradphi);
	(nl_td_ml_prob.type_elem[kelt][order_ind2]->*(nl_td_ml_prob.type_elem[kelt][order_ind2])->Jacobian_ptr)(vx_old,ig,Weight_old,phi_old,gradphi_old);
	(nl_td_ml_prob.type_elem[kelt][order_ind2]->*(nl_td_ml_prob.type_elem[kelt][order_ind2])->Jacobian_ptr)(vx_hat,ig,Weight_hat,phi_hat,gradphi_hat);
	phi1=nl_td_ml_prob.type_elem[kelt][order_ind1]->GetPhi(ig);
	if (flag_mat==2) Weight_nojac = nl_td_ml_prob.type_elem[kelt][order_ind2]->GetGaussWeight(ig);

	// ---------------------------------------------------------------------------
	// displacement and velocity
	for(int i=0; i<2*dim; i++){
	  SolVAR[i]=0.;
	  SolOldVAR[i]=0.;
	  for(int j=0; j<dim; j++) {
	    GradSolVAR[i][j]=0.;
	    GradSolOldVAR[i][j]=0.;
	    GradSolhatVAR[i][j]=0.;
	    GradSolOldhatVAR[i][j]=0.;
	  }
	    
	  for (unsigned inode=0; inode<nve; inode++) {
	    unsigned sol_dof = metis_node2[inode];
	      
	    double soli = (*mysolution->_Sol[indVAR[i]])(sol_dof);
	    SolVAR[i]+=phi[inode]*soli;
	      
	    double soli_old = (*mysolution->_SolOld[indVAR[i]])(sol_dof);
	    SolOldVAR[i]+=phi[inode]*soli_old;
	      
	    for(int j=0; j<dim; j++) {
	      GradSolVAR[i][j]+=gradphi[inode*dim+j]*soli;
	      GradSolOldVAR[i][j]+=gradphi[inode*dim+j]*soli_old;
	      if(i<dim){ 
		GradSolhatVAR[i][j]   +=gradphi_hat[inode*dim+j]*soli;
		GradSolOldhatVAR[i][j]+=gradphi_hat[inode*dim+j]*soli_old;
	      }
	    }	      
	  }
	}  
		  
	// pressure
	SolVAR[2*dim]=0.;
	for (unsigned inode=0; inode<nve1; inode++) {
	  double soli = (*mysolution->_Sol[indVAR[2*dim]])(metis_node1[inode]);
	  SolVAR[2*dim]+=phi1[inode]*soli;
	}
  
	// acceleration (solid only)  
	if(flag_mat!=2){
	  for(int i=2*dim+1; i<3*dim+1; i++){
	    SolVAR[i]=0.;
	    for (unsigned inode=0; inode<nve; inode++) {
	      double soli = (*mysolution->_Sol[indVAR[i]])(metis_node2[inode]);
	      SolVAR[i]+=phi[i]*soli;  
	    }
	  }
	}
  
	// ---------------------------------------------------------------------------
	 	  
	//BEGIN FLUID ASSEMBLY ============
	  
	if(flag_mat==2){
	  
          //divergence of the velocity
	  double div_vel=0.;
	  double div_w=0.;
	  for(int i=0; i<dim; i++) {
	    div_vel+=GradSolVAR[dim+i][i];
	    div_w+=(GradSolVAR[i][i]-GradSolOldVAR[i][i])*(1./dt);
	  }
	  
          { // Laplace operator + adection operator + Mass operator

            const double *gradfi=&gradphi[0];
            const double *fi=&phi[0];

            // *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim,fi++) {

              //BEGIN RESIDUALS A block ===========================
 	      double LapmapVAR[3] = {0., 0., 0.};
 	      for(int idim=0; idim<dim; idim++) {
 		for(int idim2=0; idim2<dim; idim2++) {
 	          LapmapVAR[idim] += dt*( _mu_ale[idim2]*gradphi_hat[i*dim+idim2]*GradSolhatVAR[idim][idim2] );
 		}
 	      }
      
	      // Residual ALE equations
 	      for(int idim=0; idim<dim; idim++) {
 	        Rhs[indexVAR[idim]][i]+=(!solidmark[i])*(-LapmapVAR[idim]*Weight_nojac);
 	      }
	      
	      //-------------------------------------------------------------------------------
	      
	      double LapvelVAR[3]={0.,0.,0.};
	      double AdvaleVAR[3]={0.,0.,0.};
	      for(int idim=0.; idim<dim; idim++) {
		for(int idim2=0.; idim2<dim; idim2++) {
		  LapvelVAR[idim]+=gradphi[i*dim+idim2]*GradSolVAR[dim+idim][idim2];
		  AdvaleVAR[idim]+=((SolVAR[dim+idim2]*dt - (SolVAR[idim2]-SolOldVAR[idim2]))*GradSolVAR[dim+idim][idim2])*phi[i];;
		}
	      }
	      // Residual Momentum equations 
	      for(int idim=0; idim<dim; idim++) {
		Rhs[indexVAR[dim+idim]][i]+= (
					      -phi[i]*_betans*SolVAR[dim+idim]*Weight
					      +phi[i]*_betans*SolOldVAR[dim+idim]*Weight_old
					      -AdvaleVAR[idim]*0.5*(Weight+Weight_old)
					      -0.5*dt*div_vel*SolVAR[dim+idim]*phi[i]*0.5*(Weight+Weight_old)
					      +dt*div_w*SolVAR[dim+idim]*phi[i]*0.5*(Weight+Weight_old)
					      -dt*_IRe*LapvelVAR[idim]*Weight
					      +dt*SolVAR[2*dim]*gradphi[i*dim+idim]*Weight
					      );
	      }
	      //END RESIDUALS A block ===========================
	      
              const double *gradfj=&gradphi[0];
              const double *fj=&phi[0];
              //  *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim,fj++) {

                //Laplacian
                double Lap=0.;
		for(int idim=0; idim<dim; idim++) {
		  Lap+=(*(gradfi+idim))*(*(gradfj+idim));
		}
                double LapXweight = Lap*Weight;

                //advection term I
		double Adv1=0.;
		for(int idim=0; idim<dim; idim++) {
		  Adv1+= ((SolVAR[dim+idim]*dt - (SolVAR[idim]-SolOldVAR[idim]))*(*(gradfj+idim))*(*(fi)))*0.5*(Weight+Weight_old);
		}
				
		double div_stab = 0.5*div_vel*((*fi))*((*fj))*0.5*(Weight+Weight_old);
		
		double div_ale = div_w*((*fi))*((*fj))*0.5*(Weight+Weight_old);

                double Mass = ((*fi))*((*fj))*Weight;

		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[dim+idim]][i*nve+j] += dt*_IRe*LapXweight + Adv1 + dt*div_stab - dt*div_ale + _betans*Mass;
		}
		
                for(int idim=0; idim<dim; idim++) {
		  for(int idim2=0; idim2<dim; idim2++) {
		    B[indexVAR[dim+idim]][indexVAR[0+idim2]][i*nve+j] += _betans*SolVAR[dim+idim]*((*fi))*(*(gradfi+idim2))*Weight;
		  }
		}

		double Lap_ale=0.;
		for(int idim=0; idim<dim; idim++) {
		  Lap_ale+=_mu_ale[idim]*(*(gradfi+idim))*(*(gradfj+idim));
		}
		  
		// Laplacian ALE map
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[0+idim]][indexVAR[0+idim]][i*nve+j] += (!solidmark[i])*dt*Lap_ale*Weight_nojac;  
		}

              } // end phi_j loop
            } // end phi loop
          } // end A + Bt

          { //Gradient of Pressure operator

            const double *gradfi=&gradphi[0];
            const double *fi=phi1;

            /// *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim,fi++) {

              const double *fj=phi1;
              /// *** phi_j loop ***
              for (unsigned j=0; j<nve1; j++,fj++) {
                for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[2*dim]][i*nve1+j] -= dt*((*(gradfi+idim))*(*fj))*Weight;
		}
              } // end phi_j loop
            } // end phi_i loop
          } // End Bt
        
          { // Divergence of the Velocity operator
            const double *fi=phi1;
            // *** phi_i loop ***
            for (unsigned i=0; i<nve1; i++,fi++) {

              //BEGIN RESIDUALS B block ===========================
              Rhs[indexVAR[2*dim]][i] += -(-((*fi))*div_vel)*Weight;
              //END RESIDUALS  B block ===========================

              const double *gradfj=&gradphi[0];
              // *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim) {
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[2*dim]][indexVAR[idim+dim]][i*nve+j] -= ((*fi)*(*(gradfj+idim)))*Weight;
		}
              }
            }
          }
	}   
	//END FLUID ASSEMBLY ============
	//*******************************************************************************************************
	//BEGIN SOLID ASSEMBLY ============
	  
	else{
	  //------------------------------------------------------------------------------------------------------------
          if (_solid_model==0) {

            //computation of the stress tensor
            e[0][0] = GradSolhatVAR[0][0];
            e[0][1] = 0.5*(GradSolhatVAR[0][1] + GradSolhatVAR[1][0]);
            e[0][2] = 0.;
            e[0][2] = 0.*0.5*(GradSolhatVAR[0][2] + GradSolhatVAR[2][0]);

            e[1][0] = 0.5*(GradSolhatVAR[1][0] + GradSolhatVAR[0][1]);
            e[1][1] = GradSolhatVAR[1][1];
            e[1][2] = 0.;
            e[1][2] = 0.*0.5*(GradSolhatVAR[1][2] + GradSolhatVAR[2][1]);

            e[2][0] = 0.*0.5*(GradSolhatVAR[0][2] + GradSolhatVAR[2][0]);
            e[2][1] = 0.*0.5*(GradSolhatVAR[1][2] + GradSolhatVAR[2][1]);
            e[2][2] = 0.*GradSolhatVAR[2][2];

            I_e = e[0][0] + e[1][1] + e[2][2];


            // Cauchy stress tensor
            for (int irow=0; irow<3; ++irow) {
              for (int jcol=0; jcol<3; ++jcol) {
                //compressible
		// 	   Cauchy[irow][jcol] = _lambda*I_e*Id2th[irow][jcol] + 2*_mus*e[irow][jcol];
                //incompressible
                Cauchy[irow][jcol] = 2*_mus*e[irow][jcol];
              }
            }


            //computation of the older stress tensor
            e_old[0][0] = GradSolOldhatVAR[0][0];
            e_old[0][1] = 0.5*(GradSolOldhatVAR[0][1] + GradSolOldhatVAR[1][0]);
            e_old[0][2] = 0.;
            e_old[0][2] = 0.*0.5*(GradSolOldhatVAR[0][2] + GradSolOldhatVAR[2][0]);

            e_old[1][0] = 0.5*(GradSolOldhatVAR[1][0] + GradSolOldhatVAR[0][1]);
            e_old[1][1] = GradSolOldhatVAR[1][1];
            e_old[1][2] = 0.;
            e_old[1][2] = 0.*0.5*(GradSolOldhatVAR[1][2] + GradSolOldhatVAR[2][1]);

            e_old[2][0] = 0.*0.5*(GradSolOldhatVAR[0][2] + GradSolOldhatVAR[2][0]);
            e_old[2][1] = 0.*0.5*(GradSolOldhatVAR[1][2] + GradSolOldhatVAR[2][1]);
            e_old[2][2] = 0.*GradSolOldhatVAR[2][2];

            //Iold_e = e_old[0][0] + e_old[1][1] + e_old[2][2];


            // Cauchy stress tensor
            for (int irow=0; irow<3; ++irow) {
              for (int jcol=0; jcol<3; ++jcol) {
                //compressible
		// 	   Cauchy[irow][jcol] = _lambda*I_e*Id2th[irow][jcol] + 2*_mus*e[irow][jcol];
                //incompressible
                Cauchy_old[irow][jcol] = 2*_mus*e_old[irow][jcol];
              }
            }
          }

          else if (_solid_model==1) {

            //deformation gradient
            F[0][0] = 1. + GradSolhatVAR[0][0];
            F[0][1] = GradSolhatVAR[0][1];
            F[0][2] = GradSolhatVAR[0][2];

            F[1][0] = GradSolhatVAR[1][0];
            F[1][1] = 1. + GradSolhatVAR[1][1];
            F[1][2] = GradSolhatVAR[1][2];

            F[2][0] = GradSolhatVAR[2][0];
            F[2][1] = GradSolhatVAR[2][1];
            F[2][2] = 1. + GradSolhatVAR[2][2];

            Jnp1_hat =   F[0][0]*F[1][1]*F[2][2] + F[0][1]*F[1][2]*F[2][0] + F[0][2]*F[1][0]*F[2][1]
	      - F[2][0]*F[1][1]*F[0][2] - F[2][1]*F[1][2]*F[0][0] - F[2][2]*F[1][0]*F[0][1];

            // computation of the the three deformation tensor b
            for (int I=0; I<3; ++I) {
              for (int J=0; J<3; ++J) {
                b_left[I][J]=0.;
                for (int K=0; K<3; ++K) {
                  //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
                  b_left[I][J] += F[I][K]*F[J][K];
                }
                Cauchy[I][J] = (_mus/Jnp1_hat)*(b_left[I][J] - Id2th[I][J]);
              }
            }

            I_bleft = b_left[0][0] + b_left[1][1] + b_left[2][2];
	    
	    //compressible case
	    //             for (int ii=0; ii<3; ++ii) {
	    //               for (int jj=0; jj<3; ++jj) {
	    //                 //Cauchy stress tensor
	    //                 Cauchy[ii][jj] = (_mus/Jnp1_hat)*(b_left[ii][jj] - Id2th[ii][jj]) + (_lambda/Jnp1_hat)*log(Jnp1_hat)*Id2th[ii][jj];
	    //                 for (int k=0; k<3; ++k) {
	    //                   for (int l=0; l<3; ++l) {
	    //                   }
	    //                 }
	    //               }
	    //             }
	    
	    //for the incompressible(nearly incompressible) case
            for (int ii=0; ii<3; ++ii) {
              for (int jj=0; jj<3; ++jj) {
                for (int kk=0; kk<3; ++kk) {
                  for (int ll=0; ll<3; ++ll) {
                    C_mat[ii][jj][kk][ll] = 2.*_mus*pow(Jnp1_hat,-1.6666666666666)*(
										    0.333333333333*I_bleft*Id2th[ii][kk]*Id2th[jj][ll]              //1/3*I_c*i
										    // 	                        +0.111111111111*I_C*Id2th[i][j]*Id2th[k][l]             //1/9*I_b*IxI
										    // 				-0.333333333333*b_left[i][j]*Id2th[k][l]                //-1/3*b*I
										    // 				-0.333333333333*Id2th[i][j]*b_left[k][l]                //-1/3*b*I
										    )
		      -SolVAR[2*dim]*(Id2th[ii][jj]*Id2th[kk][ll]-2.*Id2th[ii][kk]*Id2th[jj][ll] );  // -p(IxI-2i)

                  }
                }
              }
            }

            //Old deformation gradient
            F[0][0] = 1. + GradSolOldhatVAR[0][0];
            F[0][1] = GradSolOldhatVAR[0][1];
            F[0][2] = GradSolOldhatVAR[0][2];

            F[1][0] = GradSolOldhatVAR[1][0];
            F[1][1] = 1. + GradSolOldhatVAR[1][1];
            F[1][2] = GradSolOldhatVAR[1][2];

            F[2][0] = GradSolOldhatVAR[2][0];
            F[2][1] = GradSolOldhatVAR[2][1];
            F[2][2] = 1. + GradSolOldhatVAR[2][2];

            Jn_hat =   F[0][0]*F[1][1]*F[2][2] + F[0][1]*F[1][2]*F[2][0] + F[0][2]*F[1][0]*F[2][1]
	      - F[2][0]*F[1][1]*F[0][2] - F[2][1]*F[1][2]*F[0][0] - F[2][2]*F[1][0]*F[0][1] ;

            // computation of the the three deformation tensor b
            for (int I=0; I<3; ++I) {
              for (int J=0; J<3; ++J) {
                b_left[I][J]=0.;
                for (int k=0; k<3; ++k) {
                  //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
                  b_left[I][J] += F[I][k]*F[J][k];
                }
                Cauchy_old[I][J] = (_mus/Jn_hat)*(b_left[I][J] - Id2th[I][J]);
              }
            }



          }
          //----------------------------------------------------------------------------------------------------------------------------

          /////////////

          ///Mass + Stiffness operator
          {
            const double *gradfi=&gradphi[0];
            const double *fi=&phi[0];

            /// *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim,fi++) {


              //BEGIN RESIDUALS A + Bt block ===========================
	      
	      // Residual ALE equations
 	      for(int idim=0; idim<dim; idim++) {
	        Rhs[indexVAR[idim]][i] += (
					   -phi[i]*eps_pen*(
							    +SolVAR[idim] - SolOldVAR[idim]
							    -dt*gamma*SolVAR[dim+idim]
							    -dt*(1.-gamma)*SolOldVAR[dim+idim]
							    )
					   )*Weight_hat;
              }

              
              double CauchyDIR[3]={0.,0.,0.};
	      for(int idim=0.; idim<dim; idim++) {
		for(int idim2=0.; idim2<dim; idim2++) {
		  CauchyDIR[idim]+= gradphi[i*dim+idim2]*Cauchy[idim][idim2];
		}
	      }

              // Residual Momentum equations
              for(int idim=0; idim<dim; idim++) {
	        Rhs[indexVAR[dim+idim]][i] += (
					       phi[i]*dt*_gravity[idim]*Weight_hat
					       -phi[i]*_betafsi*(SolVAR[dim+idim] - SolOldVAR[dim+idim])*Weight_hat*(1./gamma)
					       +phi[i]*_betafsi*SolVAR[2*dim+1+idim]*Weight_hat*gammaratio*dt
					       -dt*CauchyDIR[idim]*Weight
					       +dt*SolVAR[2*dim]*gradphi[i*dim+idim]*Weight
					       );

              }
              
              //---------------------------------------------------------------------------------------------------------------------------------

              //END RESIDUALS A + Bt block ===========================

              const double *gradfj=&gradphi[0];
              const double *fj=&phi[0];
              // *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim,fj++) {

                /// Mass term
                // (v_n+1,psi)
                //gamma = 1 Backward Euler
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[dim+idim]][i*nve+j] += _betafsi*(*(fi))*(*(fj))*Weight_hat*(1./gamma);
		}

                //Da collaudare il 3D
                for (int icount=0; icount<dim; ++icount) {
                  for (int jcount=0; jcount<dim; ++jcount) {
                    tg_stiff_matrix[icount][jcount] = 0.;
                    for (int kcount=0; kcount<dim; ++kcount) {
                      for (int lcount=0; lcount<dim; ++lcount) {
                        tg_stiff_matrix[icount][jcount] += (*(gradfi+kcount))*0.25*(
										    C_mat[icount][kcount][jcount][lcount]+C_mat[icount][kcount][lcount][jcount]
										    +C_mat[kcount][icount][jcount][lcount]+C_mat[kcount][icount][lcount][jcount]
										    )*(*(gradfj+lcount));
                      }
                    }
                  }
                }
                
                //geometric tangent stiffness matrix
                double geom_tg_stiff_matrx = 0.;
                for(int kcount=0; kcount<dim; ++kcount) {
                  for(int lcount=0; lcount<dim; ++lcount) {
                    geom_tg_stiff_matrx += (*(gradfi+kcount))*Cauchy[kcount][lcount]*(*(gradfj+lcount));
                  }
                }


                /// Stiffness operator -- Elasticity equation (Linear or not)
                for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[0+idim]][i*nve+j] += dt*geom_tg_stiff_matrx*Weight;
 		  for(int idim2=0; idim2<dim; idim2++) {
 		    B[indexVAR[dim+idim]][indexVAR[0+idim2]][i*nve+j] += dt*tg_stiff_matrix[0+idim][0+idim2]*Weight;
 		  }
 		}

                /// Kinematic equation v = du/dt
                //   -_theta*dt*(v_n+1,eta)
		for(int idim=0; idim<dim; idim++) {
		  //   -_theta*dt*(v_n+1,eta)
		  B[indexVAR[0+idim]][indexVAR[dim+idim]][i*nve+j] -= eps_pen*gamma*dt*(*(fi))*(*(fj))*Weight_hat;
		  // (u_n+1,eta)
		  B[indexVAR[0+idim]][indexVAR[0+idim]][i*nve+j] += eps_pen*(*(fi))*(*(fj))*Weight_hat;
		}
	      }
            }
          }
          ////////////
          { ///Gradient of Pressure
            const double *gradfi=&gradphi[0];
            // *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim) {
              const double *fj=phi1;
              // *** phi_j loop ***
              for (unsigned j=0; j<nve1; j++,fj++) {
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[2*dim]][i*nve1+j] -= dt*((*(gradfi+idim))*(*fj))*Weight;
		}
              }
            }
          }
          ////////////
          { ///Divergence of the Displacement
            const double *fi=phi1;
            // *** phi_i loop ***
            for (unsigned i=0; i<nve1; i++,fi++) {

              //BEGIN RESIDUALS B block ===========================

              if (_solid_model==0) {
                Rhs[indexVAR[2*dim]][i] += -(-((*fi))*(I_e + (1./_lambda)*SolVAR[2*dim] ) )*Weight_hat;
              }
              else if (_solid_model==1) {
                Rhs[indexVAR[2*dim]][i] += -(-((*fi))*( log(Jnp1_hat)/Jnp1_hat + (1./_lambda)*SolVAR[2*dim] ) )*Weight_hat;
              }

              //END RESIDUALS B block ===========================

              const double *gradfj=&gradphi[0];
              // *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim) {
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[2*dim]][indexVAR[idim]][i*nve+j] -= ((*fi)*(*(gradfj+idim)))*Weight;
		}
              }
            }
          }
          //  /////////////
          {  ///Pressure Mass term
            const double *fi=phi1;
            // *** phi_i loop ***
            for (unsigned i=0; i<nve1; i++,fi++) {
              const double *fj=phi1;
              // *** phi_j loop ***
              for (unsigned j=0; j<nve1; j++,fj++) {
                B[indexVAR[2*dim]][indexVAR[2*dim]][i*nve1+j] -= (1./_lambda)*((*fi)*(*fj))*Weight_hat;
              }
            }
          }  //end pressure mass term
	  //---------------------------------------------------------------------------------------------------------------------------------
	}
	  
	//END SOLID ASSEMBLY ============
      }
    }

    //BEGIN local to global assembly 
    // ALE mapping
    for(int i=0; i<dim; i++) {
      myRES->add_vector_blocked(Rhs[indexVAR[i]],dofsVAR[i]);
      myKK ->add_matrix_blocked(B[indexVAR[i]][indexVAR[i]],dofsVAR[i],dofsVAR[i]);  
      if(flag_mat!=2){ //Solid only
	myKK->add_matrix_blocked(B[indexVAR[i]][indexVAR[i+dim]],dofsVAR[i],dofsVAR[i+dim]);
      }
    }
    
    // Momentum equation
    for(int i=0; i<dim; i++) {
      myRES->add_vector_blocked(Rhs[indexVAR[dim+i]],dofsVAR[dim+i]);
      myKK->add_matrix_blocked(B[indexVAR[dim+i]][indexVAR[dim+i]],dofsVAR[dim+i],dofsVAR[dim+i]);
      myKK->add_matrix_blocked(B[indexVAR[dim+i]][indexVAR[2*dim]],dofsVAR[dim+i],dofsVAR[2*dim]);
      for(int j=0; j<dim; j++) {
	myKK->add_matrix_blocked(B[indexVAR[dim+i]][indexVAR[j]],dofsVAR[dim+i],dofsVAR[j]);
      }
    }
   
    //P-continuity equation
    myRES->add_vector_blocked(Rhs[indexVAR[2*dim]],dofsVAR[2*dim]);
    for(int i=0; i<dim; i++) {
      if(flag_mat==2){ //Fluid only
	myKK->add_matrix_blocked(B[indexVAR[2*dim]][indexVAR[dim+i]],dofsVAR[2*dim],dofsVAR[dim+i]);
      }
      else{ //Solid only
	myKK->add_matrix_blocked(B[indexVAR[2*dim]][indexVAR[i]],dofsVAR[2*dim],dofsVAR[i]);  
      } 
    }
    myKK->add_matrix_blocked(B[indexVAR[2*dim]][indexVAR[2*dim]],dofsVAR[2*dim],dofsVAR[2*dim]);
    //END local to global assembly
   
  } //end list of elements loop


  //----------------------------------------------------------------------------------------------------------
  myKK->close();
  
  myRES->close();
  
  // *************************************
  end_time=clock();
  AssemblyTime+=(end_time-start_time);
  // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

  //----------------------------------------------------------------------------------------------


  return 0;

}

