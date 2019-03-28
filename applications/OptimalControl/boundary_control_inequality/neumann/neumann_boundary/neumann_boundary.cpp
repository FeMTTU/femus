#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "NumericVector.hpp"

#include "ElemType.hpp"

#define FACE_FOR_CONTROL 3  //we do control on the right (=2) face
#define AXIS_DIRECTION_CONTROL_SIDE  0  //change this accordingly to the other variable above
#include "../../param.hpp"

using namespace femus;

double InitialValueActFlag(const std::vector < double >& x) {
  return 0.;
}

double InitialValueContReg(const std::vector < double >& x) {
  return ControlDomainFlag_bdry(x);
}

double InitialValueTargReg(const std::vector < double >& x) {
  return ElementTargetFlag(x);
}

double InitialValueState(const std::vector < double >& x) {
  return 0.;
}

double InitialValueAdjoint(const std::vector < double >& x) {
  return 0.;
}

double InitialValueMu(const std::vector < double >& x) {
  return 0.;
}

double InitialValueControl(const std::vector < double >& x) {
  return 0.;
}

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet = true; //dirichlet
  value = 0;

  if(!strcmp(name,"control")) { //To enforce Neumann condition of control on Gamma_c
  if (faceName == 3)
    dirichlet = false;
  }

  if(!strcmp(name,"state")) {  //"state" corresponds to the first block row (u=q)
  if (faceName == 3)
    dirichlet = false;
  }

  if(!strcmp(name,"adjoint")) { //To enforce Neumann condition of adjoint on Gamma_c
  if (faceName == 3)
    dirichlet = false;
  }

  if(!strcmp(name,"mu")) {
//       value = 0.;
//   if (faceName == 3)
    dirichlet = false;
  }  
  
//     if(!strcmp(name,"adjoint")) {  //"adjoint" corresponds to the third block row
//   if (faceName == 3)    value = 1.;
//   }

  
  return dirichlet;
}


void ComputeIntegral(const MultiLevelProblem& ml_prob);

void AssembleOptSys(MultiLevelProblem& ml_prob);


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
  // ======= Files ========================
  Files files; 
  files.CheckIODirectories();
  files.RedirectCout();

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;

  mlMsh.GenerateCoarseBoxMesh(NSUB_X,NSUB_Y,0,0.,1.,0.,1.,0.,0.,QUAD9,"fifth");
   //1: bottom  //2: right  //3: top  //4: left
  
 /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("state", LAGRANGE, FIRST);
  mlSol.AddSolution("control", LAGRANGE, FIRST);
  mlSol.AddSolution("adjoint", LAGRANGE, FIRST);
  mlSol.AddSolution("mu", LAGRANGE, FIRST);  
  mlSol.AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  mlSol.AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  const unsigned int fake_time_dep_flag = 2;
  const std::string act_set_flag_name = "act_flag";
  mlSol.AddSolution(act_set_flag_name.c_str(), LAGRANGE, FIRST,fake_time_dep_flag);               //this variable is not solution of any eqn, it's just a given field
  
  mlSol.Initialize("All");    // initialize all varaibles to zero

  mlSol.Initialize("state", InitialValueState);
  mlSol.Initialize("control", InitialValueControl);
  mlSol.Initialize("adjoint", InitialValueAdjoint);
  mlSol.Initialize("mu", InitialValueMu);
  mlSol.Initialize("TargReg", InitialValueTargReg);
  mlSol.Initialize("ContReg", InitialValueContReg);
  mlSol.Initialize(act_set_flag_name.c_str(), InitialValueActFlag);

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("state");
  mlSol.GenerateBdc("control");
  mlSol.GenerateBdc("adjoint");
  mlSol.GenerateBdc("mu");  //we need add this to make the matrix iterations work...

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);
  
  mlProb.SetFilesHandler(&files);

 // add system  in mlProb as a Linear Implicit System
  NonLinearImplicitSystemWithPrimalDualActiveSetMethod& system = mlProb.add_system < NonLinearImplicitSystemWithPrimalDualActiveSetMethod > ("LiftRestr");
  
  system.SetActiveSetFlagName(act_set_flag_name);
//   system.SetMaxNumberOfNonLinearIterations(50);

  system.AddSolutionToSystemPDE("state");  
  system.AddSolutionToSystemPDE("control");  
  system.AddSolutionToSystemPDE("adjoint");  
  system.AddSolutionToSystemPDE("mu");  
  
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleOptSys);
  
  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->SetDebugOutput(true);

  system.SetDebugNonlinear(true);
  system.SetDebugFunction(ComputeIntegral);  //weird error if I comment this line, I expect nothing to happen but something in the assembly gets screwed up in memory I guess
   
  // initialize and solve the system
  system.init();
  system.MGsolve();
 
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
  mlSol.GetWriter()->Write(files.GetOutputPath()/*DEFAULT_OUTPUTDIR*/, "biquadratic", variablesToBePrinted);

  return 0;
}


void AssembleOptSys(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

 //*************************************************** 
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  vector < vector < double > > x(dim);
  vector < vector < double> >  x_bdry(dim);
  for (unsigned i = 0; i < dim; i++) {
         x[i].reserve(maxSize);
	 x_bdry[i].reserve(maxSize);
  }

 //***************************************************  

 //***************************************************  
  double weight = 0.; // gauss point weight
  double weight_bdry = 0.; // gauss point weight on the boundary


 //******************** state ************************ 
 //*************************************************** 
  vector <double> phi_u;  // local test function
  vector <double> phi_u_x; // local test function first order partial derivatives
  vector <double> phi_u_xx; // local test function second order partial derivatives

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * dim);
  phi_u_xx.reserve(maxSize * dim2);
  
  unsigned solIndex_u    = mlSol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u     = mlSol->GetSolutionType(solIndex_u);    // get the finite element type for "state"
  unsigned solPdeIndex_u = mlPdeSys->GetSolPdeIndex("state");    // get the position of "state" in the pdeSys object

  vector < double >  sol_u;     sol_u.reserve(maxSize);
  vector< int > l2GMap_u;    l2GMap_u.reserve(maxSize);
 //***************************************************  
 //***************************************************  

  
 //********************* adjoint *********************
 //*************************************************** 
  vector <double> phi_adj;  // local test function
  vector <double> phi_adj_x; // local test function first order partial derivatives
  vector <double> phi_adj_xx; // local test function second order partial derivatives

  phi_adj.reserve(maxSize);
  phi_adj_x.reserve(maxSize * dim);
  phi_adj_xx.reserve(maxSize * dim2);
 
  unsigned solIndex_adj    = mlSol->GetIndex("adjoint");    // get the position of "state" in the ml_sol object
  unsigned solType_adj     = mlSol->GetSolutionType(solIndex_adj);    // get the finite element type for "state"
  unsigned solPdeIndex_adj = mlPdeSys->GetSolPdeIndex("adjoint");    // get the position of "state" in the pdeSys object

  vector < double >  sol_adj;   sol_adj.reserve(maxSize);
  vector < int > l2GMap_adj; l2GMap_adj.reserve(maxSize);

  //boundary adjoint shape functions  
  vector <double> phi_adj_bdry;  
  vector <double> phi_adj_x_bdry; 

  phi_adj_bdry.reserve(maxSize);
  phi_adj_x_bdry.reserve(maxSize * dim);

 //*************************************************** 
 //***************************************************  

  
 //********************* bdry cont *******************
 //*************************************************** 
  vector <double> phi_ctrl_bdry;  
  vector <double> phi_ctrl_x_bdry; 

  phi_ctrl_bdry.reserve(maxSize);
  phi_ctrl_x_bdry.reserve(maxSize * dim);
  
  unsigned solIndex_ctrl = mlSol->GetIndex("control");
  unsigned solType_ctrl = mlSol->GetSolutionType(solIndex_ctrl);

  unsigned solPdeIndex_ctrl = mlPdeSys->GetSolPdeIndex("control");

  std::string ctrl_name = "control";
  
 vector < double >  sol_ctrl;   sol_ctrl.reserve(maxSize);
 vector< int > l2GMap_ctrl;   l2GMap_ctrl.reserve(maxSize);
 //*************************************************** 
 //*************************************************** 
  
 //****************** mu ******************************  
 //***************************************************  
  unsigned solIndex_mu;
  solIndex_mu = mlSol->GetIndex("mu");    // get the position of "mu" in the ml_sol object
   
  unsigned solPdeIndex_mu;
  solPdeIndex_mu = mlPdeSys->GetSolPdeIndex("mu");
  
  unsigned solType_mu = mlSol->GetSolutionType(solIndex_mu);    // get the finite element type for "mu"
  vector < double >  sol_mu;   sol_mu.reserve(maxSize);
  vector < int > l2GMap_mu;   l2GMap_mu.reserve(maxSize);

  
  //************** act flag **************************** 
  std::string act_flag_name = "act_flag";
  unsigned int solIndex_act_flag = mlSol->GetIndex(act_flag_name.c_str());
  unsigned int solFEType_act_flag = mlSol->GetSolutionType(solIndex_act_flag); 
     if(sol->GetSolutionTimeOrder(solIndex_act_flag) == 2) {
       *(sol->_SolOld[solIndex_act_flag]) = *(sol->_Sol[solIndex_act_flag]);
     }
  
  //********* variables for ineq constraints *****************
  const int ineq_flag = INEQ_FLAG;
  const double c_compl = C_COMPL;
  vector < double/*int*/ >  sol_actflag;   sol_actflag.reserve(maxSize); //flag for active set
  vector < double >  ctrl_lower;   ctrl_lower.reserve(maxSize);
  vector < double >  ctrl_upper;   ctrl_upper.reserve(maxSize);
  //***************************************************  

 
 //***************************************************  
  //********* WHOLE SET OF VARIABLES *****************
  const int solType_max = 2;  //biquadratic

  const int n_unknowns = 4;
 
  vector< int > l2GMap_AllVars; // local to global mapping
  l2GMap_AllVars.reserve(n_unknowns*maxSize);
  
  vector< double > Res; // local redidual vector
  Res.reserve(n_unknowns*maxSize);

  vector < double > Jac;
  Jac.reserve( n_unknowns*maxSize * n_unknowns*maxSize);
 //*************************************************** 

  
 //********************* DATA ************************ 
  double u_des = DesiredTarget();
  double alpha = ALPHA_CTRL_BDRY;
  double beta  = BETA_CTRL_BDRY;
  double penalty_outside_control_boundary = 1.e50;       // penalty for zero control outside Gamma_c
  double penalty_strong_bdry = 1.e20;  // penalty for boundary equation on Gamma_c
  double penalty_ctrl = 1.e10;         //penalty for u=q
 //*************************************************** 
  
  
  if (assembleMatrix)  KK->zero();

    
 // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned kelGeom = msh->GetElementType(iel);    // element geometry type

 //******************** GEOMETRY **********************
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    for (int i = 0; i < dim; i++)  x[i].resize(nDofx);
    
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof // via local to global solution node

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

   // elem average point 
    vector < double > elem_center(dim);   
    for (unsigned j = 0; j < dim; j++) {  elem_center[j] = 0.;  }
    for (unsigned j = 0; j < dim; j++) {  
      for (unsigned i = 0; i < nDofx; i++) {
         elem_center[j] += x[j][i];
       }
    }
    
   for (unsigned j = 0; j < dim; j++) { elem_center[j] = elem_center[j]/nDofx; }
 //*************************************************** 
  
 //************* set target domain flag **************
   int target_flag = 0;
   target_flag = ElementTargetFlag(elem_center);
 //*************************************************** 
   
 //************** set control flag *******************
  int control_el_flag = 0;
        control_el_flag = ControlDomainFlag_bdry(elem_center);
  std::vector<int> control_node_flag(nDofx,0);
//   if (control_el_flag == 0) std::fill(control_node_flag.begin(), control_node_flag.end(), 0);
 //*************************************************** 
    
 //********************* state ***********************
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);    // number of solution element dofs
    sol_u   .resize(nDof_u);
    l2GMap_u.resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
      unsigned solDof_u = msh->GetSolutionDof(i,iel, solType_u);    // global to global mapping between solution node and solution dof
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);      // global extraction and local storage for the solution
      l2GMap_u[i] = pdeSys->GetSystemDof(solIndex_u, solPdeIndex_u, i, iel);    // global to global mapping between solution node and pdeSys dof
    }
 //********************* state ***********************

 //******************* bdry cont *********************
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);    // number of solution element dofs
    sol_ctrl    .resize(nDof_ctrl);
    l2GMap_ctrl.resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);    // global to global mapping between solution node and solution dof
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);      // global extraction and local storage for the solution
      l2GMap_ctrl[i] = pdeSys->GetSystemDof(solIndex_ctrl, solPdeIndex_ctrl, i, iel);    // global to global mapping between solution node and pdeSys dof
    } 
 //***************************************************  

 //**************** adjoint **************************
    unsigned nDof_adj  = msh->GetElementDofNumber(iel, solType_adj);    
    sol_adj    .resize(nDof_adj);
    l2GMap_adj.resize(nDof_adj);
    for (unsigned i = 0; i < sol_adj.size(); i++) {
      unsigned solDof_adj = msh->GetSolutionDof(i, iel, solType_adj);
      sol_adj[i] = (*sol->_Sol[solIndex_adj])(solDof_adj);
      l2GMap_adj[i] = pdeSys->GetSystemDof(solIndex_adj, solPdeIndex_adj, i, iel);
    } 
 //*************** adjoint **************************** 

 //*************** mu ********************************* 
    unsigned nDof_mu  = msh->GetElementDofNumber(iel, solType_mu);    // number of solution element dofs
    sol_mu   .resize(nDof_mu);
    l2GMap_mu.resize(nDof_mu);
    for (unsigned i = 0; i < sol_mu.size(); i++) {
      unsigned solDof_mu = msh->GetSolutionDof(i, iel, solType_mu);   // global to global mapping between solution node and solution dof
      sol_mu[i] = (*sol->_Sol[solIndex_mu])(solDof_mu);      // global extraction and local storage for the solution 
      l2GMap_mu[i] = pdeSys->GetSystemDof(solIndex_mu, solPdeIndex_mu, i, iel);   // global to global mapping between solution node and pdeSys dof
    }
 
 //********************* ALL VARS ********************* 
    unsigned nDof_AllVars = nDof_u + nDof_ctrl + nDof_adj + nDof_mu; 
    int nDof_max    =  nDof_u;   // TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
    if(nDof_adj > nDof_max) 
    {
      nDof_max = nDof_adj;
      }
    
    if(nDof_ctrl > nDof_max)
    {
      nDof_max = nDof_ctrl;
    }
    
    
    Res.resize(nDof_AllVars);
    std::fill(Res.begin(), Res.end(), 0.);

    Jac.resize(nDof_AllVars * nDof_AllVars);
    std::fill(Jac.begin(), Jac.end(), 0.);
    
    l2GMap_AllVars.resize(0);
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_u.begin(),l2GMap_u.end());
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_ctrl.begin(),l2GMap_ctrl.end());
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_adj.begin(),l2GMap_adj.end());
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_mu.begin(),l2GMap_mu.end());
 //***************************************************

    
 //===================================================   

	// Perform face loop over elements that contain some control face
	if (control_el_flag == 1) {
	  
	  double tau=0.;
	  vector<double> normal(dim,0);
	       
	  // loop on faces of the current element

	  for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {
            std::vector < double > xyz_bdc(3,0.);  //not being used, because the boundaries are identified by the face numbers
	    // look for boundary faces
	    if(el->GetFaceElementIndex(iel,jface) < 0) {
	      unsigned int face = -( msh->el->GetFaceElementIndex(iel,jface)+1);
	      
		
// 	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"U",tau,face,0.) && tau!=0.){
	      if(  face == 3) { //control face

 //===================================================   
		//we use the dirichlet flag to say: if dirichlet = true, we set 1 on the diagonal. if dirichlet = false, we put the boundary equation
	      bool  dir_bool = mlSol->GetBdcFunction()(xyz_bdc,ctrl_name.c_str(),tau,face,0.);

 //===================================================   

		
		unsigned nve_bdry = msh->GetElementFaceDofNumber(iel,jface,solType_ctrl);
            for (unsigned idim = 0; idim < dim; idim++) {  x_bdry[idim].resize(nve_bdry); }
		const unsigned felt_bdry = msh->GetElementFaceType(iel, jface);    
		for(unsigned i=0; i < nve_bdry; i++) {
		  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i);
                  unsigned iDof = msh->GetSolutionDof(i_vol, iel, xType);
		  for(unsigned idim=0; idim<dim; idim++) {
		      x_bdry[idim][i]=(*msh->_topology->_Sol[idim])(iDof);
		  }
		}

//========================================================================================================== 

         //************** update active set flag for current nonlinear iteration **************************** 
         // 0: inactive; 1: active_a; 2: active_b
            assert(nDof_mu == nDof_ctrl);
            sol_actflag.resize(nve_bdry/*nDof_mu*/);
            ctrl_lower.resize(nve_bdry/*nDof_mu*/);
            ctrl_upper.resize(nve_bdry/*nDof_mu*/);
           std::fill(sol_actflag.begin(), sol_actflag.end(), 0);
           std::fill(ctrl_lower.begin(), ctrl_lower.end(), 0.);
           std::fill(ctrl_upper.begin(), ctrl_upper.end(), 0.);
   
		      for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
//             for (unsigned i = 0; i < sol_actflag.size(); i++) {
        std::vector<double> node_coords_i(dim,0.);
        for (unsigned d = 0; d < dim; d++) node_coords_i[d] = x_bdry[d][i_bdry];
        ctrl_lower[i_bdry] = InequalityConstraint(node_coords_i,false);
        ctrl_upper[i_bdry] = InequalityConstraint(node_coords_i,true);

        if      ( (sol_mu[i_vol] + c_compl * (sol_ctrl[i_vol] - ctrl_lower[i_bdry] )) < 0 )  sol_actflag[i_bdry] = 1;
        else if ( (sol_mu[i_vol] + c_compl * (sol_ctrl[i_vol] - ctrl_upper[i_bdry] )) > 0 )  sol_actflag[i_bdry] = 2;
            }
            
        //************** act flag **************************** 
      for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
      unsigned solDof_actflag = msh->GetSolutionDof(i_vol, iel, solFEType_act_flag); 
      (sol->_Sol[solIndex_act_flag])->set(solDof_actflag,sol_actflag[i_bdry]);     
    }
    

 // ===================================================
 //node-based insertion on the boundary ===============
 // ===================================================
    
 //============= delta_mu row ===============================
      std::vector<double> Res_mu (nDof_mu); std::fill(Res_mu.begin(),Res_mu.end(), 0.);
      
      for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
//     for (unsigned i = 0; i < sol_actflag.size(); i++) {
      if (sol_actflag[i_bdry] == 0){  //inactive
         Res_mu [i_vol] = - ineq_flag * ( 1. * sol_mu[i_vol] - 0. ); 
// 	 Res_mu [i] = Res[nDof_u + nDof_ctrl + nDof_adj + i]; 
      }
      else if (sol_actflag[i_bdry] == 1){  //active_a 
	 Res_mu [i_vol] = - ineq_flag * ( c_compl *  sol_ctrl[i_vol] - c_compl * ctrl_lower[i_bdry]);
      }
      else if (sol_actflag[i_bdry] == 2){  //active_b 
	Res_mu [i_vol]  =  - ineq_flag * ( c_compl *  sol_ctrl[i_vol] - c_compl * ctrl_upper[i_bdry]);
      }
    }

    
    RES->insert(Res_mu, l2GMap_mu);    
 //============= delta_mu row - end ===============================
    
 //============= delta_mu-delta_ctrl row ===============================
 //auxiliary volume vector for act flag
 unsigned nDof_actflag_vol  = msh->GetElementDofNumber(iel, solFEType_act_flag);
 std::vector<double> sol_actflag_vol(nDof_actflag_vol); 


 for (unsigned i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++) if (sol_actflag[i_bdry] != 0 ) sol_actflag[i_bdry] = ineq_flag * c_compl;    
 
 std::fill(sol_actflag_vol.begin(), sol_actflag_vol.end(), 0.);
    for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
       sol_actflag_vol[i_vol] = sol_actflag[i_bdry];
    }
 
 KK->matrix_set_off_diagonal_values_blocked(l2GMap_mu, l2GMap_ctrl, sol_actflag_vol);
 //============= delta_mu-delta_ctrl row - end ===============================

 //============= delta_mu-delta_mu row ===============================
  for (unsigned i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++) sol_actflag[i_bdry] =  ineq_flag * (1 - sol_actflag[i_bdry]/c_compl)  + (1-ineq_flag) * 1.;  //can do better to avoid division, maybe use modulo operator 

 std::fill(sol_actflag_vol.begin(), sol_actflag_vol.end(), 0.);
    for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
       sol_actflag_vol[i_vol] = sol_actflag[i_bdry];
    }
  
  KK->matrix_set_off_diagonal_values_blocked(l2GMap_mu, l2GMap_mu, sol_actflag_vol );
 //============= delta_mu-delta_mu row - end ===============================
    
 // =========================================================
 //node-based insertion on the boundary - end ===============
 // =========================================================
		
//========= initialize gauss quantities on the boundary =============================================
                double sol_ctrl_bdry_gss = 0.;
                double sol_adj_bdry_gss = 0.;
                std::vector<double> sol_ctrl_x_bdry_gss(dim);   std::fill(sol_ctrl_x_bdry_gss.begin(), sol_ctrl_x_bdry_gss.end(), 0.);

//========= initialize gauss quantities on the boundary ============================================
		
		for(unsigned ig_bdry=0; ig_bdry < msh->_finiteElement[felt_bdry][solType_ctrl]->GetGaussPointNumber(); ig_bdry++) {
		  
		  msh->_finiteElement[felt_bdry][solType_ctrl]->JacobianSur(x_bdry,ig_bdry,weight_bdry,phi_ctrl_bdry,phi_ctrl_x_bdry,normal);
		  msh->_finiteElement[felt_bdry][solType_adj]->JacobianSur(x_bdry,ig_bdry,weight_bdry,phi_adj_bdry,phi_adj_x_bdry,normal);

//========== temporary soln for surface gradient on a face parallel to the X axis ===================
		  double dx_dxi = 0.;
		 const elem_type_1D * myeltype = static_cast<const elem_type_1D*>(msh->_finiteElement[felt_bdry][solType_ctrl]);
		 const double * myptr = myeltype->GetDPhiDXi(ig_bdry);
		      for (int inode = 0; inode < nve_bdry/*_nc*/; inode++) dx_dxi += myptr[inode] * x_bdry[0][inode];
  
		      for (int inode = 0; inode < nve_bdry/*_nc*/; inode++) {
                            for (int d = 0; d < dim; d++) {
                              if (d==0 ) phi_ctrl_x_bdry[inode + d*nve_bdry/*_nc*/] = myptr[inode]* (1./ dx_dxi);
                              else  phi_ctrl_x_bdry[inode + d*nve_bdry/*_nc*/] = 0.;
                         }
                     }
//========== temporary soln for surface gradient on a face parallel to the X axis ====================
		  
//========== compute gauss quantities on the boundary ================================================
		  sol_ctrl_bdry_gss = 0.;
		  sol_adj_bdry_gss = 0.;
                  std::fill(sol_ctrl_x_bdry_gss.begin(), sol_ctrl_x_bdry_gss.end(), 0.);
		      for (int i_bdry = 0; i_bdry < nve_bdry/*_nc*/; i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
			
			sol_adj_bdry_gss  +=  sol_adj[i_vol] * phi_adj_bdry[i_bdry];
			sol_ctrl_bdry_gss +=  sol_ctrl[i_vol] * phi_ctrl_bdry[i_bdry];
                            for (int d = 0; d < dim; d++) {
			      sol_ctrl_x_bdry_gss[d] += sol_ctrl[i_vol] * phi_ctrl_x_bdry[i_bdry + d*nve_bdry];
			    }
		      }  

//========== compute gauss quantities on the boundary =================================================

		  // *** phi_i loop ***
		  for(unsigned i_bdry=0; i_bdry < nve_bdry; i_bdry++) {
		    
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);

                 double lap_rhs_dctrl_ctrl_bdry_gss_i = 0.;
                 for (unsigned d = 0; d < dim; d++) {
                       if ( i_vol < nDof_ctrl )  lap_rhs_dctrl_ctrl_bdry_gss_i +=  phi_ctrl_x_bdry[i_bdry + d*nve_bdry] * sol_ctrl_x_bdry_gss[d];
                 }
                 
//=============== construct control node flag field on the go  =========================================  
	      /* (control_node_flag[i])       picks nodes on \Gamma_c
	         (1 - control_node_flag[i])   picks nodes on \Omega \setminus \Gamma_c
	       */
	      if (dir_bool == false) { 
		std::cout << " found boundary control nodes ==== " << std::endl;
			for(unsigned k=0; k<control_node_flag.size(); k++) {
				  control_node_flag[i_vol] = 1;
			}
              }
//=============== construct control node flag field on the go  =========================================

//============ Bdry Residuals ==================	

        if (i_vol < nDof_u)     Res[ (0 + i_vol) ]                    += 0.;

        if (i_vol < nDof_ctrl)  Res[ (nDof_u + i_vol) ]               +=  - control_node_flag[i_vol] *  weight_bdry *
                                                              (     alpha * phi_ctrl_bdry[i_bdry] * sol_ctrl_bdry_gss
							         +  beta * lap_rhs_dctrl_ctrl_bdry_gss_i 
							         +  phi_ctrl_bdry[i_bdry]*sol_adj_bdry_gss
							        );  //boundary optimality condition

        if (i_vol < nDof_adj)   Res[ (nDof_u + nDof_ctrl +  i_vol) ]  += - control_node_flag[i_vol] *  weight_bdry * phi_adj_bdry[i_bdry]*sol_ctrl_bdry_gss; 
//============ Bdry Residuals ==================    
		    
		    for(unsigned j_bdry=0; j_bdry < nve_bdry; j_bdry ++) {
		    unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, jface, j_bdry);

		    
//============ boundary control eqns ===========
		    
// SECOND BLOCK ROW
//========block delta_control / control ========
   	           if ( i_vol < nDof_ctrl && j_vol < nDof_ctrl ) {
                     Jac[  
		        (nDof_u + i_vol) * nDof_AllVars +
	                (nDof_u + j_vol)             ]  += control_node_flag[i_vol] *  weight_bdry * (alpha * phi_ctrl_bdry[i_bdry] * phi_ctrl_bdry[j_bdry]);

                     double  lap_mat_dctrl_ctrl_bdry_gss = 0.;
		     for (unsigned d = 0; d < dim; d++) {  lap_mat_dctrl_ctrl_bdry_gss += phi_ctrl_x_bdry[i_bdry + d*nve_bdry] * phi_ctrl_x_bdry[j_bdry + d*nve_bdry];   }
		     
	             Jac[
		        (nDof_u + i_vol) * nDof_AllVars +
	                (nDof_u + j_vol)             ]  += control_node_flag[i_vol] * weight_bdry * beta *  lap_mat_dctrl_ctrl_bdry_gss;
                }
//========== block delta_control/adjoint ========
		   if ( i_vol < nDof_ctrl    && j_vol < nDof_adj )   
		     Jac[ 
			(nDof_u + i_vol) * nDof_AllVars +
		        (nDof_u + nDof_ctrl + j_vol) ]  += control_node_flag[i_vol] * (weight_bdry * phi_adj_bdry[j_bdry] * phi_ctrl_bdry[i_bdry]);      
		      

// THIRD BLOCK ROW
//========= block delta_adjoint/control =========
		   if ( i_vol < nDof_adj    && j_vol < nDof_ctrl )   
		     Jac[ 
			(nDof_u + nDof_ctrl + i_vol) * nDof_AllVars  +
		        (nDof_u + j_vol)             ]  += control_node_flag[i_vol] * (weight_bdry * phi_ctrl_bdry[j_bdry] * phi_adj_bdry[i_bdry]);      
			
//============ boundary control eqn =============			 
		   				
		    }  //end j loop
		  } //end i loop
		}
	      }
	    }
	  }    
	  
	} //end if control element flag
	
	else { //here we set the diagonal to 1 and the rhs to 0
	  
	  
	}
    
//========= gauss value quantities ==============   
	double sol_u_gss = 0.;
	double sol_adj_gss = 0.;
	std::vector<double> sol_u_x_gss(dim);     std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::vector<double> sol_adj_x_gss(dim);   std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
//===============================================   
 
 
      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][solType_max]->GetGaussPointNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
        //  ==== state 
	msh->_finiteElement[kelGeom][solType_u]  ->Jacobian(x, ig, weight, phi_u, phi_u_x, phi_u_xx);
        //  ==== adj 
        msh->_finiteElement[kelGeom][solType_adj]->Jacobian(x, ig, weight, phi_adj, phi_adj_x, phi_adj_xx);
          
	sol_u_gss = 0.;
	sol_adj_gss = 0.;
	std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	
	for (unsigned i = 0; i < nDof_u; i++) {
	                                                sol_u_gss      += sol_u[i] * phi_u[i];
                   for (unsigned d = 0; d < dim; d++)   sol_u_x_gss[d] += sol_u[i] * phi_u_x[i * dim + d];
          }
	
	for (unsigned i = 0; i < nDof_adj; i++) {
	                                                sol_adj_gss      += sol_adj[i] * phi_adj[i];
                   for (unsigned d = 0; d < dim; d++)   sol_adj_x_gss[d] += sol_adj[i] * phi_adj_x[i * dim + d];
        }

//==========FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {
	  
              double laplace_rhs_du_adj_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_u )         laplace_rhs_du_adj_i      +=  (phi_u_x   [i * dim + kdim] * sol_adj_x_gss[kdim]);
	      }
	      
              double laplace_rhs_dadj_u_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_adj )         laplace_rhs_dadj_u_i      +=  (phi_adj_x   [i * dim + kdim] * sol_u_x_gss[kdim]);
	      }
	      
//============ Volume residuals ==================	    
        if (i < nDof_u)     Res[0                  + i] += - weight * ( target_flag * phi_u[i] * ( sol_u_gss - u_des)
	                                                                  - laplace_rhs_du_adj_i); 
        if (i < nDof_ctrl)  Res[nDof_u             + i] += - penalty_outside_control_boundary * ( (1 - control_node_flag[i]) * (  sol_ctrl[i] - 0.)  );
	      
        if (i < nDof_adj)   Res[nDof_u + nDof_ctrl + i] += - weight * (-1) * (laplace_rhs_dadj_u_i);
            
        if (i < nDof_mu)    Res[nDof_u + nDof_ctrl + nDof_adj + i] += - penalty_outside_control_boundary * ( (1 - control_node_flag[i]) * (  sol_mu[i] - 0.)  );
//============  Volume Residuals ==================	    
	      
	      
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {
              double laplace_mat_dadj_u = 0.;
              double laplace_mat_du_adj = 0.;

              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_adj && j < nDof_u )     laplace_mat_dadj_u        +=  (phi_adj_x [i * dim + kdim] * phi_u_x   [j * dim + kdim]);
              if ( i < nDof_u   && j < nDof_adj )   laplace_mat_du_adj        +=  (phi_u_x   [i * dim + kdim] * phi_adj_x [j * dim + kdim]);
		
	      }

              //============ delta_state row ============================
              // BLOCK delta_state / state	      
              if ( i < nDof_u && j < nDof_u )   
		Jac[ i * nDof_AllVars +
		   (0 + j)                                ]  += weight * target_flag *  phi_u[i] * phi_u[j];   
		   
	      //BLOCK delta_state / adjoint
              if ( i < nDof_u && j < nDof_adj )   
		Jac[ i * nDof_AllVars +
		   (nDof_u + nDof_ctrl + j)               ]  += weight * (-1) * laplace_mat_du_adj;
	      
	      
              //=========== delta_control row ===========================
              //enforce control zero outside the control boundary
	      if ( i < nDof_ctrl && j < nDof_ctrl && i==j)
		Jac[    
		   (nDof_u + i) * nDof_AllVars  +
		   (nDof_u + j)                           ]  += penalty_outside_control_boundary * ( (1 - control_node_flag[i]));

		
	      //=========== delta_adjoint row ===========================
	      // BLOCK delta_adjoint / state
	      if ( i < nDof_adj && j < nDof_u ) {
		Jac[    
		   (nDof_u + nDof_ctrl + i) * nDof_AllVars  +
		   (0 + j)                                ]  += weight * (-1) * laplace_mat_dadj_u;
           
           //============= delta_mu row ===============================
	        if ( i < nDof_mu && j < nDof_mu && i==j )   
		  Jac[ (nDof_u + nDof_ctrl + nDof_adj + i) * nDof_AllVars +
		       (nDof_u + nDof_ctrl + nDof_adj + j)]  += penalty_outside_control_boundary * ( (1 - control_node_flag[i]));  

	      }
	    
	    } // end phi_j loop
	    
          } // endif assemble_matrix

        } // end phi_i loop
        
      } // end gauss point loop

/*	if (control_el_flag == 1) {
	  
         for (unsigned i = 0; i < nDof_max; i++) {
            for (unsigned j = 0; j < nDof_max; j++) {
	      std::cout << " " << std::setfill(' ') << std::setw(10) << Jac[ (nDof_u + i) * nDof_AllVars + (nDof_u + j) ];
// 	      std::cout << " " << std::setfill(' ') << std::setw(10) << Jac[ (nDof_u + i) * nDof_AllVars + (nDof_u + nDof_ctrl + j) ];
	    }
	      std::cout << std::endl;
	 }

    std::cout << " ========== " << iel << " ================== " << std::endl;      
         for (unsigned i = 0; i < nDof_max; i++) {
            for (unsigned j = 0; j < nDof_max; j++) {
// 	      std::cout << " " << std::setfill(' ') << std::setw(10) << Jac[ (nDof_u + i) * nDof_AllVars + (nDof_u + j) ];
	      std::cout << " " << std::setfill(' ') << std::setw(10) << Jac[ (nDof_u + i) * nDof_AllVars + (nDof_u + nDof_ctrl + j) ];
	    }
	      std::cout << std::endl;
	 }
	 
	 
	}
    std::cout << " ========== " << iel << " ================== " << std::endl;   */   
    
    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    RES->add_vector_blocked(Res, l2GMap_AllVars);

    if (assembleMatrix) {
      //store K in the global matrix KK
      KK->add_matrix_blocked(Jac, l2GMap_AllVars, l2GMap_AllVars);
    }
    
    
    //========== dof-based part, without summation
 


    
 //============= delta_ctrl-delta_mu row ===============================
 KK->matrix_set_off_diagonal_values_blocked(l2GMap_ctrl, l2GMap_mu, ineq_flag * 1.);//------------------------------->>>>>>
  
  
  } //end element loop for each process
  

  
  unsigned int ctrl_index = mlPdeSys->GetSolPdeIndex("control");
  unsigned int mu_index = mlPdeSys->GetSolPdeIndex("mu");

  unsigned int global_ctrl_size = pdeSys->KKoffset[ctrl_index+1][iproc] - pdeSys->KKoffset[ctrl_index][iproc];
  
  std::vector<double>  one_times_mu(global_ctrl_size, 0.);
  std::vector<int>    positions(global_ctrl_size);

  for (unsigned i = 0; i < positions.size(); i++) {
    positions[i] = pdeSys->KKoffset[ctrl_index][iproc] + i;
    one_times_mu[i] = ineq_flag * 1. * (*sol->_Sol[solIndex_mu])(i/*position_mu_i*/) ;
  }
    RES->add_vector_blocked(one_times_mu, positions);
    
  // ***************** END ASSEMBLY *******************

    if (assembleMatrix) KK->close();
    std::ostringstream mat_out; mat_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "matrix_" << mlPdeSys->GetNonlinearIt()  << ".txt";
    KK->print_matlab(mat_out.str(),"ascii"); //  KK->print();

    RES->close();
    std::ostringstream res_out; res_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "res_" << mlPdeSys->GetNonlinearIt()  << ".txt";
    std::filebuf res_fb;
    res_fb.open (res_out.str().c_str(),std::ios::out);
    std::ostream  res_file_stream(&res_fb);
    RES->print(res_file_stream);

  return;
}


 
  
void ComputeIntegral(const MultiLevelProblem& ml_prob)    {
  
  
  const NonLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystem> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

 //*************************************************** 
   const int xType = 2;
  vector < vector < double > > x(dim);
  vector < vector < double> >  x_bdry(dim);
  for (unsigned i = 0; i < dim; i++) {
         x[i].reserve(maxSize);
	 x_bdry[i].reserve(maxSize);
  }
 //***************************************************

 //***************************************************
  double weight; // gauss point weight
  double weight_bdry = 0.; // gauss point weight on the boundary

 //***************************************************
  double alpha = ALPHA_CTRL_BDRY;
  double beta  = BETA_CTRL_BDRY;
  
 //******** state ************************************ 
 //***************************************************
  vector <double> phi_u;     phi_u.reserve(maxSize);             // local test function
  vector <double> phi_u_x;   phi_u_x.reserve(maxSize * dim);     // local test function first order partial derivatives
  vector <double> phi_u_xx;  phi_u_xx.reserve(maxSize * dim2);   // local test function second order partial derivatives

 
  unsigned solIndex_u;
  solIndex_u = mlSol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u = mlSol->GetSolutionType(solIndex_u);    // get the finite element type for "state"

  vector < double >  sol_u; // local solution
  sol_u.reserve(maxSize);
  
  double u_gss = 0.;
 //***************************************************
 //***************************************************

  
 //************ desired ******************************
 //*************************************************** 
  vector <double> phi_udes;  // local test function
  vector <double> phi_udes_x; // local test function first order partial derivatives
  vector <double> phi_udes_xx; // local test function second order partial derivatives

    phi_udes.reserve(maxSize);
    phi_udes_x.reserve(maxSize * dim);
    phi_udes_xx.reserve(maxSize * dim2);
 
  
//  unsigned solIndexTdes;
//   solIndexTdes = mlSol->GetIndex("Tdes");    // get the position of "state" in the ml_sol object
//   unsigned solTypeTdes = mlSol->GetSolutionType(solIndexTdes);    // get the finite element type for "state"

  vector < double >  sol_udes; // local solution
  sol_udes.reserve(maxSize);

  double udes_gss = 0.;
 //***************************************************
 //***************************************************

 //************ cont *********************************
 //***************************************************
  vector <double> phi_ctrl_bdry;  
  vector <double> phi_ctrl_x_bdry; 

  phi_ctrl_bdry.reserve(maxSize);
  phi_ctrl_x_bdry.reserve(maxSize * dim);

  unsigned solIndex_ctrl = mlSol->GetIndex("control");
  unsigned solType_ctrl = mlSol->GetSolutionType(solIndex_ctrl);

   vector < double >  sol_ctrl;   sol_ctrl.reserve(maxSize);
 //***************************************************
 //***************************************************
  

 //*************************************************** 
 //********* WHOLE SET OF VARIABLES ****************** 
  const int solType_max = 2;  //biquadratic
 //***************************************************

  
 //********** DATA *********************************** 
  double u_des = DesiredTarget();
 //***************************************************
  
  double integral_target = 0.;
  double integral_alpha  = 0.;
  double integral_beta   = 0.;

    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned kelGeom = msh->GetElementType(iel);    // element geometry type
    
 //******************* GEOMETRY ********************** 
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    for (int i = 0; i < dim; i++)  x[i].resize(nDofx);
    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

   // elem average point 
    vector < double > elem_center(dim);   
    for (unsigned j = 0; j < dim; j++) {  elem_center[j] = 0.;  }
    for (unsigned j = 0; j < dim; j++) {  
      for (unsigned i = 0; i < nDofx; i++) {
         elem_center[j] += x[j][i];
       }
    }
    
   for (unsigned j = 0; j < dim; j++) { elem_center[j] = elem_center[j]/nDofx; }
 //*************************************************** 
  
 //***** set target domain flag ********************** 
   int target_flag = 0;
   target_flag = ElementTargetFlag(elem_center);
 //***************************************************

   
 //*********** state *********************************  
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);    // number of solution element dofs
    sol_u    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
      unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);    // global to global mapping between solution node and solution dof
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);      // global extraction and local storage for the solution
    }
 //*********** state ********************************* 


 //*********** cont ********************************** 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);    // number of solution element dofs
    sol_ctrl    .resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);    // global to global mapping between solution node and solution dof
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);      // global extraction and local storage for the solution
    } 

 //*********** cont ********************************** 
 
 
 //*********** udes ********************************** 
    unsigned nDof_udes  = msh->GetElementDofNumber(iel, solType_u);    // number of solution element dofs
    sol_udes    .resize(nDof_udes);
    for (unsigned i = 0; i < sol_udes.size(); i++) {
            sol_udes[i] = u_des;  //dof value
    } 
 //*********** udes ********************************** 

 
 //********** ALL VARS *******************************  
    int nDof_max    =  nDof_u;   // TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
    if(nDof_udes > nDof_max) 
    {
      nDof_max = nDof_udes;
      }
    
    if(nDof_ctrl > nDof_max)
    {
      nDof_max = nDof_ctrl;
    }
    
 //***************************************************

 // ==================================================
 //***** set control flag *************************** 
  int control_el_flag = 0;
        control_el_flag = ControlDomainFlag_bdry(elem_center);
  std::vector<int> control_node_flag(nDofx,0);
//   if (control_el_flag == 0) std::fill(control_node_flag.begin(), control_node_flag.end(), 0);
 //***************************************************

  
  	if (control_el_flag == 1) {
	  
	  double tau=0.;
	  vector<double> normal(dim,0);
	       
	  // loop on faces of the current element

	  for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {
            std::vector < double > xyz_bdc(3,0.);  //not being used, because the boundaries are identified by the face numbers
	    // look for boundary faces
	    if(el->GetFaceElementIndex(iel,jface) < 0) {
	      unsigned int face = -( msh->el->GetFaceElementIndex(iel,jface)+1);
	      
		
// 	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"U",tau,face,0.) && tau!=0.){
	      if(  face == 3) { //control face

		unsigned nve_bdry = msh->GetElementFaceDofNumber(iel,jface,solType_ctrl);
	        for (unsigned idim = 0; idim < dim; idim++) {  x_bdry[idim].resize(nve_bdry); }
		const unsigned felt_bdry = msh->GetElementFaceType(iel, jface);    
		for(unsigned i=0; i < nve_bdry; i++) {
		  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i);
                  unsigned iDof = msh->GetSolutionDof(i_vol, iel, xType);
		  for(unsigned idim=0; idim<dim; idim++) {
		      x_bdry[idim][i]=(*msh->_topology->_Sol[idim])(iDof);
		  }
		}
		
		//============ initialize gauss quantities on the boundary ==========================================
                double sol_ctrl_bdry_gss = 0.;
                std::vector<double> sol_ctrl_x_bdry_gss;  sol_ctrl_x_bdry_gss.reserve(dim);
		//============ initialize gauss quantities on the boundary ==========================================
		
		for(unsigned ig_bdry=0; ig_bdry < msh->_finiteElement[felt_bdry][solType_ctrl]->GetGaussPointNumber(); ig_bdry++) {
		  
		  msh->_finiteElement[felt_bdry][solType_ctrl]->JacobianSur(x_bdry,ig_bdry,weight_bdry,phi_ctrl_bdry,phi_ctrl_x_bdry,normal);

		 //========== temporary soln for surface gradient on a face parallel to the X axis ===================
		 double dx_dxi = 0.;
		 const elem_type_1D * myeltype = static_cast<const elem_type_1D*>(msh->_finiteElement[felt_bdry][solType_ctrl]);
		 const double * myptr = myeltype->GetDPhiDXi(ig_bdry);
		      for (int inode = 0; inode < nve_bdry/*_nc*/; inode++) dx_dxi += myptr[inode] * x_bdry[0][inode];
  
		      for (int inode = 0; inode < nve_bdry/*_nc*/; inode++) {
                            for (int d = 0; d < dim; d++) {
                              if (d==0 ) phi_ctrl_x_bdry[inode + d*nve_bdry/*_nc*/] = myptr[inode]* (1./ dx_dxi);
                              else  phi_ctrl_x_bdry[inode + d*nve_bdry/*_nc*/] = 0.;
                         }
                     }
		 //========== temporary soln for surface gradient on a face parallel to the X axis ===================
		  
		 //========= compute gauss quantities on the boundary ================================================
		  sol_ctrl_bdry_gss = 0.;
                  std::fill(sol_ctrl_x_bdry_gss.begin(), sol_ctrl_x_bdry_gss.end(), 0.);
		      for (int i_bdry = 0; i_bdry < nve_bdry/*_nc*/; i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
			
			sol_ctrl_bdry_gss +=  sol_ctrl[i_vol] * phi_ctrl_bdry[i_bdry];
                            for (int d = 0; d < dim; d++) {
			      sol_ctrl_x_bdry_gss[d] += sol_ctrl[i_vol] * phi_ctrl_x_bdry[i_bdry + d*nve_bdry];
			    }
		      }  

                 //========= compute gauss quantities on the boundary ================================================
                  integral_alpha +=  weight * sol_ctrl_bdry_gss * sol_ctrl_bdry_gss; 
                  integral_beta  +=  weight * (sol_ctrl_x_bdry_gss[0] * sol_ctrl_x_bdry_gss[0] /*+ sol_ctrl_x_bdry_gss[1] * sol_ctrl_x_bdry_gss[1]*/);
                 
		}
	      } //end face == 3
	      
	    } //end if boundary faces
	  }  // loop over element faces   
	  
	} //end if control element flag

//=====================================================================================================================  
//=====================================================================================================================  
//=====================================================================================================================  
  
  
   
      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][solType_max]->GetGaussPointNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
        //  ==== state 
	msh->_finiteElement[kelGeom][solType_u]   ->Jacobian(x, ig, weight, phi_u, phi_u_x, phi_u_xx);
        //  ==== adjoint
        msh->_finiteElement[kelGeom][solType_u/*solTypeTdes*/]->Jacobian(x, ig, weight, phi_udes, phi_udes_x, phi_udes_xx);

	u_gss = 0.;  for (unsigned i = 0; i < nDof_u; i++) u_gss += sol_u[i] * phi_u[i];		
	udes_gss  = 0.; for (unsigned i = 0; i < nDof_udes; i++)  udes_gss  += sol_udes[i]  * phi_udes[i];  

               integral_target += target_flag * weight * (u_gss  - udes_gss) * (u_gss - udes_gss);
	  
      } // end gauss point loop
      
  } //end element loop

  double total_integral = 0.5 * integral_target + 0.5 * alpha * integral_alpha + 0.5 * beta * integral_beta;
  
  std::cout << "The value of the integral_target is " << std::setw(11) << std::setprecision(10) << integral_target << std::endl;
  std::cout << "The value of the integral_alpha  is " << std::setw(11) << std::setprecision(10) << integral_alpha << std::endl;
  std::cout << "The value of the integral_beta   is " << std::setw(11) << std::setprecision(10) << integral_beta << std::endl;
  std::cout << "The value of the total integral  is " << std::setw(11) << std::setprecision(10) << total_integral << std::endl;
 
return;
  
}
  
  

