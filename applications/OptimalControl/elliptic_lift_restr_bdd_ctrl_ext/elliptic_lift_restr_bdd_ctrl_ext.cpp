#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "NumericVector.hpp"

#include "../elliptic_param.hpp"

using namespace femus;

double InitialValueContReg(const std::vector < double >& x) {
  return ControlDomainFlag_internal_restriction(x);
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
  value = 0.;
  
  if(!strcmp(name,"control")) {
      value = 0.;
  if (faceName == 3)
    dirichlet = false;
  
  }
  
  if(!strcmp(name,"mu")) {
      value = 0.;
  if (faceName == 3)
    dirichlet = false;
  }
  
  return dirichlet;
}


double ComputeIntegral(MultiLevelProblem& ml_prob);

void AssembleLiftExternalProblem(MultiLevelProblem& ml_prob);


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
    // ======= Files ========================
  Files files; 
        files.CheckIODirectories();
	files.RedirectCout();

  double scalingFactor = 1.;
  
  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  mlMsh.ReadCoarseMesh("./input/ext_box.neu", "seventh", scalingFactor);


  //mlMsh.GenerateCoarseBoxMesh(NSUB_X,NSUB_Y,0,0.,1.,0.,1.,0.,0.,QUAD9,"seventh");
 /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfUniformLevels = 2;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  mlMsh.PrintInfo();
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("state", LAGRANGE, FIRST);
  mlSol.AddSolution("control", LAGRANGE, FIRST);
  mlSol.AddSolution("adjoint", LAGRANGE, FIRST);
  mlSol.AddSolution("mu", LAGRANGE, FIRST);  
  mlSol.AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  mlSol.AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field

  
  mlSol.Initialize("All");    // initialize all varaibles to zero

  mlSol.Initialize("state", InitialValueState);
  mlSol.Initialize("control", InitialValueControl);
  mlSol.Initialize("adjoint", InitialValueAdjoint);
  mlSol.Initialize("mu", InitialValueMu);
  mlSol.Initialize("TargReg", InitialValueTargReg);
  mlSol.Initialize("ContReg", InitialValueContReg);

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
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("LiftRestr");

  system.AddSolutionToSystemPDE("state");  
  system.AddSolutionToSystemPDE("control");  
  system.AddSolutionToSystemPDE("adjoint");  
  system.AddSolutionToSystemPDE("mu");  
  
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleLiftExternalProblem);
  
  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->SetDebugOutput(true);
  
  system.SetDebugNonlinear(true);
//   system.SetMaxNumberOfNonLinearIterations(4);

  // initilaize and solve the system
  system.init();
  system.MGsolve();
  
  ComputeIntegral(mlProb);
 
  // print solutions
  std::vector < std::string > variablesToBePrinted;
     variablesToBePrinted.push_back("all");

    // ******* Print solution *******
  mlSol.GetWriter()->Write(files.GetOutputPath()/*DEFAULT_OUTPUTDIR*/, "biquadratic", variablesToBePrinted);

  return 0;
}


void AssembleLiftExternalProblem(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  NonLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystem> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object

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
  vector < vector < double > > x(dim);    // local coordinates
  vector < vector < double > > x_bdry(dim);    // local coordinates
  for (unsigned idim = 0; idim < dim; idim++) {
    x[idim].reserve(maxSize);
    x_bdry[idim].reserve(maxSize);
  }
 //***************************************************   

 //***************************************************  
  double weight; // gauss point weight
  double weight_bdry = 0.; // gauss point weight on the boundary

 //********************* state *********************** 
 //***************************************************  
  vector <double> phi_u;  // local test function
  vector <double> phi_u_x; // local test function first order partial derivatives
  vector <double> phi_u_xx; // local test function second order partial derivatives

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * dim);
  phi_u_xx.reserve(maxSize * dim2);
  
 
  unsigned solIndex_u;
  solIndex_u = mlSol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u = mlSol->GetSolutionType(solIndex_u);    // get the finite element type for "state"

  unsigned solPdeIndex_u;
  solPdeIndex_u = mlPdeSys->GetSolPdeIndex("state");    // get the position of "state" in the pdeSys object

  vector < double >  sol_u; // local solution
  sol_u.reserve(maxSize);
  vector< int > l2GMap_u;
  l2GMap_u.reserve(maxSize);
 //***************************************************  
 //***************************************************  

  
 //******************** control ********************** 
 //***************************************************   
  vector <double> phi_ctrl;  // local test function
  vector <double> phi_ctrl_x; // local test function first order partial derivatives
  vector <double> phi_ctrl_xx; // local test function second order partial derivatives

  phi_ctrl.reserve(maxSize);
  phi_ctrl_x.reserve(maxSize * dim);
  phi_ctrl_xx.reserve(maxSize * dim2);
  
  unsigned solIndex_ctrl;
  solIndex_ctrl = mlSol->GetIndex("control");
  unsigned solType_ctrl = mlSol->GetSolutionType(solIndex_ctrl);

  unsigned solPdeIndex_ctrl;
  solPdeIndex_ctrl = mlPdeSys->GetSolPdeIndex("control");

  vector < double >  sol_ctrl; // local solution
  sol_ctrl.reserve(maxSize);
  vector< int > l2GMap_ctrl;
  l2GMap_ctrl.reserve(maxSize);
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
 
  
  unsigned solIndex_adj;
  solIndex_adj = mlSol->GetIndex("adjoint");    // get the position of "adjoint" in the ml_sol object
  unsigned solType_adj = mlSol->GetSolutionType(solIndex_adj);    // get the finite element type for "adjoint"

  unsigned solPdeIndex_adj;
  solPdeIndex_adj = mlPdeSys->GetSolPdeIndex("adjoint");    // get the position of "adjoint" in the pdeSys object

  vector < double >  sol_adj; // local solution
    sol_adj.reserve(maxSize);
  vector< int > l2GMap_adj;
    l2GMap_adj.reserve(maxSize);
 //***************************************************  
 //***************************************************  
 //boundary adjoint shape functions  
  vector <double> phi_adj_bdry;  
  vector <double> phi_adj_x_bdry; 

  phi_adj_bdry.reserve(maxSize);
  phi_adj_x_bdry.reserve(maxSize * dim);
  
  vector <double> phi_adj_vol_at_bdry;  // local test function
  vector <double> phi_adj_x_vol_at_bdry; // local test function first order partial derivatives
  phi_adj_vol_at_bdry.reserve(maxSize);
  phi_adj_x_vol_at_bdry.reserve(maxSize * dim);
  vector <double> sol_adj_x_vol_at_bdry_gss;
  sol_adj_x_vol_at_bdry_gss.reserve(dim);
 //*************************************************** 
 //*************************************************** 

  
 //********************* bdry cont *******************
 //*************************************************** 
  vector <double> phi_ctrl_bdry;  
  vector <double> phi_ctrl_x_bdry; 

  phi_ctrl_bdry.reserve(maxSize);
  phi_ctrl_x_bdry.reserve(maxSize * dim);
  
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

  //********* variables for ineq constraints *****************
  double ctrl_lower = CTRL_BOX_LOWER;
  double ctrl_upper = CTRL_BOX_UPPER;
  assert(ctrl_lower < ctrl_upper);
  double c_compl = 1.;
  vector < double/*int*/ >  sol_actflag;   sol_actflag.reserve(maxSize); //flag for active set
  //***************************************************  

 //***************************************************  
 //********* WHOLE SET OF VARIABLES ****************** 
  const int solType_max = 2;  //biquadratic

  const int n_unknowns = 4;
 
  vector< int > l2GMap_AllVars; // local to global mapping
  l2GMap_AllVars.reserve(n_unknowns*maxSize);
  
  vector< double > Res; // local redidual vector
  Res.reserve(n_unknowns*maxSize);

  vector < double > Jac;
  Jac.reserve( n_unknowns*maxSize * n_unknowns*maxSize);
  

  vector < std::string > Solname(n_unknowns);
  Solname[0] = "state";
  Solname[1] = "control";
  Solname[2] = "adjoint";
  Solname[3] = "mu";
  
  vector < unsigned > SolPdeIndex(n_unknowns);
  vector < unsigned > SolIndex(n_unknowns);  
  vector < unsigned > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolPdeIndex[ivar]	= mlPdeSys->GetSolPdeIndex(Solname[ivar].c_str());
    SolIndex[ivar]	= mlSol->GetIndex        (Solname[ivar].c_str());
    SolFEType[ivar]	= mlSol->GetSolutionType(SolIndex[ivar]);
  }

    vector < unsigned > Sol_n_el_dofs(n_unknowns);
 //***************************************************  

  
 //********************* DATA ************************ 
  double u_des = DesiredTarget();
  double alpha = ALPHA_CTRL_VOL;
  double beta  = BETA_CTRL_VOL;
  double penalty_strong = 10e+14;
  double penalty_interface = 1.e+10;         //penalty for u=q
 //***************************************************  

  RES->zero();
  if (assembleMatrix)  KK->zero();

    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    
    int group_flag         = msh->GetElementGroup(iel);
    short unsigned kelGeom = msh->GetElementType(iel);    // element geometry type
    std::cout << " ======= grp_flag === " << group_flag << " ================== " << std::endl; 
//     int face_no         = msh->GetElementFaceNumber(iel);
//     std::cout << " ======= face# === " << face_no << " ================== " << std::endl; 
 //******************** GEOMETRY ********************* 
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    for (int i = 0; i < dim; i++)  x[i].resize(nDofx);
    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);  // global to global mapping between coordinates node and coordinate dof

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
  
 //****** set target domain flag ********************* 
   int target_flag = 0;
   target_flag = ElementTargetFlag(elem_center);
 //*************************************************** 
   
    
 //**************** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);
    sol_u    .resize(nDof_u);
    l2GMap_u.resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);
      l2GMap_u[i] = pdeSys->GetSystemDof(solIndex_u, solPdeIndex_u, i, iel);
    }
 //***************************************************  
 
 //***************** control ************************* 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);    // number of solution element dofs
    sol_ctrl    .resize(nDof_ctrl);
    l2GMap_ctrl.resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);    // global to global mapping between solution node and solution dof
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);      // global extraction and local storage for the solution
      l2GMap_ctrl[i] = pdeSys->GetSystemDof(solIndex_ctrl, solPdeIndex_ctrl, i, iel);    // global to global mapping between solution node and pdeSys dof
    } 
 //*************************************************** 
 

 //************** adjoint **************************** 
    unsigned nDof_adj  = msh->GetElementDofNumber(iel, solType_adj);    // number of solution element dofs
        sol_adj    .resize(nDof_adj);
        l2GMap_adj.resize(nDof_adj);
    for (unsigned i = 0; i < sol_adj.size(); i++) {
      unsigned solDof_adj = msh->GetSolutionDof(i, iel, solType_adj);   // global to global mapping between solution node and solution dof
      sol_adj[i] = (*sol->_Sol[solIndex_adj])(solDof_adj);      // global extraction and local storage for the solution
      l2GMap_adj[i] = pdeSys->GetSystemDof(solIndex_adj, solPdeIndex_adj, i, iel);   // global to global mapping between solution node and pdeSys dof
    } 
 //***************************************************  
 
 //************** mu **************************** 
    unsigned nDof_mu  = msh->GetElementDofNumber(iel, solType_mu);    // number of solution element dofs
    sol_mu   .resize(nDof_mu);
    l2GMap_mu.resize(nDof_mu);
    for (unsigned i = 0; i < sol_mu.size(); i++) {
      unsigned solDof_mu = msh->GetSolutionDof(i, iel, solType_mu);   // global to global mapping between solution node and solution dof
      sol_mu[i] = (*sol->_Sol[solIndex_mu])(solDof_mu);      // global extraction and local storage for the solution 
      l2GMap_mu[i] = pdeSys->GetSystemDof(solIndex_mu, solPdeIndex_mu, i, iel);   // global to global mapping between solution node and pdeSys dof
    }
    
    
 //************** update active set flag for current nonlinear iteration **************************** 
 // 0: inactive; 1: active_a; 2: active_b
   assert(nDof_mu == nDof_ctrl);
   sol_actflag.resize(nDof_mu);
     std::fill(sol_actflag.begin(), sol_actflag.end(), 0);
   
    for (unsigned i = 0; i < sol_actflag.size(); i++) {  
    if      ( (sol_mu[i] + c_compl * (sol_ctrl[i] - ctrl_lower )) < 0 )  sol_actflag[i] = 1;
    else if ( (sol_mu[i] + c_compl * (sol_ctrl[i] - ctrl_upper )) > 0 )  sol_actflag[i] = 2;
    }
 
 //******************** ALL VARS ********************* 
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
    
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
	Sol_n_el_dofs[k] = ndofs_unk;
  }
  
   vector <double> interface_flag;   interface_flag.reserve(nDof_u); //flag for boundary interface
   std::fill(interface_flag.begin(), interface_flag.end(), 0.);
   
  //************ Boundary loops *************************************** 
  	  double tau=0.;
	  vector<double> normal(dim,0);

  	  for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {

 //========= compute coordinates of boundary nodes on each element ========================================== 
		unsigned nDofx_bdry    = msh->GetElementFaceDofNumber(iel,jface,xType);
		unsigned nDofu_bdry    = msh->GetElementFaceDofNumber(iel,jface,solType_u);
		unsigned nDofctrl_bdry = msh->GetElementFaceDofNumber(iel,jface,solType_ctrl);
		if (nDofu_bdry != nDofctrl_bdry) { std::cout << "State and control need to have the same FE space" << std::endl; abort(); }
	        for (unsigned idim = 0; idim < dim; idim++) {  x_bdry[idim].resize(nDofx_bdry); }
		const unsigned felt_bdry = msh->GetElementFaceType(iel, jface);    
		for(unsigned i=0; i < nDofx_bdry; i++) {
		  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i);
                  unsigned iDof = msh->GetSolutionDof(i_vol, iel, xType);
		  for(unsigned idim=0; idim<dim; idim++) {
		      x_bdry[idim][i]=(*msh->_topology->_Sol[idim])(iDof);
		  }
		}
 //===================================================   
 
 //============== face elem average point =====================================   
    vector < double > elem_center_bdry(dim);   
    for (unsigned j = 0; j < dim; j++) {  elem_center_bdry[j] = 0.;  }
    for (unsigned j = 0; j < dim; j++) {  
      for (unsigned i = 0; i < nDofx_bdry; i++) {
         elem_center_bdry[j] += x_bdry[j][i];
       }
    }
   for (unsigned j = 0; j < dim; j++) { elem_center_bdry[j] = elem_center_bdry[j]/nDofx_bdry; }
 //===================================================   
    
    
 //============ find interface nodes (now we do with coordinates, later we can do also with flag in gambit) ======================================= 
 double my_eps = 1.e-6;
    if (elem_center_bdry[0] > 1. - my_eps   && elem_center_bdry[0] < 1.   + my_eps  && 
        elem_center_bdry[1] > 0.25 - my_eps && elem_center_bdry[1] < 0.75 + my_eps
    ) {
       std::cout << " bdry_interface " <<"("<<elem_center_bdry[0]<<","<< elem_center_bdry[1]<<")"<< std::endl; 
      		      for (int i_bdry = 0; i_bdry < nDofu_bdry; i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
		    interface_flag[i_vol] = 1.;
		    
//========= initialize gauss quantities on the boundary ============================================
		
		for(unsigned ig_bdry=0; ig_bdry < msh->_finiteElement[felt_bdry][solType_ctrl]->GetGaussPointNumber(); ig_bdry++) {
		  
		  msh->_finiteElement[felt_bdry][solType_ctrl]->JacobianSur(x_bdry,ig_bdry,weight_bdry,phi_ctrl_bdry,phi_ctrl_x_bdry,normal);
		  msh->_finiteElement[felt_bdry][solType_adj]->JacobianSur(x_bdry,ig_bdry,weight_bdry,phi_adj_bdry,phi_adj_x_bdry,normal);
		  msh->_finiteElement[kelGeom][solType_adj]->ShapeAtBoundary(x,ig_bdry,phi_adj_vol_at_bdry,phi_adj_x_vol_at_bdry);

//========== temporary soln for surface gradient on a face parallel to the X axis ===================
		  double dx_dxi = 0.;
		 const elem_type_1D * myeltype = static_cast<const elem_type_1D*>(msh->_finiteElement[felt_bdry][solType_ctrl]);
		 const double * myptr = myeltype->GetDPhiDXi(ig_bdry);
		      for (int inode = 0; inode < nDofu_bdry/*_nc*/; inode++) dx_dxi += myptr[inode] * x_bdry[0][inode];
  
		      for (int inode = 0; inode < nDofu_bdry/*_nc*/; inode++) {
                            for (int d = 0; d < dim; d++) {
                              if (d==0 ) phi_ctrl_x_bdry[inode + d*nDofu_bdry/*_nc*/] = myptr[inode]* (1./ dx_dxi);
                              else  phi_ctrl_x_bdry[inode + d*nDofu_bdry/*_nc*/] = 0.;
                         }
                     }
//========== temporary soln for surface gradient on a face parallel to the X axis ===================
		  

//=============== grad dot n for residual ========================================= 
//     compute gauss quantities on the boundary through VOLUME interpolation
           std::fill(sol_adj_x_vol_at_bdry_gss.begin(), sol_adj_x_vol_at_bdry_gss.end(), 0.);
		      for (int iv = 0; iv < nDof_adj; iv++)  {
			
                            for (int d = 0; d < dim; d++) {
//    std::cout << " ivol " << iv << std::endl;
//    std::cout << " adj dofs " << sol_adj[iv] << std::endl;
			      sol_adj_x_vol_at_bdry_gss[d] += sol_adj[iv] * phi_adj_x_vol_at_bdry[iv * dim + d];//notice that the convention of the orders x y z is different from vol to bdry
			    }
		      }  
		      
    double grad_dot_n_adj_res = 0.;
        for(unsigned d=0; d<dim; d++) {
	  grad_dot_n_adj_res += sol_adj_x_vol_at_bdry_gss[d]*normal[d];  
	}
//=============== grad dot n  for residual =========================================      
		    
//============ Bdry Residuals ==================	
                if (i_vol < nDof_u)     Res[ (0 + i_vol) ]                    +=  -  penalty_interface * ( sol_u[i_vol] - sol_ctrl[i_vol] );   // u = q
		
		if (i_vol < nDof_ctrl)  Res[ (nDof_u + i_vol) ]               +=  -  weight_bdry * ( - grad_dot_n_adj_res * phi_ctrl_bdry[i_bdry] );  //boundary optimality condition
//============ Bdry Residuals ==================	
		    for(unsigned j_bdry=0; j_bdry < nDofu_bdry; j_bdry ++) {
		         unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, jface, j_bdry);
//============ Bdry Jacobians ==================	


// FIRST BLOCK ROW
//============ u = q ===========================	    
// block delta_state/state =====================
		if (i_vol < nDof_u && j_vol < nDof_u && i_vol == j_vol)  {
		  Jac[    
		(0 + i_vol) * nDof_AllVars  +
		(0 + j_vol)                                ]  += penalty_interface * ( 1.);
		  
		         }

// block delta_state/control ===================
	      if ( i_vol < nDof_u && j_vol < nDof_ctrl && i_vol == j_vol) {
		Jac[    
		(0     + i_vol) * nDof_AllVars  +
		(nDof_u + j_vol)                           ]  += penalty_interface  * (-1.);
	
	                 }
//============ u = q ===========================		    
		    
		    } //end j_vol 
		    
//===================loop over j in the VOLUME (while i is in the boundary)	      
	for(unsigned j=0; j < nDof_max; j ++) {
  
//=============== grad dot n  =========================================    
    double grad_adj_dot_n_mat = 0.;
        for(unsigned d=0; d<dim; d++) {
	  grad_adj_dot_n_mat += phi_adj_x_vol_at_bdry[j * dim + d]*normal[d];  //notice that the convention of the orders x y z is different from vol to bdry
	}
//=============== grad dot n  =========================================    

  std::cout << " gradadjdotn " << grad_adj_dot_n_mat << std::endl;
  
		      
//==========block delta_control/adjoint ========
		   if ( i_vol < nDof_ctrl    && j < nDof_adj)   
		     Jac[ 
			(nDof_u + i_vol) * nDof_AllVars  +
		        (nDof_u + nDof_ctrl + j)             ]  += /*control_node_flag[i_vol] **/ (-1) *
								  ( weight_bdry * grad_adj_dot_n_mat * phi_ctrl_bdry[i_bdry] );    		      
		}   //end loop i_bdry // j_vol
	      }  //end ig_bdry loop
	    }
      
        }
    }  //end face loop
	  
  //************ Boundary loops *************************************** 
	    
 //*************************************************** 
    
//  //***** set control flag ****************************
//   int control_el_flag = 0;
//   control_el_flag = ControlDomainFlag(elem_center);
//   std::vector<int> control_node_flag(nDof_ctrl,0);
//   if (control_el_flag == 1) std::fill(control_node_flag.begin(), control_node_flag.end(), 1);
//  //*************************************************** 
//   
//  //***** set state flag ****************************
//   int state_el_flag = 0;
//   state_el_flag = StateDomainFlag(elem_center);
//   std::vector<int> state_node_flag(nDof_ctrl,0);
//   if (state_el_flag == 1) std::fill(state_node_flag.begin(), state_node_flag.end(), 1);
//  //*************************************************** 
  
 //========= gauss value quantities ==================   
	double sol_u_gss = 0.;
	double sol_adj_gss = 0.;
	double sol_ctrl_gss = 0.;
	std::vector<double> sol_u_x_gss(dim);       std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::vector<double> sol_adj_x_gss(dim);     std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	std::vector<double> sol_ctrl_x_gss(dim);    std::fill(sol_ctrl_x_gss.begin(), sol_ctrl_x_gss.end(), 0.);
 //===================================================   

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][solType_max]->GetGaussPointNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
	msh->_finiteElement[kelGeom][solType_u]   ->Jacobian(x, ig, weight, phi_u, phi_u_x, phi_u_xx);
        msh->_finiteElement[kelGeom][solType_ctrl]->Jacobian(x, ig, weight, phi_ctrl, phi_ctrl_x, phi_ctrl_xx);
        msh->_finiteElement[kelGeom][solType_adj] ->Jacobian(x, ig, weight, phi_adj, phi_adj_x, phi_adj_xx);
	
	std::fill(sol_u_x_gss.begin(),sol_u_x_gss.end(), 0.);
	std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	std::fill(sol_ctrl_x_gss.begin(), sol_ctrl_x_gss.end(), 0.);
	
	for (unsigned i = 0; i < nDof_u; i++) {
	                                                sol_u_gss      += sol_u[i] * phi_u[i];
                   for (unsigned d = 0; d < dim; d++)   sol_u_x_gss[d] += sol_u[i] * phi_u_x[i * dim + d];
          }
	
	for (unsigned i = 0; i < nDof_adj; i++) {
	                                                sol_adj_gss      += sol_adj[i] * phi_adj[i];
                   for (unsigned d = 0; d < dim; d++)   sol_adj_x_gss[d] += sol_adj[i] * phi_adj_x[i * dim + d];
        }
	
	for (unsigned i = 0; i < nDof_ctrl; i++) {
	                                                sol_ctrl_gss      += sol_ctrl[i] * phi_ctrl[i];
                   for (unsigned d = 0; d < dim; d++)   sol_ctrl_x_gss[d] += sol_ctrl[i] * phi_ctrl_x[i * dim + d];
        }
        
//==========FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {
	  
	      double laplace_rhs_du_adj_i = 0.; //
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_u )         laplace_rhs_du_adj_i             +=  (phi_u_x   [i * dim + kdim] * sol_adj_x_gss[kdim]);
	      }
	      
              double laplace_rhs_dctrl_ctrl_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_ctrl )         laplace_rhs_dctrl_ctrl_i      +=  (phi_ctrl_x   [i * dim + kdim] * sol_ctrl_x_gss[kdim]);
	      }
	      
	      double laplace_rhs_dctrl_adj_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_ctrl )         laplace_rhs_dctrl_adj_i       +=  (phi_ctrl_x   [i * dim + kdim] * sol_adj_x_gss[kdim]);
	      }
	      
	      double laplace_rhs_dadj_u_i = 0.;  //
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_adj )         laplace_rhs_dadj_u_i           +=  (phi_adj_x   [i * dim + kdim] * sol_u_x_gss[kdim]);
	      }
	      
	      double laplace_rhs_dadj_ctrl_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_adj )         laplace_rhs_dadj_ctrl_i        +=  (phi_adj_x   [i * dim + kdim] * sol_ctrl_x_gss[kdim]);
	      }
//======================Residuals=======================
          // FIRST ROW
	  if (i < nDof_u)  {
	     if ( group_flag == 12 )            Res[0      + i] += - weight * (target_flag * phi_u[i] * ( sol_u_gss /*+ sol_ctrl_gss*/ - u_des) - laplace_rhs_du_adj_i - 0.);
	  
	     else if ( group_flag == 13 )       Res[0      + i] +=  (1-interface_flag[i]) * (- penalty_strong) * (sol_u[i] - 0.);
	  }
          // SECOND ROW
	  if (i < nDof_ctrl)  {
	     if ( group_flag == 13 )            Res[nDof_u + i] +=  /*(control_node_flag[i]) **/ - weight * (/*target_flag * phi_ctrl[i] * ( sol_u_gss + sol_ctrl_gss - u_des) */
													      + alpha * phi_ctrl[i] * sol_ctrl_gss
		                                                                                              - laplace_rhs_dctrl_adj_i 
		                                                                                              + beta * laplace_rhs_dctrl_ctrl_i
													      /*+ 1. * sol_mu[i]*/ - 0.);
	     else if ( group_flag == 12 )       Res[nDof_u + i] +=  (1-interface_flag[i]) * (- penalty_strong) * (sol_ctrl[i] - 0.);
	  }
          // THIRD ROW
          if (i < nDof_adj) {  
	     if ( group_flag == 12 )      Res[nDof_u + nDof_ctrl + i] += - weight *  ( - laplace_rhs_dadj_u_i /*- laplace_rhs_dadj_ctrl_i*/ - 0.) ;
	     
	     else if ( group_flag == 13 ) Res[nDof_u + nDof_ctrl + i] += - weight *  (- laplace_rhs_dadj_ctrl_i - 0.) ;
	  }
//======================Residuals=======================
	      
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {
              //double laplace_mat_du_u = 0.;
              double laplace_mat_du_adj = 0.;
              double laplace_mat_dadj_u = 0.;
              double laplace_mat_dctrl_adj = 0.;
              double laplace_mat_dadj_ctrl = 0.;
              double laplace_mat_dctrl_ctrl = 0.;

              for (unsigned kdim = 0; kdim < dim; kdim++) {
              //if ( i < nDof_u && j < nDof_u )           laplace_mat_du_u           += (phi_u_x   [i * dim + kdim] * phi_u_x   [j * dim + kdim]);
              if ( i < nDof_u && j < nDof_adj )         laplace_mat_du_adj         += (phi_u_x   [i * dim + kdim] * phi_adj_x [j * dim + kdim]);
              if ( i < nDof_adj && j < nDof_u )         laplace_mat_dadj_u         += (phi_adj_x [i * dim + kdim] * phi_u_x   [j * dim + kdim]);  //equal to the previous
              if ( i < nDof_ctrl && j < nDof_adj )      laplace_mat_dctrl_adj      += (phi_ctrl_x[i * dim + kdim] * phi_adj_x [j * dim + kdim]);
              if ( i < nDof_adj && j < nDof_ctrl )      laplace_mat_dadj_ctrl      += (phi_adj_x [i * dim + kdim] * phi_ctrl_x[j * dim + kdim]);  //equal to the previous
              if ( i < nDof_ctrl && j < nDof_ctrl )     laplace_mat_dctrl_ctrl     += (phi_ctrl_x  [i * dim + kdim] * phi_ctrl_x  [j * dim + kdim]);
	      }

              //============ delta_state row ============================
              if ( group_flag == 12 ) { 
		
              //DIAG BLOCK delta_state - state
	      if ( i < nDof_u && j < nDof_u )       
		Jac[ (0 + i) * nDof_AllVars   +
		     (0 + j)                            ]  += weight * target_flag * phi_u[j] *  phi_u[i];
              
// 	      // BLOCK  delta_state - control
//               if ( i < nDof_u && j < nDof_ctrl )   
// 		Jac[ (0 + i) * nDof_AllVars   +
//                      (nDof_u + j)                       ]  += weight * target_flag  * phi_ctrl[j] *  phi_u[i];
	      
              // BLOCK  delta_state - adjoint
              if ( i < nDof_u && j < nDof_adj )  
		Jac[  (0 + i) * nDof_AllVars  +
                      (nDof_u + nDof_ctrl + j)          ]  += weight * (-1) * laplace_mat_du_adj;
	      }
		      
	      else if ( group_flag == 13 ) {  
		
              //BLOCK delta_state - state
              if ( i < nDof_u   && j < nDof_u  &&  i==j ) {
		 Jac[ (0 + i) * nDof_AllVars +
		      (0 + j)                           ]  += (1-interface_flag[i]) * 1. * penalty_strong;
	        }
	        
	      }
              
	      //=========== delta_control row ===========================     
	      if ( group_flag == 13 )  {

// 	      //BLOCK delta_control - state
//               if ( i < nDof_ctrl   && j < nDof_u   ) 
// 		Jac[ (nDof_u + i) * nDof_AllVars  +
// 		     (0 + j)                            ]  += ( control_node_flag[i]) * weight * target_flag * phi_u[j] * phi_ctrl[i];
		
	      //BLOCK delta_control - control
              if ( i < nDof_ctrl   && j < nDof_ctrl   )
		Jac[ (nDof_u + i) * nDof_AllVars +
		     (nDof_u + j)                       ]  += /*( control_node_flag[i]) **/ weight * ( beta /** control_el_flag*/  * laplace_mat_dctrl_ctrl 
		                                           + alpha /** control_el_flag*/ * phi_ctrl[i] * phi_ctrl[j] /*+ target_flag  * phi_ctrl[i] * phi_ctrl[j]*/ );
              
	      //BLOCK delta_control - adjoint
              if ( i < nDof_ctrl   && j < nDof_adj  ) 
		Jac[ (nDof_u + i) * nDof_AllVars  + 
		     (nDof_u + nDof_ctrl + j)           ]  += /*( control_node_flag[i]) **/ weight * (-1) * laplace_mat_dctrl_adj;
	      	      
	        }
	      
	      else if ( group_flag == 12 )  {  
		
              //BLOCK delta_control - control
              if ( i < nDof_ctrl   && j < nDof_ctrl &&  i==j ) {
		 Jac[ (nDof_u + i) * nDof_AllVars +
		      (nDof_u + j)                      ]  += (1-interface_flag[i]) * 1. * penalty_strong;
		}
	      
	      }
	      
// 	      //BLOCK delta_control - mu
//            if ( i < nDof_ctrl   && j < nDof_mu && i==j ) 
// 		Jac[ (nDof_u + i) * nDof_AllVars  + 
// 		     (nDof_u + nDof_ctrl + nDof_adj + j)]   =  /*control_node_flag[i] **/ 1.;
		     
		     
	      //=========== delta_adjoint row ===========================
	      if ( group_flag == 12 ){
		
              // BLOCK delta_adjoint - state	      
              if ( i < nDof_adj && j < nDof_u )   
		Jac[ (nDof_u + nDof_ctrl + i) * nDof_AllVars +
		     (0 + j)                            ]  += weight * (-1) * laplace_mat_dadj_u;   
	      }
	      
	      else if ( group_flag == 13 ){
		
              // BLOCK delta_adjoint - control   
              if ( i < nDof_adj && j < nDof_ctrl )  
		Jac[ (nDof_u + nDof_ctrl + i)  * nDof_AllVars +
		     (nDof_u  + j)                      ]  += weight * (-1) * laplace_mat_dadj_ctrl; 
	      }
		     
// 	      // BLOCK delta_adjoint - adjoint   
//               if ( i < nDof_adj && j < nDof_adj )  
// 		Jac[ (nDof_u + nDof_ctrl + i)  * nDof_AllVars +
// 		     (nDof_u + nDof_ctrl + j)                      ]  += weight * phi_adj[j] *  phi_adj[i]; 
    
	      
	      //============= delta_mu row ===============================
//	      if (sol_actflag[i] == 0) //inactive
//	      { // BLOCK delta_mu - mu	      
// 	        if ( i < nDof_mu && j < nDof_mu && i==j )   
// 		  Jac[ (nDof_u + nDof_ctrl + nDof_adj + i) * nDof_AllVars +
// 		       (nDof_u + nDof_ctrl + nDof_adj + j)]  = 1. ;  
// 	     // }
// 	      else //active
// 	      { // BLOCK delta_mu - ctrl	      
//                 if ( i < nDof_mu && j < nDof_ctrl && i==j )   
// 		  Jac[ (nDof_u + nDof_ctrl + nDof_adj + i) * nDof_AllVars +
// 		       (nDof_u + j)                       ]  = c_compl * 1. ; 
	     // }
	      
            } // end phi_j loop
          } // endif assemble_matrix

        } // end phi_i loop
        
      } // end gauss point loop

      
      
    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
   // std::cout << " ************* Element ************** " << iel << " **************************************** " << std::endl;     

// // //     if (control_el_flag == 0) {  //elements that should have zero control
// //          for (unsigned i_unk = 0; i_unk < n_unknowns; i_unk++) {
// //     std::cout << " ======= Row === " << i_unk << " =================================================== " << std::endl;     
// //         unsigned int row_block_offset = 0;
// // 	         for (unsigned k = 0; k < i_unk; k++) row_block_offset += Sol_n_el_dofs[k];
// //          for (unsigned j_unk = 0; j_unk < n_unknowns; j_unk++) {
// //     std::cout << " ======= Column === " << j_unk << " ================== " << std::endl;     
// //         unsigned int column_block_offset = 0;
// // 	         for (unsigned k = 0; k < j_unk; k++) column_block_offset += Sol_n_el_dofs[k];
// // 	  
// //          for (unsigned i = 0; i < Sol_n_el_dofs[i_unk]; i++) {
// // // 	      std::cout << Res[nDof_u + nDof_ctrl + nDof_adj + i ] << " " << std::endl;
// // 	   for (unsigned j = 0; j < Sol_n_el_dofs[j_unk]; j++) {
// // 	      std::cout <<  " " << std::setfill(' ') << std::setw(10) << Jac[ (row_block_offset + i) * nDof_AllVars + ( column_block_offset + j) ] << " ";
// // 	    }
// // 	      std::cout << std::endl;
// // 	 }
// // 
// // 	 } //j_unk
// // 	} //i_unk
	 
	 
// // // 	}
    std::vector<double> Res_ctrl (nDof_ctrl); std::fill(Res_ctrl.begin(),Res_ctrl.end(), 0.);
    for (unsigned i = 0; i < sol_ctrl.size(); i++){
     if (  group_flag == 13 ){
	Res[nDof_u + i] = - ( - Res[nDof_u + i] + sol_mu[i] /*- ( 0.4 + sin(M_PI * x[0][i]) * sin(M_PI * x[1][i]) )*/ );
	Res_ctrl[i] = Res[nDof_u + i];
      }
    }
    
//     std::vector<double> Res_u (nDof_u); std::fill(Res_u.begin(),Res_u.end(), 0.);
//     for (unsigned i = 0; i < sol_u.size(); i++){
// 	Res[0 + i] = - ( sol_u[i] - 8. );
// 	Res_u[i] = Res[0 + i];
//     }
    
//     std::vector<double> Res_adj (nDof_adj); std::fill(Res_adj.begin(),Res_adj.end(), 0.);
//     for (unsigned i = 0; i < sol_adj.size(); i++){
// 	Res[nDof_u + nDof_ctrl + i] = - (sol_adj[i] - 7.);
// 	Res_adj[i] = Res[nDof_u + nDof_ctrl + i];
//     }
    

 //========== sum-based part

    //copy the value of the adept::adoube aRes in double Res and store
    RES->add_vector_blocked(Res, l2GMap_AllVars);
      if (assembleMatrix)  KK->add_matrix_blocked(Jac, l2GMap_AllVars, l2GMap_AllVars);
    
    
 //========== dof-based part, without summation
 
 //============= delta_mu row ===============================
      std::vector<double> Res_mu (nDof_mu); std::fill(Res_mu.begin(),Res_mu.end(), 0.);
    for (unsigned i = 0; i < sol_actflag.size(); i++){
      if (sol_actflag[i] == 0){  //inactive
         Res[nDof_u + nDof_ctrl + nDof_adj + i]  = - ( 1. * sol_mu[i] - 0. ); 
	 Res_mu [i] = Res[nDof_u + nDof_ctrl + nDof_adj + i]; 
      }
      else if (sol_actflag[i] == 1){  //active_a 
	 Res[nDof_u + nDof_ctrl + nDof_adj + i]  = - ( c_compl *  sol_ctrl[i] - c_compl * ctrl_lower);
         Res_mu [i] = Res[nDof_u + nDof_ctrl + nDof_adj + i] ;
      }
      else if (sol_actflag[i] == 2){  //active_b 
	Res[nDof_u + nDof_ctrl + nDof_adj + i]  =  - ( c_compl *  sol_ctrl[i] - c_compl * ctrl_upper);
	Res_mu [i] = Res[nDof_u + nDof_ctrl + nDof_adj + i] ;
      }
    }
//          Res[nDof_u + nDof_ctrl + nDof_adj + i]  = c_compl * (  (2 - sol_actflag[i]) * (ctrl_lower - sol_ctrl[i]) + ( sol_actflag[i] - 1 ) * (ctrl_upper - sol_ctrl[i])  ) ;
//          Res_mu [i] = Res[nDof_u + nDof_ctrl + nDof_adj + i] ;

    
    RES->insert(Res_mu, l2GMap_mu);
    RES->insert(Res_ctrl, l2GMap_ctrl);
//     RES->insert(Res_u, l2GMap_u);
//     RES->insert(Res_adj, l2GMap_adj);
    
//  //============= delta_state-delta_state row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_u, l2GMap_u, 1.);

//  //============= delta_ctrl-delta_ctrl row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_ctrl, l2GMap_ctrl, 1.);
 
//  //============= delta_adj-delta_adj row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_adj, l2GMap_adj, 1.);
  
 //============= delta_ctrl-delta_mu row ===============================
 KK->matrix_set_off_diagonal_values_blocked(l2GMap_ctrl, l2GMap_mu, 1.);
  
 //============= delta_mu-delta_ctrl row ===============================
 for (unsigned i = 0; i < sol_actflag.size(); i++) if (sol_actflag[i] != 0 ) sol_actflag[i] = c_compl;    
  
 KK->matrix_set_off_diagonal_values_blocked(l2GMap_mu, l2GMap_ctrl, sol_actflag);

 //============= delta_mu-delta_mu row ===============================
  for (unsigned i = 0; i < sol_actflag.size(); i++) sol_actflag[i] = 1 - sol_actflag[i]/c_compl;  //can do better to avoid division, maybe use modulo operator 

  KK->matrix_set_off_diagonal_values_blocked(l2GMap_mu, l2GMap_mu, sol_actflag);
  
  } //end element loop for each process
  
  RES->close();

  if (assembleMatrix) KK->close();
  //KK->print();
  //RES->print();
  
  // ***************** END ASSEMBLY *******************

  return;
}



double ComputeIntegral(MultiLevelProblem& ml_prob)    {
  
  
  NonLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystem> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

 //*************************************************** 
  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }
 //*************************************************** 

 //*************************************************** 
  double weight; // gauss point weight
  
 //***************************************************  
  double alpha = ALPHA_CTRL_VOL;
  double beta  = BETA_CTRL_VOL;

 //******************** state ************************ 
 //*************************************************** 
  vector <double> phi_u;  // local test function
  vector <double> phi_u_x; // local test function first order partial derivatives
  vector <double> phi_u_xx; // local test function second order partial derivatives

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * dim);
  phi_u_xx.reserve(maxSize * dim2);
  
 
  unsigned solIndex_u;
  solIndex_u = mlSol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u = mlSol->GetSolutionType(solIndex_u);    // get the finite element type for "state"

  vector < double >  sol_u; // local solution
  sol_u.reserve(maxSize);
  
  double u_gss = 0.;
 //*************************************************** 
 //*************************************************** 

 //******************** control ********************** 
 //*************************************************** 
  vector <double> phi_ctrl;  // local test function
  vector <double> phi_ctrl_x; // local test function first order partial derivatives
  vector <double> phi_ctrl_xx; // local test function second order partial derivatives

  phi_ctrl.reserve(maxSize);
  phi_ctrl_x.reserve(maxSize * dim);
  phi_ctrl_xx.reserve(maxSize * dim2);
  
  unsigned solIndex_ctrl;
  solIndex_ctrl = mlSol->GetIndex("control");
  unsigned solType_ctrl = mlSol->GetSolutionType(solIndex_ctrl);

  vector < double >  sol_ctrl; // local solution
  sol_ctrl.reserve(maxSize);
//   vector< int > l2GMap_ctrl;
//   l2GMap_ctrl.reserve(maxSize);
  
  double ctrl_gss = 0.;
  double ctrl_x_gss = 0.;
 //*************************************************** 
 //***************************************************  

  
 //********************* desired ********************* 
 //*************************************************** 
  vector <double> phi_udes;  // local test function
  vector <double> phi_udes_x; // local test function first order partial derivatives
  vector <double> phi_udes_xx; // local test function second order partial derivatives

    phi_udes.reserve(maxSize);
    phi_udes_x.reserve(maxSize * dim);
    phi_udes_xx.reserve(maxSize * dim2);
 
  
//  unsigned solIndex_udes;
//   solIndex_udes = mlSol->GetIndex("Tdes");    // get the position of "state" in the ml_sol object
//   unsigned solType_udes = mlSol->GetSolutionType(solIndex_udes);    // get the finite element type for "state"

  vector < double >  sol_udes; // local solution
  sol_udes.reserve(maxSize);
//   vector< int > l2GMap_udes;
//   l2GMap_udes.reserve(maxSize);
   double udes_gss = 0.;
 //*************************************************** 
 //*************************************************** 


 //*************************************************** 
 //********* WHOLE SET OF VARIABLES ****************** 
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
 //*************************************************** 
  
  double integral_target = 0.;
  double integral_alpha  = 0.;
  double integral_beta   = 0.;

    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned kelGeom = msh->GetElementType(iel);    // element geometry type
 
 //***************** GEOMETRY ************************ 
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    for (int i = 0; i < dim; i++)  x[i].resize(nDofx);
    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

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
  
 //************** set target domain flag *************
   int target_flag = 0;
   target_flag = ElementTargetFlag(elem_center);
 //*************************************************** 

   
 //**************** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);    // number of solution element dofs
        sol_u    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);  // global to global mapping between solution node and solution dof
            sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);            // global extraction and local storage for the solution
    }
 //*************************************************** 


 //************** control **************************** 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);    // number of solution element dofs
    sol_ctrl    .resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);    // global to global mapping between solution node and solution dof
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);      // global extraction and local storage for the solution
    } 
 //***************************************************  
 
 
 //**************** u_des **************************** 
    unsigned nDof_udes  = msh->GetElementDofNumber(iel, solType_u);    // number of solution element dofs
    sol_udes    .resize(nDof_udes);
    for (unsigned i = 0; i < sol_udes.size(); i++) {
      sol_udes[i] = u_des;  //dof value
    } 
 //*************************************************** 

 
 //******************* ALL VARS ********************** 
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
   
      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][solType_max]->GetGaussPointNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
        //  ==== State 
	msh->_finiteElement[kelGeom][solType_u]   ->Jacobian(x, ig, weight, phi_u, phi_u_x, phi_u_xx);
        //  ==== Adjoint 
        msh->_finiteElement[kelGeom][solType_u/*solTypeTdes*/]->Jacobian(x, ig, weight, phi_udes, phi_udes_x, phi_udes_xx);
        //  ==== Control 
        msh->_finiteElement[kelGeom][solType_ctrl]  ->Jacobian(x, ig, weight, phi_ctrl, phi_ctrl_x, phi_ctrl_xx);

	u_gss = 0.;  for (unsigned i = 0; i < nDof_u; i++) u_gss += sol_u[i] * phi_u[i];		
	ctrl_gss = 0.; for (unsigned i = 0; i < nDof_ctrl; i++) ctrl_gss += sol_ctrl[i] * phi_ctrl[i];  
	udes_gss  = 0.; for (unsigned i = 0; i < nDof_udes; i++)  udes_gss  += sol_udes[i]  * phi_udes[i]; 
        ctrl_x_gss  = 0.; for (unsigned i = 0; i < nDof_ctrl; i++)  
        {
          for (unsigned idim = 0; idim < dim; idim ++) ctrl_x_gss  += sol_ctrl[i] * phi_ctrl_x[i + idim * nDof_ctrl];
        }

               integral_target += target_flag * weight * (u_gss +  ctrl_gss - udes_gss) * (u_gss +  ctrl_gss - udes_gss);
               integral_alpha  += target_flag * alpha * weight * ctrl_gss * ctrl_gss;
               integral_beta   += target_flag * beta * weight * ctrl_x_gss * ctrl_x_gss;
	  
      } // end gauss point loop
  } //end element loop

  std::cout << "The value of the integral_target is " << std::setw(11) << std::setprecision(10) << integral_target << std::endl;
  std::cout << "The value of the integral_alpha  is " << std::setw(11) << std::setprecision(10) << integral_alpha << std::endl;
  std::cout << "The value of the integral_beta   is " << std::setw(11) << std::setprecision(10) << integral_beta << std::endl;
  std::cout << "The value of the total integral  is " << std::setw(11) << std::setprecision(10) << integral_target + integral_alpha + integral_beta << std::endl;
  
return integral_target + integral_alpha + integral_beta;
  
}
  
  

