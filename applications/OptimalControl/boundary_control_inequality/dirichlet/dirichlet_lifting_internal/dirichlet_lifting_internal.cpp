#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "NumericVector.hpp"
#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"

#define FACE_FOR_CONTROL   3

#include "../../param.hpp"

#define FE_DOMAIN  2

using namespace femus;



double Solution_set_initial_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;

    if(!strcmp(name,"state")) {
        value = 0.;
    }
    else if(!strcmp(name,"control")) {
        value = 0.;
    }
    else if(!strcmp(name,"adjoint")) {
        value = 0.;
    }
    else if(!strcmp(name,"mu")) {
        value = 0.;
    }
    else if(!strcmp(name,"TargReg")) {
        value = ElementTargetFlag(x);
    }
    else if(!strcmp(name,"ContReg")) {
        value = ControlDomainFlag_internal_restriction(x);
    }
    else if(!strcmp(name,"act_flag")) {
        value = 0.;
    }


    return value;
}



bool Solution_set_boundary_conditions(const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet = true; //dirichlet
  value = 0.;
  
  if(!strcmp(name,"control")) {
      value = 0.;
    if (faceName == FACE_FOR_CONTROL) {
        if (x[ axis_direction_Gamma_control(faceName) ] > GAMMA_CONTROL_LOWER - 1.e-5 && x[ axis_direction_Gamma_control(faceName) ] < GAMMA_CONTROL_UPPER + 1.e-5)    
            dirichlet = false;
    }
  }
  
  if(!strcmp(name,"mu")) {
//       value = 0.;
//   if (faceName == FACE_FOR_CONTROL)
    dirichlet = false;
  }
  
  return dirichlet;
}


void ComputeIntegral(const MultiLevelProblem& ml_prob);

void AssembleLiftRestrProblem(MultiLevelProblem& ml_prob);



int main(int argc, char** args) {

  // ======= Init ========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
  // ======= Files ========================
  Files files; 
        files.CheckIODirectories();
        files.RedirectCout();

  // ======= Quad Rule ========================
  std::string fe_quad_rule("seventh");
  
  // ======= Mesh  ==================
  MultiLevelMesh ml_mesh;
   
  std::string input_file = "square_parametric.med";
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_file;
  const std::string infile = mystream.str();
  const double Lref = 1.;
  ml_mesh.ReadCoarseMesh(infile.c_str(),fe_quad_rule.c_str(),Lref);

  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfUniformLevels = 6;
  unsigned numberOfSelectiveLevels = 0;
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  ml_mesh.EraseCoarseLevels(numberOfUniformLevels - 1);
  ml_mesh.PrintInfo();

    // ======= Solution ========================
  MultiLevelSolution ml_sol(&ml_mesh);

  ml_sol.AddSolution("state", LAGRANGE, FIRST);
  ml_sol.AddSolution("control", LAGRANGE, FIRST);
  ml_sol.AddSolution("adjoint", LAGRANGE, FIRST);
  ml_sol.AddSolution("mu", LAGRANGE, FIRST);  
  ml_sol.AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  ml_sol.AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  const unsigned int fake_time_dep_flag = 2;  //this is needed to be able to use _SolOld
  const std::string act_set_flag_name = "act_flag";
  ml_sol.AddSolution(act_set_flag_name.c_str(), LAGRANGE, FIRST,fake_time_dep_flag);               //this variable is not solution of any eqn, it's just a given field
  
  // ======= Problem ========================
  MultiLevelProblem ml_prob(&ml_sol);
  
  ml_prob.SetFilesHandler(&files);
  ml_prob.SetQuadratureRuleAllGeomElems(fe_quad_rule);
  ml_prob.set_all_abstract_fe();

  // ======= Solution: Initial Conditions ==================
  ml_sol.Initialize("All");    // initialize all variables to zero

  ml_sol.Initialize("state",       Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize("control",     Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize("adjoint",     Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize("mu",          Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize("TargReg",     Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize("ContReg",     Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize(act_set_flag_name.c_str(), Solution_set_initial_conditions, & ml_prob);

  // ======= Solution: Boundary Conditions ==================
  ml_sol.AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions);
  ml_sol.GenerateBdc("state");
  ml_sol.GenerateBdc("control");
  ml_sol.GenerateBdc("adjoint");
  ml_sol.GenerateBdc("mu");


  // ======= System ==================
  NonLinearImplicitSystemWithPrimalDualActiveSetMethod& system = ml_prob.add_system < NonLinearImplicitSystemWithPrimalDualActiveSetMethod > ("LiftRestr");
  
  system.SetActiveSetFlagName(act_set_flag_name);

  system.AddSolutionToSystemPDE("state");  
  system.AddSolutionToSystemPDE("control");  
  system.AddSolutionToSystemPDE("adjoint");  
  system.AddSolutionToSystemPDE("mu");  
  
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleLiftRestrProblem);
  
  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);
  
  system.SetDebugNonlinear(true);
  system.SetDebugFunction(ComputeIntegral);
//   system.SetMaxNumberOfNonLinearIterations(2);

  // initialize and solve the system
  system.init();
  system.MGsolve();
  
  // ======= Print ==================
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
  ml_sol.GetWriter()->Write(files.GetOutputPath()/*DEFAULT_OUTPUTDIR*/, "biquadratic", variablesToBePrinted);

  return 0;
}


void AssembleLiftRestrProblem(MultiLevelProblem& ml_prob) {

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  
  //=============== Geometry ========================================
  unsigned solType_coords = FE_DOMAIN; //we do linear FE this time
 
  CurrentElem < double > geom_element(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
   //***************************************************  
  
 //***************************************************  
  double weight; // gauss point weight
  

 //********************* state *********************** 
 //***************************************************  
  vector <double> phi_u;  // local test function
  vector <double> phi_u_x; // local test function first order partial derivatives
  vector <double> phi_u_xx; // local test function second order partial derivatives

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * dim);
  phi_u_xx.reserve(maxSize * dim2);
  
 
  unsigned solIndex_u;
  solIndex_u = ml_sol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u = ml_sol->GetSolutionType(solIndex_u);    // get the finite element type for "state"

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
  solIndex_ctrl = ml_sol->GetIndex("control");
  unsigned solType_ctrl = ml_sol->GetSolutionType(solIndex_ctrl);

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
  solIndex_adj = ml_sol->GetIndex("adjoint");    // get the position of "adjoint" in the ml_sol object
  unsigned solType_adj = ml_sol->GetSolutionType(solIndex_adj);    // get the finite element type for "adjoint"

  unsigned solPdeIndex_adj;
  solPdeIndex_adj = mlPdeSys->GetSolPdeIndex("adjoint");    // get the position of "adjoint" in the pdeSys object

  vector < double >  sol_adj; // local solution
    sol_adj.reserve(maxSize);
  vector< int > l2GMap_adj;
    l2GMap_adj.reserve(maxSize);
 //***************************************************  
 //***************************************************  

 //****************** mu ******************************  
 //***************************************************  
  vector <double> phi_mu;  // local test function
  vector <double> phi_mu_x; // local test function first order partial derivatives
  vector <double> phi_mu_xx; // local test function second order partial derivatives

  phi_mu.reserve(maxSize);
  phi_mu_x.reserve(maxSize * dim);
  phi_mu_xx.reserve(maxSize * dim2);
    
  unsigned solIndex_mu;
  solIndex_mu = ml_sol->GetIndex("mu");    // get the position of "mu" in the ml_sol object
   
  unsigned solPdeIndex_mu;
  solPdeIndex_mu = mlPdeSys->GetSolPdeIndex("mu");
  
  unsigned solType_mu = ml_sol->GetSolutionType(solIndex_mu);    // get the finite element type for "mu"
  vector < double >  sol_mu;   sol_mu.reserve(maxSize);
  vector < int > l2GMap_mu;   l2GMap_mu.reserve(maxSize);

  
 //************** act flag **************************** 
   std::string act_flag_name = "act_flag";
   unsigned int solIndex_act_flag = ml_sol->GetIndex(act_flag_name.c_str());
   unsigned int solFEType_act_flag = ml_sol->GetSolutionType(solIndex_act_flag); 
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
 //********* WHOLE SET OF VARIABLES ****************** 

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
    SolIndex[ivar]	= ml_sol->GetIndex        (Solname[ivar].c_str());
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);
  }

    vector < unsigned > Sol_n_el_dofs(n_unknowns);
 //***************************************************  

  
 //********************* DATA ************************ 
  double u_des = DesiredTarget();
  double alpha = ALPHA_CTRL_VOL;
  double beta  = BETA_CTRL_VOL;
  double penalty_strong = 1.e50;
 //***************************************************  

  RES->zero();
  if (assembleMatrix)  KK->zero();

    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    geom_element.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    const short unsigned ielGeom = geom_element.geom_type();
    
   geom_element.set_elem_center(iel, solType_coords);

  //************* set target domain flag **************
   int target_flag = 0;
   target_flag = ElementTargetFlag(geom_element.get_elem_center());
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
    
    
//   update_active_set_flag_for_current_nonlinear_iteration
//          (msh, sol, iel, coords_at_dofs, sol_eldofs, Sol_n_el_dofs, pos_mu, pos_ctrl, c_compl, ctrl_lower, ctrl_upper, sol_actflag, solFEType_act_flag, solIndex_act_flag);
    
    
 //************** update active set flag for current nonlinear iteration **************************** 
 // 0: inactive; 1: active_a; 2: active_b
   assert(nDof_mu == nDof_ctrl);
   sol_actflag.resize(nDof_mu);
   ctrl_lower.resize(nDof_mu);
   ctrl_upper.resize(nDof_mu);
     std::fill(sol_actflag.begin(), sol_actflag.end(), 0);
     std::fill(ctrl_lower.begin(), ctrl_lower.end(), 0.);
     std::fill(ctrl_upper.begin(), ctrl_upper.end(), 0.);
   
    for (unsigned i = 0; i < sol_actflag.size(); i++) {
        std::vector<double> node_coords_i(dim,0.);
        for (unsigned d = 0; d < dim; d++) node_coords_i[d] = geom_element.get_coords_at_dofs()[d][i];
        ctrl_lower[i] = InequalityConstraint(node_coords_i,false);
        ctrl_upper[i] = InequalityConstraint(node_coords_i,true);

        if      ( (sol_mu[i] + c_compl * (sol_ctrl[i] - ctrl_lower[i] )) < 0 )  sol_actflag[i] = 1;
        else if ( (sol_mu[i] + c_compl * (sol_ctrl[i] - ctrl_upper[i] )) > 0 )  sol_actflag[i] = 2;
    }

 //************** act flag **************************** 
    unsigned nDof_act_flag  = msh->GetElementDofNumber(iel, solFEType_act_flag);    // number of solution element dofs
    
    for (unsigned i = 0; i < nDof_act_flag; i++) {
      unsigned solDof_mu = msh->GetSolutionDof(i, iel, solFEType_act_flag); 
      (sol->_Sol[solIndex_act_flag])->set(solDof_mu,sol_actflag[i]);     
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
    
 //*************************************************** 
    
 //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = ControlDomainFlag_internal_restriction(geom_element.get_elem_center());
  std::vector<int> control_node_flag(nDof_ctrl,0);
  if (control_el_flag == 1) std::fill(control_node_flag.begin(), control_node_flag.end(), 1);
 //*************************************************** 
  
 //========= gauss value quantities ==================   
	double sol_u_gss = 0.;
	double sol_adj_gss = 0.;
	double sol_ctrl_gss = 0.;
	double sol_mu_gss = 0.;
	std::vector<double> sol_u_x_gss(dim);       std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::vector<double> sol_adj_x_gss(dim);     std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	std::vector<double> sol_ctrl_x_gss(dim);    std::fill(sol_ctrl_x_gss.begin(), sol_ctrl_x_gss.end(), 0.);
	std::vector<double> sol_mu_x_gss(dim);      std::fill(sol_mu_x_gss.begin(), sol_mu_x_gss.end(), 0.);
 //===================================================   

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
	msh->_finiteElement[ielGeom][solType_u]   ->Jacobian(geom_element.get_coords_at_dofs_3d(), ig, weight, phi_u, phi_u_x, phi_u_xx);
    msh->_finiteElement[ielGeom][solType_ctrl]->Jacobian(geom_element.get_coords_at_dofs_3d(), ig, weight, phi_ctrl, phi_ctrl_x, phi_ctrl_xx);
    msh->_finiteElement[ielGeom][solType_adj] ->Jacobian(geom_element.get_coords_at_dofs_3d(), ig, weight, phi_adj, phi_adj_x, phi_adj_xx);
	msh->_finiteElement[ielGeom][solType_mu]  ->Jacobian(geom_element.get_coords_at_dofs_3d(), ig, weight, phi_mu, phi_mu_x, phi_mu_xx);

	
    sol_u_gss = 0.;
	sol_adj_gss = 0.;
	sol_ctrl_gss = 0.;
	sol_mu_gss = 0.;
    
	std::fill(sol_u_x_gss.begin(),sol_u_x_gss.end(), 0.);
	std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	std::fill(sol_ctrl_x_gss.begin(), sol_ctrl_x_gss.end(), 0.);
	std::fill(sol_mu_x_gss.begin(), sol_mu_x_gss.end(), 0.);
	
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
        
        for (unsigned i = 0; i < nDof_mu; i++) {
	                                                sol_mu_gss      += sol_mu[i] * phi_mu[i];
		   for (unsigned d = 0; d < dim; d++)   sol_mu_x_gss[d] += sol_mu[i] * phi_mu_x[i * dim + d];
          
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
	  if (i < nDof_u)                      Res[0      + i] += - weight * (target_flag * phi_u[i] * ( sol_u_gss + sol_ctrl_gss - u_des) - laplace_rhs_du_adj_i );
          // SECOND ROW
	  if (i < nDof_ctrl)  {
	     if ( control_el_flag == 1)        Res[nDof_u + i] +=  /*(control_node_flag[i]) **/ - weight * (target_flag * phi_ctrl[i] * ( sol_u_gss + sol_ctrl_gss - u_des) 
													      + alpha * phi_ctrl[i] * sol_ctrl_gss
		                                                                                              - laplace_rhs_dctrl_adj_i 
		                                                                                              + beta * laplace_rhs_dctrl_ctrl_i
													      /*+ ineq_flag * sol_mu_gss*/ );
	      else if ( control_el_flag == 0)  Res[nDof_u + i] +=  /*(1 - control_node_flag[i]) **/ (- penalty_strong) * (sol_ctrl[i]);
	  }
          // THIRD ROW
          if (i < nDof_adj)        Res[nDof_u + nDof_ctrl + i] += /*-weight * phi_adj[i] * sol_adj_gss - 6.;*/- weight *  ( - laplace_rhs_dadj_u_i - laplace_rhs_dadj_ctrl_i ) ;

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
              //DIAG BLOCK delta_state - state
	      if ( i < nDof_u && j < nDof_u )       
		Jac[ (0 + i) * nDof_AllVars   +
		     (0 + j)                            ]  += weight * target_flag * phi_u[j] *  phi_u[i];
              
	      // BLOCK  delta_state - control
              if ( i < nDof_u && j < nDof_ctrl )   
		Jac[ (0 + i) * nDof_AllVars   +
                     (nDof_u + j)                       ]  += weight * target_flag  * phi_ctrl[j] *  phi_u[i];
	      
              // BLOCK  delta_state - adjoint
              if ( i < nDof_u && j < nDof_adj )  
		Jac[  (0 + i) * nDof_AllVars  +
                      (nDof_u + nDof_ctrl + j)          ]  += weight * (-1) * laplace_mat_du_adj;
              
	      //=========== delta_control row ===========================     
	      if ( control_el_flag == 1)  {

	      //BLOCK delta_control - state
              if ( i < nDof_ctrl   && j < nDof_u   ) 
		Jac[ (nDof_u + i) * nDof_AllVars  +
		     (0 + j)                            ]  += ( control_node_flag[i]) * weight * target_flag * phi_u[j] * phi_ctrl[i];
		
	      //BLOCK delta_control - control
              if ( i < nDof_ctrl   && j < nDof_ctrl   )
		Jac[ (nDof_u + i) * nDof_AllVars +
		     (nDof_u + j)                       ]  += ( control_node_flag[i]) * weight * ( beta * control_el_flag  * laplace_mat_dctrl_ctrl 
		                                                                                + alpha * control_el_flag * phi_ctrl[i] * phi_ctrl[j] 
		                                                                                            + target_flag * phi_ctrl[i] * phi_ctrl[j] );
              
	      //BLOCK delta_control - adjoint
              if ( i < nDof_ctrl   && j < nDof_adj  ) 
		Jac[ (nDof_u + i) * nDof_AllVars  + 
		     (nDof_u + nDof_ctrl + j)           ]  += ( control_node_flag[i]) * weight * (-1) * laplace_mat_dctrl_adj;
	      	      
	        }
	      
	      else if ( control_el_flag == 0)  {  
		
              //BLOCK delta_control - control
              if ( i < nDof_ctrl   && j < nDof_ctrl &&  i==j ) {
		 Jac[ (nDof_u + i) * nDof_AllVars +
		      (nDof_u + j)                      ]  += (1-control_node_flag[i]) * penalty_strong;
		}
	      
	      }
	      
// 	      //BLOCK delta_control - mu
//              if ( i < nDof_ctrl   && j < nDof_mu && i==j ) 
// 		Jac[ (nDof_u + i) * nDof_AllVars  + 
// 		     (nDof_u + nDof_ctrl + nDof_adj + j)]   =  ineq_flag * control_node_flag[i] *  phi_mu[j] * phi_ctrl[i] ;
		     
		     
	      //=========== delta_adjoint row ===========================
              // BLOCK delta_adjoint - state	      
              if ( i < nDof_adj && j < nDof_u )   
		Jac[ (nDof_u + nDof_ctrl + i) * nDof_AllVars +
		     (0 + j)                            ]  += weight * (-1) * laplace_mat_dadj_u;   
	      
              // BLOCK delta_adjoint - control   
              if ( i < nDof_adj && j < nDof_ctrl )  
		Jac[ (nDof_u + nDof_ctrl + i)  * nDof_AllVars +
		     (nDof_u  + j)                      ]  += weight * (-1) * laplace_mat_dadj_ctrl; 
		     
		     
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
    //std::vector<double> Res_ctrl (nDof_ctrl); std::fill(Res_ctrl.begin(),Res_ctrl.end(), 0.);
    for (unsigned i = 0; i < sol_ctrl.size(); i++){
       unsigned n_els_that_node = 1;
     if ( control_el_flag == 1) {
// 	Res[nDof_u + i] += - ( + n_els_that_node * ineq_flag * sol_mu[i] /*- ( 0.4 + sin(M_PI * x[0][i]) * sin(M_PI * x[1][i]) )*/ );
// 	Res_ctrl[i] =  Res[nDof_u + i];
      }
    }//------------------------------->>>>>>
    
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
    

 //========== end of integral-based part

    //copy the value of the adept::adoube aRes in double Res and store
    RES->add_vector_blocked(Res, l2GMap_AllVars);
      if (assembleMatrix)  KK->add_matrix_blocked(Jac, l2GMap_AllVars, l2GMap_AllVars);
    
    
 //========== dof-based part, without summation
 
 //============= delta_mu row ===============================
      std::vector<double> Res_mu (nDof_mu); std::fill(Res_mu.begin(),Res_mu.end(), 0.);
      
    for (unsigned i = 0; i < sol_actflag.size(); i++) {
      if (sol_actflag[i] == 0){  //inactive
         Res_mu [i] = - ineq_flag * ( 1. * sol_mu[i] - 0. ); 
// 	 Res_mu [i] = Res[nDof_u + nDof_ctrl + nDof_adj + i]; 
      }
      else if (sol_actflag[i] == 1){  //active_a 
	 Res_mu [i] = - ineq_flag * ( c_compl *  sol_ctrl[i] - c_compl * ctrl_lower[i]);
      }
      else if (sol_actflag[i] == 2){  //active_b 
	Res_mu [i]  =  - ineq_flag * ( c_compl *  sol_ctrl[i] - c_compl * ctrl_upper[i]);
      }
    }
//          Res[nDof_u + nDof_ctrl + nDof_adj + i]  = c_compl * (  (2 - sol_actflag[i]) * (ctrl_lower[i] - sol_ctrl[i]) + ( sol_actflag[i] - 1 ) * (ctrl_upper[i] - sol_ctrl[i])  ) ;
//          Res_mu [i] = Res[nDof_u + nDof_ctrl + nDof_adj + i] ;

    
    RES->insert(Res_mu, l2GMap_mu);
//     RES->insert(Res_ctrl, l2GMap_ctrl);
//     RES->insert(Res_u, l2GMap_u);
//     RES->insert(Res_adj, l2GMap_adj);
    
//  //============= delta_state-delta_state row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_u, l2GMap_u, 1.);

//  //============= delta_ctrl-delta_ctrl row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_ctrl, l2GMap_ctrl, 1.);
 
//  //============= delta_adj-delta_adj row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_adj, l2GMap_adj, 1.);
  
 //============= delta_ctrl-delta_mu row ===============================
 KK->matrix_set_off_diagonal_values_blocked(l2GMap_ctrl, l2GMap_mu, ineq_flag * 1.);//------------------------------->>>>>>
  
 //============= delta_mu-delta_ctrl row ===============================
 for (unsigned i = 0; i < sol_actflag.size(); i++) if (sol_actflag[i] != 0 ) sol_actflag[i] = ineq_flag * c_compl;    
  
 KK->matrix_set_off_diagonal_values_blocked(l2GMap_mu, l2GMap_ctrl, sol_actflag);

 //============= delta_mu-delta_mu row ===============================
  for (unsigned i = 0; i < sol_actflag.size(); i++) sol_actflag[i] =  ineq_flag * (1 - sol_actflag[i]/c_compl)  + (1-ineq_flag) * 1.;  //can do better to avoid division, maybe use modulo operator 

  KK->matrix_set_off_diagonal_values_blocked(l2GMap_mu, l2GMap_mu, sol_actflag );
  
//     assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs, 10, 5);
//     assemble_jacobian<double,double>::print_element_jacobian(iel, Jac, Sol_n_el_dofs, 10, 5);
  
  } //end element loop for each process
  
  RES->close();

  if (assembleMatrix) KK->close();
 std::ostringstream mat_out; mat_out << "matrix" << mlPdeSys->GetNonlinearIt()  << ".txt";
  KK->print_matlab(mat_out.str(),"ascii"); //  KK->print();
  
  // ***************** END ASSEMBLY *******************
  unsigned int ctrl_index = mlPdeSys->GetSolPdeIndex("control");
  unsigned int mu_index = mlPdeSys->GetSolPdeIndex("mu");

  unsigned int global_ctrl_size = pdeSys->KKoffset[ctrl_index+1][iproc] - pdeSys->KKoffset[ctrl_index][iproc];
  
  std::vector<double>  one_times_mu(global_ctrl_size, 0.);
  std::vector<int>    positions(global_ctrl_size);
//  double position_mu_i;
  for (unsigned i = 0; i < positions.size(); i++) {
    positions[i] = pdeSys->KKoffset[ctrl_index][iproc] + i;
//     position_mu_i = pdeSys->KKoffset[mu_index][iproc] + i;
//     std::cout << position_mu_i << std::endl;
    one_times_mu[i] = ineq_flag * 1. * (*sol->_Sol[solIndex_mu])(i/*position_mu_i*/) ;
  }
    RES->add_vector_blocked(one_times_mu, positions);
    RES->print();
    

  return;
}



void ComputeIntegral(const MultiLevelProblem& ml_prob)    {
  
  
  const NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("LiftRestr");
  const unsigned          level      = mlPdeSys->GetLevelToAssemble();

  Mesh*                          msh = ml_prob._ml_msh->GetLevel(level);

  MultiLevelSolution*          ml_sol = ml_prob._ml_sol;
  Solution*                      sol = ml_prob._ml_sol->GetSolutionLevel(level);

  const unsigned     dim = msh->GetDimension();                                 // get the domain dimension of the problem
  const unsigned    dim2 = (3 * (dim - 1) + !(dim - 1));                        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  const unsigned   iproc = msh->processor_id(); 

 //********** Geometry ***************************************** 
 unsigned solType_coords = FE_DOMAIN; //we do linear FE this time
 
  CurrentElem < double > geom_element(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
 //*************************************************** 

 //*************************************************** 
  double weight; // gauss point weight
  
 //***************************************************  
  double alpha = ALPHA_CTRL_VOL;
  double beta  = BETA_CTRL_VOL;

 //******************** state ************************ 
 //*************************************************** 
  vector <double> phi_u;
  vector <double> phi_u_x;
  vector <double> phi_u_xx;

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * dim);
  phi_u_xx.reserve(maxSize * dim2);
  
 
  unsigned solIndex_u;
  solIndex_u = ml_sol->GetIndex("state");
  unsigned solType_u = ml_sol->GetSolutionType(solIndex_u);

  vector < double >  sol_u; // local solution
  sol_u.reserve(maxSize);
  
  double u_gss = 0.;
 //*************************************************** 
 //*************************************************** 

 //******************** control ********************** 
 //*************************************************** 
  vector <double> phi_ctrl;
  vector <double> phi_ctrl_x;
  vector <double> phi_ctrl_xx;

  phi_ctrl.reserve(maxSize);
  phi_ctrl_x.reserve(maxSize * dim);
  phi_ctrl_xx.reserve(maxSize * dim2);
  
  unsigned solIndex_ctrl;
  solIndex_ctrl = ml_sol->GetIndex("control");
  unsigned solType_ctrl = ml_sol->GetSolutionType(solIndex_ctrl);

  vector < double >  sol_ctrl; // local solution
  sol_ctrl.reserve(maxSize);
  
  double ctrl_gss = 0.;
  double ctrl_x_gss = 0.;
 //*************************************************** 
 //***************************************************  

  
 //********************* desired ********************* 
 //*************************************************** 
  vector <double> phi_udes;
  vector <double> phi_udes_x;
  vector <double> phi_udes_xx;

  phi_udes.reserve(maxSize);
  phi_udes_x.reserve(maxSize * dim);
  phi_udes_xx.reserve(maxSize * dim2);
 
  
//  unsigned solIndex_udes;
//   solIndex_udes = ml_sol->GetIndex("Tdes");    // get the position of "state" in the ml_sol object
//   unsigned solType_udes = ml_sol->GetSolutionType(solIndex_udes);    // get the finite element type for "state"

  vector < double >  sol_udes; // local solution
  sol_udes.reserve(maxSize);

  double udes_gss = 0.;
 //*************************************************** 
 //*************************************************** 

 //********************* DATA ************************ 
  double u_des = DesiredTarget();
 //*************************************************** 
  
  double integral_target = 0.;
  double integral_alpha  = 0.;
  double integral_beta   = 0.;

    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
 
    geom_element.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    const short unsigned ielGeom = geom_element.geom_type();

  //************* set target domain flag **************
   geom_element.set_elem_center(iel, solType_coords);

   int target_flag = 0;
   target_flag = ElementTargetFlag(geom_element.get_elem_center());
 //*************************************************** 

 //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = ControlDomainFlag_internal_restriction(geom_element.get_elem_center());
 //*************************************************** 
   
 //**************** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);
        sol_u    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);
            sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);
    }
 //*************************************************** 


 //************** control **************************** 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);
    sol_ctrl    .resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);
    } 
 //***************************************************  
 
 
 //**************** u_des **************************** 
    unsigned nDof_udes  = msh->GetElementDofNumber(iel, solType_u);
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
      for (unsigned ig = 0; ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
	    msh->_finiteElement[ielGeom][solType_u]               ->Jacobian(geom_element.get_coords_at_dofs(), ig, weight, phi_u, phi_u_x, phi_u_xx);
        msh->_finiteElement[ielGeom][solType_u/*solTypeudes*/]->Jacobian(geom_element.get_coords_at_dofs(), ig, weight, phi_udes, phi_udes_x, phi_udes_xx);
        msh->_finiteElement[ielGeom][solType_ctrl]            ->Jacobian(geom_element.get_coords_at_dofs(), ig, weight, phi_ctrl, phi_ctrl_x, phi_ctrl_xx);

	u_gss = 0.;  for (unsigned i = 0; i < nDof_u; i++) u_gss += sol_u[i] * phi_u[i];		
	ctrl_gss = 0.; for (unsigned i = 0; i < nDof_ctrl; i++) ctrl_gss += sol_ctrl[i] * phi_ctrl[i];  
	udes_gss  = 0.; for (unsigned i = 0; i < nDof_udes; i++)  udes_gss  += sol_udes[i]  * phi_udes[i]; 
        ctrl_x_gss  = 0.; for (unsigned i = 0; i < nDof_ctrl; i++)  
        {
          for (unsigned idim = 0; idim < dim; idim ++) ctrl_x_gss  += sol_ctrl[i] * phi_ctrl_x[i + idim * nDof_ctrl];
        }

               integral_target += target_flag * weight * (u_gss +  ctrl_gss - udes_gss) * (u_gss +  ctrl_gss - udes_gss);
               integral_alpha  += control_el_flag * weight * ctrl_gss * ctrl_gss;
               integral_beta   += control_el_flag * weight * ctrl_x_gss * ctrl_x_gss;
	  
      } // end gauss point loop
  } //end element loop

  std::cout << "The value of the integral_target is " << std::setw(11) << std::setprecision(10) << integral_target << std::endl;
  std::cout << "The value of the integral_alpha  is " << std::setw(11) << std::setprecision(10) << integral_alpha << std::endl;
  std::cout << "The value of the integral_beta   is " << std::setw(11) << std::setprecision(10) << integral_beta << std::endl;
  std::cout << "The value of the total integral  is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_target + 0.5 * alpha * integral_alpha + 0.5 * beta * integral_beta << std::endl;

return;
  
}
  
  

