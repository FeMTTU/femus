#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "NumericVector.hpp"
#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"

#define FACE_FOR_CONTROL   2
#define FACE_FOR_TARGET    2

#include "../../../param.hpp"

#define FE_DOMAIN  2

using namespace femus;


/// @todo here I have to do the separation of add and insert parts for the parallel, like I did in the other app



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
  const bool use_output_time_folder = false;
  const bool redirect_cout_to_file = true;
  Files files; 
        files.CheckIODirectories(use_output_time_folder);
        files.RedirectCout(redirect_cout_to_file);

  // ======= Quad Rule ========================
  std::string fe_quad_rule("seventh");
  
  // ======= Mesh  ==================
  MultiLevelMesh ml_mesh;
   
//   std::string input_file = "square_4x5.med";
  std::string input_file = "parametric_square_2x2.med";
//   std::string input_file = "square_parametric.med";
//   std::string input_file = "Mesh_3_groups_with_bdry_nodes.med";
//   std::string input_file = "Mesh_3_groups_with_bdry_nodes_coarser.med";
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_file;
  const std::string infile = mystream.str();
  const double Lref = 1.;
  ml_mesh.ReadCoarseMesh(infile.c_str(),fe_quad_rule.c_str(),Lref);

  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfUniformLevels = N_UNIFORM_LEVELS;
  unsigned numberOfSelectiveLevels = 0;
  const unsigned erased_levels = N_ERASED_LEVELS;
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  ml_mesh.EraseCoarseLevels(erased_levels/*numberOfUniformLevels - 1*/);
  ml_mesh.PrintInfo();

    // ======= Solution ========================
  MultiLevelSolution ml_sol(&ml_mesh);

  ml_sol.AddSolution("state", LAGRANGE, /*SECOND*/FIRST);
  ml_sol.AddSolution("control", LAGRANGE, /*SECOND*/FIRST);
  ml_sol.AddSolution("adjoint", LAGRANGE, /*SECOND*/FIRST);
  ml_sol.AddSolution("mu", LAGRANGE, /*SECOND*/FIRST);  
  ml_sol.AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  ml_sol.AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  const unsigned int fake_time_dep_flag = 2;  //this is needed to be able to use _SolOld
  const std::string act_set_flag_name = "act_flag";
  ml_sol.AddSolution(act_set_flag_name.c_str(), LAGRANGE, /*SECOND*/FIRST, fake_time_dep_flag);               //this variable is not solution of any eqn, it's just a given field
  
  // ======= Problem ========================
  MultiLevelProblem ml_prob(&ml_sol);
  
  ml_prob.SetFilesHandler(&files);
  ml_prob.SetQuadratureRuleAllGeomElems(fe_quad_rule);
  ml_prob.set_all_abstract_fe_multiple();

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
//   system.assemble_call_before_boundary_conditions(1);
  
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

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             KK = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

  unsigned    iproc = msh->processor_id();

  constexpr bool print_algebra_global = true;
  constexpr bool print_algebra_local = true;
  
  //=============== Geometry ========================================
  unsigned solType_coords = FE_DOMAIN;
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
   //***************************************************  
  
 //***************************************************  
  double weight; // gauss point weight
  

 //********************* state *********************** 
 //***************************************************  
  vector <double> phi_u; 
  vector <double> phi_u_x;
  vector <double> phi_u_xx;

  phi_u.reserve(max_size);
  phi_u_x.reserve(max_size * dim);
  phi_u_xx.reserve(max_size * dim2);
  
 
  unsigned solIndex_u;
  solIndex_u = ml_sol->GetIndex("state");
  unsigned solType_u = ml_sol->GetSolutionType(solIndex_u);

  unsigned solPdeIndex_u;
  solPdeIndex_u = mlPdeSys->GetSolPdeIndex("state");

  vector < double >  sol_u;
  sol_u.reserve(max_size);
  vector< int > l2GMap_u;
  l2GMap_u.reserve(max_size);
 //***************************************************  
 //***************************************************  

  
 //******************** control ********************** 
 //***************************************************   
  vector <double> phi_ctrl;
  vector <double> phi_ctrl_x;
  vector <double> phi_ctrl_xx;

  phi_ctrl.reserve(max_size);
  phi_ctrl_x.reserve(max_size * dim);
  phi_ctrl_xx.reserve(max_size * dim2);
  
  unsigned solIndex_ctrl;
  solIndex_ctrl = ml_sol->GetIndex("control");
  unsigned solType_ctrl = ml_sol->GetSolutionType(solIndex_ctrl);

  unsigned solPdeIndex_ctrl;
  solPdeIndex_ctrl = mlPdeSys->GetSolPdeIndex("control");

  vector < double >  sol_ctrl;
  sol_ctrl.reserve(max_size);
  vector< int > l2GMap_ctrl;
  l2GMap_ctrl.reserve(max_size);
 //***************************************************  
 //***************************************************  
  
  
 //********************* adjoint ********************* 
 //***************************************************  
  vector <double> phi_adj;
  vector <double> phi_adj_x;
  vector <double> phi_adj_xx;

  phi_adj.reserve(max_size);
  phi_adj_x.reserve(max_size * dim);
  phi_adj_xx.reserve(max_size * dim2);
 
  
  unsigned solIndex_adj;
  solIndex_adj = ml_sol->GetIndex("adjoint");
  unsigned solType_adj = ml_sol->GetSolutionType(solIndex_adj);

  unsigned solPdeIndex_adj;
  solPdeIndex_adj = mlPdeSys->GetSolPdeIndex("adjoint");

  vector < double >  sol_adj;
    sol_adj.reserve(max_size);
  vector< int > l2GMap_adj;
    l2GMap_adj.reserve(max_size);
 //***************************************************  
 //***************************************************  

 //****************** mu ******************************  
 //***************************************************  
  vector <double> phi_mu;
  vector <double> phi_mu_x;
  vector <double> phi_mu_xx;

  phi_mu.reserve(max_size);
  phi_mu_x.reserve(max_size * dim);
  phi_mu_xx.reserve(max_size * dim2);
    
  unsigned solIndex_mu;
  solIndex_mu = ml_sol->GetIndex("mu");
   
  unsigned solPdeIndex_mu;
  solPdeIndex_mu = mlPdeSys->GetSolPdeIndex("mu");
  
  unsigned solType_mu = ml_sol->GetSolutionType(solIndex_mu);
  vector < double >  sol_mu;   sol_mu.reserve(max_size);
  vector < int > l2GMap_mu;   l2GMap_mu.reserve(max_size);

  
  //************** act flag ****************************   
  unsigned int solIndex_act_flag_sol; 
  unsigned int solFEType_act_flag_sol;
  store_act_flag_in_old(mlPdeSys, ml_sol, sol,
                        solIndex_act_flag_sol, //this becomes a vector
                        solFEType_act_flag_sol //remove this one, only Index
                       );
    
  
  
  //********* variables for ineq constraints *****************
  const int ineq_flag = INEQ_FLAG;
  const double c_compl = C_COMPL;
  vector < double/*int*/ >  sol_actflag;   sol_actflag.reserve(max_size);
  vector < double >  ctrl_lower;   ctrl_lower.reserve(max_size);
  vector < double >  ctrl_upper;   ctrl_upper.reserve(max_size);
  //***************************************************  

 //***************************************************  
 //********* WHOLE SET OF VARIABLES ****************** 

  const int n_unknowns =  mlPdeSys->GetSolPdeIndex().size();
 
  vector< int > l2GMap_AllVars; // local to global mapping
  l2GMap_AllVars.reserve(n_unknowns*max_size);
  
  vector< double > Res;
  Res.reserve(n_unknowns*max_size);

  vector < double > Jac;
  Jac.reserve( n_unknowns*max_size * n_unknowns*max_size);
  

//***************************************************
    enum Pos_in_matrix {pos_mat_state = 0, pos_mat_ctrl, pos_mat_adj, pos_mat_mu}; //these are known at compile-time 
                    ///@todo these are the positions in the MlSol object or in the Matrix? I'd say the matrix, but we have to check where we use it...

                    
    assert(pos_mat_state   == mlPdeSys->GetSolPdeIndex("state"));
    assert(pos_mat_ctrl    == mlPdeSys->GetSolPdeIndex("control"));
    assert(pos_mat_adj     == mlPdeSys->GetSolPdeIndex("adjoint"));
    assert(pos_mat_mu      == mlPdeSys->GetSolPdeIndex("mu"));
//***************************************************

  vector < std::string > Solname(n_unknowns);
  Solname[0] = "state";
  Solname[1] = "control";
  Solname[2] = "adjoint";
  Solname[3] = "mu";
  

 //***************************************************
   enum Pos_in_Sol {pos_sol_state = 0, pos_sol_ctrl, pos_sol_adj, pos_sol_mu, pos_sol_targreg, pos_sol_contreg, pos_sol_actflag}; //these are known at compile-time 

        assert(pos_sol_actflag == solIndex_act_flag_sol);
//***************************************************
    
  


  vector < unsigned > SolPdeIndex(n_unknowns);
  vector < unsigned > SolIndex(n_unknowns);  
  vector < unsigned > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolPdeIndex[ivar]	= mlPdeSys->GetSolPdeIndex(Solname[ivar].c_str());
    SolIndex[ivar]	= ml_sol->GetIndex        (Solname[ivar].c_str());
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);
  }

    vector < unsigned > Sol_n_el_dofs_Mat_vol(n_unknowns);
 //***************************************************  

//***************************************************
    vector < vector < double > >  sol_eldofs_Mat(n_unknowns);  //should have Mat order
    for(int k = 0; k < n_unknowns; k++) {        sol_eldofs_Mat[k].reserve(max_size);    }


    //----------- quantities (at dof objects) ------------------------------
    std::vector< int >       L2G_dofmap_Mat_AllVars;
    L2G_dofmap_Mat_AllVars.reserve( n_unknowns * max_size );
    vector < vector < int > >     L2G_dofmap_Mat(n_unknowns);     //should have Mat order
    for(int i = 0; i < n_unknowns; i++) {
        L2G_dofmap_Mat[i].reserve(max_size);
    }

    
    
 //********************* DATA ************************ 
  double u_des = DesiredTarget();
  double alpha = ALPHA_CTRL_VOL;
  double beta  = BETA_CTRL_VOL;
  double penalty_outside_control_domain = 1.e50;         // penalty for zero control outside
 //***************************************************  

  RES->zero();
  if (assembleMatrix)  KK->zero();

    

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    const short unsigned ielGeom = geom_element_iel.geom_type();
    
   geom_element_iel.set_elem_center_3d(iel, solType_coords);

  //************* set target domain flag **************
   int target_flag = 0;
   target_flag = ElementTargetFlag(geom_element_iel.get_elem_center_3d());
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
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);
    sol_ctrl    .resize(nDof_ctrl);
    l2GMap_ctrl.resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);
      l2GMap_ctrl[i] = pdeSys->GetSystemDof(solIndex_ctrl, solPdeIndex_ctrl, i, iel);
    } 
 //*************************************************** 
 

 //************** adjoint **************************** 
    unsigned nDof_adj  = msh->GetElementDofNumber(iel, solType_adj);
        sol_adj    .resize(nDof_adj);
        l2GMap_adj.resize(nDof_adj);
    for (unsigned i = 0; i < sol_adj.size(); i++) {
      unsigned solDof_adj = msh->GetSolutionDof(i, iel, solType_adj);
      sol_adj[i] = (*sol->_Sol[solIndex_adj])(solDof_adj);
      l2GMap_adj[i] = pdeSys->GetSystemDof(solIndex_adj, solPdeIndex_adj, i, iel);
    } 
 //***************************************************  
 
 //************** mu **************************** 
    unsigned nDof_mu  = msh->GetElementDofNumber(iel, solType_mu);
    sol_mu   .resize(nDof_mu);
    l2GMap_mu.resize(nDof_mu);
    for (unsigned i = 0; i < sol_mu.size(); i++) {
      unsigned solDof_mu = msh->GetSolutionDof(i, iel, solType_mu);
      sol_mu[i] = (*sol->_Sol[solIndex_mu])(solDof_mu);
      l2GMap_mu[i] = pdeSys->GetSystemDof(solIndex_mu, solPdeIndex_mu, i, iel);
    }
 //***************************************************  
    
 
    
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
	Sol_n_el_dofs_Mat_vol[k] = ndofs_unk;
  }
    
 //*************************************************** 
    
 //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = ControlDomainFlag_internal_restriction(geom_element_iel.get_elem_center_3d());
  std::vector<int> control_node_flag(nDof_ctrl, 0);
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
	msh->_finiteElement[ielGeom][solType_u]   ->Jacobian(geom_element_iel.get_coords_at_dofs_3d(), ig, weight, phi_u, phi_u_x, phi_u_xx);
    msh->_finiteElement[ielGeom][solType_ctrl]->Jacobian(geom_element_iel.get_coords_at_dofs_3d(), ig, weight, phi_ctrl, phi_ctrl_x, phi_ctrl_xx);
    msh->_finiteElement[ielGeom][solType_adj] ->Jacobian(geom_element_iel.get_coords_at_dofs_3d(), ig, weight, phi_adj, phi_adj_x, phi_adj_xx);
	msh->_finiteElement[ielGeom][solType_mu]  ->Jacobian(geom_element_iel.get_coords_at_dofs_3d(), ig, weight, phi_mu, phi_mu_x, phi_mu_xx);

	
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
	     if ( control_el_flag == 1)    {    Res[nDof_u + i] +=  /*(control_node_flag[i]) **/ - weight * (target_flag * phi_ctrl[i] * ( sol_u_gss + sol_ctrl_gss - u_des) 
													      + alpha * phi_ctrl[i] * sol_ctrl_gss
		                                                                                              - laplace_rhs_dctrl_adj_i 
		                                                                                              + beta * laplace_rhs_dctrl_ctrl_i
													      /*+ ineq_flag * sol_mu_gss*/ ); }
	      else if ( control_el_flag == 0) { Res[nDof_u + i] +=   (- penalty_outside_control_domain)  *  (1 - control_node_flag[i]) * (sol_ctrl[i] - 0.); }
	  }
          // THIRD ROW
          if (i < nDof_adj)        Res[nDof_u + nDof_ctrl + i] += /*-weight * phi_adj[i] * sol_adj_gss - 6.;*/- weight *  ( - laplace_rhs_dadj_u_i - laplace_rhs_dadj_ctrl_i ) ;

       if (i < nDof_mu)  {
            if ( control_el_flag == 0) {  Res[nDof_u + nDof_ctrl + nDof_adj + i] +=   (- penalty_outside_control_domain) *  (1 - control_node_flag[i]) * (sol_mu[i] - 0.); }
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
                      (nDof_u + nDof_ctrl + j)          ]  += weight * (-1.) * laplace_mat_du_adj;
              
	      //=========== delta_control row ===========================     
	      if ( control_el_flag == 1)  {

	      //BLOCK delta_control - state
              if ( i < nDof_ctrl   && j < nDof_u   ) 
		Jac[ (nDof_u + i) * nDof_AllVars  +
		     (0 + j)                            ]  += ( control_node_flag[i]) * weight * target_flag * phi_u[j] * phi_ctrl[i];
		
	      //BLOCK delta_control - control
              if ( i < nDof_ctrl   && j < nDof_ctrl   )
		Jac[ (nDof_u + i) * nDof_AllVars +
		     (nDof_u + j)                       ]  += ( control_node_flag[i]) * weight * ( beta * /*control_el_flag  **/ laplace_mat_dctrl_ctrl 
		                                                                                + alpha * /*control_el_flag **/ phi_ctrl[i] * phi_ctrl[j] 
		                                                                                            + target_flag * phi_ctrl[i] * phi_ctrl[j] );
              
	      //BLOCK delta_control - adjoint
              if ( i < nDof_ctrl   && j < nDof_adj  ) 
		Jac[ (nDof_u + i) * nDof_AllVars  + 
		     (nDof_u + nDof_ctrl + j)           ]  += ( control_node_flag[i]) * weight * (-1.) * laplace_mat_dctrl_adj;
	      	      
	        }
	      
	      else if ( control_el_flag == 0)  {  
		
              //BLOCK delta_control - control
              if ( i < nDof_ctrl   && j < nDof_ctrl &&  i==j ) {
		 Jac[ (nDof_u + i) * nDof_AllVars +
		      (nDof_u + j)                      ]  +=  penalty_outside_control_domain * (1 - control_node_flag[i]);
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
                if ( i < nDof_mu && j < nDof_mu && i==j )  {   
	        if ( control_el_flag == 0)  {  
		  Jac[   (nDof_u + nDof_ctrl + nDof_adj + i) * nDof_AllVars +  (nDof_u + nDof_ctrl + nDof_adj + j) ]  += penalty_outside_control_domain * (1 - control_node_flag[i]);    //MU
                }
	      }
	      
//           if (sol_actflag[i] == 0) //inactive
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
// // 	         for (unsigned k = 0; k < i_unk; k++) row_block_offset += Sol_n_el_dofs_Mat_vol[k];
// //          for (unsigned j_unk = 0; j_unk < n_unknowns; j_unk++) {
// //     std::cout << " ======= Column === " << j_unk << " ================== " << std::endl;     
// //         unsigned int column_block_offset = 0;
// // 	         for (unsigned k = 0; k < j_unk; k++) column_block_offset += Sol_n_el_dofs_Mat_vol[k];
// // 	  
// //          for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[i_unk]; i++) {
// // // 	      std::cout << Res[nDof_u + nDof_ctrl + nDof_adj + i ] << " " << std::endl;
// // 	   for (unsigned j = 0; j < Sol_n_el_dofs_Mat_vol[j_unk]; j++) {
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

      
     if (print_algebra_local) {
         assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat_vol, 10, 5);
         assemble_jacobian<double,double>::print_element_jacobian(iel, Jac, Sol_n_el_dofs_Mat_vol, 10, 5);
     }

       
  } //end element loop for each process
      
  //MU
  
  
// // // 
// // //   unsigned int global_ctrl_size = pdeSys->KKoffset[ctrl_index+1][iproc] - pdeSys->KKoffset[ctrl_index][iproc];
// // //   
// // //   std::vector<double>  one_times_mu(global_ctrl_size, 0.);
// // //   std::vector<int>    positions(global_ctrl_size);
// // // //  double position_mu_i;
// // //   for (unsigned i = 0; i < positions.size(); i++) {
// // //     positions[i] = pdeSys->KKoffset[ctrl_index][iproc] + i;
// // // //     position_mu_i = pdeSys->KKoffset[mu_index][iproc] + i;
// // // //     std::cout << position_mu_i << std::endl;
// // //     one_times_mu[i] = ineq_flag * 1. * (*sol->_Sol[solIndex_mu])(i/*position_mu_i*/) ;
// // //   }
// // //     RES->add_vector_blocked(one_times_mu, positions);
  unsigned int ctrl_index = mlPdeSys->GetSolPdeIndex("control");
  unsigned int mu_index = mlPdeSys->GetSolPdeIndex("mu");
    
add_one_times_mu_res_ctrl(iproc,
                               ineq_flag,
                               ctrl_index,
                               mu_index,
                               SolIndex,
                               sol,
                               mlPdeSys,
                               pdeSys,
                               RES);
    
    
    
     
RES->close();
if (assembleMatrix) KK->close();  /// This is needed for the parallel, when splitting the add part from the insert part!!!
      // ***************** ADD PART - END  *******************
    
    
 


//   ***************** INSERT PART - BEGIN (must go AFTER the sum, clearly) *******************
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

// -------
   geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
      
   geom_element_iel.set_elem_center_3d(iel, solType_coords);
// -------
   
// -------
    el_dofs_unknowns_vol(sol, msh, pdeSys, iel,
                        SolFEType,
                        SolIndex,
                        SolPdeIndex,
                        Sol_n_el_dofs_Mat_vol,
                        sol_eldofs_Mat,
                        L2G_dofmap_Mat);
      
    //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = ControlDomainFlag_internal_restriction(geom_element_iel.get_elem_center_3d());
 
    
    if (control_el_flag == 1) {

        
  update_active_set_flag_for_current_nonlinear_iteration
  (msh,
   sol,
   iel,
   geom_element_iel.get_coords_at_dofs/*_3d*/(),
   sol_eldofs_Mat,
   Sol_n_el_dofs_Mat_vol,
   pos_mat_mu,
   pos_mat_ctrl,
   c_compl,
   ctrl_lower,
   ctrl_upper,
   sol_actflag,
   solFEType_act_flag_sol,
   solIndex_act_flag_sol);
  
      


    node_insertion(iel,
                   msh,
                   L2G_dofmap_Mat,
                   pos_mat_mu,
                   pos_mat_ctrl,
                   sol_eldofs_Mat,
                   Sol_n_el_dofs_Mat_vol,
                   sol_actflag,
                   solFEType_act_flag_sol,
                   ineq_flag,
                   c_compl,
                   ctrl_lower,
                   ctrl_upper,
                   KK,
                   RES,
                   assembleMatrix
                   );
   

       }

 
    
    
 //============= delta_ctrl-delta_mu row ===============================
  if (assembleMatrix) { KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap_Mat[pos_mat_ctrl], L2G_dofmap_Mat[pos_mat_mu], ineq_flag * 1.); }
  
  
  }
//   ***************** INSERT PART - END (must go AFTER the sum, clearly) *******************
  
  RES->close();

  if (assembleMatrix) KK->close();
  
  // ***************** END ASSEMBLY *******************
  
  
    const unsigned nonlin_iter = mlPdeSys->GetNonlinearIt();
  if (print_algebra_global) {
    assemble_jacobian< double, double >::print_global_jacobian(assembleMatrix, ml_prob, KK, nonlin_iter);
//     assemble_jacobian< double, double >::print_global_residual(ml_prob, RES,  mlPdeSys->GetNonlinearIt());

    RES->close();
    std::ostringstream res_out; res_out << ml_prob.GetFilesHandler()->GetOutputPath() << "./" << "res_" << nonlin_iter  << ".txt";
    pdeSys->print_with_structure_matlab_friendly(iproc, res_out.str().c_str(), RES);

  }
  

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
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

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

  phi_u.reserve(max_size);
  phi_u_x.reserve(max_size * dim);
  phi_u_xx.reserve(max_size * dim2);
  
 
  unsigned solIndex_u;
  solIndex_u = ml_sol->GetIndex("state");
  unsigned solType_u = ml_sol->GetSolutionType(solIndex_u);

  vector < double >  sol_u; // local solution
  sol_u.reserve(max_size);
  
  double u_gss = 0.;
 //*************************************************** 
 //*************************************************** 

 //******************** control ********************** 
 //*************************************************** 
  vector <double> phi_ctrl;
  vector <double> phi_ctrl_x;
  vector <double> phi_ctrl_xx;

  phi_ctrl.reserve(max_size);
  phi_ctrl_x.reserve(max_size * dim);
  phi_ctrl_xx.reserve(max_size * dim2);
  
  unsigned solIndex_ctrl;
  solIndex_ctrl = ml_sol->GetIndex("control");
  unsigned solType_ctrl = ml_sol->GetSolutionType(solIndex_ctrl);

  vector < double >  sol_ctrl; // local solution
  sol_ctrl.reserve(max_size);
  
  double ctrl_gss = 0.;
  double ctrl_x_gss = 0.;
 //*************************************************** 
 //***************************************************  

  
 //********************* desired ********************* 
 //*************************************************** 
  vector <double> phi_udes;
  vector <double> phi_udes_x;
  vector <double> phi_udes_xx;

  phi_udes.reserve(max_size);
  phi_udes_x.reserve(max_size * dim);
  phi_udes_xx.reserve(max_size * dim2);
 
  
//  unsigned solIndex_udes;
//   solIndex_udes = ml_sol->GetIndex("Tdes");    // get the position of "state" in the ml_sol object
//   unsigned solType_udes = ml_sol->GetSolutionType(solIndex_udes);    // get the finite element type for "state"

  vector < double >  sol_udes; // local solution
  sol_udes.reserve(max_size);

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
   geom_element.set_elem_center_3d(iel, solType_coords);

   int target_flag = 0;
   target_flag = ElementTargetFlag(geom_element.get_elem_center_3d());
 //*************************************************** 

 //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = ControlDomainFlag_internal_restriction(geom_element.get_elem_center_3d());
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

  ////////////////////////////////////////
  double total_integral = 0.5 * integral_target + 0.5 * alpha * integral_alpha + 0.5 * beta * integral_beta;
  
  std::cout << "total integral on processor " << iproc << ": " << total_integral << std::endl;

  double integral_target_parallel = 0.; MPI_Allreduce( &integral_target, &integral_target_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  double integral_alpha_parallel = 0.; MPI_Allreduce( &integral_alpha, &integral_alpha_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  double integral_beta_parallel = 0.;  MPI_Allreduce( &integral_beta, &integral_beta_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  double total_integral_parallel = 0.; MPI_Allreduce( &total_integral, &total_integral_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  
  std::cout << "The value of the integral_target is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_target_parallel << std::endl;
  std::cout << "The value of the integral_alpha  is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_alpha_parallel << std::endl;
  std::cout << "The value of the integral_beta   is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_beta_parallel << std::endl;
  std::cout << "The value of the total integral  is " << std::setw(11) << std::setprecision(10) << total_integral_parallel << std::endl;

return;
  
}
  
  

